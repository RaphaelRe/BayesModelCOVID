import numpy as np
import copy
import scipy.stats as stats
import pandas as pd

import basics
import update


class LatentVariables:
    '''
    The Latent variable class holds the true number of infections.
    It also holds the derived quantities.

    '''
    def __init__(self, model, data, parameter_values,
                 proposal_distribution='Uniform', proposal_range=-1,
                 start=-99999, nb_blocks=10, fix_latent_variable=False,
                 adapt_reporting_weekend=True
                 ):
        self.model = model
        self.target = 0.2
        self.start = start
        self.data = basics.dictionarize_data(data)
        self.countries = self.data.keys()

        if "hospitalizations" in data.columns:
            self.hospitalization_exist = {country: not pd.isna(self.data[country]["hospitalizations"].iloc[0]) for country in self.countries}
        else:
            self.hospitalization_exist = {country: False for country in self.countries}

        if "intensiveCare" in data.columns:
            self.intensiveCare_exist = {country: not pd.isna(self.data[country]["intensiveCare"].iloc[0]) for country in self.countries}
        else:
            self.intensiveCare_exist = {country: False for country in self.countries}
        self.i = 0

        if self.model == 'infections':
            self.infections = self.initialize_values(data, fix_latent_variable=fix_latent_variable, parameter_values=parameter_values)
            self.acceptance_infections = {country: np.empty(30000) for country in self.countries}
            self.infections_view = {country: self.infections[np.min(np.where(data['country'].values == country)):(np.max(np.where(data['country'].values == country)) + 1)] for country in self.countries}
            self.cases = np.empty_like(self.infections)
            self.cases[:] = -1
            self.cases_view = {country: self.cases[np.min(np.where(data['country'].values == country)):\
                    (np.max(np.where(data['country'].values == country)) + 1)] for country in self.countries}
            for country in self.countries:
                self.cases_view[country][:] = basics.calculate_cases(self.infections_view[country], parameter_values)

        elif self.model == 'cases':
            self.cases = self.initialize_values(data, fix_latent_variable=fix_latent_variable)
            self.cases_view = {country: self.cases[np.min(np.where(data['country'].values == country)):(np.max(np.where(data['country'].values == country)) + 1)] for country in self.countries}


        if adapt_reporting_weekend:
            # import pudb;pu.db
            self.Xi_R = {}
            self.Xi_D = {}
            for country_iter in self.countries:
                self.Xi_R[country_iter] = basics.initialize_country_specific_Xi_R(parameter_values, country_iter)
                self.Xi_D[country_iter] = basics.initialize_country_specific_Xi_D(parameter_values, country_iter)



        # block view is a dict which stores pointers on cases for mor cases
        # model or infections for infections-model
        self.block_view = {country: self.initialize_blocks(country, nb_blocks) for country in self.countries}

        self.proposal_distribution = proposal_distribution  # brauchen wir das noch??

        self.proposal_range = {}
        for country in self.countries:
            self.proposal_range[country] = {block: proposal_range for block in self.block_view[country].keys()}

        self.acceptance_cases = {}  # actual this are infections
        for country in self.countries:
            self.acceptance_cases[country] = {block: np.empty(30000) for block in self.block_view[country].keys()}

        # initialize Rt
        # Rt contains R0, VOCs and alpha, NO correction factors
        self.Rt = basics.calculate_Rt(parameter_values, self.data)
        self.update_sumut(parameter_values, start=start)

        # initialize correction factors
        if 'first_vaccination' in data.columns:
            self.correction_factor1 = basics.calculate_correction_factor1(self.infections_view, parameter_values)
            self.correction_factor2 = basics.calculate_correction_factor2(self.data, parameter_values)
        else:
            self.correction_factor1 = {country_iter: 0 for country_iter in self.data.keys()}
            self.correction_factor2 = {country_iter: 0 for country_iter in self.data.keys()}

        self.likelihood_ratio = lambda \
                data, values_parameters, current_values, candidate_values,\
                country, hospitalization_model, intensiveCare_model:\
                update.ratio_likelihood_nbinom(model, data, values_parameters,\
                values_parameters, current_values, candidate_values,\
                country=country, start=start,\
                hospitalization_model=hospitalization_model,
                intensiveCare_model=intensiveCare_model
                )



    def update_Xi_R(self, Xi_R, country):
        self.Xi_R[country] = Xi_R[country]

    def update_Xi_D(self, Xi_D, country):
        self.Xi_D[country] = Xi_D[country]




    def update_sumut(self, parameter_values, start):
        if self.model == 'cases':
            self.sumut = basics.calculate_sumut(self.cases_view, parameter_values, start)
        elif self.model == 'infections':
            self.sumut = basics.calculate_sumut(self.infections_view, parameter_values, start)


    def update_Rt(self, Rt):
        for country_iter in self.countries:
            self.Rt[country_iter] = Rt[country_iter]


    def get_values(self):
        values = {'cases': self.cases,
                  'sumut': self.sumut,
                  'Rt': self.Rt,
                  'correction_factor1': self.correction_factor1,
                  'correction_factor2': self.correction_factor2,
                  }

        if self.model == 'cases':
            values.update(self.cases_view)

        if self.model == 'infections':
            values['infections'] = self.infections
            values['cases_view'] = self.cases_view
            values.update(self.infections_view)

        # if self.adapt_reporting_weekend:
        if hasattr(self, 'Xi_R'):  # von sabine checken lassen
            values['Xi_R'] = self.Xi_R
            values['Xi_D'] = self.Xi_D

        return(values)


    def update(self, parameter_values, start):
        for country in self.countries:
            self.update_cases(parameter_values, country, start)
        self.i += 1


    def update_cases(self, parameter_values, country, start):  # this should be update_infections in case of model == 'infections'
        current_values = self.get_values()
        for block in self.block_view[country].keys():
            block_selection = self.block_view[country][block]
            candidate_values = {country: copy.deepcopy(current_values[country])}
            candidate_values[country][block_selection] = self.propose_values(current_values[country][block_selection], country, block)
            candidate_values['sumut'] = basics.calculate_sumut(candidate_values,
                    parameter_values, start=start)

            candidate_values['Rt'] = current_values['Rt']
            candidate_values['correction_factor1'] = current_values['correction_factor1']
            candidate_values['correction_factor2'] = current_values['correction_factor2']

            if not isinstance(candidate_values['correction_factor1'][country], int):
                candidate_values['correction_factor1'][country] = basics.calculate_correction_factor1_country(candidate_values[country], parameter_values['N_'+country], parameter_values['probability_reinfection'])

            if ('Xi_R' in current_values.keys()): # necessary sinc we need Xi_R in the likelihood
                candidate_values['Xi_R'] = current_values['Xi_R']
                candidate_values['Xi_D'] = current_values['Xi_D']

            if self.model == 'infections':
                candidate_values['cases_view'] =\
                {country: basics.calculate_cases(candidate_values[country],\
                    parameter_values)}

            ratio_likelihood = self.likelihood_ratio(self.data[country],
                    parameter_values, current_values, candidate_values, 
                    country=country,
                    hospitalization_model=self.hospitalization_exist[country],
                    intensiveCare_model=self.intensiveCare_exist[country]
                    )

            ratio_proposals = self.proposal_ratio(current_values[country][block_selection],
                    candidate_values[country][block_selection],country, block)
            ratio = ratio_likelihood * ratio_proposals
            accept = (np.random.uniform(0, 1) < ratio)
            # print(accept)
            values_new = current_values[country][block_selection] + (candidate_values[country][block_selection]-current_values[country][block_selection]) * accept
#            self.cases_view[country][0:(len(current_values[country])+1)] = values_new

            if self.model == 'cases':
                self.cases_view[country][block_selection] = values_new
            elif self.model == 'infections':
                self.infections_view[country][block_selection] = values_new
                self.cases_view[country][:] =\
                basics.calculate_cases(self.infections_view[country], parameter_values)

            if accept:
                self.correction_factor1[country] = candidate_values['correction_factor1'][country]
                self.sumut[country] = candidate_values['sumut'][country]

            self.acceptance_cases[country][block][self.i] = accept



    def propose_values(self, current_values, country, block):
        if self.proposal_distribution == 'Poisson':
            values = stats.poisson.rvs(current_values)
        elif self.proposal_distribution == 'Uniform':
            mean = current_values
            percent = self.proposal_range[country][block]
            lower = np.maximum(0, np.floor((1 - percent) * mean))
            upper = np.maximum(2, np.ceil((1 + percent) * mean) + 1)  # +1 is required because python counts x-1 at upper bound
            values = stats.randint.rvs(low=lower, high=upper)
        return(values)

    def proposal_ratio(self, current_values, candidate_values, country, block):
        if self.proposal_distribution == 'Poisson':
            ratio = np.prod(stats.poisson.pmf(current_values, candidate_values)/stats.poisson.pmf(candidate_values, current_values))

        elif self.proposal_distribution == 'Uniform':
            mean_curr = current_values
            mean_cand = candidate_values

            percent = self.proposal_range[country][block]
            lower_curr = np.floor((1 - percent) * mean_curr) 
            lower_cand = np.floor((1 - percent) * mean_cand) 
            upper_curr = np.maximum(2, np.ceil((1 + percent) * mean_curr) + 1)  # +1 is required because python counts x-1 at upper bound
            upper_cand = np.maximum(2, np.ceil((1 + percent) * mean_cand) + 1)  # +1 is required because python counts x-1 at upper bound
            selection = (lower_curr < 0) | (lower_cand < 0) | (mean_curr == 0) | (mean_cand == 0)
            if (not np.any(selection)):
                ratio = 1
            else:
                pmf_cand = basics.calculate_pmf_uniform_trunc(lower_curr[selection],
                        upper_curr[selection]) # probavility to move from current value to candidate value (q(theta_cand|theta_curr))
                pmf_curr = basics.calculate_pmf_uniform_trunc(lower_cand[selection],
                        upper_cand[selection]) # probability to move from candidate value to current value (q(theta_curr|theta_cand))
                ratio = np.prod(pmf_curr/pmf_cand)
        return(ratio)



    def adapt_proposal(self, nb_iterations, phase):
        for country in self.countries:
            for block in self.block_view[country].keys():
                acceptance_rate = np.mean(self.acceptance_cases[country][block][phase * nb_iterations:(phase + 1) * nb_iterations - 1])
                diff = acceptance_rate - self.target
                change = abs(diff) > 0.02
                sign = np.sign(diff)
                self.proposal_range[country][block] *= (1 + 0.1 * sign * change)
                # with open('proposal_sds_lv' ,'a') as ff:
                    # ff.write(country + ':' + str(self.proposal_range) + '\n')

                print("Acceptance rate Latent_Variable " + country + ':' + str(round(acceptance_rate, 4)))


    def reset_values(self, iterations):
        for country in self.countries:
            self.acceptance_cases[country] = {block: np.empty(iterations) for block in self.block_view[country].keys()}
        self.i = 0


    def get_proposal_sd(self):
        return self.proposal_range


    def set_proposal_sd(self, proposal_sd):
        for country in self.countries:
            for block in self.block_view[country].keys():
                self.proposal_range[country][block] = proposal_sd[country][str(block)]


    def set_state(self, parameter_values: dict, path: str, chain: str, thin: int) -> None:
        """
        Set the state of the Latent variable
        """
        # set state for infections
        for country in self.countries:
            len_lv = self.infections_view[country].shape[0]
            x = np.loadtxt(f'{path}/{chain}_predictions_infections_{country}.txt')[:len_lv, -1]

            self.infections_view[country][:] = x

            self.cases_view[country][:] = basics.calculate_cases(self.infections_view[country], parameter_values)

            self.Xi_R[country] = basics.initialize_country_specific_Xi_R(parameter_values, country)
            self.Xi_D[country] = basics.initialize_country_specific_Xi_D(parameter_values, country)

            self.Rt = basics.calculate_Rt(parameter_values, self.data)
            self.update_sumut(parameter_values, start=self.start)
            self.correction_factor1 = basics.calculate_correction_factor1(self.infections_view, parameter_values)
            self.correction_factor2 = basics.calculate_correction_factor2(self.data, parameter_values)

            self.i = -9999999  # harakiri. should not need this anymore (only for adaptive phase.)
        # set state for cases and other stuff like Rt, cfs and whatever is there Oo




    def save_chains(self, path_results):
        for key in self.statistics:
            np.savetxt(path_results + '_' + key + ".txt",
                           self.statistics[key])
        for key in self.samples:
            np.savetxt(path_results  + '_' + key + ".txt",
                               self.samples[key])


    def initialize_values(self, data, fix_latent_variable=False, parameter_values=None):
        # import pudb;pu.db
        if self.model == 'cases':
            # values = stats.poisson.rvs(data['reported_cases'].to_numpy()*10)
            values = data['cases'].to_numpy()
        elif self.model == 'infections':
            # if simulated data, the real data are available....therefore best initializations

            # multiply deaths data with inverse pi_D and shift it
            # problem are the very small values: e.g. 1: then we force that we observe at least pi_D**-1 cases. this is especially for the start a problem
            values = (data['deaths'] + 1)  # +1 to get rid of multiplication with 0
            values[values <= 2] = 0.2
            # values *= (parameter_values['pi_D']**-1) 
            values *= (0.005**-1) 
            countries = data.country.unique()
            dd = values[data.country == countries[0]]
            dd[:] = np.convolve(dd, np.ones(7) / 7, mode='same')
            dd = dd.shift(-18-7).fillna(method='bfill').fillna(method = 'ffill') # shift with 18, mean of generation time distribution and 5 for incubation time
            dd = np.convolve(dd, np.ones(7) / 7, mode='same').round()
            for cc in countries[1:]:
                d = values[data.country == cc]
                d[:] = np.convolve(d, np.ones(7) / 7, mode='same')
                d = d.shift(-18-5).fillna(method = 'bfill')
                d.fillna(method='ffill', inplace=True)
                d_smooth = np.convolve(d, np.ones(7)/7, mode='same')
                dd = np.concatenate((dd, d_smooth)).round() + 1
            values = dd


        if not fix_latent_variable:
            values = stats.poisson.rvs(values) + 1
        elif fix_latent_variable:
            values = data['infections'].to_numpy()
        return(values)
         # muss noch intialisiert werden


    def initialize_blocks(self, country, nb_blocks = 10):
        np.random.seed(0)
        if self.model == 'cases':
            cases_country = self.cases_view[country]
        elif self.model == 'infections':
            cases_country = self.infections_view[country]
        block_length = np.int(np.ceil(cases_country.shape[0] / nb_blocks))
        blocks = {}
        for block in range(nb_blocks):
            maximum = np.minimum(cases_country.shape[0], block_length*(block+1))
            blocks[block] = list(range(block_length * block, maximum))
        return(blocks)

