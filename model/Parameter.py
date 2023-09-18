import numpy as np
import scipy.stats as stats
import copy
import update


class Parameter:

    def __init__(self, name, model, precision, proposal_sd, prior_parameters,
                 start_values, data, start=-9999, name_parent=None, country='all',
                 informative_priors=True,
                 reporting_model=False, renewal_model=True,
                 death_model=False, hospitalization_model=False,
                 intensiveCare_model=False):

        if not name_parent:
            name_parent = name
        self.name = name
        self.target = 0.4
        self.precision = precision[name_parent]
        self.samples = np.empty(30000)  # 30k empty samples - gets overwritten for actual sampling. Should only be large enough for adaptive phases + burnin

        if name_parent == 'rho':
            if 'rho_period' in data.columns:
                self.samples[0] = start_values['rho_' + country][name]
            else:
                self.samples[0] = start_values['rho_' + country]

            reporting_model = True
            renewal_model = False

            self.country = country

        elif name_parent == 'alpha' and len(name.split('_')) == 3:
            m = start_values['alpha_' + name.split('_')[1]]
            sd = start_values['alpha_sd'] * 1.1
            self.samples[0] = stats.norm.rvs(loc=m, scale=sd)  # Gaussian initial values

        else:
            self.samples[0] = start_values[name]

        self.i = 0
        self.acceptance = np.empty(30000, dtype=bool)
        self.proposal_sd = proposal_sd[name_parent]

        self.prior_parameters = prior_parameters[name_parent]
        self.prior_ratio = lambda theta_cand: update.priors[name_parent](self.samples[self.i], theta_cand, self.prior_parameters, informative_priors)

        self.propose_value = lambda: update.proposals[name_parent](self.samples[self.i], self.proposal_sd)

        self.calculate_proposal_ratio = update.proposals[name_parent].__name__ in ['proposal_trunc', 'proposal_double_trunc', 'propsal_trunc_uniform', 'proposal_double_trunc_rho']

        self.likelihood_ratio = lambda \
                values_t, values_cand, latent_variable_curr, latent_variable_cand:\
                update.ratio_likelihood_nbinom(model, data, values_t, values_cand,\
                        latent_variable_curr, latent_variable_cand, \
                        country=country, death_model=death_model,\
                        reporting_model=reporting_model,\
                        renewal_model_lv=renewal_model,\
                        # renewal_model_parameter could allow an accelterated version of the renewal model but is currently not impleented
                        renewal_model_parameter=False,\
                        hospitalization_model=hospitalization_model,\
                        intensiveCare_model=intensiveCare_model,
                        start=start)

        if name_parent == 'phi':
            name_parent = name

        self.update = lambda current_values, latent_variable: update.parameter_update(self, current_values, latent_variable, name_parent)


    def adapt_proposal(self, nb_iterations, phase):
        acceptance_rate = np.mean(self.acceptance[phase * nb_iterations:(phase + 1) * nb_iterations - 1])
        diff = acceptance_rate - self.target
        change = abs(diff) > 0.02
        sign = np.sign(diff)
        self.proposal_sd *= (1 + 0.1 * sign * change)
        print(self.name)
        print("Acceptance rate: " + str(round(acceptance_rate, 4)))


    def write_samples(self, path_results, thin):
        if (isinstance(self.name, np.integer)):
            file_name = path_results + 'results_' + 'rho_' + self.country + str(self.name) + ".txt"

        elif (isinstance(self, PriorParameter)):
            file_name = path_results + 'results_' + self.name_parent + '_' + self.name + '.txt'
        else:
            file_name = path_results + 'results_' + self.name + ".txt"

        index = np.arange(0, len(self.samples[:self.i + 1]), thin)
        np.savetxt(file_name, self.samples[:self.i + 1][index])

    def get_current_value(self):
        current_value = self.samples[self.i]
        return(current_value)


    def proposal_ratio(self, theta_cand, name_parent):
        theta_curr = self.get_current_value()
        proposal_sd = self.proposal_sd
        a_trunc_curr = (0 - theta_cand) / proposal_sd
        a_trunc_cand = (0 - theta_curr) / proposal_sd

        if update.proposals[name_parent].__name__ == 'proposal_trunc':
            ratio = stats.truncnorm.pdf(theta_curr, a=a_trunc_curr, b=np.inf, loc=theta_cand, scale=proposal_sd) / \
            stats.truncnorm.pdf(theta_cand, a=a_trunc_cand, b=np.inf, loc=theta_curr, scale=proposal_sd)
        elif update.proposals[name_parent].__name__ == 'proposal_double_trunc':
            b_trunc_curr = (1 - theta_cand) / proposal_sd
            b_trunc_cand = (1 - theta_curr) / proposal_sd
            ratio = stats.truncnorm.pdf(theta_curr, a=a_trunc_curr, b=b_trunc_curr, loc=theta_cand, scale=proposal_sd) / \
            stats.truncnorm.pdf(theta_cand, a=a_trunc_cand, b=b_trunc_cand, loc=theta_curr, scale=proposal_sd)
        elif update.proposals[name_parent].__name__ == 'proposal_double_trunc_rho':
            b_trunc_curr = (10 - theta_cand) / proposal_sd
            b_trunc_cand = (10 - theta_curr) / proposal_sd
            ratio = stats.truncnorm.pdf(theta_curr, a=a_trunc_curr, b=b_trunc_curr, loc=theta_cand, scale=proposal_sd) / \
            stats.truncnorm.pdf(theta_cand, a=a_trunc_cand, b= b_trunc_cand, loc=theta_curr, scale=proposal_sd)

        return(ratio)


    def update_priors(self, current_prior_values):
        for prior_parm in self.prior_parameters:
            self.prior_parameters[prior_parm] = current_prior_values[prior_parm]


    def reset_values(self, iterations):
        current_value = self.get_current_value()
        self.samples = np.empty(iterations + 2)
        self.samples[0] = current_value
        self.acceptance = np.empty(iterations + 2, dtype=np.bool)
        self.i = 0


    def get_statistics(self):
        statistics = {}
        i = self.i
        statistics['median'] = round(np.median(self.samples[:i]), self.precision)
        statistics['mean'] = round(np.mean(self.samples[:i]), self.precision)
        statistics['IC_2.5'] = round(np.percentile(self.samples[:i], 2.5), self.precision)
        statistics['IC_97.5'] = round(np.percentile(self.samples[:i], 97.5), self.precision)
        statistics['acceptance'] = round(np.mean(self.acceptance[self.i]), 4)
        return(statistics)


    def get_proposal_sd(self):
        return self.proposal_sd


    def set_proposal_sd(self, proposal_sd):
        self.proposal_sd = proposal_sd


    def set_state(self, path: str, chain: str, thin: int) -> None:
        """
        Function to set the current state of the parameter. It loads the full
        trjectory of the chain (assumes np.array as txt) from path and stes the
        parameter at this state

        :param path: The path where to search for the np.array
        :param chain: Info which chain to load
        :param thin: Info about the thinning to use. A thinning of 10 implies a
            fill trajectory length of the parameter x*(10-1)+1
        """
        if (isinstance(self.name, np.integer)):
            # rho is special since its deepest level has not the full name, only an integer as name
            name_parameter = 'rho_' + self.country + str(self.name)
        else:
            name_parameter = self.name

        x = np.loadtxt(f'{path}/{chain}_results_{name_parameter}.txt')
        # recycle values thin times since we only have the thinned values (exception is the last one)
        x_long = np.append(np.repeat(x[:-1], thin), x[-1])
        self.i = x_long.shape[0] - 1  # -1 here because i points always to the current value and shape is always 1 longer sisnce python starts as 0
        self.samples[:self.i + 1] = x_long




class FixedParameter:

    def __init__(self, name, value):
        self.name = name
        self.value = value
        self.i = 0
        self.proposal_sd = 0

    def adapt_proposal(self, nb_iterations, phase):
        pass

    def update(self, current_values, latent_variable):
        pass

    def update_priors(self, current_prior_values):
        pass

    def write_samples(self, path_results, thin=None):
        pass

    def reset_values(self, iterations):
        pass

    def get_current_value(self):
        current_value = self.value
        return(current_value)

    def get_proposal_sd(self):
        return self.proposal_sd

    def set_proposal_sd(self, proposal_sd):
        pass

    def set_state(self, path: str, chain: str, thin: int):
        pass

    def get_statistics(self):
        statistics = {}
        statistics['median'] = 'fixed'
        statistics['mean'] = 'fixed'
        statistics['IC_2.5'] = 'fixed'
        statistics['IC_97.5'] = 'fixed'
        statistics['acceptance'] = 'fixed'
        return(statistics)


class ParameterVector:

    def __init__(self, name, keys, model, precision, proposal_sd, prior_parameters,
                 start_values, data, start=-9999, informative_priors=True):
        """
        This class is a wrapper for a single Parameter from outside
        it behaves like a single parameter but is doing the functionality
        for many of them. If there is a column rho_preiod in the data, the
        ParameterVector is used recursively (therefore vector of vectors of
        parameters). The same logic logic counts for NPIs. There check for a
        nested dict to identify the hierarhcical order
        """

        name_parent = name.split('_')[0]

        if name_parent in ['beta', 'betaD']:
            name_parent = name

        hierarchical_alpha = isinstance(keys, dict)  # checkt whether the algo is in the first level of the nested structure of alpha
        name_intervention = None
        if hierarchical_alpha:
            keys_country = keys['countries']
            keys = keys['interventions']
        if name == 'alpha':
            data = {alpha_key: data for alpha_key in keys}  # necessary as parameterVector expects data dict - not the best solution but atm fine
            country = {alpha_key: 'all' for alpha_key in keys}
        else:
            country = {country_key: country_key for country_key in keys}

        self.name = name
        self.parameters = {}
        self.prior_parameters = {}

        # this condition is only for parameter vectors with priors, else the prior parameters are FixedParameters
        if (name in ['R0', 'tau']) or name[0:6] == 'alpha_':
            # init hierarchical alpha prior
            name_intervention = name.split('_')[-1]
            for parm in prior_parameters[name_parent].keys():
                self.prior_parameters[parm] = PriorParameter(parm,
                                                             precision,
                                                             proposal_sd,
                                                             prior_parameters,
                                                             start_values,
                                                             name_parent=name_parent,
                                                             name_intervention=name_intervention)
                prior_parameters = copy.deepcopy(prior_parameters)
                prior_parameters[name_parent][parm] = self.prior_parameters[parm].get_current_value()

        else:
            prior_name = name.split('_')[0]  # necessary for the recursion function call for rho

            if prior_name in ['beta', 'betaD']:
                prior_name = name

            for parm in prior_parameters[prior_name].keys():
                self.prior_parameters[parm] = FixedParameter(parm, prior_parameters[prior_name][parm])

            # init alpha, rho priors
        if name == 'rho':
            # ifelse condition necessary because the should be flexible concerning the rho_period column
            if 'rho_period' in data[list(keys)[0]].columns:
                for parm in keys:
                    self.parameters[name + '_' + parm] = ParameterVector(name + '_' + parm,
                                                                         keys=data[parm]['rho_period'].unique(),
                                                                         model=model,
                                                                         precision=precision,
                                                                         proposal_sd=proposal_sd,
                                                                         prior_parameters=prior_parameters,
                                                                         start_values=start_values,
                                                                         data=data[parm],
                                                                         informative_priors=informative_priors)
            else:
                for parm in keys:
                    self.parameters[name + '_' + parm] = Parameter(name + '_' + parm, model,
                                                                   precision,
                                                                   proposal_sd,
                                                                   prior_parameters,
                                                                   start_values,
                                                                   data[parm], start=start,
                                                                   name_parent=name_parent,
                                                                   country=country[parm],
                                                                   informative_priors=informative_priors,
                                                                   reporting_model=True, renewal_model=False)


        elif hierarchical_alpha:
            for parm in keys:
                self.parameters[name + '_' + parm] = ParameterVector(name + '_' + parm,
                                                                     keys=keys_country,
                                                                     model=model,
                                                                     precision=precision,
                                                                     proposal_sd=proposal_sd,
                                                                     prior_parameters=prior_parameters,
                                                                     start_values=start_values,
                                                                     data=data[parm],
                                                                     start=start,
                                                                     informative_priors=informative_priors)


        # The following loop is called by a recursive call of Parameter vector which is done for rho (parameter vector of parameter vectors)
        elif name[0:4] == 'rho_':
            # keys is here a list with values defined by data["rho_period"].unique()
            for period in keys:
                self.parameters[period] = Parameter(period, model,
                                                    precision,
                                                    proposal_sd,
                                                    prior_parameters,
                                                    start_values,
                                                    data,  # gets all data of a country
                                                    start=start,
                                                    name_parent=name_parent,
                                                    country=self.name[4:],
                                                    informative_priors=informative_priors,
                                                    reporting_model=True, renewal_model=False
                                                    )

        else:
            for parm in keys:
                if name in ['beta_sat', 'beta_sun', 'beta_mon', 'beta_tue', 'beta_wed', 'beta_fri']:
                    death_model = False
                    reporting_model = True
                    renewal_model = False
                    hospitalization_model = False
                    intensiveCare_model = False
                elif name in ['betaD_sat', 'betaD_sun', 'betaD_mon', 'betaD_tue', 'betaD_wed', 'betaD_fri']:
                    death_model = True
                    reporting_model = False
                    renewal_model = False
                    hospitalization_model = False
                    intensiveCare_model = False
                elif name in ['piH']:
                    # condition if hospitalizations are not avalable at all for this country, then no parameter should be initialized
                    if not parm in list(map(lambda x: x.split("_")[1], start_values.keys())):
                        continue
                    death_model = False
                    reporting_model = False
                    renewal_model = False
                    hospitalization_model = True
                    intensiveCare_model = False
                elif name in ['piHicu']:
                    # condition if hospitalizations are not avalable at all for this country, then no parameter should be initialized
                    if not parm in list(map(lambda x: x.split("_")[1], start_values.keys())):
                        continue
                    death_model = False
                    reporting_model = False
                    renewal_model = False
                    hospitalization_model = False
                    intensiveCare_model = True
                else:
                    # the case for R0, tau, alpha (and more?)
                    death_model = False
                    reporting_model = False
                    renewal_model = True
                    hospitalization_model = False
                    intensiveCare_model = False

                self.parameters[name + '_' + parm] = Parameter(name + '_' + parm, model,
                                                               precision,
                                                               proposal_sd,
                                                               prior_parameters,
                                                               start_values,
                                                               data[parm], 
                                                               start=start,
                                                               name_parent=name_parent,
                                                               country=country[parm],
                                                               informative_priors=informative_priors,
                                                               reporting_model=reporting_model,
                                                               renewal_model=renewal_model,
                                                               death_model=death_model,
                                                               hospitalization_model=hospitalization_model,
                                                               intensiveCare_model=intensiveCare_model
                                                               )




    def update(self, current_values, latent_variable):
        for p in self.parameters.keys():
            current_values.update(current_values[self.name])
            self.parameters[p].update(current_values, latent_variable)
            current_values[self.parameters[p].name] = self.parameters[p].get_current_value()
        current_prior_values = self.get_current_prior_value()
        for prior_p in self.prior_parameters.keys():
            self.prior_parameters[prior_p].update(current_values, current_prior_values)
            current_prior_values[prior_p] = self.prior_parameters[prior_p].get_current_value()
        if not isinstance(list(self.parameters.values())[0], ParameterVector):
            for p in self.parameters.keys():
                self.parameters[p].update_priors(current_prior_values)

    def adapt_proposal(self, nb_iterations, phase):
        for p in self.parameters.keys():
            self.parameters[p].adapt_proposal(nb_iterations, phase)
        for prior_p in self.prior_parameters.keys():
            self.prior_parameters[prior_p].adapt_proposal(nb_iterations, phase)


    def write_samples(self, path_results, thin=1):
        for parameter in self.parameters.keys():
            self.parameters[parameter].write_samples(path_results=path_results, thin=thin)
        for prior_p in self.prior_parameters.keys():
            self.prior_parameters[prior_p].write_samples(path_results=path_results, thin=thin)



    def reset_values(self, iterations):
        for parameter in self.parameters.keys():
            self.parameters[parameter].reset_values(iterations)
        for prior_parameter in self.prior_parameters.keys():
            self.prior_parameters[prior_parameter].reset_values(iterations)




    def update_priors(self, current_prior_values):
        for parameter in self.parameters.keys():
            self.parameters[parameter].update_priors(current_prior_values)



    def get_current_value(self):
        current_values = {parm: self.parameters[parm].get_current_value() for parm in self.parameters.keys()}
        return(current_values)


    def get_current_prior_value(self):
        current_prior_values = {prior_parm: self.prior_parameters[prior_parm].get_current_value() for prior_parm in self.prior_parameters.keys()}
        return(current_prior_values)


    def get_statistics(self):
        statistics = {}
        for parameter in self.parameters.keys():
            statistics[parameter] = self.parameters[parameter].get_statistics()
        for prior_parameter in self.prior_parameters.keys():
            statistics[prior_parameter] = self.prior_parameters[prior_parameter].get_statistics()
        return(statistics)


    def get_proposal_sd(self):
        proposal_sds = {}
        for parm in self.parameters:
            # specific case for rho: keys are np.int64. json can only dump normal ints
            if isinstance(parm, np.int64):
                proposal_sds[int(parm)] = self.parameters[parm].get_proposal_sd()
            else:
                proposal_sds[parm] = self.parameters[parm].get_proposal_sd()
        proposal_sds['priors'] = {}
        for prior in self.prior_parameters:
            proposal_sds['priors'][prior] = self.prior_parameters[prior].get_proposal_sd()
        return proposal_sds


    def set_proposal_sd(self, proposal_sd):
        for parm in self.parameters:
            if isinstance(parm, np.int64):
                # loaded json makes str instead of np.int64
                self.parameters[parm].set_proposal_sd(proposal_sd[str(parm)])
            else:
                self.parameters[parm].set_proposal_sd(proposal_sd[parm])
        for prior in self.prior_parameters:
            self.prior_parameters[prior].set_proposal_sd(proposal_sd['priors'][prior])


    def set_state(self, path: str, chain: str, thin: int):
        """
        Function to set the current state of the ParameterVector. It iterates
        over the parameters in the vector.

        :param path: The path where to search for the np.array
        :param chain: Info which chain to load
        :param thin: Info about the thinning to use. A thinning of 10 implies a
            fill trjectory length of the parameter x*(10-1)+1
        """
        for parm in self.parameters:
            self.parameters[parm].set_state(path, chain, thin)

        for prior in self.prior_parameters:
            self.prior_parameters[prior].set_state(path, chain, thin)


class PriorParameter(Parameter):

    def __init__(self, name, precision, proposal_sd, prior_parameters, start_values, name_parent, name_intervention=None):
        model = 'harakiriii'
        if name_parent == 'alpha' and name == 'mean':
            start_values[name] = start_values[name_parent + '_' + name_intervention]
        else:
            start_values[name] = start_values[name_parent + '_' + name]
        super().__init__(name, model, precision, proposal_sd, prior_parameters,
                         start_values, data=None, start=-1,
                         name_parent=name_parent, country='harakiriii', informative_priors=True, renewal_model=False,
                         reporting_model=False, death_model=False, hospitalization_model=False, intensiveCare_model=False)

        if name_parent == 'alpha':
            self.name_parent = name_parent + '_' + name_intervention
        else:
            self.name_parent = name_parent

        # condition for season
        if (name_intervention[:6] == 'season') and (name == 'mean'):
            self.prior_parameters = prior_parameters[name_parent + '_season_' + name]
        else:
            self.prior_parameters = prior_parameters[name_parent + '_' + name]

        self.prior_ratio = lambda theta_cand: update.priors[name_parent + '_' + name](self.samples[self.i], theta_cand, self.prior_parameters, informative_priors=True)  # overwritten version
        self.likelihood_ratio = update.prior_likelihoods[name_parent]
        self.update = lambda current_values, current_prior_values: update.prior_update(self, current_values, current_prior_values, name_parent + '_' + name)

        self.propose_value = lambda: update.proposals[name_parent + '_' + name](self.samples[self.i], self.proposal_sd)

        self.calculate_proposal_ratio = update.proposals[name_parent + '_' + name].__name__ in ['proposal_trunc', 'proposal_double_trunc']

    def extract_values(self, current_values):
        values_list = list(current_values[self.name_parent].values())
        values = np.array(values_list)
        return(values)


    def get_proposal_sd(self):
        return self.proposal_sd


    def set_state(self, path: str, chain: str, thin: int) -> None:
        """
        Basically the same funct as for Parameter. However, here we need a custom
        method since self.name is just the name. TO load it we also need the
        self.name_parent.
        For detials see doc of this method for parameter.
        """
        name_parameter = self.name_parent + '_' + self.name

        x = np.loadtxt(f'{path}/{chain}_results_{name_parameter}.txt')
        # recycle values thin times since we only have the thinned values (exception is the last one)
        x_long = np.append(np.repeat(x[:-1], thin), x[-1])
        self.i = x_long.shape[0] - 1
        self.samples[:self.i + 1] = x_long
