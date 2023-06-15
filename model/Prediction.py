import numpy as np
import pandas as pd
import scipy.stats as stats
import basics
import copy
import datetime
import random



class Prediction:
    '''

    '''
    def __init__(self, model, data, nb_futue_values,
                 sample_size=1000, parameter_values=None, oos_data=None):
        self.model = model
        self.data = basics.dictionarize_data(data)
        if oos_data is not None:
            # oos country is assumed to be only ONE country
            self.oos_data = oos_data
            self.oos_country = oos_data.country.unique()[0]
            self.nb_days_oos_init = self.data[self.oos_country].shape[0]
            # calculate cf2 only at start
            self.oos_cf2 = basics.calculate_correction_factor2_country(
                    self.oos_data, parameter_values,
                    self.oos_country)
        else:
            self.oos_country = None

        self.countries = self.data.keys()

        self.future_indices = dict.fromkeys(self.countries)
        self.cases = dict.fromkeys(self.countries)
        for country_iter in self.countries:
            if country_iter == self.oos_country:
                self.future_indices[country_iter] = np.arange(self.data[country_iter].shape[0], self.oos_data.shape[0])
                self.cases[country_iter] = np.zeros((self.oos_data.shape[0], sample_size), dtype=int)
            else:
                self.future_indices[country_iter] = np.arange(self.data[country_iter].shape[0], self.data[country_iter].shape[0] + nb_futue_values)
                self.cases[country_iter] = np.zeros((self.data[country_iter].shape[0] + nb_futue_values, sample_size), dtype=int)

        if self.model == 'infections':
            self.infections = copy.deepcopy(self.cases)


            self.correction_factor_2 = {}
            cf2 = {}
            for country_iter in self.countries:
                if country_iter == self.oos_country:
                    cf2[country_iter] = basics.calculate_correction_factor2_country(self.oos_data, parameter_values, country_iter)
                else:
                    cf2[country_iter] = basics.calculate_correction_factor2_country(self.data[country_iter][:(self.future_indices[country_iter][0])], parameter_values, country_iter)
                    # cf2_mean_diff = cf2[country_iter].diff().values[-7:].mean()
                    cf2_mean_diff = np.diff(cf2[country_iter])[-7:].mean()

                    # self.correction_factor_2[country_iter] = np.concatenate((cf2[country_iter].values, np.zeros(nb_futue_values)))
                    self.correction_factor_2[country_iter] = np.concatenate((cf2[country_iter], np.zeros(nb_futue_values)))

                    for t in self.future_indices[country_iter]:
                        self.correction_factor_2[country_iter][t] = self.correction_factor_2[country_iter][t - 1] + cf2_mean_diff


        self.deaths = copy.deepcopy(self.cases)
        self.hospitalizations = copy.deepcopy(self.cases)
        self.intensiveCare = copy.deepcopy(self.cases)
        self.reported_cases = copy.deepcopy(self.cases)
        self.Rt = copy.deepcopy(self.cases)

        self.i = 0


        # append data
        for country_iter in self.countries:
            for ind in range(nb_futue_values):
                new_date = self.data[country_iter].date.iloc[-1] + datetime.timedelta(days=1)
                new_weekday = datetime.datetime.weekday(new_date)

                if 'improved_treatment' in self.data[country_iter].columns:
                    new_improved_treatment = self.data[country_iter].improved_treatment.iloc[-1]
                    # self.data[country_iter] = self.data[country_iter].append({
                    # 'date': new_date,
                    # 'weekday':new_weekday,
                    # 'improved_treatment': new_improved_treatment
                    # },ignore_index=True)
                    self.data[country_iter] = pd.concat([self.data[country_iter], pd.DataFrame({
                                                             'date': [new_date],
                                                             'weekday': [new_weekday],
                                                             'improved_treatment': new_improved_treatment
                                                             })],
                                                        ignore_index=True)

                else:
                    # self.data[country_iter] = self.data[country_iter].append({  # this version uses append which is depreciated in futere pandas versions
                    # 'date': new_date,
                    # 'weekday':new_weekday,
                    # },ignore_index=True)
                    self.data[country_iter] = pd.concat([self.data[country_iter], pd.DataFrame({
                                                             'date': [new_date],
                                                             'weekday': [new_weekday],
                                                             })],
                                                        ignore_index=True)




                if 'rho_period' in self.data[country_iter].columns:
                    # import pudb; pu.db
                    # self.data[country_iter].rho_period.iloc[-1] = self.data[country_iter].rho_period.iloc[-2]
                    self.data[country_iter].rho_period.iat[-1] = self.data[country_iter].rho_period.iat[-2]


    def make_predictions(self, parameter_values, latent_variables, current_prior_values):
        self.predict_cases(parameter_values, latent_variables)
        if not self.oos_country is None:
            self.make_oos_cases(parameter_values, current_prior_values, latent_variables)
        self.predict_observed_quantities(parameter_values, latent_variables)
        self.i += 1












    def predict_cases(self, parameter_values, latent_variables,
                      intervention_change={}):
        # intevention_change is a dictionary which indicates for every changeed
        # intervention wether it was added (+1) or removed (-1) 

        
        # ToDo: Aktuell wird angenommen, dass Rt in der Zukunft konstant bleibt
        # (letzter wert). Um zu berücksichtigen, dass neue interventionen in der
        # Zukunft sein könnten, müssten diese interventionen predit_cases über-
        # geben werden und Rt müsste neu berechnet werden (es müsste ein
        # dataframe übergeben werden auf das calculate_Rt angewendet werden kann)

        for country_iter in self.countries:
            if country_iter == self.oos_country:
                pass
            else:
                # calculate Rt for predictions
                Rt = latent_variables['Rt'][country_iter][-1]

                self.Rt[country_iter][:,self.i] = np.concatenate((latent_variables['Rt'][country_iter], np.repeat(Rt, self.future_indices[country_iter].shape[0])))

                # falls die interventionen sich in Zukunft ändern könnte das so geschehen:
                for intervention in intervention_change.keys():

                    # multiplication with -1 to invert the sign, since the formula
                    # for Rt is exp(-alpha)
                    Rt *= np.exp(intervention_change[intervention] * (-1) * parameter_values['alpha']['alpha_' + intervention])


                if self.model == 'infections':
                    infections = self.infections[country_iter][:,self.i] # Attention! This should be a pointer
                    infections[:(self.future_indices[country_iter][0])] = latent_variables[country_iter]
                    # cf2 = basics.calculate_correction_factor2_country(self.data[country_iter][:(self.future_indices[country_iter][0])], parameter_values, country_iter).iloc[-1]
                    cf2 = basics.calculate_correction_factor2_country(self.data[country_iter][:(self.future_indices[country_iter][0])], parameter_values, country_iter)[-1]
                    for t in self.future_indices[country_iter]:
                        cf1 = basics.calculate_correction_factor1_country(infections[:t+1], parameter_values['N_'+country_iter], parameter_values['probability_reinfection'])[t] # here t+1 in infection sis fine because the last value gets thrown away
                        future_value = self.generate_future_value(infections, parameter_values, t, Rt, cf1, cf2)
                        infections[t] = future_value

                    # predict all cases including future cases (yes this is redundant...)
                    self.cases[country_iter][:,self.i] = basics.calculate_cases(infections, parameter_values)


    def make_oos_cases(self, parameter_values, current_prior_values, latent_variables):
        list_tmp = list(latent_variables['Xi_D'].keys())
        list_tmp.remove(self.oos_country)
        random_country = random.choice(list_tmp)

        alpha_oos= {alpha_key: stats.norm.rvs(current_prior_values['alpha'][alpha_key]['mean'], current_prior_values['alpha'][alpha_key]['sd']) for alpha_key in current_prior_values['alpha'].keys()}
        parameter_values_oos = {
                'phi_infections':parameter_values['phi_infections'],
                # 'beta_voc': parameter_values['beta_voc'],
                'beta_alpha': parameter_values['beta_alpha'],
                'beta_delta': parameter_values['beta_delta'],
                'R0':{'R0_' + self.oos_country: parameter_values['R0']['R0_'+ self.oos_country]},
                'alpha': {alpha_key: {alpha_key + '_' + self.oos_country: alpha_oos[alpha_key]} for alpha_key in alpha_oos},
                'N': parameter_values['N_' + self.oos_country] # this is actually N + self.oos_country, currently single country and therefore N is sufficient
                }

        # caluclate Rt
        Rt = basics.calculate_Rt_country(parameter_values=parameter_values_oos, 
                                         data=self.oos_data, country=self.oos_country)
        self.Rt[self.oos_country][:,self.i] = Rt
        # given Rt cf2 one must calculate iteratively the infections since cf1 depends on the past
        

        prob_reinfection = parameter_values['probability_reinfection']

        infections = np.zeros((self.oos_data.shape[0]), dtype = int)
        infections[:self.nb_days_oos_init] = latent_variables[self.oos_country] 
        N_tmp = parameter_values_oos['N']

        # cf1_oos_t = (np.sum(infections[:self.nb_days_oos_init])/parameter_values_oos['N'])*(1 - prob_reinfection)
        cf1_oos_t = (np.minimum(np.sum(infections[:self.nb_days_oos_init]), N_tmp)/N_tmp)*(1 - prob_reinfection)

        for t in range(self.nb_days_oos_init+1, len(infections)):
            infections[t] = self.generate_future_value(infections, parameter_values, t, Rt[t-1], cf1 = cf1_oos_t, cf2 = self.oos_cf2[t-1])
            # cf1_oos_t = (np.sum(infections[:t])/parameter_values_oos['N'])*(1 - prob_reinfection)
            cf1_oos_t = (np.minimum(np.sum(infections[:t]), N_tmp)/N_tmp)*(1 - prob_reinfection)

        self.infections[self.oos_country][:,self.i] = infections
        self.cases[self.oos_country][:,self.i] = basics.calculate_cases(infections, parameter_values)




    def generate_future_value(self,cases, parameter_values, t, Rt, cf1=0, cf2=0):
        gamma = parameter_values['gamma']
        sumut = 0
        ##### jit this for loop?
        for u in range(t):
            sumut += cases[u] * gamma[t-u-1]

        # EC_t = sumut * Rt * (1-cf1-cf2)# this is actually EI_t in the case of model == infections
        cf_full = 1 - cf1 - cf2*(1-cf1)
        EC_t = sumut * Rt * cf_full# this is actually EI_t in the case of model == infections
        future_value = basics.nbinom_alt_rvs(EC_t, parameter_values['phi_infections'])
        return(future_value)
            
    
    


    def predict_observed_quantities(self, parameter_values, latent_variables):
        # get mean piH for countries where no piH is available 
        # mean_piH = np.array([parameter_values['piH'][cc] for cc in parameter_values['piH'].keys()]).mean()
        # now sample, see in the if condition

        # required for piH
        parameter_values_tmp = copy.deepcopy(parameter_values)

        if not self.oos_country is None:
            random_country = random.choice(list(latent_variables['Xi_D'].keys()))

            if 'Xi_R' in latent_variables.keys():
                latent_variables['Xi_D'][self.oos_country] = latent_variables['Xi_D'][random_country]
                latent_variables['Xi_R'][self.oos_country] = latent_variables['Xi_R'][random_country]
            rho_oos = parameter_values['rho']['rho_' + random_country]
            # calc oos_piH
            random_country_piH =  random.choice(list(parameter_values['piH'].keys())) # is required because not every country has its own piH
            random_country_piHicu =  random.choice(list(parameter_values['piHicu'].keys())) # is required because not every country has its own piH
            pi_H_predictions = parameter_values['piH'][random_country_piH]
            pi_Hicu_predictions = parameter_values['piHicu'][random_country_piHicu]
            pi_H_oos = 1/3*pi_H_predictions + 2/3*parameter_values['piH']['piH_' + self.oos_country]
            pi_Hicu_oos = 1/3*pi_Hicu_predictions + 2/3*parameter_values['piHicu']['piHicu_' + self.oos_country]

            parameter_values_tmp['piH']['piH_' + self.oos_country] = pi_H_oos
            parameter_values_tmp['piHicu']['piHicu_' + self.oos_country] = pi_Hicu_oos

            parameter_values_tmp['rho']['rho_' + self.oos_country] = rho_oos

            parameter_values_tmp['pi_D'][self.oos_country] = parameter_values['pi_D'][self.oos_country + '_oos']
            parameter_values_tmp['correction_hospitalization_'+self.oos_country] = parameter_values['correction_hospitalization_'+self.oos_country + '_oos']

        for country_iter in self.countries:
            cases = self.cases[country_iter][:, self.i]

            if 'Xi_R' in latent_variables.keys():
                Xi_R = latent_variables['Xi_R'][country_iter]['values']
                Xi_D = latent_variables['Xi_D'][country_iter]['values']
                simple_xid = False # flag if xid is Tx7 or not
            else:
                Xi_R = parameter_values['xi_R']
                Xi_D = parameter_values['xi_D']
                simple_xid = True # flag if xid is Tx7 or not


            # param_values_tmp = copy.deepcopy(parameter_values)

            # deaths
             # expand pi_D to the full length T + predictions
            last_val = parameter_values_tmp['pi_D'][country_iter][-1]
            fut_n = next(iter(self.future_indices.values())).shape[0]
            fut_n = self.future_indices[country_iter].shape[0]
            
            if country_iter == self.oos_country:
                data = self.oos_data    
            else:
                data = self.data[country_iter]
                parameter_values_tmp['pi_D'][country_iter] = np.concatenate((parameter_values_tmp['pi_D'][country_iter], np.repeat(last_val, fut_n)))

                 # expand correction for piH to the full length T + predictions
                last_val = parameter_values_tmp['correction_hospitalization_' + country_iter][-1]
                parameter_values_tmp['correction_hospitalization_' + country_iter] = np.concatenate((parameter_values_tmp['correction_hospitalization_' + country_iter], np.repeat(last_val, fut_n)))

            if simple_xid:
                ED_t = basics.calculate_EDt(cases, parameter_values_tmp, data, country=country_iter,Xi_D=None)
            else:
                ED_t = basics.calculate_EDt(cases, parameter_values_tmp, data, country=country_iter,Xi_D=Xi_D)
            # deaths = stats.poisson.rvs(ED_t)
            deaths = basics.nbinom_alt_rvs(ED_t, parameter_values['phi_deaths'])

            # hospitalizations
            if not 'piH_' + country_iter in  parameter_values_tmp['piH'].keys():
                # parameter_values_tmp['piH']['piH_' + country_iter] = mean_piH
                random_country_piH = random.choice(list(parameter_values_tmp['piH'].keys()))
                parameter_values_tmp['piH']['piH_' + country_iter] = parameter_values_tmp['piH'][random_country_piH]
            if not 'piHicu_' + country_iter in  parameter_values_tmp['piHicu'].keys():
                # parameter_values_tmp['piH']['piH_' + country_iter] = mean_piH
                random_country_piHicu = random.choice(list(parameter_values_tmp['piHicu'].keys()))
                parameter_values_tmp['piHicu']['piHicu_' + country_iter] = parameter_values_tmp['piHicu'][random_country_piHicu]

            EH_t = basics.calculate_EHt(cases, parameter_values_tmp, country_iter)
            EHicu_t = basics.calculate_EHicut(cases, parameter_values_tmp, country_iter)
            # EH_t = basics.calculate_EHt(cases, parameter_values, country_iter)
            # hospitalizations = stats.poisson.rvs(EH_t)
            hospitalizations = basics.nbinom_alt_rvs(EH_t, parameter_values['phi_hospitalizations'])
            intensiveCare = basics.nbinom_alt_rvs(EHicu_t, parameter_values['phi_intensiveCare'])

            # reported_cases
            if country_iter == self.oos_country: # in the short vector there are not all periods available, !!assumes rho_period in data - will break when not
                # sample aus parameter_values['rho'] ein dict und für die nicht vorhandenen prerioden
                rho_dict = parameter_values['rho']['rho_' + country_iter] # to short - only has periods for the first K obs
                # expand with samples from the rest:
                rho_countries_diff = list(parameter_values['rho'].keys()) # get all countries
                rho_countries_diff.remove('rho_' + country_iter)  # remove oos country
                for period_iter in range(max(rho_dict.keys()) + 1, self.oos_data.rho_period.values.max() + 1):
                    rho_dict[period_iter] = parameter_values['rho'][random.choice(rho_countries_diff)][period_iter]  # sample from rest of ccs and assign its period value        
                rho = rho_dict
            else:
                rho = parameter_values['rho']['rho_' + country_iter]
            if 'rho_period' in self.data[country_iter].columns:
                rho = np.array([rho[period] for period in data['rho_period'].values])
                # rho = np.append(rho_vector, np.ones_like(self.future_indices[country_iter]) * rho_vector[-1])
            
            # E_Crt = basics.calculate_E_Crt(cases, parameter_values, self.data[country_iter], rho)
            E_Crt = basics.calculate_E_Crt(cases, Xi_R, data, rho)
            
            # reported_cases =  stats.binom.rvs(cases, rho)
            reported_cases = basics.nbinom_alt_rvs(E_Crt, parameter_values['phi_rep_cases'])
           
            self.deaths[country_iter][:, self.i] = deaths
            self.reported_cases[country_iter][:, self.i] = reported_cases
            self.hospitalizations[country_iter][:, self.i] = hospitalizations
            self.intensiveCare[country_iter][:, self.i] = intensiveCare

            

    
    def write_predictions(self, path_results):
        # saving:
            # deaths
            # hospitalizations
            # intensiveCare
            # reported_cases
            # cases
        for country_iter in self.countries:
            np.savetxt(path_results + 'predictions' + "_deaths_" + country_iter + ".txt", self.deaths[country_iter][:, :(self.i + 1)])
            np.savetxt(path_results + 'predictions' + "_reported_cases_" + country_iter + ".txt", self.reported_cases[country_iter][:, :(self.i + 1)])
            np.savetxt(path_results + 'predictions' + "_cases_" + country_iter + ".txt", self.cases[country_iter][:, :(self.i + 1)])
            np.savetxt(path_results + 'predictions' + "_hospitalizations_" + country_iter + ".txt", self.hospitalizations[country_iter][:, :(self.i + 1)])
            np.savetxt(path_results + 'predictions' + "_intensiveCare_" + country_iter + ".txt", self.intensiveCare[country_iter][:, :(self.i + 1)])
            np.savetxt(path_results + 'predictions' + "_Rt_" + country_iter + ".txt", self.Rt[country_iter][:, :(self.i + 1)])

            if self.model == 'infections':
                np.savetxt(path_results + 'predictions' + "_infections_" + country_iter + ".txt", self.infections[country_iter][:, :(self.i + 1)])


    def set_state(self, path, chain):
        # if self.oos_data is not None:
        if hasattr(self, "oos_data"):
            raise NotImplementedError("Not yet implemented for the case where there is an OOS country")
        for country_iter in self.countries:
            infections = np.loadtxt(f'{path}/{chain}_predictions_infections_{country_iter}.txt')
            deaths = np.loadtxt(f'{path}/{chain}_predictions_deaths_{country_iter}.txt')
            repCases = np.loadtxt(f'{path}/{chain}_predictions_reported_cases_{country_iter}.txt')
            cases = np.loadtxt(f'{path}/{chain}_predictions_cases_{country_iter}.txt')
            hospitalizations = np.loadtxt(f'{path}/{chain}_predictions_hospitalizations_{country_iter}.txt')
            intensiveCare = np.loadtxt(f'{path}/{chain}_predictions_intensiveCare_{country_iter}.txt')
            Rt = np.loadtxt(f'{path}/{chain}_predictions_Rt_{country_iter}.txt')

            # shape should be the same for all, if here something breaks this should not be used
            if not (deaths.shape == repCases.shape == cases.shape == hospitalizations.shape == intensiveCare.shape == Rt.shape == infections.shape):
                raise ValueError("Shapes of the read posterior predictions is not the same. Chains should be at the same state for all!")

            # atm we assume that the number of rows of the predictions is the same as in the attributes of self here
            # if self.deaths.shape[0] != deaths.shape[0]:
                # raise ValueError("Length of the times series is not the same as the provided one in the data!")

            n_predicitons_made = infections.shape[1]

            self.infections[country_iter][:, :n_predicitons_made] = infections
            self.deaths[country_iter][:, :n_predicitons_made] = deaths
            self.reported_cases[country_iter][:, :n_predicitons_made] = repCases
            self.cases[country_iter][:, :n_predicitons_made] = cases
            self.hospitalizations[country_iter][:, :n_predicitons_made] = hospitalizations
            self.intensiveCare[country_iter][:, :n_predicitons_made] = intensiveCare
            self.Rt[country_iter][:, :n_predicitons_made] = Rt

            self.i += n_predicitons_made - 1  # -1 here correct???
            print("Other attributes to set?????")
