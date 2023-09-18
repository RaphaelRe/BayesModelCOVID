import time
import json

from basics import (split_oos_data, dictionarize_data)
from LatentVariables import LatentVariables
from Parameter import Parameter
from Parameter import FixedParameter
from Parameter import ParameterVector
from Prediction import Prediction


class MCMC(object):
    """
    Main MCMC class

    Parameters are:
    :data: pandas DataFrame
    :model: str, should either be cases or infections
    :informative_priors: bool - currently always True
    :path_results: str, the path where the results get stored
    :proposal_sd: dict, a dictionary with information about all
                       standard deviations for proposals
    :prior_parameters: dict, a dictionary with information about all
                            the prior parameters
    :start_values: dict, a dictionary with info where all parameterts
                        get initialized
    :chain: str, name of the chain, useful to run several chains in parallel
    :nb_future_values: number of predictions
    :epidemic_start: day when the epidemic starts. The cases before will be used as a seed. Note that counting starts at day 1, not zero
    """

    def __init__(self, data, model_specification,
                 informative_priors=False,
                 path_results='/home/rrehms/Corona_Modelling/Corona_modeling/corona_modeling/python/results/',
                 proposal_sd={},
                 prior_parameters={},
                 start_values={},
                 chain='chain1', nb_future_values=10, write_predictions=True,
                 fix_latent_variable=False, fixed_parameters=['gamma', 'pi_D', 'xi_D', 'xi_C', 'xi_H', 'xi_Hicu', 'xi_R'],
                 epidemic_start=2, exceptions_intervention=[],
                 oos_country=None, nb_days_oos_init=42):


        if oos_country is not None:
            self.data, self.oos_data = split_oos_data(data, oos_country, nb_days_oos_init)
            start_values[chain]['pi_D'][oos_country + '_oos'] = start_values[chain]['pi_D'][oos_country]
            start_values[chain]['pi_D'][oos_country] = start_values[chain]['pi_D'][oos_country][:self.data[self.data.country == oos_country].shape[0]]
            start_values[chain]['correction_hospitalization_' + oos_country + '_oos'] = start_values[chain]['correction_hospitalization_' + oos_country]
            start_values[chain]['correction_hospitalization_' + oos_country] = start_values[chain]['correction_hospitalization_' + oos_country][:self.data[self.data.country == oos_country].shape[0]]
            fixed_parameters.append('correction_hospitalization_' + oos_country + '_oos')
        else:
            self.data = data
            self.oos_data = None

        self.epidemic_start = epidemic_start
        self.model = model_specification['model']
        self.interventions = model_specification['interventions']
        self.hierarchical_interventions = model_specification['hierarchical_interventions']  # logical - definiert whether hierarchisches alpha is used 

        self.adapt_reporting_weekend = model_specification['adapt_reporting_weekend']

        self.path = path_results + chain + '_'
        
        self.fix_latent_variable = fix_latent_variable 
        self.fixed_parameters = fixed_parameters

        self.filename = self.path + 'statistics' + '.txt'
        self.write_predictions = write_predictions 
        self.initialize_parameters(proposal_sd,
                                   prior_parameters,
                                   start_values,
                                   chain,
                                   informative_priors)

        # some countries dont have all interventions
        for parm_iter in exceptions_intervention:
            split = parm_iter.split('_')
            intervention_tmp = split[1]
            country_tmp = split[2]
            del self.parameters['alpha'].parameters['alpha_' + intervention_tmp].parameters['alpha_' + intervention_tmp + '_' + country_tmp]

        self.latent_variables =  LatentVariables(self.model, self.data,
                parameter_values = self.get_current_values(),
                proposal_distribution = 'Uniform',
                proposal_range = proposal_sd['cases'], start = epidemic_start, 
                fix_latent_variable = fix_latent_variable , adapt_reporting_weekend = self.adapt_reporting_weekend
                )

        self.predictions = Prediction(self.model, self.data, nb_future_values, parameter_values=self.get_current_values(),oos_data=self.oos_data)


    def get_current_values(self):
        current_values = {key: self.parameters[key].get_current_value() for key in self.parameters}
        return(current_values)

    def reset_values(self, iterations):
        for key in self.parameters.keys():
            self.parameters[key].reset_values(iterations)
        self.latent_variables.reset_values(iterations)


    def update_chains(self):
        for key in self.parameters.keys():
            self.parameters[key].update(self.get_current_values(),
                                    self.latent_variables)
        if not self.fix_latent_variable:
            self.update_latentVariables()

    def update_latentVariables(self):
        current_parameter_values = self.get_current_values()
        self.latent_variables.update(current_parameter_values, start=self.epidemic_start)

    def adapt_proposals(self, nb_iterations, phase):
        for key in self.parameters:
            self.parameters[key].adapt_proposal(nb_iterations, phase)
        if not self.fix_latent_variable:
            self.latent_variables.adapt_proposal(nb_iterations, phase)


    def run_adaptive_phase(self, nb_iterations, nb_phases):
        print("starting adaptive phase......")
        for phase in range(nb_phases):
            print('Adaptive phase '+ str(phase))
            for _ in range(nb_iterations):
                self.update_chains()
            self.adapt_proposals(nb_iterations, phase)


    def run_burnin(self, nb_burnin):
        print('Start burnin-phase')
        self.reset_values(nb_burnin)
        for i in range(nb_burnin):
            self.update_chains()
            if i > 1 and i % 1000 == 0:
                print(f'Iteration {i}/{nb_burnin} [BURNIN] finished')
        print('End Burnin-phase')


    def run_algorithm(self, nb_iterations, thin=1, prediction_interval=300, save_chains=True):
        # check whether the Prediciton object has enough space
        planned_predictions = nb_iterations / prediction_interval
        if any(map(lambda x: x.shape[1] < planned_predictions, self.predictions.cases.values())):
            raise ValueError("The algorithm will sample more posterior \n \
                              predictions than planned. Raise the number the \n \
                              sample_size parameter for the Prediction object \n \
                              or the prediction_interval parameter!")

        print('Start algorithm')
        self.reset_values(nb_iterations)
        for i in range(nb_iterations):
            self.update_chains()

            if i > 0 and i % prediction_interval == 0:
                current_prior_values = {}
                current_prior_values['R0'] = self.parameters['R0'].get_current_prior_value()
                current_prior_values['alpha'] = {alpha_key: self.parameters['alpha'].parameters[alpha_key].get_current_prior_value() for alpha_key in self.parameters['alpha'].parameters.keys()}
                self.predictions.make_predictions(self.get_current_values(),
                                                  self.latent_variables.get_values(),
                                                  current_prior_values
                                                  )
                

            if i > 0 and i % 1000 == 0:
                print(f'Iteration {i}/{nb_iterations} [SAMPLING] finished')
                self.write_statistics()
                if save_chains:
                    self.write_chains(thin)
        # final write
        if save_chains:
            self.write_chains(thin)


    def run_adaptive_algorithm(self, iterations, burnin, adaptive_phases=100,
                               thin=1, prediction_interval=300, save_chains=True):
        """
        Run the full chain with adaptive phase, burnin and sampling

        This function calls sequentially run_adaptive_phase, run_burnin, run_algorithm

        :param iterations: The number of sampling steps of the chain
        :param burnin: The number of burnin iterations
        :param adaptive_phases: The number of adaptive phases, each phase runs 59 iterations
        :param thin: Specifies the thinning, for example thin=10 keeps only each 10th iteration
        :param prediction_interval: After how many samples should a posterior prediction be done? WARNING: This is coputationally expensive. Value should be > 100
        """
        start = time.perf_counter()
        self.run_adaptive_phase(50, adaptive_phases)
        print(time.perf_counter() - start)

        start = time.perf_counter()
        self.run_burnin(burnin)
        print(time.perf_counter() - start)

        start = time.perf_counter()
        self.run_algorithm(iterations, thin, prediction_interval, save_chains)
        print(time.perf_counter() - start)

    def write_chains(self, thin):
        for key in self.parameters:
            self.parameters[key].write_samples(self.path, thin=thin)

        if self.write_predictions:
            self.predictions.write_predictions(self.path)



    def write_statistics(self):
        with open(self.filename, 'w') as result:
            result.write('The final results of the algorithm are:')
            parameters = self.parameters
            for key in parameters:
                stats = parameters[key].get_statistics()
                result.write(key + ':')
                result.write(str(stats))
                result.write('\n \n')

            result.write('\n \n')



    def write_proposal_sds(self, name='proposal_sds.json'):
        proposal_sds = {}
        for key in self.parameters:
            proposal_sds[key] = self.parameters[key].get_proposal_sd()
        proposal_sds['lv'] = self.latent_variables.get_proposal_sd()

        with open(self.path + name, 'w') as results:
            results.write(json.dumps(proposal_sds))



    def set_proposal_sds(self, json_file):
        with open(json_file) as f:
            proposal_sds = json.load(f)
        for key in self.parameters:
            self.parameters[key].set_proposal_sd(proposal_sds[key])
        self.latent_variables.set_proposal_sd(proposal_sds['lv'])



    def set_state_algorithm(self, path: str, chain: str, nb_iterations: int, thin: int = 100):
        """
        Function to set the current state of the algorithm. It iterates
        over the parameters to set them to the state provided in the path
        and the chain accordingly. It also reuqires the info abot the length
        of the full samples to do. For example: the states provided are
        [1,2,1,2] and you want 10 samples in the end it sets the first 4
        samples on the values. and 6 samples are left to make (thin is 1 here).

        :param path: The path where to search for the np.array
        :param chain: Info which chain to load
        :param nb_iterations: The number of iterations to make for the full chain
        :param thin: Info about the thinning to use. A thinning of 10 implies a
            fill trjectory length of the parameter x*(10-1)+1
        """
        self.reset_values(nb_iterations)  # sets samples and acceptance in all params nd lv to initial state
        for key in self.parameters:
            self.parameters[key].set_state(path, chain, thin)

        self.latent_variables.set_state(self.get_current_values(), path, chain,
                                        thin)

        self.predictions.set_state(path, chain)



    def initialize_parameters(self, proposal_sd, prior_parameters, start_values, chain,
                              informative_priors,
                              precision={'sigma_R': 2,
                                         'rho': 2,
                                         'alpha': 2,
                                         'tau': 2,
                                         'R0': 2,
                                         'phi_rep_cases': 2,
                                         'phi_deaths': 2,
                                         'phi_hospitalizations': 2,
                                         'phi_intensiveCare': 2,
                                         'phi_infections': 2,
                                         'piH': 2,
                                         'piHicu': 2,
                                         'beta_D': 2,
                                         'beta_alpha': 2,
                                         'beta_delta': 2,
                                         'beta_sat': 2,
                                         'beta_sun': 2,
                                         'beta_mon': 2,
                                         'beta_tue': 2,
                                         'beta_wed': 2,
                                         'beta_fri': 2,
                                         'betaD_sat': 2,
                                         'betaD_sun': 2,
                                         'betaD_mon': 2,
                                         'betaD_tue': 2,
                                         'betaD_wed': 2,
                                         'betaD_fri': 2,
                                         }
                              ):
        parameter_names = start_values[chain].keys()
        data_dict = dictionarize_data(self.data)
        alpha_keys = self.interventions
        alpha_data = self.data
        if self.hierarchical_interventions:
            alpha_keys = {'interventions': alpha_keys,
                          'countries': data_dict.keys()
                          }
            alpha_data = data_dict

        model = self.model

        self.parameters = {}
        for parameter in parameter_names:
            if parameter in self.fixed_parameters:
                self.parameters[parameter] = FixedParameter(parameter,
                                                            start_values[chain][parameter])
            else:
                if isinstance(start_values[chain][parameter], dict):  # checks if parameter should be a parameterVector
                    if parameter == 'alpha':
                        self.parameters[parameter] = ParameterVector('alpha',
                                                                     alpha_keys,
                                                                     model,
                                                                     precision,
                                                                     proposal_sd,
                                                                     prior_parameters,
                                                                     start_values[chain][parameter],
                                                                     alpha_data,
                                                                     start=self.epidemic_start)
                    else:
                        self.parameters[parameter] = ParameterVector(parameter,
                                                                     data_dict.keys(),
                                                                     model,
                                                                     precision,
                                                                     proposal_sd,
                                                                     prior_parameters,
                                                                     start_values[chain][parameter],
                                                                     data_dict,
                                                                     start=self.epidemic_start)
                else:
                    data = self.data
                    if parameter == 'phi_infections':
                        renewal_model = True
                        death_model = False
                        hospitalization_model = False
                        intensiveCare_model = False
                        reporting_model = False
                    elif parameter in ['phi_hospitalizations']:
                        renewal_model = False
                        death_model = False
                        hospitalization_model = True
                        intensiveCare_model = False
                        reporting_model = False
                        data = dictionarize_data(data)
                    elif parameter in ['phi_intensiveCare']:
                        renewal_model = False
                        death_model = False
                        hospitalization_model = False
                        intensiveCare_model = True
                        reporting_model = False
                        data = dictionarize_data(data)
                    elif parameter == 'phi_deaths':
                        renewal_model = False
                        death_model = True
                        hospitalization_model = False
                        intensiveCare_model = False
                        reporting_model = False
                        data = dictionarize_data(data)
                    elif parameter == 'phi_rep_cases':
                        renewal_model = False
                        death_model = False
                        hospitalization_model = False
                        intensiveCare_model = False
                        reporting_model = True
                        data = dictionarize_data(data)
                    elif parameter in ['tau', 'beta_alpha', 'beta_delta']:
                        if parameter == 'tau':
                            raise RuntimeError("tau should be a ParameterVector")
                        renewal_model = True
                        death_model = False
                        hospitalization_model = False
                        intensiveCare_model = False
                        reporting_model = False
                        data = self.data
                    elif parameter == 'beta_D':
                        renewal_model = False
                        death_model = True
                        hospitalization_model = False
                        intensiveCare_model = False
                        reporting_model = False
                        data = dictionarize_data(data)

                    self.parameters[parameter] = Parameter(parameter,
                                                           model,
                                                           precision,
                                                           proposal_sd,
                                                           prior_parameters,
                                                           start_values[chain],
                                                           data,
                                                           start=self.epidemic_start,
                                                           renewal_model=renewal_model,
                                                           death_model=death_model,
                                                           hospitalization_model=hospitalization_model,
                                                           intensiveCare_model=intensiveCare_model,
                                                           reporting_model=reporting_model
                                                           )
