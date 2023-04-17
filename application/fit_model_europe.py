import os
import numpy as np
import pandas as pd
from toolz.dicttoolz import valmap
import time
import scipy.stats as stats
import multiprocessing

import sys
sys.path.append('../model')
from MCMC import MCMC
from basics import load_data, load_time_distribution, load_reporting_delay_distributions


PATH = os.getcwd()
NB_CHAINS = 8
CORES = NB_CHAINS

np.random.seed(seed=123)


def define_prior_values():
    """
    Function just returns a dict with values of prior distributions
    """
    prior_parameters = {'alpha': {'mean': None, 'sd': None},
                        'alpha_mean': {'mean': 0, 'sd': 0.3},
                        'alpha_season_mean': {'mean': 0, 'sd': 0.3},
                        'alpha_sd': {'mean': 0, 'sd': 0.015},
                        'tau': {'mean': None, 'sd': None},
                        'tau_mean': {'shape': 10 / 1, 'scale': 1},
                        'tau_sd': {'mean': 0, 'sd': 10},
                        'R0': {'mean': None, 'sd': None},
                        'R0_mean': {'mean': 3.25, 'sd': 0.05},
                        # 'R0_sd': {'mean': 0, 'sd': 0.035},
                        'R0_sd': {'mean': 0, 'sd': 0.01},
                        'rho': {'a': 1, 'b': 1, 'scale': 3},
                        'phi_infections': {'inv_mean': 0, 'inv_sd': 0.015},
                        'phi_deaths': {'inv_mean': 0, 'inv_sd': 0.015},
                        'phi_hospitalizations': {'inv_mean': 0, 'inv_sd': 0.015},
                        'phi_intensiveCare': {'inv_mean': 0, 'inv_sd': 0.015},
                        'phi_rep_cases': {'inv_mean': 0, 'inv_sd': 0.015},
                        'beta_alpha': {'mean': 0.6, 'sd': 0.01},
                        'beta_delta': {'mean': 1.5, 'sd': 0.01},
                        'beta_sat': {'min': 0, 'max': 10},
                        'beta_sun': {'min': 0, 'max': 10},
                        'beta_mon': {'min': 0, 'max': 10},
                        'beta_tue': {'min': 0, 'max': 10},
                        'beta_wed': {'min': 0, 'max': 10},
                        'beta_fri': {'min': 0, 'max': 10},
                        'betaD_sat': {'min': 0, 'max': 10},
                        'betaD_sun': {'min': 0, 'max': 10},
                        'betaD_mon': {'min': 0, 'max': 10},
                        'betaD_tue': {'min': 0, 'max': 10},
                        'betaD_wed': {'min': 0, 'max': 10},
                        'betaD_fri': {'min': 0, 'max': 10},
                        'piH': {'a': 1, 'b': 1, 'scale': 10},
                        'piHicu': {'a': 1, 'b': 1, 'scale': 10}
                        }
    return prior_parameters


def define_proposal_sds():
    """
    Function just returns dict with initial proposal_sds
    """
    proposal_sds = {'alpha': 0.05,
                    'tau': 3,
                    'rho': 0.05,
                    'R0': 0.06,
                    'cases': 0.2,
                    'phi_infections': 0.005,
                    'phi_deaths': 0.01,
                    'phi_hospitalizations': 0.001,
                    'phi_intensiveCare': 0.001,
                    'phi_rep_cases': 0.002,
                    'beta_alpha': 0.05,
                    'beta_delta': 0.05,
                    'beta_sat': 0.15,
                    'beta_sun': 0.15,
                    'beta_mon': 0.15,
                    'beta_tue': 0.05,
                    'beta_wed': 0.05,
                    'beta_fri': 0.05,
                    'betaD_sat': 0.05,
                    'betaD_sun': 0.05,
                    'betaD_mon': 0.05,
                    'betaD_tue': 0.05,
                    'betaD_wed': 0.05,
                    'betaD_fri': 0.05,
                    'piH': 0.01,
                    'piHicu': 0.01
                    }
    return proposal_sds


def define_start_values(file_name='data_europe.csv', rnd_seed=-1):
    """
    Function defines start values
    """

    np.random.seed(seed=rnd_seed)

    # load data, will be used afterwards to init start values
    data_path = f'{PATH}/../data/'
    # read one dataset, required to iterate over contries etc
    data = load_data(f'{data_path}/{file_name}')
    # data.head().T

    # make first eho_period longer to get better power
    # rhotmp = data.rho_period.values
    # rhotmp[(rhotmp == 1) | (rhotmp == 2)] = 2
    # data.rho_period = rhotmp
    # data.rho_period = data.rho_period - 1


    # drop summer since it is the reference and rename seasons to correct names
    data.drop('summer', axis=1, inplace=True)
    data.rename(columns={'winter': 'seasonWinter',
                         'spring': 'seasonSpring',
                         'autumn': 'seasonAutumn'},
                inplace=True)

    # define intervention names
    used_interventions = ['closingSchools',
                          'restrictGatherings',
                          'lockdown',
                          'subsequentLockdown',
                          'generalBehavioralChanges']
    used_interventions.extend(['seasonWinter', 'seasonSpring', 'seasonAutumn'])

    print('Running the model with the following interventions:')
    print(used_interventions)

    # # define start values as dict
    # define random functions to initialize differenz start values
    rnd = lambda: stats.uniform(loc=0.9, scale=0.2).rvs(1)[0]

    alpha_starts = {f'alpha_{intervention}': 0.2 * rnd() for intervention in used_interventions}
    alpha_starts['alpha_sd'] = 0.02

    # # start vals
    start_vals = {'alpha': alpha_starts,
                  'phi_infections': 5 * rnd(),
                  'phi_deaths': 5 * rnd(),
                  'phi_hospitalizations': 5,
                  'phi_intensiveCare': 5,
                  'phi_rep_cases': 5 * rnd(),
                  'beta_alpha': 0.6,
                  'beta_delta': 1.5,
                  'tau': {'tau_' + cc: 5 * rnd() for cc in data.country.unique()},
                  'piH': {'piH_' + cc: 0.9 * rnd() for cc in data.country.unique()},
                  'piHicu': {'piHicu_' + cc: 0.25 * rnd() for cc in data.country.unique()},
                  'R0': {'R0_' + country_key: 3.25 * rnd() for country_key in data['country'].unique()}
                  }

    # init mean effect for R0 and tau
    start_vals['R0'].update({'R0_mean': 3.28 * rnd(), 'R0_sd': 0.1 * rnd()})
    start_vals['tau'].update({'tau_mean': 5 * rnd(), 'tau_sd': 5 * rnd()})

    # init rho
    rhos = {}
    countries = data.country.unique()
    df_rhos = data[['rho_period', 'country']].drop_duplicates()
    for cc in countries:
        rhos['rho_' + cc] = {pp: 0.5 for pp in df_rhos[df_rhos.country == cc].rho_period}
    start_vals['rho'] = rhos

    # # init parameters for weekday effects
    beta_Xi_sats = {'beta_sat_' + country_key: 1.0 * rnd() for _, country_key in enumerate(data['country'].unique())}
    beta_Xi_suns = {'beta_sun_' + country_key: 1.0 * rnd() for _, country_key in enumerate(data['country'].unique())}
    beta_Xi_mons = {'beta_mon_' + country_key: 1.0 * rnd() for _, country_key in enumerate(data['country'].unique())}
    beta_Xi_tues = {'beta_tue_' + country_key: 1.0 * rnd() for _, country_key in enumerate(data['country'].unique())}
    beta_Xi_weds = {'beta_wed_' + country_key: 1.0 * rnd() for _, country_key in enumerate(data['country'].unique())}
    beta_Xi_fris = {'beta_fri_' + country_key: 1.0 * rnd() for _, country_key in enumerate(data['country'].unique())}
    betaD_Xi_sats = {'betaD_sat_' + country_key: 1.0 * rnd() for _, country_key in enumerate(data['country'].unique())}
    betaD_Xi_suns = {'betaD_sun_' + country_key: 1.0 * rnd() for _, country_key in enumerate(data['country'].unique())}
    betaD_Xi_mons = {'betaD_mon_' + country_key: 1.0 * rnd() for _, country_key in enumerate(data['country'].unique())}
    betaD_Xi_tues = {'betaD_tue_' + country_key: 1.0 * rnd() for _, country_key in enumerate(data['country'].unique())}
    betaD_Xi_weds = {'betaD_wed_' + country_key: 1.0 * rnd() for _, country_key in enumerate(data['country'].unique())}
    betaD_Xi_fris = {'betaD_fri_' + country_key: 1.0 * rnd() for _, country_key in enumerate(data['country'].unique())}

    # for france and spain specify better start values sice we have issues with numerical stability
    beta_Xi_sats['beta_sat_France'] = 0.3
    beta_Xi_suns['beta_sun_France'] = 0.3
    beta_Xi_sats['beta_sat_Spain'] = 0.5
    beta_Xi_suns['beta_sun_Spain'] = 0.1

    start_vals['beta_sat'] = beta_Xi_sats
    start_vals['beta_sun'] = beta_Xi_suns
    start_vals['beta_mon'] = beta_Xi_mons
    start_vals['beta_tue'] = beta_Xi_tues
    start_vals['beta_wed'] = beta_Xi_weds
    start_vals['beta_fri'] = beta_Xi_fris
    start_vals['betaD_sat'] = betaD_Xi_sats
    start_vals['betaD_sun'] = betaD_Xi_suns
    start_vals['betaD_mon'] = betaD_Xi_mons
    start_vals['betaD_tue'] = betaD_Xi_tues
    start_vals['betaD_wed'] = betaD_Xi_weds
    start_vals['betaD_fri'] = betaD_Xi_fris

    # # add fixed parameters
    # time shifting distributions
    fixed_parameters = {'gamma': load_time_distribution(f'{data_path}/gamma_generation_time.csv'),
                        'xi_C': load_time_distribution(f'{data_path}/Xi_C_incubation_period.csv'),
                        'xi_R': load_reporting_delay_distributions(f'{data_path}/Xi_R_reporting_times_weekdays_estimated_lgl.csv'),
                        'xi_D': load_reporting_delay_distributions(f'{data_path}/Xi_D_symptoms_to_death_weekdays.csv')
                        }
    XIHs = pd.read_csv(f'{data_path}/XiH_all.csv')
    fixed_parameters['xi_H'] = XIHs.wards.values
    fixed_parameters['xi_Hicu'] = XIHs.icu.values

    # test if length is enough
    d_lengths = data.groupby("country").apply(lambda df: df.shape[0]).to_dict()
    xis_lengths = valmap(lambda d: d.shape[0], fixed_parameters)

    for cc in d_lengths:
        for xx in xis_lengths:
            if d_lengths[cc] + 14 > xis_lengths[xx]:
                raise ValueError(f'length of data (assumes additional 14 days for prediction) in country {cc} is shorter than {xx}')

    # define population
    start_vals.update({'N_Austria': 8859000,
                       'N_Belgium': 11460000,
                       'N_Czechia': 10650000,
                       'N_Denmark': 5806000,
                       'N_Finland': 5518000,
                       'N_France': 67060000,
                       'N_Germany': 83020000,
                       'N_Greece': 10720000,
                       'N_Hungary': 9773000,
                       'N_Ireland': 4904000,
                       'N_Italy': 60360000,
                       'N_Netherlands': 17280000,
                       'N_Norway': 5328000,
                       'N_Poland': 37970000,
                       'N_Portugal': 10280000,
                       'N_Slovenia': 2081000,
                       'N_Spain': 46940000,
                       'N_Sweden': 10230000,
                       'N_Switzerland': 8545000,
                       'N_UnitedKingdom': 66650000
                       })

    fixed_parameters['probability_reinfection'] = 0.159
    fixed_parameters['correction_first_vaccination'] = 0.5
    fixed_parameters['correction_second_vaccination'] = 0.35

    # add pi_D
    pi_D_all = data[['country', 'ifr_t_m']]
    pi_D = {cc: pi_D_all[pi_D_all.country == cc].ifr_t_m.values for cc in pi_D_all.country.unique()}
    fixed_parameters['pi_D'] = pi_D

    # add correction factor hospitalizations
    cf_hospitalization = valmap(lambda x: x / x[0], pi_D)
    cfnames = cf_hospitalization.keys()
    correction_hosp = {'correction_hospitalization_' + name: cf_hospitalization[name] for name in cfnames}
    fixed_parameters.update(correction_hosp)

    start_vals.update(fixed_parameters)
    return start_vals


def run_chain(chain, path, file_name, rnd_seed):
    """
    Run one chain

    Function initializes a MCMC object and lets the chain run

    :param chain: An string indicating the j-th chain, e.g. 'chain1'
    :param path: A string defining the path where the reults are written
    :param file_name: A string defining the used data
    :param rnd_seed: An integer which is used to initialize a random state for RNG
    """

    # set random state for the sampling scheme
    np.random.seed(rnd_seed)

    # load data
    data = load_data(f'{PATH}/../data/{file_name}')
    data.rename(columns={'winter': 'seasonWinter',
                         'spring': 'seasonSpring',
                         'autumn': 'seasonAutumn'},
                inplace=True)

    # make first eho_period longer to get better power
    # rhotmp = data.rho_period.values
    # rhotmp[(rhotmp == 1) | (rhotmp == 2)] = 2
    # data.rho_period = rhotmp
    # data.rho_period = data.rho_period - 1

    # nowcasting not necessary, therefore 1
    data['pi_nd'] = 1
    data['pi_nc'] = 1

    # set dtypes
    cols = data.columns[2:]
    for col in cols:
        data[col] = pd.to_numeric(data[col])
    # data.dtypes

    # reset index
    data = data.reset_index(drop=True)  # not necessary but I am superstitious

    # get the used interventions
    used_interventions = ['closingSchools',
                          'restrictGatherings',
                          'lockdown',
                          'subsequentLockdown',
                          'generalBehavioralChanges'
                          ]
    used_interventions.extend(['seasonWinter', 'seasonSpring', 'seasonAutumn'])

    prior_parameters = define_prior_values()
    proposal_sds = define_proposal_sds()

    seed_generation_start_vals = np.random.randint(1, 100, 1)
    start_values = {chain: define_start_values(file_name, seed_generation_start_vals)}

    model_specification = {'model': 'infections',
                           'interventions': used_interventions,
                           'hierarchical_interventions': True,
                           'adapt_reporting_weekend': True,
                           }

    # Define fixed parameters
    # ATTENTION! fixed_params is a string list defining the fixed parameters
    fixed_params = ['gamma', 'pi_D', 'xi_D', 'xi_C', 'xi_H', 'xi_Hicu', 'xi_R']
    fixed_params.append('probability_reinfection')
    fixed_params.append('correction_first_vaccination')
    fixed_params.append('correction_second_vaccination')
    fixed_params.extend(['correction_hospitalization_' + cc for cc in data.country.unique()])

    for country_iter in data.country.unique():
        fixed_params.append('N_' + country_iter)

    path_results = f'{path}/'
    if os.path.exists(path_results):
        print('Attention! Result path exists. Results may get overwritten!')
    else:
        os.makedirs(path_results)

    print('==================================================================')
    print('Writing results to:')
    print(path_results)
    print('==================================================================')

    # check if countries do not have implemented all interventions,
    # if not delete them, this is actually not necessary for simulated data
    exceptions_intervention = []
    for inter in used_interventions:
        for cc in data.country.unique():
            if np.all(data[data.country == cc][inter] == 0):
                exceptions_intervention.append('alpha_' + inter + '_' + cc)

    # run chains
    algo = MCMC(data=data,
                model_specification=model_specification,
                informative_priors=True,
                path_results=path_results,
                proposal_sd=proposal_sds,
                prior_parameters=prior_parameters,
                start_values=start_values,
                chain=chain,
                nb_future_values=1,
                fix_latent_variable=False,
                fixed_parameters=fixed_params,
                exceptions_intervention=exceptions_intervention,
                oos_country=None
                )

    # set good proposals (optional)
    algo.set_proposal_sds("../data/proposal_sds_real_data.json")

    algo.run_adaptive_algorithm(iterations=50000,
                                burnin=20000,
                                adaptive_phases=10,
                                thin=100,
                                prediction_interval=300)

    # short version. Can be used check whether the model is working and writes
    # results
    # algo.run_adaptive_algorithm(iterations=11,
                                # burnin=1,
                                # adaptive_phases=0,
                                # thin=1,
                                # prediction_interval=3)


if __name__ == '__main__':
    seeds = np.random.randint(1, 999999, size=NB_CHAINS)
    path = f'{PATH}/results/main_results/'
    data_file = '/real_data.csv'

    t = time.time()
    # one chain for debug
    # import pudb; pu.db
    # run_chain('chain1', path, data_file, 123)

    np.seterr('ignore')
    with multiprocessing.Pool(CORES) as pool:
        res = pool.starmap(run_chain, [(f'chain{i+1}', path, data_file, rnd_seed) for i, rnd_seed in enumerate(seeds)])

    print("Full calculation time: " + str(time.time() - t))
