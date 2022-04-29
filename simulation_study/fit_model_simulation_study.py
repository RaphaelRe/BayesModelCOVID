import os
import numpy as np
import pandas as pd
from toolz.dicttoolz import valmap
import time
import scipy.stats as stats
from string import ascii_uppercase

import sys
sys.path.append('../model')
from MCMC import MCMC
from basics import load_data, load_time_distribution, load_reporting_delay_distributions


PATH = os.getcwd()
NB_DATASETS = 2
CORES = 2

np.random.seed(seed=123)


def define_prior_values():
    """
    Function just returns a dict with values of prior distributions
    """
    prior_parameters = {'alpha': {'mean': None, 'sd': None},
                        'alpha_mean': {'mean': 0, 'sd': 1},
                        'alpha_sd': {'mean': 0, 'sd': 0.1},
                        'tau': {'mean': None, 'sd': None},
                        'tau_mean': {'shape': 10 / 1, 'scale': 1},
                        'tau_sd': {'mean': 0, 'sd': 4},
                        'R0': {'mean': None, 'sd': None},
                        'R0_mean': {'mean': 3.25, 'sd': 1},
                        'R0_sd': {'mean': 0, 'sd': 0.1},
                        'rho': {'a': 1, 'b': 1, 'scale': 3},
                        'phi_infections': {'inv_mean': 0, 'inv_sd': 0.015},
                        'phi_deaths': {'inv_mean': 0, 'inv_sd': 0.015},
                        'phi_hospitalizations': {'inv_mean': 0, 'inv_sd': 0.015},
                        'phi_intensiveCare': {'inv_mean': 0, 'inv_sd': 0.015},
                        'phi_rep_cases': {'inv_mean': 0, 'inv_sd': 0.015},
                        'beta_alpha': {'mean': 0.4, 'sd': 0.01},
                        'beta_delta': {'mean': 1.4, 'sd': 0.02},
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
                    'tau': 5,
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


def define_start_values(file_name='data_sim_5NPIs_1.csv'):
    """
    Function defines start values
    """

    # load data, will be used afterwards to init start values (e.g. country iter)
    data_path = f'{PATH}/../data/'

    data = load_data(f'{data_path}/simulated_data/{file_name}')

    # not required in simulated data, but here for completeness
    # data.rename(columns={'winter': 'seasonWinter',
                         # 'spring': 'seasonSpring',
                         # 'autumn': 'seasonAutumn'},
                # inplace=True)

    # get the used interventions
    used_interventions = [f'NPI{i}' for i in range(1, 6)]
    # used_interventions.extend(['seasonWinter', 'seasonSpring', 'seasonAutumn'])

    # define random functions to initialize differenz start values
    rnd = lambda: stats.uniform(loc=0.5, scale=1.5).rvs(1)[0]
    rnd_small = lambda: stats.uniform(loc=0.9, scale=0.2).rvs(1)[0]

    alpha_starts = {f'alpha_{intervention}': 0.2 * rnd() for intervention in used_interventions}
    alpha_starts['alpha_sd'] = 0.02

    # # start vals
    start_vals = {'alpha': alpha_starts,
                  'phi_infections': 50 * rnd_small(),
                  'phi_deaths': 50 * rnd_small(),
                  'phi_hospitalizations': 50 * rnd_small(),
                  'phi_intensiveCare': 50 * rnd_small(),
                  'phi_rep_cases': 50 * rnd_small(),
                  'beta_alpha': 0.4,
                  'beta_delta': 1.4,
                  'piH': {'piH_' + cc: 0.9 * rnd_small() for cc in data.country.unique()},
                  'piHicu': {'piHicu_' + cc: 0.25 * rnd_small() for cc in data.country.unique()},
                  'R0': {'R0_' + country_key: 3.25 * rnd_small() for country_key in data['country'].unique()},
                  'tau': {'tau_' + country_key: 30 * rnd_small() for country_key in data['country'].unique()}
                  }

    # # init mean effect for R0
    start_vals['R0'].update({'R0_mean': 3.25, 'R0_sd': 0.01})
    start_vals['tau'].update({'tau_mean': 30 * rnd(), 'tau_sd': 5})

    # # init rho
    rhos = {}
    countries = data.country.unique()
    df_rhos = data[['rho_period', 'country']].drop_duplicates()
    for cc in countries:
        rhos['rho_' + cc] = {pp: 0.5 for pp in df_rhos[df_rhos.country == cc].rho_period}
    start_vals['rho'] = rhos

    # # init parameters for weekday effects
    beta_Xi_sats = {'beta_sat_' + country_key: 1.0 * rnd_small() for _, country_key in enumerate(data['country'].unique())}
    beta_Xi_suns = {'beta_sun_' + country_key: 1.0 * rnd_small() for _, country_key in enumerate(data['country'].unique())}
    beta_Xi_mons = {'beta_mon_' + country_key: 1.0 * rnd_small() for _, country_key in enumerate(data['country'].unique())}
    beta_Xi_tues = {'beta_tue_' + country_key: 1.0 * rnd_small() for _, country_key in enumerate(data['country'].unique())}
    beta_Xi_weds = {'beta_wed_' + country_key: 1.0 * rnd_small() for _, country_key in enumerate(data['country'].unique())}
    beta_Xi_fris = {'beta_fri_' + country_key: 1.0 * rnd_small() for _, country_key in enumerate(data['country'].unique())}
    betaD_Xi_sats = {'betaD_sat_' + country_key: 1.0 * rnd_small() for _, country_key in enumerate(data['country'].unique())}
    betaD_Xi_suns = {'betaD_sun_' + country_key: 1.0 * rnd_small() for _, country_key in enumerate(data['country'].unique())}
    betaD_Xi_mons = {'betaD_mon_' + country_key: 1.0 * rnd_small() for _, country_key in enumerate(data['country'].unique())}
    betaD_Xi_tues = {'betaD_tue_' + country_key: 1.0 * rnd_small() for _, country_key in enumerate(data['country'].unique())}
    betaD_Xi_weds = {'betaD_wed_' + country_key: 1.0 * rnd_small() for _, country_key in enumerate(data['country'].unique())}
    betaD_Xi_fris = {'betaD_fri_' + country_key: 1.0 * rnd_small() for _, country_key in enumerate(data['country'].unique())}
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

    # test if length of loaded Xis is is enough
    d_lengths = data.groupby("country").apply(lambda df: df.shape[0]).to_dict()
    xis_lengths = valmap(lambda d: d.shape[0], fixed_parameters)

    for cc in d_lengths:
        for xx in xis_lengths:
            if d_lengths[cc] + 14 > xis_lengths[xx]:  # 14 here arbitrary, 1 would be enough
                raise ValueError(f'length of data (assumes additional 14 days for prediction) in country {cc} is shorter than {xx}')

    # define population
    fixed_parameters.update({f'N_country{ascii_uppercase[i]}': 1e8 for i in range(10)})  # same for all countries in simulation
    fixed_parameters['probability_reinfection'] = 0.16
    fixed_parameters['correction_first_vaccination'] = 0.6
    fixed_parameters['correction_second_vaccination'] = 0.3

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


def run_chains(path, file_name='data_sim_5NPIs_1.csv', nb_chains=1):
    """
    Functions serves to run chains for one given dataset
    """
    # load data
    data = load_data(f'{PATH}/../data/simulated_data/{file_name}')
    # not required, no seasons
    # data.rename(columns={'winter': 'seasonWinter',
                         # 'spring': 'seasonSpring',
                         # 'autumn': 'seasonAutumn'},
                # inplace=True)

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
    used_interventions = [f'NPI{i}' for i in range(1, 6)]
    # used_interventions.extend(['seasonWinter', 'seasonSpring', 'seasonAutumn'])

    prior_parameters = define_prior_values()
    proposal_sds = define_proposal_sds()
    start_values = {f'chain{i}': define_start_values(file_name) for i in range(1, nb_chains + 1)}

    # check if pi_d and the data have the correct shape
    len_piD = start_values['chain1']['pi_D'][next(iter(start_values['chain1']['pi_D']))].shape[0]
    len_data = data[data.country == next(iter(start_values['chain1']['pi_D']))].shape[0]

    if len_piD != len_data:
        raise ValueError('Length of pi_D must be the same as the data for a country')

    model_specification = {'model': 'infections',
                           'interventions': used_interventions,
                           'hierarchical_interventions': True,
                           'adapt_reporting_weekend': True,
                           }

    # Define fixed parameters
    # fixed_params here defines which parameters are fixed! (In the function before it is a dict with the used values)
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

    print('writing results to:')
    print(path_results)

    # check if countries do not have implemented all interventions,
    # this is actually not necessary for simulated data
    exceptions_intervention = []
    for inter in used_interventions:
        for cc in data.country.unique():
            if np.all(data[data.country == cc][inter] == 0):
                exceptions_intervention.append('alpha_' + inter + '_' + cc)

    # run chains
    for i in range(1, nb_chains + 1):
        algo = MCMC(data=data,
                    model_specification=model_specification,
                    informative_priors=True,
                    path_results=path_results,
                    proposal_sd=proposal_sds,
                    prior_parameters=prior_parameters,
                    start_values=start_values,
                    chain=f'chain{i}',
                    nb_future_values=1,
                    fix_latent_variable=False,
                    fixed_parameters=fixed_params,
                    exceptions_intervention=exceptions_intervention,
                    oos_country=None
                    # oos_country=oos_country
                    )

        algo.run_adaptive_algorithm(iterations=50001, burnin=20000,
                                    adaptive_phases=10, thin=100,
                                    prediction_interval=300)

        # very short period to test, whether algo is sampling everything
        # algo.run_adaptive_algorithm(iterations=11, burnin=1,
                                    # adaptive_phases=0, thin=1,
                                    # prediction_interval=3)


if __name__ == '__main__':
    import multiprocessing
    t = time.time()

    # for debugging
    # import pudb; pu.db
    # run_chains('results')

    np.seterr('ignore')
    with multiprocessing.Pool(CORES) as pool:
        res = pool.starmap(run_chains, [(f'results/res_dataset_{i}', f'data_sim_5NPIs_{i}.csv', 2) for i in range(1, NB_DATASETS + 1)])
    print("Full calculation time: " + str(time.time() - t))
