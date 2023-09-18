from numba import jit
import numpy as np
import pandas as pd
import copy
import datetime
import scipy.stats as stats
from scipy.special import gamma

def is_not_empty(x):
    if x:
        not_empty = True
    else:
        not_empty = False

    return(not_empty)


def load_time_distribution(filename):
    variables = ['index', 'description', 'probability_mass']
    data = pd.read_csv(filename, sep=',', quotechar='"', skiprows=1,
                       names=variables)
    result = data['probability_mass'].values
    return(result)


def load_reporting_delay_distributions(filename):
    variables = ['index', 'description', 1, 2, 3, 4, 5, 6, 7]

    data = pd.read_csv(filename, sep=',', quotechar='"', skiprows=1,
                       names=variables)
    result = data.iloc[:, 2:].values
    return(result)


def get_wday(date):
    weekday = datetime.datetime.weekday(date)
    return(weekday)


def load_data(name, data_type='simulated', model='infections', dictionarize=False):
    """ The function load_data takes the name of a csv file and loads it as a
    pandas data frame in Python.
    """
    if data_type in ['real', 'simulated']:
        data = pd.read_csv(name, sep=',', quotechar='"', header=0)
    else:
        return None

    data = data.rename(columns={"repC_t": "reported_cases",
                                "D_t": "deaths",
                                "C_t": "cases",
                                "I_t": "infections",
                                "H_t": "hospitalizations",
                                "Hicu_t": "intensiveCare"
                                }
                       )
    data.date = pd.to_datetime(data.date)
    data['weekday'] = data.date.apply(get_wday)
    if 'pi_nc' not in data:
        data['pi_nd'] = 1
        data['pi_nc'] = 1

    if dictionarize:
        data = dictionarize_data(data)

    return(data)


def dictionarize_data(data):
    ind = data.country.unique()
    res = {country: data[data['country'].values == country].copy(deep=True) for country in ind}
    return(res)


def split_data(data, country):
    """
    split the data in 2 datasets:
    one with all country except the specified one
    the second only with the specified one
    """
    data_2 = data[data.country == country]
    data_1 = data.drop(data.loc[data.country == country].index, inplace=False)

    if data_2.shape[0] == 0:
        raise Exception("OOS dataset has 0 rows! Probably the defined country is not in the full dataset")

    return data_1, data_2


def split_oos_data(data, oos_country, nb_days_oos_init):
    data_oos = data[data.country == oos_country]
    i = data[data.country == oos_country].index
    i_remove= i[nb_days_oos_init:].to_list()

    full_index = data.index.to_list()
    item_list = [e for e in full_index if e not in i_remove]
    data_trunc = data.iloc[item_list]

    return data_trunc, data_oos



def initialize_country_specific_Xi_R(parameter_values, country):
    # assumes that the first row is Monday, second Tuesday, etc.
    if parameter_values['xi_R'].shape[0] % 7 != 1:
        curr_shape = parameter_values['xi_R'].shape
        raise ValueError(f'% 7 of number of rows in Xi_R must be 1! i.e. n % 7 == 1; currently shape is {curr_shape}')
    result = {}
    xi_new = np.copy(parameter_values['xi_R'], order='F')
    result['values'] = xi_new
    tmp = result['values'].reshape(-1, order='A')  # only used to define the view. numpy just can create view on slices...trick only works when we get nb %1
    result['view_sat'] = tmp[5::7]
    result['view_sat'][:] = result['view_sat'] * parameter_values['beta_sat']['beta_sat_' + country]
    result['view_sun'] = tmp[6::7]
    result['view_sun'][:] = result['view_sun'] * parameter_values['beta_sun']['beta_sun_' + country]
    result['view_mon'] = tmp[0::7]
    result['view_mon'][:] = result['view_mon'] * parameter_values['beta_mon']['beta_mon_' + country]
    result['view_tue'] = tmp[1::7]
    result['view_tue'][:] = result['view_tue'] * parameter_values['beta_tue']['beta_tue_' + country]
    result['view_wed'] = tmp[2::7]
    result['view_wed'][:] = result['view_wed'] * parameter_values['beta_wed']['beta_wed_' + country]
    result['view_fri'] = tmp[4::7]
    result['view_fri'][:] = result['view_fri'] * parameter_values['beta_fri']['beta_fri_' + country]
    for j in range(7):
        vals = result['values'][:, j]
        result['values'][:, j] = vals / vals.sum()

    return result


def initialize_country_specific_Xi_D(parameter_values, country):
    # assumes that the first row is Monday, second Tuesday, etc.
    if parameter_values['xi_D'].shape[0] % 7 != 1:
        raise ValueError('% 7 of number of rows in Xi_D must be 1! i.e. n % 7 == 1')
    result = {}
    xi_new = np.copy(parameter_values['xi_D'], order='F') # new version - should work better, more control about memory allocation order
    result['values'] = xi_new
    tmp = result['values'].reshape(-1, order = 'A') # only used to define the view. numpy just can create view on slices...trick only works when we get nb %1
    result['view_sat'] = tmp[5::7]
    result['view_sat'][:] = result['view_sat'] * parameter_values['betaD_sat']['betaD_sat_' + country]
    result['view_sun'] = tmp[6::7]
    result['view_sun'][:] = result['view_sun'] * parameter_values['betaD_sun']['betaD_sun_' + country]
    result['view_mon'] = tmp[0::7]
    result['view_mon'][:] = result['view_mon'] * parameter_values['betaD_mon']['betaD_mon_' + country]
    result['view_tue'] = tmp[1::7]
    result['view_tue'][:] = result['view_tue'] * parameter_values['betaD_tue']['betaD_tue_' + country]
    result['view_wed'] = tmp[2::7]
    result['view_wed'][:] = result['view_wed'] * parameter_values['betaD_wed']['betaD_wed_' + country]
    result['view_fri'] = tmp[4::7]
    result['view_fri'][:] = result['view_fri'] * parameter_values['betaD_fri']['betaD_fri_' + country]
    for j in range(7):
        vals = result['values'][:,j]
        result['values'][:,j] = vals/vals.sum()

    return result


def create_candidate_Xi_R(Xi_R):
    # expects a country specific Xi_R dict with keys: values, view_sat, view_sun, ...
    result = {}
    result['values'] = copy.deepcopy(Xi_R['values'])  # np array mit ~number_day x 7
    tmp = result['values'].reshape(-1, order='F')  # only used to define the view. numpy just can create view on slices...trick only works when we get nb % 1
    result['view_sat'] = tmp[5::7]
    result['view_sun'] = tmp[6::7]
    result['view_mon'] = tmp[0::7]
    result['view_tue'] = tmp[1::7]
    result['view_wed'] = tmp[2::7]
    result['view_fri'] = tmp[4::7]
    return result


def create_candidate_Xi_D(Xi_D):
    # expects a country specific Xi_D dict with keys: values, view_sat, view_sun,...
    result = {}
    result['values'] = copy.deepcopy(Xi_D['values'])
    tmp = result['values'].reshape(-1, order='F')
    result['view_sat'] = tmp[5::7]
    result['view_sun'] = tmp[6::7]
    result['view_mon'] = tmp[0::7]
    result['view_tue'] = tmp[1::7]
    result['view_wed'] = tmp[2::7]
    result['view_fri'] = tmp[4::7]
    return result


def adapt_Xi_R(Xi_R, theta_curr, theta_cand, weekday):
    # here Xi_R is a dict with a np.array(~days x 7), and 6 views (i.e. teh specific values on sat,sun,etc)
    if weekday == 'sat':
        Xi_R['view_sat'][:] = Xi_R['view_sat'] * theta_cand / theta_curr
    elif weekday == 'sun':
        Xi_R['view_sun'][:] = Xi_R['view_sun'] * theta_cand / theta_curr
    elif weekday == 'mon':
        Xi_R['view_mon'][:] = Xi_R['view_mon'] * theta_cand / theta_curr
    elif weekday == 'tue':
        Xi_R['view_tue'][:] = Xi_R['view_tue'] * theta_cand / theta_curr
    elif weekday == 'wed':
        Xi_R['view_wed'][:] = Xi_R['view_wed'] * theta_cand / theta_curr
    elif weekday == 'fri':
        Xi_R['view_fri'][:] = Xi_R['view_fri'] * theta_cand / theta_curr

    for j in range(7):
        vals = Xi_R['values'][:,j]
        Xi_R['values'][:,j] = vals / vals.sum()

    return Xi_R


def adapt_Xi_D(Xi_D, theta_curr, theta_cand, weekday):
    # here Xi_D is a dict with a np.array(~days x 7), and 6 views (i.e. teh specific values on sat, sun, ...)
    if weekday == 'sat':
        Xi_D['view_sat'][:] = Xi_D['view_sat']*theta_cand/theta_curr
    elif weekday == 'sun':
        Xi_D['view_sun'][:] = Xi_D['view_sun']*theta_cand/theta_curr
    elif weekday == 'mon':
        Xi_D['view_mon'][:] = Xi_D['view_mon']*theta_cand/theta_curr
    elif weekday == 'tue':
        Xi_D['view_tue'][:] = Xi_D['view_tue']*theta_cand/theta_curr
    elif weekday == 'wed':
        Xi_D['view_wed'][:] = Xi_D['view_wed']*theta_cand/theta_curr
    elif weekday == 'fri':
        Xi_D['view_fri'][:] = Xi_D['view_fri']*theta_cand/theta_curr

    for j in range(7):
        vals = Xi_D['values'][:,j]
        Xi_D['values'][:,j] = vals/vals.sum()

    return Xi_D


def calculate_correction_factor1(infections, parameter_values):
    correction_factors1 = {country_key: calculate_correction_factor1_country(
        infections[country_key], parameter_values['N_' + country_key], parameter_values['probability_reinfection']) for country_key in infections}
    return correction_factors1


# faster version - see above for basic idea
@jit(nopython=True)
def calculate_correction_factor1_country(infections_country, N, probability_reinfection):
    cf_unshifted = np.cumsum(infections_country) / N * (1 - probability_reinfection)
    # shift by one day because the sum runs only until the previous day
    correction_factor1 = np.concatenate((np.array([0]), cf_unshifted[:-1]))
    return correction_factor1


# correction factor for vaccinated 1
def calculate_correction_factor2(data_dict, parameter_values):
    correction_factor2 = {country_key: calculate_correction_factor2_country(
        data_dict[country_key], parameter_values) for country_key in data_dict.keys()}
    return correction_factor2


# correction factor for vaccinated 2
def calculate_correction_factor2_country(data_country, parameter_values):
    cf1 = parameter_values['correction_first_vaccination']
    cf2 = parameter_values['correction_second_vaccination']
    correction_factor2 = data_country['first_vaccination'] * cf1 + data_country['second_vaccination'] * cf2
    return correction_factor2.values


# only required for initialization of Rt
def calculate_Rt(parameter_values, data_dict, lv=None):
    Rt = {country_key: calculate_Rt_country(parameter_values,
          data_dict[country_key], country_key, lv=lv) for country_key in data_dict.keys()}
    return(Rt)


# only required for initialization of Rt
def calculate_Rt_country(parameter_values, data, country, lv=None):
    R0 = parameter_values['R0']['R0_' + country]
    if set(['alpha', 'delta']).issubset(data.columns):
        R0_1 = R0
        R0_alpha = R0 * (1 + parameter_values['beta_alpha'])
        R0_delta = R0 * (1 + parameter_values['beta_delta'])
        prevalence_alpha = data.alpha.values
        prevalence_delta = data.delta.values
        R0 = R0_1 * (1 - prevalence_alpha - prevalence_delta) + R0_alpha * prevalence_alpha + R0_delta * prevalence_delta

    intervention_names = [alpha_key[6:] for alpha_key in parameter_values['alpha'].keys()]
    interventions = data[intervention_names]
    alpha = parameter_values['alpha']

    hierarchical_interventions = isinstance(alpha['alpha_' + intervention_names[0]], dict)

    Rt = R0
    for intervention in intervention_names:
        if hierarchical_interventions:
            # condition necessary because some countries dont implement all interventions
            if 'alpha_' + intervention + '_' + country in alpha['alpha_' + intervention]:
                Rt *= np.exp(interventions[intervention] * (-alpha['alpha_' + intervention]['alpha_' + intervention + '_' + country]))
        else:
            Rt *= np.exp(interventions[intervention] * (-alpha['alpha_' + intervention]))
    Rt = Rt.values

    return(Rt)


def adapt_Rt(current_Rt, current_values, candidate_values, parent_name,
             parameter_name, data, lv=None):

    if parent_name == 'R0':
        Rt = copy.deepcopy(current_Rt)
        country = parameter_name[3:]
        R0_old = current_values['R0'][parameter_name]
        R0_new = candidate_values['R0'][parameter_name]

        Rt[country] = R0_new * current_Rt[country] / R0_old

    elif parent_name == 'alpha':
        name_split = parameter_name.split('_')
        intervention = name_split[1]
        hierarchical_interventions = len(name_split) == 3  # check whether alpha is hierarchical
        Rt = copy.deepcopy(current_Rt)

        if hierarchical_interventions:
            country = name_split[2]
            factor_old = np.exp(-current_values['alpha']['alpha_' + intervention][parameter_name])
            factor_new = np.exp(-candidate_values['alpha']['alpha_' + intervention][parameter_name])
            line_selection = data[country][intervention] == 1
            Rt[country][line_selection] = factor_new*current_Rt[country][line_selection]/factor_old
        else:
            factor_old = np.exp(-current_values['alpha'][parameter_name])
            factor_new = np.exp(-candidate_values['alpha'][parameter_name])
            for country_iter in data.keys():
                line_selection = data[country_iter][intervention] == 1
                Rt[country_iter][line_selection] = factor_new*current_Rt[country_iter][line_selection]/factor_old
    elif parent_name in ['beta_alpha', 'beta_delta']:
        Rt = copy.deepcopy(current_Rt)
        for country_iter in data.keys():
            prevalence_alpha = data[country_iter].alpha.values
            prevalence_delta = data[country_iter].delta.values
            factor_old = (1-prevalence_alpha-prevalence_delta) + (1+current_values['beta_alpha'])*prevalence_alpha + (1+current_values['beta_delta'])*prevalence_delta
            factor_new = (1-prevalence_alpha-prevalence_delta) + (1+candidate_values['beta_alpha'])*prevalence_alpha + (1+candidate_values['beta_delta'])*prevalence_delta
            Rt[country_iter] = factor_new*current_Rt[country_iter]/factor_old
    else: 
        Rt = current_Rt
    return(Rt)


def calculate_ECt(sumut, Rt, data, country, parameter_values, start, correction_factor1=None, correction_factor2=None):
    if correction_factor1 is None:
        correction_factor1 = {country_iter: 0 for country_iter in data.country.unique()}
    if correction_factor2 is None:
        correction_factor2 = {country_iter: 0 for country_iter in data.country.unique()}

    if country == 'all':
        EC_t = np.zeros(data.shape[0])
        countries = data['country'].unique()
        for country_iter in countries:
            sumut_country = sumut[country_iter]
            cf_full = 1 - correction_factor1[country_iter] - correction_factor2[country_iter] * (1 - correction_factor1[country_iter])
            EC_t_country = sumut_country * Rt[country_iter] * cf_full
            EC_t_country[0:(start - 1)] = parameter_values['tau']['tau_' + country_iter]
            EC_t[data['country'].values == country_iter] = EC_t_country
    else:
        sumut_country = sumut[country]
        cf_full = 1 - correction_factor1[country] - correction_factor2[country] * (1 - correction_factor1[country])
        EC_t = sumut_country * Rt[country] * cf_full
        EC_t[0:(start - 1)] = parameter_values['tau']['tau_' + country]

    return(EC_t)


def calculate_sumut(cases_view, parameter_values, start):
    # assumes country specific values
    sumut = {}
    for country in cases_view:
        sumut_country = np.zeros(cases_view[country].shape)
        sumut_country[0:(start - 1)] = None
        sumut_country = calc_sumut_fast(cases=cases_view[country],
                                        sumut=sumut_country,
                                        gamma=parameter_values['gamma'],
                                        start=start)
        sumut[country] = sumut_country
    return(sumut)


@jit(nopython=True)
def calc_sumut_fast(cases, sumut, gamma, start):
    for t in range(start - 1, len(cases)):
        element_sumut = 0
        for u in range(t):
            element_sumut += cases[u] * gamma[t - u]
        sumut[t] = element_sumut

    return(sumut)




# Version using convolve - faster but less exact
def calculate_sumut_alternative(cases_view, parameter_values, start):
    sumut = {}
    for country in cases_view:
        sumut_country = calc_sumut_faster(cases = cases_view[country], 
                gamma = parameter_values['gamma'], start = start, lim = cases_view[country].size)
        sumut_country_start = np.repeat(None, start-1)
        sumut_country = np.concatenate((sumut_country_start, sumut_country))
        sumut[country] = sumut_country

    return(sumut)


# dont jit this! it is slower when jitted
def calc_sumut_faster(cases, gamma, start, lim):
    """
    Is 25% faster but not as exact.
    Most likely reason is the transformation to spectral domain...
    """
    return np.convolve(cases[:lim], gamma)[(start-1):lim]



def calculate_EDt(cases, parameter_values, data, country=None, Xi_D=None):
    if 'improved_treatment' in data.columns:
        improved_treatment = data.improved_treatment.values
        pi_D = parameter_values['pi_D'][country] * np.exp(-parameter_values['beta_D']*improved_treatment)
    else:
        pi_D = parameter_values['pi_D'][country]
    if Xi_D is None:
        Xi_D = parameter_values['xi_D']
        ED_t = np.zeros(cases.shape)
        ED_t = calculate_EDt_fast(ED_t, cases = cases, xi_D = Xi_D, pi_D = pi_D)
    else:
        weekday = data.weekday.values
        ED_t = np.zeros(cases.shape)
        ED_t = calculate_EDt_fast_weekdays(ED_t, cases = cases, xi_D = Xi_D, pi_D = pi_D, weekday=weekday)
    return(ED_t)


@jit(nopython=True)
def calculate_EDt_fast(ED_t, cases, xi_D, pi_D):
    for t in range(len(cases)):
        sum_u_t = 0
        for u in range(t+1):
            sum_u_t += cases[u] * xi_D[t-u] 
        ED_t[t] = sum_u_t
    ED_t *= pi_D
    return(ED_t)


@jit(nopython=True)
def calculate_EDt_fast_weekdays(ED_t, cases, xi_D, pi_D, weekday):
    for t in range(len(cases)):
        sum_u_t = 0
        for u in range(t+1):
            j = weekday[u] # python starts counting at 0, weekdays start at 1
            sum_u_t += cases[u] * xi_D[t-u,j]
        ED_t[t] = sum_u_t
    ED_t *= pi_D
    return(ED_t)


def calculate_E_Crt(cases, Xi_R, data, rho):
    weekday = data.weekday.values
    E_Crt = np.zeros(cases.shape)
    E_Crt = calculate_E_Crt_fast(E_Crt, cases = cases, xi_R = Xi_R, weekday = weekday) 
    E_Crt *= rho
    return(E_Crt)


@jit(nopython=True)
def calculate_E_Crt_fast(E_Crt, cases, xi_R, weekday):
    for t in range(len(cases)):
        sum_u_t = 0
        for u in range(t+1):
            j = weekday[u] # python starts counting at 0, weekdays start at 1
            sum_u_t += cases[u] * xi_R[t-u,j]
        E_Crt[t] = sum_u_t
    return(E_Crt)


# alternative using built in convolution
def calculate_EDt_alternative(cases, parameter_values,data):
    if 'improved_treatment' in data.columns:
        improved_treatment = data.improved_treatment.values
        pi_D = parameter_values['pi_D'] * np.exp(-parameter_values['beta_D']*improved_treatment)
    else:
        pi_D = parameter_values['pi_D']

    lim = cases.size # assumes 1-dim array!
    ED_t = calc_EDt_faster(cases = cases, xi_D = parameter_values['xi_D'],
            pi_D = pi_D, lim=lim)
    return(ED_t)


# dont jit this! it is slower when jitted
def calc_EDt_faster(cases, xi_D, pi_D, lim):
    """
    Is 25% faster but not as exact.
    Most likely reason is the transformation to spectral domain...
    """
    return np.convolve(cases[:lim], xi_D)[:lim]*pi_D


# @jit(nopython=True)
# def calculate_EDt_faster(ED_t, cases, xi_D, pi_D):
    # for t in range(len(cases)):
        # sum_u_t = 0
        # for u in range(max(t-188, 0) ,t+1):
            # sum_u_t += cases[u] * xi_D[t-u]
        # ED_t[t] = sum_u_t
    # ED_t *= pi_D
    # return(ED_t)


def calculate_cases(infections, parameter_values): # infections country-wise
    xi_C = parameter_values['xi_C']
    cases = calculate_cases_fast(infections, xi_C)
    return(cases)



@jit(nopython=True)
def calculate_cases_fast(infections, xi_C): # infections country-wise
    cases = np.zeros(infections.shape)
    rest = 0
    for t in range(len(cases)):
        sum_u_t = rest 
        for u in range(t+1):
            sum_u_t += infections[u] * xi_C[t-u]
        rounded_sumut = np.round(sum_u_t)
        rest = sum_u_t - rounded_sumut
        cases[t] = rounded_sumut 
    return(cases)


def calculate_EHt(cases, parameter_values, country):
    EH_t = np.zeros(cases.shape)
    EH_t = calculate_EHt_fast(EH_t, cases = cases, xi_H = parameter_values['xi_H'],
            pi_H = parameter_values['piH']['piH_' + country]) * parameter_values['correction_hospitalization_'+country]
    return(EH_t)


def calculate_EHicut(cases, parameter_values, country):
    EH_t = np.zeros(cases.shape)
    EH_t = calculate_EHt_fast(EH_t, cases = cases, xi_H = parameter_values['xi_Hicu'],
            pi_H = parameter_values['piHicu']['piHicu_' + country]) * parameter_values['correction_hospitalization_'+country]
    return(EH_t)


@jit(nopython=True)
def calculate_EHt_fast(EH_t, cases, xi_H, pi_H):
    for t in range(len(cases)):
        sum_u_t = 0
        for u in range(t+1):
            sum_u_t += cases[u] * xi_H[t-u]
        EH_t[t] = sum_u_t*pi_H
    return(EH_t)


@jit(nopython=True)
def calculate_pmf_uniform_trunc(lower, upper):
    lower_new = np.maximum(0,lower)
    upper_new = np.maximum(2, upper)
    pm = 1/(upper - lower_new)
    return(pm)


# slower
def convert_params_nbinom(mu, phi):
    '''
    Converts mu and overdispersion phi to standard parameters of Negative
    Binomial distribution
    '''
    p = phi / (mu + phi)
    return (phi, p)


def nbinom_alt_pmf(counts, mu, phi):
    '''
    Alternative version of Negative Binomial distribution with mean and
    overdispersion. 
    '''
    return stats.nbinom.pmf(counts, *convert_params_nbinom(mu, phi))


def nbinom_alt_log_pmf(counts, mu, phi):
    '''
    Alternative version of Negative Binomial distribution with mean and
    overdispersion. 
    '''
    return stats.nbinom.logpmf(counts, *convert_params_nbinom(mu, phi))


def nbinom_alt_rvs(mu, phi):
    '''
    Alternative version of Negative Binomial distribution with mean and
    overdispersion. 
    '''
    return stats.nbinom.rvs(*convert_params_nbinom(mu, phi))



def asymmetric_laplace_log_pdf(x, scale, kappa):
    # this condition is not necessary since the proposals will always be poitive - here for completeness
    # if x < 0:
        # return 0
    return np.log(scale/ (kappa + (kappa** -1))) + (
            -x* scale* np.sign(x) * (kappa** np.sign(x))
           )


def halfT_pdf(x, nu, scale):
    # since wikipedia here should stand an additional factor of 2
    # here ther eis no 2, however, comparing the density, it suggests that
    scale=scale**2
    return (gamma((nu+1)/2)/(gamma(nu/2)*np.sqrt(nu*np.pi*scale))*
    (1+(x**2)/(nu*scale))**(-0.5*(nu+1)))
