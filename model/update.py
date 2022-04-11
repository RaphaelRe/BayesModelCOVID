import pandas as pd
import numpy as np
import scipy.stats as stats
import basics
import copy
import math # only required to test for nan-value
#import warnings
#import pdb
#import sys
#import traceback
import matplotlib.pyplot as plt

#np.seterr(all='warn')
#warnings.filterwarnings('error')



def parameter_update(parameter, current_values, latent_variable, parent_name):
    theta_t = parameter.get_current_value()
    # theta_t = current_values[parameter.name]
    theta_cand = parameter.propose_value()
    # if parameter.name == "alpha_NPI1_countryA" and theta_cand > theta_t:
        # print(f'proposed value for {parameter}: {theta_cand}')
    ratio = parameter.prior_ratio(theta_cand)

    # all values are the same except one
    candidate_values = copy.deepcopy(current_values)

    if isinstance(parameter.name, str):  # if condidtion is required due to int names of rho
        splitted_name = parameter.name.split('_')  # required for hierarchical alpha test

    if parent_name == parameter.name:
        candidate_values[parameter.name] = theta_cand
    elif (isinstance(parameter.name, np.integer)):
        candidate_values['rho']['rho_' + parameter.country][parameter.name] = theta_cand
    elif ((len(splitted_name) == 3) and (splitted_name[0] != 'beta' and splitted_name[0] != 'betaD')):  # check if alpha is in the last hierarchy
        candidate_values[parent_name][parent_name + '_' + splitted_name[1]][parameter.name] = theta_cand
    else:
        candidate_values[parent_name][parameter.name] = theta_cand

    latent_variable_curr = latent_variable.get_values()
    latent_variable_cand = {'cases': latent_variable_curr['cases'],
                            'sumut': latent_variable_curr['sumut'],
                            'correction_factor1': latent_variable_curr['correction_factor1'],
                            'correction_factor2': latent_variable_curr['correction_factor2'],
                            }
    if (isinstance(parameter.name, np.integer)):  # necessary sinc we need Xi_R in rho
        latent_variable_cand['Xi_R'] = latent_variable_curr['Xi_R']


    if parent_name in ['beta_sat', 'beta_sun', 'beta_mon', 'beta_tue', 'beta_wed', 'beta_fri', 'betaD_sat', 'betaD_sun', 'betaD_mon', 'betaD_tue', 'betaD_wed', 'betaD_fri']:
        country = splitted_name[-1]  # # # check this! sould return country
        weekday = splitted_name[-2]  # # # check this!should return the weekday
        if splitted_name[0] == 'beta':
            latent_variable_cand['Xi_R'] = {}
            candidate_Xi = basics.create_candidate_Xi_R(latent_variable_curr['Xi_R'][country])
            latent_variable_cand['Xi_R'][country] = basics.adapt_Xi_R(candidate_Xi, theta_t, theta_cand, weekday)
        elif splitted_name[0] == 'betaD':
            latent_variable_cand['Xi_D'] = {}
            candidate_Xi = basics.create_candidate_Xi_D(latent_variable_curr['Xi_D'][country])
            latent_variable_cand['Xi_D'][country] = basics.adapt_Xi_D(candidate_Xi, theta_t, theta_cand, weekday)

    if latent_variable.model == 'infections':
        latent_variable_cand['infections'] = latent_variable_curr['infections']
        latent_variable_cand['cases_view'] = latent_variable_curr['cases_view']
        latent_variable_cand.update(latent_variable.infections_view) # updates the countries i.e. all keys where a country stands is a case view
    elif latent_variable.model == 'cases':
        latent_variable_cand.update(latent_variable.cases_view) # updates the countries i.e. all keys where a country stands is a case view

    latent_variable_cand['Rt'] = basics.adapt_Rt(latent_variable_curr['Rt'],
            current_values, candidate_values, parent_name,
            parameter.name, latent_variable.data)

    ratio *= parameter.likelihood_ratio(current_values, candidate_values,
            latent_variable_curr, latent_variable_cand)

    if parameter.calculate_proposal_ratio:
        ratio *= parameter.proposal_ratio(theta_cand, parent_name)
    accept = (np.random.uniform(0, 1) < ratio)
    if accept:
        # theta_new = theta_t + (theta_cand-theta_t)*accept # schnelle version von ifelse
        theta_t = theta_cand
        # if parent_name in ['R0', 'alpha', 'beta_voc']:
        if parent_name in ['R0', 'alpha', 'beta_alpha', 'beta_delta']:
            latent_variable.update_Rt(latent_variable_cand['Rt'])
        elif parent_name in ['beta_sat', 'beta_sun', 'beta_mon', 'beta_tue', 'beta_wed', 'beta_fri']:
            latent_variable.update_Xi_R(latent_variable_cand['Xi_R'], country)
        elif parent_name in ['betaD_sat', 'betaD_sun', 'betaD_mon', 'betaD_tue', 'betaD_wed', 'betaD_fri']:
            latent_variable.update_Xi_D(latent_variable_cand['Xi_D'], country)


    parameter.samples[parameter.i + 1] = theta_t
    parameter.acceptance[parameter.i] = accept
    parameter.i += 1




def prior_update(prior_parameter, current_values, current_prior_values, name_parent):
    values = prior_parameter.extract_values(current_values)
    theta_t = prior_parameter.get_current_value()
    theta_cand = prior_parameter.propose_value()
    ratio = prior_parameter.prior_ratio(theta_cand)


    # all values are the same except one
    candidate_prior_values = copy.deepcopy(current_prior_values)

    candidate_prior_values[prior_parameter.name] = theta_cand
    ratio *= prior_parameter.likelihood_ratio(values, current_prior_values, candidate_prior_values)

    if prior_parameter.calculate_proposal_ratio:
        ratio *= prior_parameter.proposal_ratio(theta_cand, name_parent)

    accept = (np.random.uniform(0, 1) < ratio)
    if accept:
        theta_t = theta_cand

    prior_parameter.samples[prior_parameter.i + 1] = theta_t
    prior_parameter.acceptance[prior_parameter.i] = accept
    prior_parameter.i += 1



# def ratio_likelihood_normal(values, current_prior_values, candidate_prior_values):
    # ratio = np.prod(
                # stats.norm.pdf(values, loc = candidate_prior_values['mean'], scale = candidate_prior_values['sd'])/
                # stats.norm.pdf(values, loc = current_prior_values['mean'], scale = current_prior_values['sd'])
            # )
    # return(ratio)


def ratio_likelihood_normal(values, current_prior_values, candidate_prior_values):
    mu_curr = current_prior_values['mean']
    sd_curr = current_prior_values['sd']
    mu_cand = candidate_prior_values['mean']
    sd_cand = candidate_prior_values['sd']

    return np.exp((values.shape[0] * (np.log(sd_curr) - np.log(sd_cand)) +
                  1 / (2 * sd_curr**2) * np.sum((values - mu_curr)**2) -
                  1 / (2 * sd_cand**2) * np.sum((values - mu_cand)**2))
                  )

# test
# for _ in range(1000):
    # x = np.random.random(1)
    # current_prior_values = {'mean':np.random.random(1), 'sd':np.random.random(1)}
    # candidate_prior_values = {'mean':np.random.random(1), 'sd':np.random.random(1)}
    # print("__________________")
    # print(ratio_likelihood_normal(x, current_prior_values, candidate_prior_values)-
            # ratio_likelihood_normal2(x, current_prior_values, candidate_prior_values)
            # )




def prior_ratio_exponential(theta_t, theta_cand, prior_parameters_theta, informative_priors):
    ratio = stats.expon.pdf(theta_cand,loc=0,scale=1/prior_parameters_theta['lambda'])/stats.expon.pdf(
                            theta_t,loc=0, scale=1/prior_parameters_theta['lambda'])
    return(ratio)

def prior_ratio_gamma(theta_t, theta_cand, prior_parameters_theta, informative_priors):

    if informative_priors:
        ratio =  stats.gamma.pdf(theta_cand, a = prior_parameters_theta['shape'],
                                 scale = prior_parameters_theta['scale'])/stats.gamma.pdf(theta_t,
                                 a = prior_parameters_theta['shape'],
                                 scale = prior_parameters_theta['scale'])
    else:
        prior_min  = prior_parameters_theta['min']
        prior_max  = prior_parameters_theta['max']
        ratio =  stats.uniform.pdf(theta_cand, loc= prior_min,
                                   scale = prior_max-prior_min)/stats.uniform.pdf(theta_t,
                                   loc = prior_min, scale = prior_max-prior_min)
    return(ratio)


# old version
# def prior_ratio_uniform(theta_t, theta_cand, prior_parameters_theta, informative_priors = False):
    # prior_min  = prior_parameters_theta['min']
    # prior_max  = prior_parameters_theta['max']
    # ratio =  (stats.uniform.pdf(theta_cand, loc= prior_min,
                                   # scale = prior_max-prior_min) / 
             # stats.uniform.pdf(theta_t, loc = prior_min, 
                                   # scale = prior_max-prior_min)
             # )
    # return(ratio)


# faster version see explanation in the function
def prior_ratio_uniform(theta_t, theta_cand, prior_parameters_theta, informative_priors = False):
    # prior ratio is 1 as long as theta_cand lies inside the min and max. 
    # otherwise it is 0. Theta_t does not affect the ratio (see desity of 
    # Uniform distribution )
    return(1.0 if ((prior_parameters_theta['min'] < theta_cand) & (theta_cand < prior_parameters_theta['max'])) else 0)

# pp = {'min': 0, 'max': 2}
# val1 = 1
# val2 = 10

# prior_ratio_uniform(1,1,pp,"a")
# prior_ratio_uniform_fast(1,1,pp,"a")


# for i in range(100):
    # m = np.random.uniform()
    # pp = {'min': m, 'max': m*10}
    # val = np.random.uniform()
    # val2 = np.random.uniform()
    # print(f"m: {m}, val1: {val}, val2: {val2}")
    # print(prior_ratio_uniform(val,val2,pp,"a"))
    # print(prior_ratio_uniform_fast(val,val2,pp,"a"))
    # print("")





def prior_ratio_beta(theta_t, theta_cand, prior_parameters_theta, informative_priors):
    ratio = stats.beta.pdf(theta_cand, a = prior_parameters_theta['a'],
                                       b = prior_parameters_theta['b'],
                                       scale = prior_parameters_theta['scale']
#                                       loc = prior_parameters_theta['min'],
#                                       scale = prior_parameters_theta['max'] -
#                                       prior_parameters_theta['min'])
                                     ) /stats.beta.pdf(theta_t,
                                       a = prior_parameters_theta['a'],
                                       b = prior_parameters_theta['b'],
                                       scale = prior_parameters_theta['scale'])
#                                       loc = prior_parameters_theta['min'],
#                                       scale = prior_parameters_theta['max'] - 
#                                       prior_parameters_theta['min'])
    return(ratio)
    
    
def prior_ratio_inv_gamma(theta_t, theta_cand, prior_parameters_theta, informative_priors):
    
    if informative_priors:
        ratio =  stats.invgamma.pdf(theta_cand, a = prior_parameters_theta['shape'],
                                 scale = prior_parameters_theta['scale'])/stats.invgamma.pdf(theta_t,
                                 a = prior_parameters_theta['shape'], 
                                 scale = prior_parameters_theta['scale']) 
    else:
        prior_min  = prior_parameters_theta['min']
        prior_max  = prior_parameters_theta['max']
        ratio =  stats.uniform.pdf(theta_cand, loc= prior_min,
                                   scale = prior_max-prior_min)/stats.uniform.pdf(theta_t,
                                   loc = prior_min, scale = prior_max-prior_min)
    return(ratio)



def prior_ratio_lognormal(theta_t, theta_cand, prior_parameters_theta):
    #ratio = ((1/theta_cand)*np.exp(-((np.log(theta_cand)-prior_parameters_theta['mean'])**2)/(2*(prior_parameters_theta['sd']**2)))/
    #(1/theta_t)*np.exp(-((np.log(theta_t)-prior_parameters_theta['mean'])**2)/(2*(prior_parameters_theta['sd']**2))))
    ratio = (stats.lognorm.pdf(theta_cand, s = prior_parameters_theta['sd'], scale = np.exp(prior_parameters_theta['mean']))/
            stats.lognorm.pdf(theta_t, s= prior_parameters_theta['sd'],scale = np.exp(prior_parameters_theta['mean'])))
    return(ratio)
    
# def prior_ratio_normal(theta_t, theta_cand, prior_parameters_theta, informative_priors):
    # if informative_priors:
        # ratio = stats.norm.pdf(theta_cand,loc=prior_parameters_theta['mean'],
                            # scale=prior_parameters_theta['sd'])/stats.norm.pdf(theta_t,
                            # loc=prior_parameters_theta['mean'], 
                            # scale=prior_parameters_theta['sd'])
    # else:
        # prior_min  = prior_parameters_theta['min']
        # prior_max  = prior_parameters_theta['max']
        # ratio =  stats.uniform.pdf(theta_cand, loc= prior_min,
                                   # scale = prior_max-prior_min)/stats.uniform.pdf(theta_t,
                                   # loc = prior_min, scale = prior_max-prior_min)
    # return(ratio)
    

def prior_ratio_normal(theta_t, theta_cand, prior_parameters_theta, informative_priors):
    mu = prior_parameters_theta['mean']
    var = prior_parameters_theta['sd']**2
    return np.exp( 1/(2*var)*( (theta_t-mu)**2 - (theta_cand-mu)**2 ))


# for _ in range(1000):
    # p = {'mean': np.random.random(1), 'sd': np.random.random(1)*5}
    # x1 = prior_ratio_normal(1,2, p,1)
    # x2 = prior_ratio_normal2(1,2, p,1)
    # print(np.isclose(x1/x2,1))




def prior_ratio_normal_inverse(theta_t, theta_cand, prior_parameters_theta, informative_priors):
    if informative_priors:
        inv_theta_t = 1 / np.sqrt(theta_t)
        inv_theta_cand = 1 / np.sqrt(theta_cand)

        ratio = stats.norm.pdf(inv_theta_cand, loc=prior_parameters_theta['inv_mean'],  # muss man die dann auf 0 und 1 stzen?
                               scale=prior_parameters_theta['inv_sd']) / stats.norm.pdf(inv_theta_t,
                               loc=prior_parameters_theta['inv_mean'],
                               scale=prior_parameters_theta['inv_sd'])
    else:
        prior_min  = prior_parameters_theta['min']
        prior_max  = prior_parameters_theta['max']
        ratio =  stats.uniform.pdf(theta_cand, loc= prior_min,
                                   scale = prior_max-prior_min)/stats.uniform.pdf(theta_t,
                                   loc = prior_min, scale = prior_max-prior_min)
    return(ratio)



def prior_ratio_cauchy(theta_t, theta_cand, prior_parameters_theta, informative_priors):
    ratio = (stats.cauchy.pdf(theta_cand, loc = prior_parameters_theta['loc'], scale = prior_parameters_theta['scale']) /
            stats.cauchy.pdf(theta_t, loc = prior_parameters_theta['loc'], scale = prior_parameters_theta['scale']))
    return ratio




def prior_ratio_halfT(theta_t, theta_cand, prior_parameters_theta, informative_priors):
    ratio = (basics.halfT_pdf(theta_cand, prior_parameters_theta['nu'], prior_parameters_theta['scale']) /
             basics.halfT_pdf(theta_t, prior_parameters_theta['nu'], prior_parameters_theta['scale'])
             )
    return ratio


def prior_ratio_asymmetric_laplace(theta_t, theta_cand, prior_parameters_theta, informative_priors):
    ratio = np.exp(
            basics.asymmetric_laplace_log_pdf(
                theta_cand, prior_parameters_theta['scale'], prior_parameters_theta['kappa']
            ) -
            basics.asymmetric_laplace_log_pdf(
                theta_t, prior_parameters_theta['scale'], prior_parameters_theta['kappa']
            ))
    return ratio

# prior_ratio_asymmetric_laplace(1,0.2, {'scale':10., 'kappa': 0.5}, "aaa")


# x = np.linspace(-1,2, 1000)
# plt.plot(x,prior_ratio_asymmetric_laplace(0.05,x, {'scale':10., 'kappa': 0.5}, "aaa"))
# plt.show()

# proposals

# def proposal_normal(theta_t, proposal_sd):
    # theta_cand = theta_t + stats.norm.rvs(loc=0, scale=proposal_sd)
    # return(theta_cand)

def proposal_normal(theta_t, proposal_sd):
    return theta_t + np.random.normal(0, proposal_sd, 1)[0]




def proposal_trunc(theta_t, proposal_sd):
    a_trunc = (0 - theta_t) / proposal_sd
    theta_cand = stats.truncnorm.rvs(a=a_trunc, b=np.inf, loc=theta_t,
                                     scale=proposal_sd)
    return(theta_cand)


def proposal_double_trunc(theta_t, proposal_sd):
    # a and b in truncnorm muessen als relative werte zum mittelwert und sd
    # angegeben weren!
    a_trunc = (0 - theta_t) / proposal_sd
    b_trunc = (1 - theta_t) / proposal_sd
    theta_cand = stats.truncnorm.rvs(a=a_trunc, b=b_trunc, loc=theta_t,
                                     scale=proposal_sd)
    return(theta_cand)



def proposal_double_trunc_rho(theta_t, proposal_sd):
    # a and b in truncnorm muessen als relative werte zum mittelwert und sd
    # angegeben weren!
    a_trunc = (0 - theta_t) / proposal_sd
    b_trunc = (10 - theta_t) / proposal_sd
    theta_cand = stats.truncnorm.rvs(a=a_trunc, b=b_trunc, loc=theta_t,
                                     scale=proposal_sd)
    return(theta_cand)



priors = {'alpha': prior_ratio_normal,
          'alpha_mean': prior_ratio_normal,
          # 'alpha_mean': prior_ratio_asymmetric_laplace,
          'alpha_sd': prior_ratio_normal,
          # 'alpha_sd': prior_ratio_halfT,
          'tau': prior_ratio_normal,
          'tau_mean': prior_ratio_gamma,
          'tau_sd': prior_ratio_normal,
          'rho': prior_ratio_beta,
          'R0': prior_ratio_normal,
          'R0_mean': prior_ratio_normal,
          'R0_sd': prior_ratio_normal,
          'phi': prior_ratio_normal_inverse,
          'piH': prior_ratio_beta,
          'piHicu': prior_ratio_beta,
          # 'phi_infections': prior_ratio_normal,
          # 'phi_deaths': prior_ratio_normal,
          # 'phi_hospitalizations': prior_ratio_normal
          'beta_D': prior_ratio_cauchy,
          # 'beta_voc': prior_ratio_normal,
          'beta_alpha': prior_ratio_normal,
          'beta_delta': prior_ratio_normal,
          'beta_sat': prior_ratio_uniform,
          'beta_sun': prior_ratio_uniform,
          'beta_mon': prior_ratio_uniform,
          'beta_tue': prior_ratio_uniform,
          'beta_wed': prior_ratio_uniform,
          'beta_fri': prior_ratio_uniform,
          'betaD_sat': prior_ratio_uniform,
          'betaD_sun': prior_ratio_uniform,
          'betaD_mon': prior_ratio_uniform,
          'betaD_tue': prior_ratio_uniform,
          'betaD_wed': prior_ratio_uniform,
          'betaD_fri': prior_ratio_uniform,
          }

priors['phi_hospitalizations'] = priors['phi']
priors['phi_intensiveCare'] = priors['phi']
priors['phi_infections'] = priors['phi']
priors['phi_deaths'] = priors['phi']
priors['phi_rep_cases'] = priors['phi']



proposals = {'alpha': proposal_normal,
             'alpha_mean': proposal_normal,
             'alpha_sd': proposal_trunc,
             # 'rho': proposal_double_trunc,
             'rho': proposal_double_trunc_rho,
             'tau': proposal_trunc,
             'tau_mean': proposal_trunc,
             'tau_sd': proposal_trunc,
             'R0': proposal_trunc,
             'R0_mean': proposal_trunc,
             'R0_sd': proposal_trunc,
             'phi': proposal_trunc,
             'piH': proposal_double_trunc_rho,
             'piHicu': proposal_double_trunc_rho,
             # 'piH':proposal_double_trunc,
             # 'piHicu':proposal_double_trunc,
             'beta_D': proposal_normal,
             # 'beta_voc': proposal_normal,
             'beta_alpha': proposal_normal,
             'beta_delta': proposal_normal,
             # 'phi_infections': proposal_trunc,
             # 'phi_deaths': proposal_trunc,
             # 'phi_hospitalizations': proposal_trunc
             'beta_sat': proposal_trunc,
             'beta_sun': proposal_trunc,
             'beta_mon': proposal_trunc,
             'beta_tue': proposal_trunc,
             'beta_wed': proposal_trunc,
             'beta_fri': proposal_trunc,
             'betaD_sat': proposal_trunc,
             'betaD_sun': proposal_trunc,
             'betaD_mon': proposal_trunc,
             'betaD_tue': proposal_trunc,
             'betaD_wed': proposal_trunc,
             'betaD_fri': proposal_trunc,
             }

proposals['phi_hospitalizations'] = proposals['phi']
proposals['phi_intensiveCare'] = proposals['phi']
proposals['phi_infections'] = proposals['phi']
proposals['phi_deaths'] = proposals['phi']
proposals['phi_rep_cases'] = proposals['phi']


# dict for priors
prior_likelihoods = {'R0': ratio_likelihood_normal,
                     'alpha': ratio_likelihood_normal,
                     'tau': ratio_likelihood_normal,
                     # 'pi_H': ratio_likelihood_beta,
                     # 'nb_trials': ratio_likelihood_discrete_iwas # für negative binomial model
                     }



def ratio_likelihood(model, data, values_curr, values_cand, latent_variables_curr,
                     latent_variables_cand, country='all', start=7,
                     death_model=True, reporting_model=True, renewal_model_lv=True,
                     renewal_model_parameter=False, hospitalization_model=False):  # values_... ist ein dict mit den Parameterwerten: jede reintrag enthält einen Parameter (im Falle von rho ein weiteres dict mit M Werten)
    # latent_variable ist ein dict mit allen latenten Variablen: {'cases':..., 'infections':...}
    '''
        Determines the likelihood ratio for two different sets of parameters and values for the latent variables.
        '''
    print("==================================================")
    print("==================================================")
    print("Using old Likelihood")
    print("==================================================")
    print("==================================================")
    ratio = 1
    if country == 'all': # the case for tau and alpha
        C_t_cand = latent_variables_cand['cases']
        C_t_curr = latent_variables_curr['cases']

        if model == 'infections':
            I_t_cand = latent_variables_cand['infections']
            I_t_curr = latent_variables_curr['infections']
        else: # the case for rho and R0
            I_t_cand = C_t_cand 
            I_t_curr = C_t_curr

    else:
        I_t_cand = latent_variables_cand[country]
        I_t_curr = latent_variables_curr[country]

        if model == 'infections':
            C_t_cand = latent_variables_cand['cases_view'][country]
            C_t_curr = latent_variables_curr['cases_view'][country]
        else:
            C_t_cand = I_t_cand 
            C_t_curr = I_t_curr


    
    if (renewal_model_lv | renewal_model_parameter): # renewal_model_parameter is always Flase currently....just an idea for getting faster updates where !k can be dropped in calculation
        R_t_cand = latent_variables_cand['Rt']
        R_t_curr= latent_variables_curr['Rt'] 
        sumut_cand = latent_variables_cand['sumut']
        sumut_curr = latent_variables_curr['sumut']
        E_Ct_cand = basics.calculate_ECt(sumut_cand, R_t_cand, data, country,
                values_cand, start)
        E_Ct_curr = basics.calculate_ECt(sumut_curr, R_t_curr, data, country,
                values_curr, start)

        ## hier koennte eine schnellere Version geschrieben werden fuer das updaten von R0 und alpha
        ## wenn die Cts also konstant bleiben: dann koennen wir k! rauskuerzen, die untere Version funktioniert nur nicht!
        #ratio_renewal_model = np.prod(
        #                           (E_Ct_cand**C_t_cand/E_Ct_curr**C_t_curr)
        #                           )*np.exp(np.sum(E_Ct_curr)-np.sum(E_Ct_cand))
        #print(stats.poisson.pmf(C_t_cand, E_Ct_cand)/
         #                             stats.poisson.pmf(C_t_curr, E_Ct_curr))
        if renewal_model_lv:
            ratio_renewal_model = np.prod(
                                          stats.poisson.pmf(I_t_cand, E_Ct_cand)/
                                          stats.poisson.pmf(I_t_curr, E_Ct_curr)
                                          )
            if np.isnan(ratio_renewal_model):
                index = np.where(stats.poisson.pmf(I_t_curr, E_Ct_curr) == 0)
                #print(pd.DataFrame({'v1':I_t_cand[index],'v2':E_Ct_cand[index],'v3':E_Ct_curr[index]}))   

                num_zeros = np.sum(stats.poisson.pmf(I_t_curr, E_Ct_curr) == 0)
                with open("zero_dev_errors", "a") as errfile:
                    errfile.write(country + str(num_zeros) +  "\n")
                #print("Setting " + str(num_zeros) +  " ratios to 1 in " + country)
                pois_ratio1 = stats.poisson.pmf(I_t_cand, E_Ct_cand)
                pois_ratio2 = stats.poisson.pmf(I_t_curr, E_Ct_curr)
                pois_ratio = np.divide(pois_ratio1, pois_ratio2,
                                       out = np.ones_like(pois_ratio1),
                                       where = pois_ratio2!=0
                                       )

                ratio_renewal_model = np.prod(pois_ratio)
        

        elif renewal_model_parameter: # schnellere Version, die k! nicht berechnen, da dies sich im ratio rauskürzt
            print('Attention!!!!!!! Currently not working!!!')
            ratio_renewal_model = np.prod(
                                          (E_Ct_cand/E_Ct_curr)**C_t_cand
                                          )*np.exp(np.sum(E_Ct_curr)-np.sum(E_Ct_cand)
                                          )

        ratio *= ratio_renewal_model

    ## Harakiri: wenn hier ein Fehler kommt, dass man EDt fuer alle Laender braeuchte dann
    ## ist das eine falsche Anfrage weil bis jetzt das death model nur fuer das updaten von
    ## Ct benutzt wird und dieses Update ohnehin laenderspezifisch sein sollte! Falls eines Tages
    ## pi_d geschaetzt werden soll (auch anhand von externer Daten) aendert sich das vielleicht
    if death_model:
        D_t = data.deaths
        E_Dt_cand = basics.calculate_EDt(C_t_cand, values_cand, data = pd.DataFrame())
        E_Dt_curr = basics.calculate_EDt(C_t_curr, values_curr, data = pd.DataFrame())
        ratio_death_model = np.prod(
                                   (E_Dt_cand/E_Dt_curr)**D_t
                                   )*np.exp(np.sum(E_Dt_curr)-np.sum(E_Dt_cand))

        #print('ratio_death_model')
        #print(ratio_death_model)
        ratio *= ratio_death_model

    # rho_t_cand = basics.calculate_long_rho(values_cand, data) # brauchen wir im Moment nicht - siehe basics.calculate long_rho
    # rho_t_curr= basics.calculate_long_rho(values_curr, data)

    # reporting model:
    if reporting_model:
        Cr_t = data.reported_cases

        if 'rho_period' in data.columns:
            rho_t_cand = np.array([values_cand['rho']['rho_' + country][period] for period in data['rho_period'].values])
            rho_t_curr= np.array([values_curr['rho']['rho_' + country][period] for period in data['rho_period'].values])
            # Achtung: Die folgenden 2 Zeilen sind fuer den Fall dass rho fix ist und rho_period vorhanden ist - dann müssen die oberen beiden auskommentiert und die entkommentiert werden
            #rho_t_cand = values_cand['rho']['rho_' + country] # Harakiri: wenn Fehler rho == 'all' existiert nicht kommt also KeyError ist was schiefgegengen
            #rho_t_curr = values_curr['rho']['rho_' + country]
        else:
            rho_t_cand = values_cand['rho']['rho_' + country] # Harakiri: wenn Fehler rho == 'all' existiert nicht kommt also KeyError ist was schiefgegengen
            rho_t_curr = values_curr['rho']['rho_' + country]

        ratio_reporting_model = np.prod(
                                        stats.binom.pmf(Cr_t, C_t_cand, rho_t_cand)/
                                        stats.binom.pmf(Cr_t, C_t_curr, rho_t_curr)
                                       )

        if(math.isnan(ratio_reporting_model)):
            print('Attention: reporting model was nan:' + country)
       #     with open("reporting_model_nan", "a") as errfile:
       #         errfile.write(str("None") + "\n")
       # else:
       #     with open("reporting_model_nan", "a") as errfile:
       #         errfile.write(str("YAY") + "\n")




        ratio *= ratio_reporting_model
    #print('ratio_reporting_model')
    #print(ratio_reporting_model)

    ## hier koennte eine schnellere Version geschrieben werden fuer das updaten von R0 und alpha
    ## wenn die Cts also konstant bleiben: dann koennen wir k! rauskuerzen, die untere Version funktioniert nur nicht!
    #ratio_renewal_model = np.prod(
    #                           (E_Ct_cand**C_t_cand/E_Ct_curr**C_t_curr)
    #                           )*np.exp(np.sum(E_Ct_curr)-np.sum(E_Ct_cand))



    if hospitalization_model:

        H_t = data.hospitalizations
        E_Ht_cand = basics.calculate_EHt(C_t_cand, values_cand)
        E_Ht_curr = basics.calculate_EHt(C_t_curr, values_curr)
        ratio_hospitalization_model = np.prod(
                                   (E_Ht_cand/E_Ht_curr)**H_t
                                   )*np.exp(np.sum(E_Ht_curr)-np.sum(E_Ht_cand))

        

        ratio *= ratio_hospitalization_model

    return(ratio)




def ratio_likelihood_nbinom(model, data, values_curr, values_cand,
        latent_variables_curr, latent_variables_cand, country = 'all', 
        start=-99999, death_model=True, reporting_model=True, renewal_model_lv=True,
        renewal_model_parameter = False, hospitalization_model = False,
        intensiveCare_model=False):

    '''
    Determines the likelihood ratio for two different sets of parameters
    and values for the latent variables.
    '''

    ratio = 1
    if country == 'all':  # the case for tau and alpha
        C_t_cand = latent_variables_cand['cases']
        C_t_curr = latent_variables_curr['cases']

        if model == 'infections':
            I_t_cand = latent_variables_cand['infections']
            I_t_curr = latent_variables_curr['infections']
        else:  # the case for rho and R0
            I_t_cand = C_t_cand
            I_t_curr = C_t_curr

    else:
        I_t_cand = latent_variables_cand[country]
        I_t_curr = latent_variables_curr[country]

        if model == 'infections':
            C_t_cand = latent_variables_cand['cases_view'][country]
            C_t_curr = latent_variables_curr['cases_view'][country]
        else:
            C_t_cand = I_t_cand
            C_t_curr = I_t_curr

    if (renewal_model_lv):
        R_t_cand = latent_variables_cand['Rt']
        R_t_curr = latent_variables_curr['Rt']
        sumut_cand = latent_variables_cand['sumut']
        sumut_curr = latent_variables_curr['sumut']

        cf1_cand = latent_variables_cand['correction_factor1']
        cf1_curr = latent_variables_curr['correction_factor1']
        cf2 = latent_variables_curr['correction_factor2'] # cf2 ändert sich nicht

        E_Ct_cand = basics.calculate_ECt(sumut_cand, R_t_cand, data, country,
                values_cand, start, cf1_cand, cf2)
        E_Ct_curr = basics.calculate_ECt(sumut_curr, R_t_curr, data, country,
                values_curr, start, cf1_curr, cf2)

        if renewal_model_lv:
            phi_infections_curr = values_curr['phi_infections']
            phi_infections_cand = values_cand['phi_infections']

            ratio_renewal_model = np.exp(
                                    np.sum(
                                          basics.nbinom_alt_log_pmf(I_t_cand, E_Ct_cand, phi_infections_cand)-
                                          basics.nbinom_alt_log_pmf(I_t_curr, E_Ct_curr, phi_infections_curr)
                                          )
                                    )
            # print(ratio_renewal_model)
            # ratio_renewal_model = np.prod(
                                          # basics.nbinom_alt_pmf(I_t_cand, E_Ct_cand, phi_infections_cand)/
                                          # basics.nbinom_alt_pmf(I_t_curr, E_Ct_curr, phi_infections_curr)
                                          # )

            if np.isnan(ratio_renewal_model):
                denominator = (basics.nbinom_alt_pmf(I_t_curr, E_Ct_curr, phi_infections_curr))
                ii = np.where(basics.nbinom_alt_pmf(I_t_curr, E_Ct_curr, phi_infections_curr) == 0)
                I = (I_t_curr)
                # EI = (E_Ct_curr.values)
                EI = E_Ct_curr
                print("NA in renewal model")
                import pudb;pu.db

        ratio *= ratio_renewal_model


    if death_model:
        phi_deaths_curr = values_curr['phi_deaths']
        phi_deaths_cand = values_cand['phi_deaths']


        if country == 'all':
            # this part of the code can only be evaluated if ratio_likelihood
            # gets called for a parameter update (phi_deaths/hospitalizations)
            # therefore latent_variables_curr is always latent_variables_cand
            ratio_death_model = 1
            for country_iter in data.keys():
                if model == 'infections':
                    C_t_curr = latent_variables_curr['cases_view'][country_iter]
                else:
                    C_t_curr = latent_variables_curr[country_iter]

                if 'Xi_D' in latent_variables_curr.keys():
                    Xi_D = latent_variables_curr['Xi_D'][country_iter]['values']  # da country = 'all' ist lv_curr = lv_cand
                else:
                    # Xi_D = values_curr['xi_D'] geänder, evtl von Sabine checken lassen
                    Xi_D = None

                D_t = data[country_iter].deaths
                E_Dt_cand = basics.calculate_EDt(C_t_curr, values_cand, data[country_iter], country=country_iter, Xi_D=Xi_D)
                # E_Dt_cand = basics.calculate_EDt(C_t_curr, values_cand, data[country_iter]) # old
                E_Dt_cand *= data[country_iter]['pi_nd']
                E_Dt_curr = basics.calculate_EDt(C_t_curr, values_curr, data[country_iter], country=country_iter, Xi_D=Xi_D)
                # E_Dt_curr = basics.calculate_EDt(C_t_curr, values_curr,data[country_iter]) # old
                E_Dt_curr *= data[country_iter]['pi_nd']

                ratio_death_model *= np.exp(np.sum(basics.nbinom_alt_log_pmf(D_t, E_Dt_cand, phi_deaths_cand) -
                                                   basics.nbinom_alt_log_pmf(D_t, E_Dt_curr, phi_deaths_curr)
                                                   )
                                            )


        else:
            if 'Xi_D' in latent_variables_curr.keys():
                Xi_D_cand = latent_variables_cand['Xi_D'][country]['values']
                Xi_D_curr = latent_variables_curr['Xi_D'][country]['values']
            else:
                # Xi_D_cand = values_curr['xi_D'] # geänder, evtl von Sabine checken lassen
                # Xi_D_curr = Xi_D_cand # geänder, evtl von Sabine checken lassen
                Xi_D_cand = None  # geänder, evtl von Sabine checken lassen
                Xi_D_curr = None  # geänder, evtl von Sabine checken lassen

            D_t = data.deaths
            # E_Dt_cand = basics.calculate_EDt(C_t_cand, values_cand, data) # old
            # E_Dt_curr = basics.calculate_EDt(C_t_curr, values_curr, data) # old
            E_Dt_cand = basics.calculate_EDt(C_t_cand, values_cand, data, country=country, Xi_D=Xi_D_cand)
            E_Dt_curr = basics.calculate_EDt(C_t_curr, values_curr, data, country=country, Xi_D=Xi_D_curr)
    
            ratio_death_model = np.exp(
                                        np.sum(
                                          basics.nbinom_alt_log_pmf(D_t, E_Dt_cand, phi_deaths_cand)-
                                          basics.nbinom_alt_log_pmf(D_t, E_Dt_curr, phi_deaths_curr)
                                          )
                                    )
            if math.isnan(ratio_death_model):
                import pudb;pu.db


        ratio *= ratio_death_model

    if reporting_model:
        # import pudb;pu.db
        phi_rep_cases_curr = values_curr['phi_rep_cases']
        phi_rep_cases_cand = values_cand['phi_rep_cases']


        if country == 'all':
            # this part of the code can only be evaluated if ratio_likelihood
            # gets called for a parameter update (phi_deaths/hospitalizations)
            # therefore latent_variables_curr is always latent_variables_cand
            ratio_reporting_model = 1
            for country_iter in data.keys():
                Cr_t = data[country_iter].reported_cases
                if model == 'infections':
                    C_t_curr = latent_variables_curr['cases_view'][country_iter]
                else:
                    C_t_curr = latent_variables_curr[country_iter]


                if 'rho_period' in data[country_iter].columns:
                    rho_t_cand = np.array([values_cand['rho']['rho_' + country_iter][period] for period in data[country_iter]['rho_period'].values])
                    rho_t_curr= np.array([values_curr['rho']['rho_' + country_iter][period] for period in data[country_iter]['rho_period'].values])
                else:
                    rho_t_cand = values_cand['rho']['rho_' + country_iter]
                    rho_t_curr = values_curr['rho']['rho_' + country_iter]

                #### Achtung! Hier wird im Grunde viel Arbeit umsonst gemacht!
                #### rho wird ganz am schluss draufmultipliziert - der rest ist gleich!

                # ANPASSUNG FÜR XI (hier st curr = cand , also egal)
                if 'Xi_R' in latent_variables_curr.keys():
                    Xi_R = latent_variables_curr['Xi_R'][country_iter]['values'] # da country = 'all' ist lv_curr = lv_cand
                else:
                    Xi_R = values_curr['xi_R'] # falls Xi_R nicht im lv-dict, ist XI_R in den fixed parameters

                # E_Crt_cand = basics.calculate_E_Crt(C_t_curr, values_cand, data[country_iter] ,rho_t_cand)
                # E_Crt_cand *= data[country_iter]['pi_nc'] 
                # E_Crt_curr = basics.calculate_E_Crt(C_t_curr, values_curr, data[country_iter], rho_t_curr)
                # E_Crt_curr *= data[country_iter]['pi_nc'] 
                E_Crt_cand = basics.calculate_E_Crt(C_t_curr, Xi_R, data[country_iter] ,rho_t_cand)
                E_Crt_cand *= data[country_iter]['pi_nc'] 
                E_Crt_curr = basics.calculate_E_Crt(C_t_curr, Xi_R, data[country_iter], rho_t_curr)
                E_Crt_curr *= data[country_iter]['pi_nc'] 


                ratio_reporting_model_tmp = np.exp(
                                            np.sum(
                                              basics.nbinom_alt_log_pmf(Cr_t, E_Crt_cand, phi_rep_cases_cand)-
                                              basics.nbinom_alt_log_pmf(Cr_t, E_Crt_curr, phi_rep_cases_curr)
                                              )
                                            )
                if(math.isnan(ratio_reporting_model_tmp)):
                    denominator = basics.nbinom_alt_pmf(Cr_t, E_Crt_curr, phi_rep_cases_curr)
                    ii = np.where(basics.nbinom_alt_pmf(Cr_t, E_Crt_curr, phi_rep_cases_curr) == 0)
                    R = (Cr_t.values)
                    ER = (E_Crt_curr.values)
                    print("NA reporting model")
                    import pudb;pu.db

                ratio_reporting_model *= np.exp(
                                            np.sum(
                                              basics.nbinom_alt_log_pmf(Cr_t, E_Crt_cand, phi_rep_cases_cand)-
                                              basics.nbinom_alt_log_pmf(Cr_t, E_Crt_curr, phi_rep_cases_curr)
                                              )
                                            )


        else:
            Cr_t = data.reported_cases
            ########################
            ## ACHTUNG! Hier braucht man dictionarized data und eine schleife über die länder von phi_rep_cases_cand

            if 'rho_period' in data.columns:
                rho_t_cand = np.array([values_cand['rho']['rho_' + country][period] for period in data['rho_period'].values])
                rho_t_curr= np.array([values_curr['rho']['rho_' + country][period] for period in data['rho_period'].values])
            else:
                rho_t_cand = values_cand['rho']['rho_' + country]
                rho_t_curr = values_curr['rho']['rho_' + country]

            # ratio_reporting_model = np.prod(
                                            # stats.binom.pmf(Cr_t, C_t_cand, rho_t_cand)/
                                            # stats.binom.pmf(Cr_t, C_t_curr, rho_t_curr)
                                           # )
            if 'Xi_R' in latent_variables_curr.keys():
                Xi_R_cand = latent_variables_cand['Xi_R'][country]['values']
                Xi_R_curr= latent_variables_curr['Xi_R'][country]['values']
            else:
                Xi_R_cand = values_curr['xi_R']
                Xi_R_curr = Xi_R_cand
                

            # E_Crt_cand = basics.calculate_E_Crt(C_t_cand, values_cand, data ,rho_t_cand)
            # E_Crt_cand *= data['pi_nc']
            # E_Crt_curr = basics.calculate_E_Crt(C_t_curr, values_curr, data, rho_t_curr)
            # E_Crt_curr*= data['pi_nc']
            E_Crt_cand = basics.calculate_E_Crt(C_t_cand, Xi_R_cand, data ,rho_t_cand)
            E_Crt_cand *= data['pi_nc']
            E_Crt_curr = basics.calculate_E_Crt(C_t_curr, Xi_R_curr, data, rho_t_curr)
            E_Crt_curr*= data['pi_nc']

            # ratio_reporting_model = np.prod(
                                          # basics.nbinom_alt_pmf(Cr_t, E_Crt_cand, phi_rep_cases_cand)/
                                          # basics.nbinom_alt_pmf(Cr_t, E_Crt_curr, phi_rep_cases_curr)
                                          # )
            ratio_reporting_model = np.exp(
                                        np.sum(
                                          basics.nbinom_alt_log_pmf(Cr_t, E_Crt_cand, phi_rep_cases_cand)-
                                          basics.nbinom_alt_log_pmf(Cr_t, E_Crt_curr, phi_rep_cases_curr)
                                          )
                                        )
            if(math.isnan(ratio_reporting_model)):
                denominator = basics.nbinom_alt_pmf(Cr_t, E_Crt_curr, phi_rep_cases_curr)
                ii = np.where(basics.nbinom_alt_pmf(Cr_t, E_Crt_curr, phi_rep_cases_curr) == 0)
                R = (Cr_t.values)
                ER = (E_Crt_curr)
                import pudb;pu.db

                print('Attention: reporting model was nan:' + country)
        ratio *= ratio_reporting_model

    if hospitalization_model:
        phi_hospitalizations_curr = values_curr['phi_hospitalizations']
        phi_hospitalizations_cand = values_cand['phi_hospitalizations']

        if country == 'all':
            # this part of the code can only be evaluated if ratio_likelihood
            # gets called for a parameter update (phi_deaths/hospitalizations)
            # therefore latent_variables_curr is always latent_variables_cand
            ratio_hospitalization_model = 1
            for country_iter in data.keys():
                if pd.isna(data[country_iter]["hospitalizations"].iloc[-1]):
                        pass
                else:
                    if model == 'infections':
                        C_t_curr = latent_variables_curr['cases_view'][country_iter]
                    else:
                        C_t_curr = latent_variables_curr[country_iter]

                    H_t = data[country_iter].hospitalizations
                    E_Ht_cand = basics.calculate_EHt(C_t_curr, values_cand, country_iter)
                    E_Ht_curr = basics.calculate_EHt(C_t_curr, values_curr, country_iter)

                    # arr = basics.nbinom_alt_pmf(H_t, E_Ht_cand, phi_hospitalizations_cand)/basics.nbinom_alt_pmf(H_t, E_Ht_curr, phi_hospitalizations_curr)
                    # masked_arr = np.ma.masked_array(arr, mask = np.isnan(arr))

                    arr = (basics.nbinom_alt_log_pmf(H_t, E_Ht_cand, phi_hospitalizations_cand)-
                            basics.nbinom_alt_log_pmf(H_t, E_Ht_curr, phi_hospitalizations_curr)
                            )
                    masked_arr = np.ma.masked_array(arr, mask = np.isnan(arr))

                    # ratio_hospitalization_model *= np.prod(masked_arr)
                    ratio_hospitalization_model *= np.exp(np.sum(masked_arr))

        else:
            H_t = data.hospitalizations
            E_Ht_cand = basics.calculate_EHt(C_t_cand, values_cand, country)
            E_Ht_curr = basics.calculate_EHt(C_t_curr, values_curr, country)

            # arr = basics.nbinom_alt_pmf(H_t, E_Ht_cand, phi_hospitalizations_cand)/basics.nbinom_alt_pmf(H_t, E_Ht_curr, phi_hospitalizations_curr)
            # masked_arr = np.ma.masked_array(arr, mask = np.isnan(arr))

            arr = (basics.nbinom_alt_log_pmf(H_t, E_Ht_cand, phi_hospitalizations_cand)-
                    basics.nbinom_alt_log_pmf(H_t, E_Ht_curr, phi_hospitalizations_curr)
                    )
            masked_arr = np.ma.masked_array(arr, mask = np.isnan(arr))

            # ratio_hospitalization_model = np.prod(masked_arr)
            ratio_hospitalization_model = np.exp(np.sum(masked_arr))

        ratio *= ratio_hospitalization_model
        # ratio *= ratio_hospitalization_model2



    if intensiveCare_model:
        phi_intensiveCare_curr = values_curr['phi_intensiveCare']
        phi_intensiveCare_cand = values_cand['phi_intensiveCare']

        if country == 'all':
            # this part of the code can only be evaluated if ratio_likelihood
            # gets called for a parameter update (phi_deaths/hospitalizations)
            # therefore latent_variables_curr is always latent_variables_cand
            ratio_intensiveCare_model = 1
            for country_iter in data.keys():
                if pd.isna(data[country_iter]["intensiveCare"].iloc[-1]):
                        pass
                else:
                    if model == 'infections':
                        C_t_curr = latent_variables_curr['cases_view'][country_iter]
                    else:
                        C_t_curr = latent_variables_curr[country_iter]

                    Hicu_t = data[country_iter].intensiveCare
                    E_Hicut_cand = basics.calculate_EHicut(C_t_curr, values_cand, country_iter)
                    E_Hicut_curr = basics.calculate_EHicut(C_t_curr, values_curr, country_iter)
                    
                    # arr = basics.nbinom_alt_pmf(Hicu_t, E_Hicut_cand, phi_intensiveCare_cand)/basics.nbinom_alt_pmf(Hicu_t, E_Hicut_curr, phi_intensiveCare_curr)
                    # masked_arr = np.ma.masked_array(arr, mask = np.isnan(arr))

                    arr = (basics.nbinom_alt_log_pmf(Hicu_t, E_Hicut_cand, phi_intensiveCare_cand)-
                            basics.nbinom_alt_log_pmf(Hicu_t, E_Hicut_curr, phi_intensiveCare_curr)
                            )
                    masked_arr = np.ma.masked_array(arr, mask = np.isnan(arr))

                    # ratio_intensiveCare_model *= np.prod(masked_arr)
                    ratio_intensiveCare_model *= np.exp(np.sum(masked_arr))
        else:
            Hicu_t = data.intensiveCare
            E_Hicut_cand = basics.calculate_EHicut(C_t_cand, values_cand, country)
            E_Hicut_curr = basics.calculate_EHicut(C_t_curr, values_curr, country)

            # arr = basics.nbinom_alt_pmf(Hicu_t, E_Hicut_cand, phi_intensiveCare_cand)/basics.nbinom_alt_pmf(Hicu_t, E_Hicut_curr, phi_intensiveCare_curr)
            # masked_arr = np.ma.masked_array(arr, mask = np.isnan(arr))
            # ratio_intensiveCare_model = np.prod(masked_arr)

            arr = (basics.nbinom_alt_log_pmf(Hicu_t, E_Hicut_cand, phi_intensiveCare_cand)-
                    basics.nbinom_alt_log_pmf(Hicu_t, E_Hicut_curr, phi_intensiveCare_curr)
                    )
            masked_arr = np.ma.masked_array(arr, mask = np.isnan(arr))
            ratio_intensiveCare_model = np.exp(np.sum(masked_arr))

            # ratio_hospitalization_model = np.prod(
                    # basics.nbinom_alt_pmf(H_t, E_Ht_cand, phi_hospitalizations_cand)/
                    # basics.nbinom_alt_pmf(H_t, E_Ht_curr, phi_hospitalizations_curr)
                                          # )

        ratio *= ratio_intensiveCare_model


    return(ratio)




