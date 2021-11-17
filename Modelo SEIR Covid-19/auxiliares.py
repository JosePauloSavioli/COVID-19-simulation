#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 12:55:05 2020

@author: jotape42p
"""

import pictures
import SEIR
import parametros as par

import os
import csv
import json
import math

from scipy import optimize
from datetime import datetime
import skill_metrics as sm
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import numpy.polynomial.polynomial as npoly
import lmfit

#                                                                             #
##########################-POPULAÇÕES_ADICIONAIS-##############################
#                                                                             #

def pop_add(info, pop):
    pop['Novos_Casos'] = [pop['Acum'][0]] + np.diff(np.array(pop['Acum'])).tolist()
    pop['Mortes_novas'] = [0] + np.diff(np.array(pop['D'])).tolist()
    if info['tipo_seir'] == 'std':
        pop['I_SIR'] = np.array(pop['E']) + np.array(pop['I']) + np.array(pop['Cf']) + np.array(pop['A'])
        pop['I_SIR'] = pop['I_SIR'].tolist()
    elif info['tipo_seir'] == 'mod':
        pop['I_SIR'] = np.array(pop['E']) + np.array(pop['I']) + np.array(pop['Q']) + np.array(pop['A']) + np.array(pop['H']) + np.array(pop['U']) + np.array(pop['N'])   
        pop['I_SIR'] = pop['I_SIR'].tolist()
        pop['Total_Hospitalizados'] = np.array(pop['H']) + np.array(pop['U']) + np.array(pop['N'])
        pop['Total_Hospitalizados'] = pop['Total_Hospitalizados'].tolist()
        pop['H_tot'] = np.array(pop['U']) + np.array(pop['N'])
        pop['H_tot'] = pop['H_tot'].tolist()
    return pop

#                                                                             #
##############################-ÍNDICES_BOOT-###################################
#                                                                             #

#Índices do boot
def indices(rf, data, info):
    
    ind, it = [], []
    for i in range(0, len(data['lockdown'][0]), 2):
        ind.append('dia_set'+str(int(i / 2 + 1)))
        ind.append('tau_set'+str(int(i / 2 + 1)))
        it.append(rf[1].params['dia_set'+str(int(i / 2 + 1))].value)
        it.append(rf[1].params['tau_set'+str(int(i / 2 + 1))].value)
    for i in range(0, len(data['lockdown'][2]), 2):
        ind.append('dia'+str(int(i / 2 + 1)))
        ind.append('beta'+str(int(i / 2 + 1)))
        it.append(rf[1].params['dia'+str(int(i / 2 + 1))].value)
        it.append(rf[1].params['beta'+str(int(i / 2 + 1))].value)
    for i in range(0, len(data['lockdown'][3]), 2):
        ind.append('dob'+str(int(i / 2 + 1)))
        it.append(rf[1].params['dob'+str(int(i / 2 + 1))].value)
        if info['tipo_seir'] == 'std': 
            ind.append('mort_UTI'+str(int(i / 2 + 1)))
            it.append(rf[1].params['mort_UTI'+str(int(i / 2 + 1))].value)
        elif info['tipo_seir'] == 'mod':
            ind.append('obt'+str(int(i / 2 + 1)))
            it.append(rf[1].params['obt'+str(int(i / 2 + 1))].value)
    
    if info['tipo_seir'] == 'mod':
        r = [rf[1].params["beta"].value, rf[1].params["alpha"].value, rf[1].params["mu"].value, rf[1].params["xi"].value, rf[1].params["zeta"].value, rf[1].params["mort_UTI"].value, rf[1].aic, rf[1].params["t_rec"].value, rf[1].params["f"].value, rf[1].params["I_init"].value] + it
        ix = ['beta', 'alpha', 'mu', 'xi', 'zeta', 'mort_UTI', 'Akaike', 't_rec', 'f', 'I_init'] + ind
    elif info['tipo_seir'] == 'std':
        r = [rf[1].params["beta"].value, rf[1].params["alpha"].value, rf[1].params["mu"].value, rf[1].params["mort_UTI"].value, rf[1].aic, rf[1].params["t_rec"].value, rf[1].params["f"].value, rf[1].params["I_init"].value] + it
        ix = ['beta', 'alpha', 'mu', 'mort_UTI', 'Akaike', 't_rec', 'f', 'I_init'] + ind

    return r, ix

#                                                                             #
#############################-LER_NOMES_CSV-###################################
#                                                                             #

#Função de leitura dos nomes dos municípios
def ler_cod(csv_):
    with open("./" + csv_ + ".csv") as csv_file:
        
        try:
            csv_reader = csv.reader(csv_file, delimiter = ';')
            for line in csv_reader:
                assert line[1]
        except:
            csv_reader = csv.reader(csv_file, delimiter = ',')
        
        lista = []
        
        for line in csv_reader:
            try:
                lista.append(line[1])
            except:
                pass
            
    return lista

#                                                                             #
###############################-TAU_FUTURO-####################################
#                                                                             #

#Função de definição dos períodos de lockdown futuros
def dados_para_tau_futuro(data_inicial):
    n = input("\nDigite o número de ciclos de lockdown desejado na simulação futura:\n-> ")
    
    if n == 'não':
        print('Seu PIMPOLHO!')
        raise Exception('É um pimpolho mesmo!')
    else:
        n = int(n)
    
    if n != 0:
        ciclos = [0]*(4*n)
        for i in range(0,len(ciclos),4):
            ciclos[i] = (datetime.strptime(input("Digite o dia do lockdown " + str(i/4 + 1) + " [ex: 20/04/2020]: "), "%d/%m/%Y").date() - data_inicial.date()).days
            ciclos[i+1] = float(input("Digite o tau do ciclo de lockdown " + str(i/4 + 1) + " [ex: 0.5 para 50%]: "))
            ciclos[i+2] = int(input("Digite o número de dias em que o ciclo ficará em vigor: ")) + ciclos[i]
            ciclos[i+3] = 0.000000001
    
        print(ciclos)
        conf = input("Confirme se os dados estão corretos [y/n]")
        while conf != 'y' and conf != 'n':
            print("Favor digitar uma das opções ('y' ou 'n'): ")
            conf = input()

        if conf == 'n':
            ciclos = dados_para_tau_futuro(data_inicial)

    else:
        ciclos = []
    
    return ciclos

#                                                                             #
#################################-LOG_FIT-#####################################
#                                                                             #

#Função para fitting do tempo de recuperação ideal (INFORMATIVO -> não faz parte do fitting)
def best_log_fit(param, ca, ob):
    tr = int(round(param['t_recs'].value))
    t_ = int(round(20 - tr))
    d  = int(round(param['d'].value))
    v = ca[tr:-t_ or None]
    w = np.concatenate((ob[:d],  ob[d:] - param['f'].value))
    f = np.std(v-w)
    return f

#Função para ativação do fitting do tempo de recuperação ideal
def best_log_ativation(control, ca, ob, title, path):
    
    variaveis = {}
    
    ob = ob[max(20, 20 + np.where(ca >= 1)[0][0]):]
    ca = ca[np.where(ca >= 1)[0][0]:]
    
    param = lmfit.Parameters()
    param.add('t_recs', value = 15, min = 4, max = 16, vary = True) #valor 0 da erro no [:-0]
    if not 1:#control[2]:
        param.add('d', value = 1, min = 0.1, max = len(ob)-1, vary = False)
        param.add('f', value = 0, min = -10, max = 10, vary = False)
    else:
        param.add('d', value = 1, min = 0.1, max = len(ob)-1, vary = True)
        param.add('f', value = 0, min = -10, max = 10, vary = True)
    
    #correção para quando os dados estão essados no csv e há valores decrescentes nos acumulados 
    ob = [0.1 if x == 0 else x for x in ob]
    ca = [0.1 if x == 0 else x for x in ca]
    
    if len(ob) < 8:
        return -1
    
    minner = lmfit.Minimizer(best_log_fit, param, fcn_args=(np.log(np.array(ca)), np.log(np.array(ob))))
    result = minner.minimize(method = 'differential_evolution')
    
    variaveis['tr'] = int(round(result.params['t_recs'].value))
    variaveis['t_'] = int(round(20 - variaveis['tr']))
    variaveis['d']  = int(round(result.params['d'].value))
    variaveis['m']  = result.params['f'].value
    
    ob = ob[:variaveis['d']] +  [x/math.exp(result.params['f'].value) for x in ob[variaveis['d']:]]
    
    v = np.log(ca[variaveis['tr']:-variaveis['t_']])
    w = np.log(ob)
    
    variaveis['g'] = np.average(v-w)
    variaveis['f'] = np.std(v-w)
    
    pictures.log_fit(variaveis, v, w, path + "/Figuras/Adicionais/bestlog_" + title)
    pictures.log_fit(variaveis, np.array(ca[variaveis['tr']:-variaveis['t_']]), np.array(ob), path + "/Figuras/Adicionais/bestac_" + title)
    
    return variaveis['tr']
    
#                                                                             #
################################-PIECEWISE-####################################
#                                                                             #
    
def find_best_piecewise_polynomial(breakpoints, x, y):
    breakpoints = tuple(map(int, sorted(breakpoints)))
    xs = np.split(x, breakpoints)
    ys = np.split(y, breakpoints)
    result = []
    for xi, yi in zip(xs, ys):
        if len(xi) < 2: continue
        coefs = npoly.polyfit(xi, yi, 1)
        f = npoly.Polynomial(coefs)
        result.append([f, xi, yi])
    return result

def piecewise(x, ca, ob, n, pais, path):
    
    num_breakpoints = n
    
    plt.figure(figsize = (15,9))
    
    for item in (ca, ob):
        
        def f(breakpoints, x, y, fcache):
            breakpoints = tuple(map(int, sorted(breakpoints)))
            if breakpoints not in fcache:
                total_error = 0
                for f, xi, yi in find_best_piecewise_polynomial(breakpoints, x, y):
                    total_error += ((f(xi) - yi)**2).sum()
                fcache[breakpoints] = total_error
            return fcache[breakpoints]
        
        if len(item) != len(x):
            x = x[len(x) - len(item):]
        
        breakpoints = optimize.brute(
                f, [slice(1, len(x), 1)]*num_breakpoints, args=(x, item, {}), finish=None)
        
        if isinstance(breakpoints, np.float64):
            breakpoints = [breakpoints]
        
        plt.scatter(x, item, c='blue', s=50)
        for f, xi, yi in find_best_piecewise_polynomial(breakpoints, x, item):
            x_interval = np.array([xi.min(), xi.max()])
            plt.plot(x_interval, f(x_interval), 'ro-')

    plt.savefig(path + "/Figuras/Adicionais/piecewise"+str(n)+pais+".png")
    plt.close()
    
#                                                                             #
################################-CRIAR_DIR-####################################
#                                                                             #
    
def criar_dir(path):
    try: 
        os.makedirs(path, exist_ok = True)
    except OSError: 
        pass     

#                                                                             #
##############################-SENSIBILIDADE-##################################
#                                                                             #
        
def s_choice(info, data, sim):
    if info['tipo_seir'] == 'mod':
        popsl = SEIR.SEIR(data['populacao'], info['periodo'], sim['result'].params, data['lockdown'], len(data['x']), info['control'], info['validacao'])
        popsl = pop_add(info, popsl)
    elif info['tipo_seir'] == 'std':
        popsl = SEIR.SEIR_std(data['populacao'], info['periodo'], sim['result'].params, data['lockdown'], len(data['x']), info['control'], info['validacao'])
        popsl = pop_add(info, popsl)
    return popsl
    
def sensibilidade(info, data, sim):
    
    if info['tipo_seir'] == 'mod':
        var = (('alpha', 'mu', 'xi', 'zeta', 'mort_UTI', 'beta'), ('Z', 'D', 'EQ', 't_UTI'))
    elif info['tipo_seir'] == 'std':
        var = (('alpha', 'mu', 'mort_UTI', 'beta'), ('Z', 'D', 'EQ'))
    
    for key in sim['result'].params.valuesdict().keys():
                
        #utilizado para guardar o valor entre as trocas de valor na sensibilidade
        temp = sim['result'].params[key].value
        
        #datas da literatura recebem máximo e mínimo
        if key in var[1]:
            sim['result'].params[key].value = sim['result'].params[key].max
            popsl1 = s_choice(info, data, sim)
            sim['result'].params[key].value = sim['result'].params[key].min
            popsl2 = s_choice(info, data, sim)
            sim['result'].params[key].value = temp
        
        #datas da simulação recebem +-1 (parâmetros do fitting [f e I_inicial] também)
        elif key[:-1] == 'dia' or key[:-1] == 'dob' or key == 't_rec' or key == 'I_init':
            sim['result'].params[key].value = temp + 1
            popsl1 = s_choice(info, data, sim)
            sim['result'].params[key].value = max(temp - 1, 0)
            popsl2 = s_choice(info, data, sim)
            sim['result'].params[key].value = temp
        
        #parâmetros que variam recebem +- 0.1, parâmetros que não variam recebem máximo e mínimo
        elif key[:-1] == 'obt' or key in var[0]:
            if sim['result'].params[key].vary:
                sim['result'].params[key].value = temp + 0.05
            elif not sim['result'].params[key].vary:
                sim['result'].params[key].value = sim['result'].params[key].max
            
            popsl1 = s_choice(info, data, sim)
            
            if sim['result'].params[key].vary:
                sim['result'].params[key].value = max(temp - 0.05, 0)
            elif not sim['result'].params[key].vary:
                sim['result'].params[key].value = sim['result'].params[key].min
            
            popsl2 = s_choice(info, data, sim)
            
            sim['result'].params[key].value = temp
        else:
            continue
        
        pictures.grafico_simulacao(info, data, sim, (popsl1, popsl2, key))
    
#                                                                             #
##############################-TAYLOR_DIAGRAM-#################################
#                                                                             #
    
def taylor(rest, validacao, ca, ob, path):
    val_c = sm.taylor_statistics(rest['Acumulados'][-validacao:], ca[-validacao:])
    val_o = sm.taylor_statistics(rest['Obitos'][-validacao:], ob[-validacao:])
    
    stdev_c = np.array([val_c['sdev'][0],  val_c['sdev'][1]])
    crmsd_c = np.array([val_c['crmsd'][0], val_c['crmsd'][1]])
    ccoef_c = np.array([val_c['ccoef'][0], val_c['ccoef'][1]]) 
    
    stdev_o = np.array([val_o['sdev'][0],  val_o['sdev'][1]])
    crmsd_o = np.array([val_o['crmsd'][0], val_o['crmsd'][1]])
    ccoef_o = np.array([val_o['ccoef'][0], val_o['ccoef'][1]])
    
    sm.taylor_diagram(stdev_c, crmsd_c, ccoef_c)
    plt.savefig(path + '/Figuras/Adicionais/taylor_casos_acumulados.png')
    plt.close()
    sm.taylor_diagram(stdev_o, crmsd_o, ccoef_o)
    plt.savefig(path + '/Figuras/Adicionais/taylor_obitos.png')
    plt.close()
    
#                                                                             #
###########################-INICIALIZAÇÕES_SIM-################################
#                                                                             #

#Inicializações
def buscar_dados(info, data):
    
    #checagem se é um loop pelo csv para o tau futuro
    if info['tipo'] == 'n':
        var_lock = dados_para_tau_futuro(data['data_inicial'])
    else:
        var_lock = []
        
    #gráficos adicionais
    if info['add']:
        criar_dir(info['path'] + '/Figuras/Adicionais')
        piecewise(data['x'], np.log(data['casos_acumulados'])[np.where(data['casos_acumulados'] >= 1)[0][0]:], np.log(data['obitos'])[np.where(data['obitos'] >= 1)[0][0]:], 1, info['nome'], info['path'])
        piecewise(data['x'], np.log(data['casos_acumulados'])[np.where(data['casos_acumulados'] >= 1)[0][0]:], np.log(data['obitos'])[np.where(data['obitos'] >= 1)[0][0]:], 2, info['nome'], info['path'])
    
    t_rec = best_log_ativation(info['control'], data['casos_acumulados'], data['obitos'], info['nome'] + '_' + info['tipo_seir'], info['path'])
    
    if t_rec == -1:
        return False, False, True
    
    #parâmetros
    parametros = par.parametros(data['lockdown'], len(data['casos_acumulados']), info['validacao'], info['control'], info['tipo_seir'], t_rec, data, info['IFR'])
    
    return var_lock, parametros, False

#                                                                             #
##########################-CENARIOS_LOCK_FUTURO-###############################
#                                                                             #

#Lockdowns de cenários futuros
def cenarios_fut(info, data, sim):
    
    #sem lock
    data['lockdown'][1].extend([len(data['casos_ac_no_validation']), 0.0])
    sim['pop_no_lock'] = s_choice(info, data, sim)
    #com lock de 15% a mais que o atual
    var = max(list(filter(lambda x: isinstance(x, float), data['lockdown'][0]))) if data['lockdown'][0] else 0.3
    data['lockdown'][1][-1] = min(0.89, 0.15 + var)
    sim['pop_wt_lock'] = s_choice(info, data, sim)
    #reset no valor
    data['lockdown'][1].pop(), data['lockdown'][1].pop()
    
    return sim

#                                                                             #
############################-ERROS_DE_FITTING-#################################
#                                                                             #

#Função para gravar dados de erros de fitting
def e_fit(sim, data, info):
    
    #csv com diferenças totais e percentuais
    rest = pd.DataFrame({'Acumulados': sim['pop']['Acum'][:len(data['casos_acumulados'])], 
                         'Dif_acum': sim['pop']['Acum'][:len(data['casos_acumulados'])] - data['casos_acumulados'], 
                         'Dif_acum_pc': np.divide(sim['pop']['Acum'][:len(data['casos_acumulados'])] - data['casos_acumulados'],data['casos_acumulados']), 
                         'Obitos': sim['pop']['D'][:len(data['obitos'])], 
                         'Dif_obitos': sim['pop']['D'][:len(data['obitos'])] - data['obitos'], 
                         'Dif_obitos_pc': np.divide(sim['pop']['D'][:len(data['obitos'])] - data['obitos'], data['obitos'])})
    
    #gráficos de taylor
    if info['add']:
        taylor(rest, info['validacao'], data['casos_acumulados'], data['obitos'], info['path'])
    
    #salvar csv
    rest.T.to_csv(info['path'] + "/Infos/diferenças_" + info['nome'] + '_' + info['tipo_seir'] + '.csv', mode = 'a')
    
    return rest

#                                                                             #
###############################-QUARTIS_SIM-###################################
#                                                                             #

def stat_sim(info, data, boot, n_fixed, quartil_id):
    
    for tp in n_fixed:
        boot['result'].params[tp[0]].value = tp[4][quartil_id]
    pop = s_choice(info, data, boot)
    
    return pop

#                                                                             #
##################################-JSON-#######################################
#                                                                             #
    
def criar_json_site(info, data, sim):
    dict_all = {}
    
    dict_all['pop']             = sim['pop']
    dict_all['data_inicial']    = datetime.strftime(data['data_inicial'], '%d/%m/%Y')
    dict_all['result']          = sim['result'].params['f'].value
    dict_all['nome']            = info['nome']
    dict_all['x']               = data['x'].tolist()
    if info['tipo'] == 'b':
        dict_all['pop_max']         = sim['pop_max']
        dict_all['pop_min']         = sim['pop_min']
        dict_all['pop_25']          = sim['pop_25']
        dict_all['pop_75']          = sim['pop_75']
        dict_all['pop_50']          = sim['pop_50']
    dict_all['casos_acumulados']= data['casos_acumulados'].tolist()
    dict_all['casos_ac_no_validation']= data['casos_ac_no_validation'].tolist()
    dict_all['obitos_no_validation']= data['obitos_no_validation'].tolist()
    dict_all['casos_novos_no_validation']= data['casos_novos_no_validation'].tolist()
    dict_all['obitos_novos_no_validation']= data['obitos_novos_no_validation']
    dict_all['pop_no_lock']     = sim['pop_no_lock']
    dict_all['pop_wt_lock']     = sim['pop_wt_lock']
    dict_all['path']            = info['path']
    
    with open(info['path'] + "/Infos/json_infos_" + info['cod'] + ".json", 'w') as f:
        json.dump(dict_all, f, indent = 4)
    
#                                                                             #
##############################-TAU_TRIGGER-####################################
#                                                                             #

def triexp_model(param, data, n):
    func = param['a'].value * np.exp(param['b'].value * data['x'])
    return func-data['casos_acumulados']
    
def log_model(param, data, n):
    func = param['c'].value / (1 + np.exp( - param['d'].value * (data['x']-param['x0'].value)))
    return func-data['casos_acumulados']

def triexp(data, path):
    
    params = lmfit.Parameters()
    params.add('a', value = 1, min = 0, max = 10000, vary = True)
    params.add('b', value = 1, min = 0, max = 1, vary = True)
    
    minner = lmfit.Minimizer(triexp_model, params, fcn_args=(data,0))
    result = minner.minimize(method = 'differential_evolution')
    
    params2 = lmfit.Parameters()
    params2.add('c', value = max(data['casos_acumulados']), min = 0, max = 10000000000, vary = False)
    params2.add('d', value = 1, min = 0, max = 100, vary = True)
    params2.add('x0', value = 1, min = 10, max = len(data['casos_acumulados'])-10, vary = True)
    
    minner2 = lmfit.Minimizer(log_model, params2, fcn_args=(data,0))
    result2 = minner2.minimize(method = 'differential_evolution')
    
    plt.figure(figsize = (15,9))
    plt.title(f"{result.params['a'].value}*exp({result.params['b'].value:.3f}*x) and {result2.params['c'].value}/(1+exp(-{result2.params['d'].value:.3f}*(x-{result2.params['x0'].value:.3f})))")
    plt.plot(data['casos_acumulados'])
    plt.plot(result.params['a'].value * np.exp(result.params['b'].value * data['x']))
    plt.plot(result2.params['c'].value / (1 + np.exp( - result2.params['d'].value * (data['x']-result2.params['x0'].value))))
    
    plt.savefig(path + '/Figuras/Adicionais/exp.png')
    plt.close()
    
#                                                                             #
##################################-DIFF-#######################################
#                                                                             #  
    
def diffs(data, path):
    d = data['casos_confirmados'].tolist()
    difs = [1/280*d[d.index(x)-4]-4/105*d[d.index(x)-3]+1/5*d[d.index(x)-2]-4/5*d[d.index(x)-1]+4/5*d[d.index(x)+1]-1/5*d[d.index(x)+2]+4/105*d[d.index(x)+3]-1/280*d[d.index(x)+4] for x in d[4:-4]]
    plt.plot(difs)
    plt.savefig(path + '/Figuras/Adicionais/diff.png')
    plt.close()