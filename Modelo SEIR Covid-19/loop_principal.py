#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 10:37:00 2020

@author: jotape42p

Arquivo que contém o corpo da simulação e as chamadas para fitting e
resultados das simulações. Neste arquivo pode encontrar:
    
    - Funções major_boot para chamadas do bootstrap e major_loop para
        chamadas -n, -s e -sl:
        -> Chamada para criação dos parâmetros
        -> Chamada para o fitting ['fitting.py']
        -> Chamada para simulação do período escolhido em dias ['SEIR.py']
        -> Chamadas para criação de outputs ['output.py' e internas]
        -> Chamadas para criação de figuras ['pictures.py']
        -> Chamadas para cenários de lockdown
        
    - Funções específicas do major_loop:
        -> Chamada para sensibilidade
        
    - Funções específicas para bootstrap:
        -> Função interna de repetição do bootstrap ['bootstrap']
        -> Função de estimativa de incerteza e teste de normalidade
        -> Criação de arquivos dos samples
        
    -/-> Para saber sobre gráficos e outputs adicionais no modo 'add' = 1,
            procure por "if info['add']:"
        
"""

import dados
import auxiliares
import fitting
import pictures
import output
import inicializacao
import scipy.stats as stats
import time

from statsmodels.stats.diagnostic import lilliefors
from statistics import mean, stdev
from p_tqdm import p_map
from functools import partial
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import lmfit

#Seed dos números randômicos
np.random.seed(200)

#Função de bootstrap
def bootstrap(info, data, boot, parametros, boot_list):
    
    #random id
    t = 1000 * time.time()
    np.random.seed(int(t) % 2**32)
    id_ = np.random.choice(np.random.ranf(int(sum(boot_list))))
    np.random.seed(200)
    
    #separação das listas para bootstrap
    data['casos_acumulados']  = boot_list[:int(len(boot_list)/2)]
    data['obitos']            = boot_list[int(len(boot_list)/2):]
    
    info['nome'] = info['nome'] + '_' + str(id_)
    
    #fitting com dados do sample
    start_fit = datetime.now()
    boot['result'] = fitting.fit_me(info, data, parametros)
    boot['time_fit'] = datetime.now() - start_fit
    
    #simulação pós-fitting
    boot['pop'] = auxiliares.s_choice(info, data, boot)
    
    #csv do boot
    boot['check'] = True
    output.organizar_dados(info, data, boot)
    
    #gráfico do fitting
    if info['add']:
        info['boot'] = True
        pictures.grafico_simulacao(info, data, boot)
    
    return (id_, boot['result'])

#Função pré-bootstrap
def major_boot(info, data, boot):
    
    #processos iniciais
    data['lockdown'][1], parametros_boot, check = auxiliares.buscar_dados(info, data)
    
    if check:
        data['check'] = False
        boot['check'] = True
        boot['nome'] = info['nome']
        boot['error_msg'] = "Número insuficiente de dias com óbitos"
        return boot
    
    print(f"\n---> Simulação para {info['nome']} iniciada")
    
    for _ in boot['par']:
        
        if boot['v_boot'] == 'all':
            info['__nome'] = info['nome']
            basename = '_(' + '_'.join(_) + ')'
            info['nome'] = info['nome'] + basename
        else:
            basename = ''
        
        #setar para todos os passos da análise a variação dos parâmetros para True
        parametros_boot['beta'].set(vary = True)
        #parametros_boot['alpha'].set(vary = True)
        parametros_boot['mu'].set(vary = True)
        parametros_boot['mort_UTI'].set(vary = True)
        #parametros_boot['t_rec'].set(vary = True)
        
        if info['tipo_seir'] == 'mod':
            parametros_boot['xi'].set(vary = True)
            parametros_boot['zeta'].set(vary = True)    
        
        #travar parâmetros escolhidos
        if _[0] != '0':
            for l in _:
                parametros_boot[l].set(vary = False)
                
        #fitting primário para coletar erros
        start_fit = datetime.now()
        boot['result'] = fitting.fit_me(info, data, parametros_boot)
        boot['time_fit'] = str(datetime.now() - start_fit)
        
        print(f"\n---> Fitting para {info['nome']} terminado em {boot['time_fit']} minutos\n")
        
        #correção conforme parâmetro 'f' do dia inicial da simulação
        for key in ('x', 'obitos', 'casos_acumulados', 'casos_confirmados', 'obitos_novos', 'obitos_no_validation', 'casos_ac_no_validation', 'obitos_novos_no_validation', 'casos_novos_no_validation'):
            data[key] = data[key][int(round(boot['result'].params['f'].value)):]
        
        #fixando f [dias] para o boot
        for i in range(0,len(data['lockdown'][0]),2):
            parametros_boot['dia_set'+ str(int(i/2) + 1)].value -= int(round(boot['result'].params['f'].value))
        parametros_boot['f'].value = 0
        parametros_boot['f'].vary  = False
        data['check_I_init'] = False
        parametros_boot['I_init'].value = boot['result'].params['I_init'].value
        parametros_boot['I_init'].vary = False
        
        #simulação pós-fitting
        boot['pop'] = auxiliares.s_choice(info, data, boot)
        
        pictures.grafico_simulacao(info, data, boot)
        
        #erros
        boot['erros_ca'] = data['casos_acumulados'] - boot['pop']['Acum'][:len(data['x'])-info['validacao']]
        boot['erros_ob'] = data['obitos'] - boot['pop']['D'][:len(data['x'])-info['validacao']]

        #dados da simulação para csv
        boot['lista_csv'].append(output.organizar_cabecalho(info['tipo_seir'], info['control'], data['lockdown']))
        boot['check'] = False
        boot['lista_csv'].append(output.organizar_dados(info, data, boot))
        output.csv_dados_fitting(boot['lista_csv'], 'boot_valores', info['path'])
        
        #cenários de lockdown
        boot = auxiliares.cenarios_fut(info, data, boot)
        
        #arquivo de saída do fitting
        with open(info['path'] + "/Infos/fit_" + info['nome'] + '_' + info['method'] + '_' + info['tipo_seir'] + ".txt", 'w') as f:
            f.write(lmfit.fit_report(boot['result']))
        
        #amostragem (samples) para bootstrap [valores preditos são do FITTING]
        boot['ca_ob'] = np.concatenate((boot['pop']['Acum'][:len(data['x'])-info['validacao']], boot['pop']['D'][:len(data['x'])-info['validacao']]))
        boot['samples'] = []
        
        for i in range(boot['n_iteracao']):
            bs_sample = np.concatenate((np.random.choice(boot['erros_ca'], size = len(boot['erros_ca']), replace = True), 
                                        np.random.choice(boot['erros_ob'], size = len(boot['erros_ob']), replace = True)))
            
            boot['samples'].append(np.add(boot['ca_ob'], bs_sample))
            
        boot['samples'] = np.array(boot['samples'])
        
        if info['add']:
            output.csv_dados_fitting(boot['samples'], 'dados_brutos_boot', info['path'])
            b_mean = np.mean(boot['samples'], axis = 0)
            d_b = data.copy()
            d_b['casos_acumulados'], d_b['obitos'] = b_mean[:len(d_b['casos_acumulados'])], b_mean[len(d_b['obitos']):]
            pictures.grafico_simulacao(info, d_b, boot)
            
        print(f"\nIniciando bootstrap {info['nome']}\n")
        
        boot['boot'] = None
        
        #multiprocessamento e execução do fitting para cada sample do bootstrap [p_map]
        for rf in p_map(partial(bootstrap, info, data, boot, parametros_boot), boot['samples'], num_cpus = 6):
            r, ind = auxiliares.indices(rf, data, info)
            if boot.get('boot') is None:
                boot['boot'] = pd.DataFrame({str(rf[0]):r}, index=ind)
            else:
                boot['boot'][str(rf[0])] = r
        
        #criação do csv para armazenar os dados
        boot['boot'].T.to_csv(info['path'] + '/Infos/boot_'+ info['nome'] + '.csv', index = False)
        
        print("\nIniciando incerteza\n")
        
        fixed = []
        n_fixed = []
        for name, col in boot['boot'].T.iteritems():
            if len(set(col.values)) != 1 and name not in {"Akaike"}:
                me   = mean(col.values)
                qt_min = np.quantile(col.values, 0.5) - 1.5*(np.quantile(col.values, 0.75)-np.quantile(col.values, 0.25))
                qt_max = np.quantile(col.values, 0.5) + 1.5*(np.quantile(col.values, 0.75)-np.quantile(col.values, 0.25))
                qt   = (max(qt_min, min(col.values)),
                        np.quantile(col.values, 0.25), np.quantile(col.values, 0.5), np.quantile(col.values, 0.75),
                        min(qt_max, max(col.values)))
                stdh = me + stdev(col.values)
                stdl = max(me - stdev(col.values),0)
                cf90 = (me + 1.645 * stdev(col.values), max(me - 1.645 * stdev(col.values),0))
                cf95 = (me + 1.960 * stdev(col.values), max(me - 1.960 * stdev(col.values),0))
                cf99 = (me + 2.576 * stdev(col.values), max(me - 2.576 * stdev(col.values),0))
                print(f"----->{name}\toriginal mean: {boot['result'].params[name].value:.5f}\tboot mean: {me:.5f}\tmax: {stdh:.5f}\tmin: {stdl:.5f}")
                n_fixed.append((name, me, stdh, stdl, qt, cf90, cf95, cf99))
            elif len(set(col.values)) == 1:
                fixed.append(name)
        
        #outputs
        with open(info['path'] + "/Infos/fit_" + info['nome'] + '_incerteza_' + info['tipo_seir'] + ".txt", 'w') as f:
            for tupl in n_fixed:
                print('\n', tupl[0], file = f)
                print('\toriginal mean: \t', boot['result'].params[tupl[0]].value, file = f)
                print('\tmean:          \t', tupl[1], file = f)
                print('\tstd_max:       \t', tupl[2], file = f)
                print('\tstd_min:       \t', tupl[3], file = f)
                print('\tquartiles:     \t', tupl[4], file = f)
                print('\tcf 90:         \t', tupl[5], file = f)
                print('\tcf 95:         \t', tupl[6], file = f)
                print('\tcf 99:         \t', tupl[7], file = f)
                print("\tnormal [liliefors]  (p-value > 0.05):\t", lilliefors(boot['boot'].T[tupl[0]].values)[1], file = f)
                print("\tnormal [kolmogorov] (p-value > 0.05):\t", stats.kstest(boot['boot'].T[tupl[0]].values, 'norm', stats.norm.fit(boot['boot'].T[tupl[0]].values))[1], file = f)
                print("\tnormal [pearson's]  (p-value > 0.05):\t", stats.normaltest(boot['boot'].T[tupl[0]].values)[1], file = f)
                print("\tnormal [shapiro]    (p-value > 0.05):\t", stats.shapiro(boot['boot'].T[tupl[0]].values)[1], file = f)
                pictures.fig_boot(boot['boot'].T, tupl[0], info, basename)
        
        #Quartis
        boot['pop_min'] = auxiliares.stat_sim(info, data, boot, n_fixed, 0)
        boot['pop_25'] = auxiliares.stat_sim(info, data, boot, n_fixed, 1)
        boot['pop_50'] = auxiliares.stat_sim(info, data, boot, n_fixed, 2)
        boot['pop_75'] = auxiliares.stat_sim(info, data, boot, n_fixed, 3)
        boot['pop_max'] = auxiliares.stat_sim(info, data, boot, n_fixed, 4)
        
        rest = auxiliares.e_fit(boot, data, info)
        
        #voltar lockdown do csv ao valor original para os gráficos
        for i in range(0,len(data['lockdown'][0]),2):
            data['lockdown'][0][i] = int(round(boot['result'].params['dia_set'+ str(int(i/2) + 1)].value))
        pictures.grafico_simulacao(info, data, boot)
        pictures.plot_log(info, data, boot)
        pictures.grafico_populacao(info, data, boot)
        pictures.grafico_populacao(info, data, boot, 1)
        pictures.grafico_populacao(info, data, boot, 2)
        pictures.picture(info, data, boot)
        
        #figuras de saída dos casos
        pictures.plot_dif(rest, info, data)
        pictures.grafico_simulacao(info, data, boot, (0,0,'boot'))
        
        if boot['v_boot'] == 'all':
            info['nome'] = info['__nome']
        
    return boot
        
        
#Função de execução do loop pelas cidades em um csv ou para uma cidade apenas
def major_loop(info, cod = None):
    
    sim = {}
    
    #leitura de dados
    if cod is not None:
        info['cod'] = cod
        info['nome'] = inicializacao.ler_loc(info['cod'], info['csv_ca'])
    
    data = dados.initial(info)
    info['IFR'] = data['IFR']
    #checagem se há dados
    if not data['check']:
        sim['check'] = False
        sim['nome'] = info['nome']
        sim['error_msg'] = data['error_msg']
        return sim, data
    else:
        
        #leitura de dados
        if cod is not None:
            info['path'] += f"/{info['nome']}"
            auxiliares.criar_dir(info['path'])
            auxiliares.criar_dir(info['path'] + '/Figuras')
            auxiliares.criar_dir(info['path'] + '/Infos')
        
        #processos iniciais
        data['lockdown'][1], parametros_loop, check = auxiliares.buscar_dados(info, data)
        if check:
            data['check'] = False
            sim['check'] = False
            sim['nome'] = info['nome']
            sim['error_msg'] = "Número insuficiente de dias com óbitos"
            return sim, data
         
        if info['add']:
            auxiliares.triexp(data, info['path'])
            auxiliares.diffs(data, info['path'])
        
        #fitting
        start_fit = datetime.now()
        sim['result'] = fitting.fit_me(info, data, parametros_loop)
        sim['time_fit'] = str(datetime.now() - start_fit)

        if info['tipo'] != 's':
            print(f"\n---> Fitting para {info['nome']} terminado em {sim['time_fit']} minutos\n")
        
        #arquivo de saída do fitting
        with open(info['path'] + "/Infos/fit_" + info['nome'] + '_' + info['method'] + '_' + info['tipo_seir'] + ".txt", 'w') as f:
            f.write(lmfit.fit_report(sim['result']))
        
        #correção conforme parâmetro 'f' do dia inicial da simulação
        for key in ('x', 'obitos', 'casos_acumulados', 'casos_confirmados', 'obitos_novos', 'obitos_no_validation', 'casos_ac_no_validation', 'obitos_novos_no_validation', 'casos_novos_no_validation'):
            data[key] = data[key][int(round(sim['result'].params['f'].value)):]
        
        #simulação pós-fitting
        sim['pop'] = auxiliares.s_choice(info, data, sim)
       
        #sensibilidade
        if info['tipo'] == 'sl':
            auxiliares.criar_dir(info['path'] + '/Figuras/Sensibilidade')
            auxiliares.sensibilidade(info, data, sim)
            
        #cenários de lockdown
        sim = auxiliares.cenarios_fut(info, data, sim)
        
        #csv com a população dos lockdowns 
        writer = pd.ExcelWriter(info['path'] + '/Infos/Cenários_' + info['nome'] + '.xlsx', engine = 'xlsxwriter')
        
        for n, item in (('current', sim['pop']), ('no_lockdown', sim['pop_no_lock']), ('full_lockdown', sim['pop_wt_lock'])):
            x_pp = np.linspace(0, len(sim['pop']['S']) - 1, len(sim['pop']['S']))
            x_p = [data['data_inicial'] + timedelta(days = int(round(sim['result'].params['f'].value)) + x) for x in x_pp]
            cen = pd.DataFrame(item, columns = item.keys(), index = x_p)
            cen.to_excel(writer, sheet_name = n)
    
        writer.save()
        
        rest = auxiliares.e_fit(sim, data, info)
        #voltar lockdown do csv ao valor original para os gráficos
        for i in range(0,len(data['lockdown'][0]),2):
            data['lockdown'][0][i] = int(round(sim['result'].params['dia_set'+ str(int(i/2) + 1)].value))

        #criação das figuras
        pictures.plot_log(info, data, sim)
        pictures.grafico_simulacao(info, data, sim)
        pictures.plot_dif(rest, info, data)
        pictures.grafico_populacao(info, data, sim)
        pictures.picture(info, data, sim)
        
        #dados da simulação para csv
        sim['lista_csv'] = output.organizar_dados(info, data, sim)
        sim['check'] = data['check']
        sim['lock_down'] = data['lockdown']
        sim['nome'] = data['nome']
        sim['error_msg'] = None
        
        return sim,data
