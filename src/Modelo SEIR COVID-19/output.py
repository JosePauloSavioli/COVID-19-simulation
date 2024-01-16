#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 03:22:23 2020

@author: jotape42p
"""

import csv
import datetime
import numpy as np

#Função que escreve os dados de saída da simulação em csv
def csv_dados_fitting(entrada, nome, path):
    
    with open(path + '/Infos/' + nome + '.csv', 'a+') as csv_file:
        writer = csv.writer(csv_file, delimiter = ';')
        writer.writerows(entrada)

#Função para escrever o cabeçalho do csv
def organizar_cabecalho(tipo_seir, control, lockdown):
    
    col_names = ['Nome do país', 
                 'Código', 
                 'Tempo de execução',
                 'Dias de fitting', 
                 'Akaike',
                 'Total de óbitos', 
                 'Total de casos acumulados', 
                 'Total de casos recuperados', 
                 'Pico de casos diários novos', 
                 'Dia do pico de casos diários novos', 
                 'Pico de infectados notificados (sintomáticos - I)', 
                 'Dia do pico de infectados notificados (sintomáticos - I)', 
                 'Pico de infectados não notificados (assintomáticos - A)', 
                 'Dia do infectados não notificados (assintomáticos - A)']
    
    if tipo_seir == 'mod':
        
        col_names.extend(['Pico de casos leves (quarentenados - Q)',
                          'Dia do pico de casos leves (quarentenados - Q)',
                          'Pico de casos severos (hospitalizações - H)',
                          'Dia do pico de casos severos (hospitalizações - H)',
                          'Pico de casos críticos (UTI - U)',
                          'Dia do pico de casos críticos (UTI - U)',
                          'Data de passagem de leitos', 
                          'alpha', 
                          'mu', 
                          'zeta', 
                          'mortalidade_UTI', 
                          'beta [0]', 
                          'R [0]', 
                          'xi [0]'])
        
    elif tipo_seir == 'std':
        
        col_names.extend(['alpha', 
                          'mu', 
                          'beta [0]', 
                          'R [0]', 
                          'mortalidade [0]'])
    
    if lockdown[0]:
        for i in range(0, len(lockdown[0]), 2):
            col_names.extend(['Data_l [' + str(i+1) + ']', 'tau [' + str(i+1) + ']', 'beta [' + str(i+1) + ']', 'R [' + str(i+1) + ']'])
        
    if control[0] == 1:
        for i in range(control[1]):
            col_names.extend(['Data_b [' + str(i+1) + ']', 'beta [' + str(i+1) + ']', 'R [' + str(i+1) + ']'])
    
    if control[2] == 1:
        for i in range(control[3]):
            if tipo_seir == 'mod':
                col_names.extend(['Data_m [' + str(i+1) + ']', 'obt [' + str(i+1) + ']', 'xi [' + str(i+1) + ']'])
            elif tipo_seir == 'std':
                col_names.extend(['Data_m [' + str(i+1) + ']', 'mortalidade [' + str(i+1) + ']'])
    
    return col_names
            
#Função para organizar os dados resultantes do fitting para o csv
def organizar_dados(info, data, sim):
    
    #Organizar dados de lockdown
    
    l = []
    
    for item in (data['lockdown'][0], data['lockdown'][1]):
        for i, value in enumerate(item):
            if not i%2:
                l.append(item[i])
            else:
                l.append(item[i])
                l.append(sim['result'].params['beta'].value * item[i])
                l.append(sim['result'].params['beta'].value * item[i] * (sim['result'].params['mu'].value * (1 - sim['result'].params['alpha'].value) * sim['result'].params['D'].value + sim['result'].params['alpha'].value * sim['result'].params['EQ'].value))
        
    for i, value in enumerate(data['lockdown'][2]):
        if not i%2:
            l.append(data['lockdown'][2][i])
        else:
            l.append(data['lockdown'][2][i])
            l.append(sim['result'].params['beta'].value * (sim['result'].params['mu'].value * (1 - sim['result'].params['alpha'].value) * sim['result'].params['D'].value + sim['result'].params['alpha'].value * sim['result'].params['EQ'].value))
    
    for i, value in enumerate(data['lockdown'][3]):
        if not i%2:
            l.append(data['lockdown'][3][i])
        else:
            l.append(data['lockdown'][3][i])
            if info['tipo_seir'] == 'mod':
                l.append(sim['result'].params['xi'].value * data['lockdown'][3][i] ** (1/3))
    
    #variável indicativa da passagem do número de leitos
    if sim['pop'].get('U') is not None:
        if len(np.where(np.array(sim['pop']["U"]) >= data['leitos'])[0]) == 0:
            dia_uti = "Não passa"
        elif data['leitos'] < 10:
            dia_uti = "Número de leitos nulo ou menor que 10 para o país"
        else:
            dia_uti = (data['data_inicial'] + datetime.timedelta(days = int(np.where(np.array(sim['pop']["U"]) >= data['leitos'])[0][0]))).strftime("%d/%m/%Y")

    #lista de informações para a localidade
    
    dados = [info['nome'], 
             info['cod'],
             sim['time_fit'],
             len(data['x']),
             sim['result'].aic,
             int(max(max(sim['pop']["D"]), max(data['obitos']))),
             max(sim['pop']["Acum"]), 
             max(sim['pop']["R"]), 
             max(sim['pop']["Novos_Casos"]), 
             (data['data_inicial'] + datetime.timedelta(days = int(np.where(np.array(sim['pop']["Novos_Casos"]) == max(sim['pop']["Novos_Casos"]))[0][0]))).strftime("%d/%m/%Y"), 
             max(sim['pop']["I"]), 
             (data['data_inicial'] + datetime.timedelta(days = int(np.where(np.array(sim['pop']["I"]) == max(sim['pop']["I"]))[0][0]))).strftime("%d/%m/%Y"), 
             max(sim['pop']["A"]), 
             (data['data_inicial'] + datetime.timedelta(days = int(np.where(np.array(sim['pop']["A"]) == max(sim['pop']["A"]))[0][0]))).strftime("%d/%m/%Y")]

    if info['tipo_seir'] == 'mod':
        dados.extend([max(sim['pop']['Q']),
                      (data['data_inicial'] + datetime.timedelta(days = int(np.where(np.array(sim['pop']["Q"]) == max(sim['pop']["Q"]))[0][0]))).strftime("%d/%m/%Y"), 
                      max(sim['pop']['H']),
                      (data['data_inicial'] + datetime.timedelta(days = int(np.where(np.array(sim['pop']["H"]) == max(sim['pop']["H"]))[0][0]))).strftime("%d/%m/%Y"), 
                      max(sim['pop']['U']),
                      (data['data_inicial'] + datetime.timedelta(days = int(np.where(np.array(sim['pop']["U"]) == max(sim['pop']["U"]))[0][0]))).strftime("%d/%m/%Y"), 
                      dia_uti, 
                      sim['result'].params['alpha'].value,
                      sim['result'].params['mu'].value,
                      sim['result'].params['zeta'].value,
                      sim['result'].params['mort_UTI'].value,
                      sim['result'].params['beta'].value,
                      sim['result'].params['R'].value,
                      sim['result'].params['xi'].value])
   
    elif info['tipo_seir'] == 'std':
        dados.extend([sim['result'].params['alpha'].value,
                      sim['result'].params['mu'].value,
                      sim['result'].params['beta'].value,
                      sim['result'].params['R'].value,
                      sim['result'].params['mort_UTI'].value])
        
    dados.extend(l)
    dados.extend(data['lockdown'])
    
    if info['tipo'] == 'b' and sim['check']:
        csv_dados_fitting([dados], "boot_valores", info['path'])
    else:
        return dados
    