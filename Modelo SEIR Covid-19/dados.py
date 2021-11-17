#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 09:51:41 2020

@author: jotape42p
"""

import csv
import numpy as np
from datetime import datetime, timedelta
        
def conferir(csv_file):
    try:
        csv_reader = csv.reader(csv_file, delimiter = ';')
        for line in csv_reader:
            assert line[1]
    except:
        csv_reader = csv.reader(csv_file, delimiter = ',')
        
    return csv_reader
        
#Função principal de leitura do csv
def ler_csv(cod, csv_ca, csv_ob, val):
    
    data = {}

    #leitura dos dados do csv
    with open("./" + csv_ca + ".csv") as csv_file:
        
        csv_reader = conferir(csv_file)
        
        #loop pelo arquivo
        for line in csv_reader:
            if line[1] == cod:
                data['nome'] = line[0]
                #população
                data['populacao'] = int(line[2])
                #lockdown
                data['lockdown'] = list(filter(None, line[5:16]))
                #número de leitos
                if line[3] != '':
                    data['leitos'] = int(line[3])
                else:
                    data['leitos'] = 0
                #dia inicial dos dados    
                day_std_csv = line[4]
                #IFR
                data['IFR'] = float(line[17])
                #casos confirmados
                casos_confirmados_csv = np.array(list(map(int, line[18:])))
                
                #checagem se há casos na localidade, corte dos casos confirmados (novos) no dia do primeiro caso e atribuição dos casos acumulados
                if not (casos_confirmados_csv >= 1).any():
                    data['check'] = False
                    data['error_msg'] = "Não há casos confirmados na localidade"
                    break
                else:
                    casos_confirmados = casos_confirmados_csv[np.where(casos_confirmados_csv >= 1)[0][0]:]
                    data['casos_confirmados'] = casos_confirmados
                    data['casos_acumulados'] = np.cumsum(casos_confirmados)
                    data['check'] = True
                    
                #limitação para o número mínimo de casos e verificação dos dados
                if (max(data['casos_acumulados']) < 400 or len(data['casos_acumulados']) < 20):
                    data['check'] = False
                    data['error_msg'] = "Número de casos ou dias insuficiente para a simulação (<400 casos acumulados ou <20 dias)"
                    break

                #correção da data inicial
                data['data_inicial'] = datetime.strptime(day_std_csv, '%d/%m/%Y') + timedelta(days = len(casos_confirmados_csv) - len(casos_confirmados))

                break

    with open("./" + csv_ob + ".csv") as csv_file:
        
        csv_reader = conferir(csv_file)
        
        #loop pelo arquivo
        for line in csv_reader:
            if line[1] == cod:
            
                #mortes acumuladas
                mortes_csv = np.array(list(map(int, line[18:])))
                
                #checagem se há mortes
                if not (mortes_csv >= 1).any():
                    data['check'] = False
                    data['error_msg'] = "Não há mortes confirmadas na localidade"
                    break
                elif (len(mortes_csv) - np.where(mortes_csv >= 1)[0][0]) <= val:
                    data['check'] = False
                    data['error_msg'] = "Dias para validação maiores que os com mortes confirmadas na localidade"
                    break
                elif (len(mortes_csv) - np.where(mortes_csv >= 1)[0][0]) < 8:
                    data['check'] = False
                    data['error_msg'] = "Número insuficiente de dias com óbitos"
                    break
                
                #checagem de comprimento e corte no vetor de mortes do csv onde houve o primeiro caso
                if len(casos_confirmados_csv) != len(mortes_csv):
                    data['check'] = False
                    data['error_msg'] = "Número de dias de dados diferente nos csvs"
                    break
                else:
                    data['obitos'] = mortes_csv[len(casos_confirmados_csv) - len(casos_confirmados):]
                    data['obitos_novos'] = [0] + np.diff(np.array(data['obitos'])).tolist()
                
                break

    return data
    
#Função de pré processamento dos dados do csv
def initial(info):
    
    #recolhimento dos dados
    data = ler_csv(info['cod'], info['csv_ca'], info['csv_ob'], info['validacao'])
    
    #checagem inicial para localidades sem dados e comprimento das listas de casos e mortes acumuladas do csv
    if not data['check']:
        return data
    else:
        
        #dia inicial como 20 dias antes da primeira morte
        dia_0 = int(np.where(data['obitos'] >= 1)[0][0])
        
        if dia_0 <= 20:
            data['obitos'] = np.concatenate((np.zeros(20 - dia_0), data['obitos']))
            data['casos_confirmados'] = np.concatenate((np.zeros(20 - dia_0), data['casos_confirmados']))
            data['obitos_novos'] = np.concatenate((np.zeros(20 - dia_0), data['obitos_novos']))
            data['casos_acumulados'] = np.concatenate((np.zeros(20 - dia_0), data['casos_acumulados']))
        elif dia_0 > 20:
            data['obitos'] = data['obitos'][dia_0 - 20:]
            data['casos_confirmados'] = data['casos_confirmados'][dia_0 - 20:]
            data['obitos_novos'] = data['obitos_novos'][dia_0 - 20:]
            data['casos_acumulados'] = data['casos_acumulados'][dia_0 - 20:]
        
        #vetor eixo x
        data['x'] = np.linspace(0, len(data['casos_acumulados']) - 1, len(data['casos_acumulados']))
        
        data['obitos_no_validation'] = data['obitos']
        data['casos_ac_no_validation'] = data['casos_acumulados']
        data['obitos_novos_no_validation'] = data['obitos_novos']
        data['casos_novos_no_validation'] = data['casos_confirmados']
        #vetor de dias do fitting
        if info['validacao'] != 0:
            for key in ('obitos', 'casos_acumulados', 'obitos_novos', 'casos_confirmados'):
                 data[key] = data[key][:-info['validacao']]
        
        #recorte dos dados de lockdown
        if data['lockdown']:
            for i in range(0, len(data['lockdown']), 2):
                data['lockdown'][i]     = max((datetime.strptime(data['lockdown'][i], '%d/%m/%Y').date() - data['data_inicial'].date()).days, 1)
                data['lockdown'][i + 1] = float(data['lockdown'][i + 1].replace(',','.'))
        
        data['lockdown'] = [data['lockdown'], [], [], []]

        return data
        
