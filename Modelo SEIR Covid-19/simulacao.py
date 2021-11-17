#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 22:09:12 2020

@author: jotape42p

Arquivo principal da simulação SEIR e SEIR com hospitalização feita para
o laboratório GYRO. Neste arquivo pode encontrar:
    
    - Chamadas principais das funções:
        -> Inicialização: dados recolhidos do terminal
        -> Escolha entre simulações [-n, -s, -sl, -b] e funções específicas
        -> Criação dos diretórios padrões [específicos em 'loop_principal.py']
        -> Chamada para função principal em 'loop_principal.py'
        -> Outputs finais para csv com informações gerais ('lista_csv')
        
"""


#Remoção das mensagens de warning das matrizes de covariância
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

#Importações de funções locais
import inicializacao
import loop_principal
import dados
import auxiliares
import output

#Importação de funções adicionais
import json
import numpy as np
from datetime import datetime
from functools import partial
from p_tqdm import p_map



#Função principal (main)
if __name__ == '__main__':

    np.random.seed(200)
    
    #Tempo inicial da simulação
    start_all = datetime.now()
    
    #leitura das variáveis dadas pelo sistema
    info = inicializacao.sistema()
    
    #variáveis iniciais
    time_dict = {}
    lista_csv = []
    
    print("\nIniciando simulação")
    
    #simulação única
    if info['tipo'] in ('n', 'sl'):
        
        #criar diretórios
        info['path'] = f"./Simulações/Sim_{info['nome']}_{info['tipo']}_{info['tipo_seir']}_{info['control'][0]}{info['control'][1]}{info['control'][2]}{info['control'][3]}_{info['method']}"
        auxiliares.criar_dir(info['path'])
        auxiliares.criar_dir(info['path'] + '/Figuras')
        auxiliares.criar_dir(info['path'] + '/Infos')
        
        print(f"\n---> Simulação para {info['nome']} iniciada")
        
        #simulação
        sim, data = loop_principal.major_loop(info)
        
        if not sim['check']:
            print(f"--> {info['nome']}: {sim['error_msg']}")
        else:
            lista_csv.append(output.organizar_cabecalho(info['tipo_seir'], info['control'], sim['lock_down']))
            #auxiliares.criar_json_site(info, data, sim)
            lista_csv.append(sim['lista_csv'])
            output.csv_dados_fitting(lista_csv, f"{info['nome']}_dados_principais", info['path'])

            print(f"\nSimulação para {info['nome']} terminada em {datetime.now() - start_all} minutos")
    
    #simulação em loop por csv
    elif info['tipo'] == 's':
        
        #criar diretórios
        info['path'] = f"./Simulações/DB_{info['csv_ca']}_{info['tipo']}_{info['tipo_seir']}_{info['control'][0]}{info['control'][1]}{info['control'][2]}{info['control'][3]}_{info['method']}"
        auxiliares.criar_dir(info['path'])
        auxiliares.criar_dir(info['path'] + '/Infos')
        
        print(f"\n---> Simulação para {info['csv_ca']} iniciada")
        
        lista = auxiliares.ler_cod(info['csv_ca'])
        
        #simulação com multiprocessamento
        for sim, data in p_map(partial(loop_principal.major_loop, info), lista, num_cpus = 7):
            if not sim['check']:
                print(f"--> {sim['nome']}: {sim['error_msg']}")
            else:
                time_dict[info['nome']] = str(sim['time_fit'])
                lista_csv.insert(0, output.organizar_cabecalho(info['tipo_seir'], info['control'], sim['lock_down']))        
                lista_csv.append(sim['lista_csv'])

        time_dict['0_total'] = str(datetime.now() - start_all)
        output.csv_dados_fitting(lista_csv, f"{info['csv_ca']}_dados_principais", info['path'])
        
        #salvar dados de tempo de execução
        with open(info['path'] + "/Infos/tempos_" + info['csv_ca'] + ".json", 'w') as fp:
            json.dump(time_dict, fp, indent = 4)
        
        print(f"\nSimulação para {info['csv_ca']}.csv terminada em {time_dict['0_total']} minutos")
        
    #bootstrap
    elif info['tipo'] == 'b':
        
        #criar diretórios
        info['path'] = f"./Simulações/Sim_{info['nome']}_{info['tipo']}_{info['tipo_seir']}_{info['control'][0]}{info['control'][1]}{info['control'][2]}{info['control'][3]}_{info['method']}"
        auxiliares.criar_dir(info['path'])
        auxiliares.criar_dir(info['path'] + '/Figuras')
        auxiliares.criar_dir(info['path'] + '/Infos')
        auxiliares.criar_dir(info['path'] + '/Figuras/Boot')
        auxiliares.criar_dir(info['path'] + '/Figuras/Adicionais')
        
        boot = {}
        
        boot['par'], boot['n_iteracao'], boot['v_boot'] = inicializacao.lista_boot(info['tipo_seir'])
        boot['lista_csv'] = lista_csv
        
        #leitura dos dados da localidade
        data = dados.initial(info)
        info['IFR'] = data['IFR']
        #checagem se há dados
        if not data['check']:
            print(f"--> {info['nome']}: {sim['error_msg']}")
        else:
            
            #simulação
            boot = loop_principal.major_boot(info, data, boot)
            
            if boot['check']:
                print(f"--> {info['nome']}: {boot['error_msg']}")
            else:
                print(f"\nSimulação para {info['nome']} terminada em: {datetime.now() - start_all} minutos\n")
