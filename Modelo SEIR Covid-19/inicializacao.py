#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 09:39:36 2020

@author: jotape42p
"""

import sys
import datetime
import csv

#Função de atribuição dos parâmetros variantes no bootstrap
def lista_boot(tipo_seir):
    
    boot_iter = int(input("\nDigite o número de iterações do bootstrap [ex: 1000]:\n-> "))
    
    p_b = []
    p_b.append(input("\nDigite as variáveis a serem travadas na análise [ex: beta, alpha, mu, xi, zeta, mort_UTI]: \n\t('all': todas as permutações; \n\t 0: simulação sem variáveis fixas; \n\t 1: simulação default [mod: beta e xi variaveis, std: beta e mort_UTI [alpha] variaveis])\n-> ").replace(" ", "").split(","))
    
    if p_b[0][0] == 'all':
        
        if tipo_seir == 'mod':
            par_list = [['alpha', 'mu', 'zeta'],
                        ['alpha', 'mu', 'xi'],
                        ['alpha', 'mu', 'mort_UTI'],
                        ['alpha', 'mu', 'zeta', 'xi'],
                        ['alpha', 'mu', 'zeta', 'mort_UTI'],
                        ['alpha', 'mu', 'xi', 'mort_UTI'],
                        ['mu', 'zeta'],
                        ['mu', 'xi'],
                        ['mu', 'mort_UTI'],
                        ['mu', 'zeta', 'xi'],
                        ['mu', 'zeta', 'mort_UTI'],
                        ['mu', 'xi', 'mort_UTI'],
                        ['mu', 'zeta', 'xi', 'mort_UTI'],
                        ['alpha', 'zeta'],
                        ['alpha', 'xi'],
                        ['alpha', 'mort_UTI'],
                        ['alpha', 'zeta', 'xi'],
                        ['alpha', 'zeta', 'mort_UTI'],
                        ['alpha', 'xi', 'mort_UTI'],
                        ['alpha', 'zeta', 'xi', 'mort_UTI'],
                        ['zeta'],
                        ['xi'],
                        ['mort_UTI'],
                        ['zeta', 'xi'],
                        ['zeta', 'mort_UTI'],
                        ['xi', 'mort_UTI'],
                        ['zeta', 'xi', 'mort_UTI']]
        elif tipo_seir == 'std':
            par_list = [['alpha', 'mu'],
                        ['alpha'],
                        ['mu'],
                        ['alpha', 'mu', 't_rec'],
                        ['alpha', 't_rec'],
                        ['mu', 't_rec'],
                        ['alpha', 'mu', 'mort_UTI'],
                        ['alpha', 'mort_UTI'],
                        ['mu', 'mort_UTI']]
        return par_list, boot_iter, p_b[0][0]
    
    elif p_b[0][0] == '1':
        if tipo_seir == 'mod':
            par_list = [["alpha", "mu", "zeta", "mort_UTI"]]
        elif tipo_seir =='std':
            par_list = [["mu"]]
        return par_list, boot_iter, p_b[0][0]
    
    elif p_b[0][0] == '0':
        par_list = [['0']]
        return par_list, boot_iter, p_b[0][0]
    
    else:
        return p_b, boot_iter, p_b[0][0]

#Função de erro para o caso de digitação errada
def erro(loop, nn):
    if loop not in ('s', 'n', 'b', 'sl') and nn==0:
        loop = input("Digite 's' (loop por localidades no csv), 'b' (bootstrap) ou 'n' (localidade única): ")
        print()
        loop = erro(loop, 0)
    elif loop not in ('std', 'mod') and nn==1:
        loop = input("Digite 'std' (seir original) ou 'mod' (seir modificado com hospitalizados): ")
        print()
        loop = erro(loop, 1)
    return loop

#Função de leitura do nome da cidade
def ler_loc(loc, csv_):
    with open("./" + csv_ + ".csv") as csv_file:
        try:
            csv_reader = csv.reader(csv_file, delimiter = ';')
            for line in csv_reader:
                assert line[1]
        except:
            csv_reader = csv.reader(csv_file, delimiter = ',')
            
        for line in csv_reader:
            if line[1] == loc:
                return line[0]
        
#Função de chamada
def sistema():

    info = {}
    info['control'] = [0]*4
    
    try:
        if sys.argv[1][0] == '-':
            if len(sys.argv[1]) > 1:
                info['cod'] = sys.argv[1][1:]
            else:
                info['cod'] = None
        if sys.argv[2][0] == '-':
            info['csv_ca'] = sys.argv[2][1:]
        if sys.argv[3][0] == '-':
            info['csv_ob'] = sys.argv[3][1:]
        if sys.argv[4][0] == '-':
            mt = int(sys.argv[4][1:])
        if sys.argv[5][0] == '-':
            info['control'][0] = int(sys.argv[5][1:2])
            info['control'][1] = int(sys.argv[5][3:4])
            info['control'][2] = int(sys.argv[5][5:6])
            info['control'][3] = int(sys.argv[5][7:8])
        if sys.argv[6][0] == '-':
            info['tipo'] = sys.argv[6][1:]
            info['tipo'] = erro(info['tipo'], 0)
        if sys.argv[7][0] == '-':
            info['periodo'] = int(sys.argv[7][1:])
        if sys.argv[8][0] == '-':
            info['validacao'] = int(sys.argv[8][1:])
        if sys.argv[9][0] == '-':
            info['tipo_seir'] = sys.argv[9][1:]
            info['tipo_seir'] = erro(info['tipo_seir'], 1)
        #if sys.argv[10][0] == '-':
        #    info['IFR'] = float(sys.argv[10][1:])/100
        if len(sys.argv) == 11:
            if sys.argv[10][0] == '-':
                info['add'] = int(sys.argv[10][1:])
        else:
            info['add'] = 0
        
        info['nome'] = ler_loc(info['cod'], info['csv_ca'])
        
        lp = {"n": "simulação específica para localidade", "s": "simulação para todas as localidades do arquivo csv", "b": "simulação bootstrap para localidade", "sl": "simulação de sensibilidade dos fatores"}
        ss = {"std": "SAEIR original sem hospitalizados", "mod": "SAEIR modificado com hospitalizados (original do trabalho)"}
        sn = ["não", "sim"]
        methods = ['basinhopping', 'differential_evolution', 'powell', 'cobyla', 'nelder', 'dual_annealing']
        
        info['method'] = methods[mt]
        
        print(f"Simulation data:\n----> tipo de simulação: {lp[info['tipo']]}")
        print(f"--------> {ss[info['tipo_seir']]} {'com gráficos adicionais' if info['add'] == 1 else ''}")
        if info['tipo'] != 's':
            print(f"----> localidade: {info['nome']} -> código: {info['cod']}")
        print(f"----> arquivos: {info['csv_ca']}.csv e {info['csv_ob']}.csv\n----> método: {info['method']}\n----> discretizações de beta no tempo: {sn[info['control'][0]]}")
        if info['control'][0] == 1:
            print(f"--------> número de conjuntos dia/beta na simulação: {info['control'][1]}")
        print(f"----> discretizações da mortalidade no tempo: {sn[info['control'][2]]}")
        if info['control'][2] == 1:
            print(f"--------> número de conjuntos dob/{'obt' if info['tipo_seir'] == 'mod' else 'mort_UTI'} na simulação: {info['control'][3]}")
        print(f"----> período de simulação: {info['periodo']} dias com validação de {info['validacao']} dias no fitting")
        print(f"----> simulação iniciada em {datetime.datetime.now()}")

    except:
        raise Exception("Declaração mal escrita")
    
    return(info)
