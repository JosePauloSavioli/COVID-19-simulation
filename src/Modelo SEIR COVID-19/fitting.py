#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 09:50:56 2020

@author: jotape42p
"""

import SEIR
import lmfit

def funcao(parametros, info, data):
    
    mortes = data['obitos'][int(round(parametros['f'].value)):]
    casos_acumulados = data['casos_acumulados'][int(round(parametros['f'].value)):]
    
    #if data.get('check_I_init') is None:
    #    parametros["I_init"].value = casos_acumulados[0] #05/06 mudança para I_inicial existente nos casos acumulados
    
    if info['control'][0]:
        for i in range(info['control'][1]):
            parametros['dia' + str(i + 1)].max = len(casos_acumulados) - max(14 - info['validacao'], 0)
            parametros['dia' + str(i + 1)].min = 1 + int(round(parametros['f'].value))
    if info['control'][2]:
        for i in range(info['control'][3]):
            parametros['dob' + str(i + 1)].max = len(casos_acumulados) - max(14 - info['validacao'], 0)
    
    if info['tipo_seir'] == 'mod':
        pop = SEIR.SEIR(data['populacao'], len(casos_acumulados), parametros, data['lockdown'], len(casos_acumulados), info['control'], info['validacao'])
    if info['tipo_seir'] == 'std':
        pop = SEIR.SEIR_std(data['populacao'], len(casos_acumulados), parametros, data['lockdown'], len(casos_acumulados), info['control'], info['validacao'])
    
    model1 = pop['Acum']    #modelagem dos casos acumulados [Acum]   
    model3 = pop['D']       #modelagem dos casos acumulados de morte [D]

    l1 = (model1 - casos_acumulados)/(max(casos_acumulados)) *100000
    l3 = (model3 - mortes)/max(mortes) *100000
    
    return (l1, l3)

def fit_me(info, data, parametros):
    
    minner = lmfit.Minimizer(funcao, parametros, fcn_args=(info, data))

    result = minner.minimize(method=info['method'])
    
    return result
