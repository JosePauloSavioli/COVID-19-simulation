#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 11:57:43 2020

@author: jotape42p
"""

import lmfit
import numpy as np

#
#Função específica para inicialização dos parâmetros default (pode ser modificada remotamente no bootstrap ou no uso de std [que não usa xi nem zeta])
#
#Fontes:
#       alpha      --> 0.69  [0.66; 0.71]   Primeiras cidades China https://science.sciencemag.org/content/368/6490/489/tab-pdf
#                   -> 0.667 [0.417; 0.717] Diamond Princess        https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.10.2000180#r13
#       mu         --> 0.42  [0.43; 0.61]   Primeiras cidades China https://science.sciencemag.org/content/368/6490/489/tab-pdf
#                   -> 0.5                  Grã-Bretanha and USA    https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf
#       zeta       --> 0.36                 Austrália               https://www.health.gov.au/resources/current-covid-19-cases-in-hospitals-and-intensive-care-units-icus
#                   -> 0.2528               China                   http://weekly.chinacdc.cn/en/article/id/e53946e2-c6c4-41e9-9a9b-fea8db1a8f51
#                   -> 0.16                 Itália                  https://www.health.gov.au/resources/current-covid-19-cases-in-hospitals-and-intensive-care-units-icus; http://weekly.chinacdc.cn/en/article/id/e53946e2-c6c4-41e9-9a9b-fea8db1a8f51 / https://jamanetwork.com/journals/jama/fullarticle/2763188
#       mort_UTI   --> 0.49                 China                   http://weekly.chinacdc.cn/en/article/id/e53946e2-c6c4-41e9-9a9b-fea8db1a8f51
#                   -> 0.62 and 0.67        Estudos em UTI únicas   https://www.thelancet.com/pdfs/journals/lanres/PIIS2213-2600(20)30161-2.pdf
#       Z          --> 5.84  [4.83; 6.85]   Compilação de estudos   https://www.medrxiv.org/content/10.1101/2020.04.01.20050138v1.full.pdf
#       D          --> 3.32  [2.92; 4.04]   Primeiras cidades China https://science.sciencemag.org/content/368/6490/489/tab-pdf
#       EQ         --> 4.82  [3.48; 6.15]   Compilação de estudos   https://www.medrxiv.org/content/10.1101/2020.04.01.20050138v1.full.pdf
#       t_UTI      --> 5.66  [4.98; 6.34]   Compilação de estudos   https://www.medrxiv.org/content/10.1101/2020.04.01.20050138v1.full.pdf
#

def parametros(lockdown, dias_fit, validacao, control, tipo_seir, t_rec, data, ifr):
    
    t_rec = 6.76
    
    #cálculo alpha
    #cfr         = max(data['obitos'])/max(data['casos_acumulados'])
    #a_i_ratio   = (cfr/ifr - 1)/(cfr/ifr)
    
    #Inicialização da classe parâmetros do fitting
    parametros = lmfit.Parameters()
    
    #
    #parâmetros para fitting de casos acumulados
    #
    
    #beta: fator de transmissibilidade (para quantas pessoas um infectado [notificado {sintomático}/não notificado {assintomático}] passa a doença em 1 dia)
    parametros.add('beta',      value = 1,      min = 0.0001,   max = 3,        vary = True  )
    #alpha: porcentagem de casos notificados {sintomáticos} em relação ao total de infectados
    #parametros.add('alpha',     value = 1 - (a_i_ratio),   min = 0,     max = 1,     vary = False ) #0.9% é o IFR do Reino Unido pelo Imperial College
    #mu: fator de redução da transmissibilidade para casos não notificados {assintomáticos}
    parametros.add('mu',        value = 0.5,   min = 0.34,     max = 0.61,     vary = False )
    
    #
    #parâmetros para fitting de óbitos
    #
    
    if tipo_seir == 'mod':
        #xi: porcentagem de casos severos {hospitalizados} em relação ao total de casos notificados no modelo 'mod'
        parametros.add('xi',        value = 0.735,  min = 0.01,     max = 1,        vary = True  )
        #zeta: porcentagem de casos de críticos {UTI} em relação ao total de casos hospitalizados no modelo 'mod'
        parametros.add('zeta',      value = 0.2528, min = 0.16,     max = 0.36,     vary = False )
        #mort_UTI: porcentagem de óbitos em relação ao total de casos críticos {UTI} no modelo 'mod'
        parametros.add('mort_UTI',  value = 0.555,  min = 0.49,     max = 0.67,     vary = False )
    elif tipo_seir == 'std':
        #mort_UTI: porcentagem de óbitos em relação ao total da casos confirmados no modelo 'std'
        parametros.add('mort_UTI',  value = 0.055,   min = 0.001,   max = 0.4,      vary = True  )
    
    #
    #parâmetros para início do fitting (ajuste variável por país do delay inicial em dias {f} a dos casos iniciais {I_init})
    #
    
    parametros.add('f',         value = 7,      min = 1,           max = 20,                           vary = False  )# DAY OF ACRE #min(np.where(data['casos_acumulados'] > 0)[0][0], min(20, int(dias_fit/3))-1),        max = min(20, int(dias_fit/3))
    parametros.add("I_init",    value = 1,      min = 0.0001,      max = data['casos_acumulados'][20], vary = True  )#data['casos_acumulados'][min(20, int(dias_fit/3))],       vary = True  ) #5/6 mudança para variação por valor atribuido no fitting
    
    #
    #parâmetros de tempo entre estágios da doença
    #
    
    #Z: dias entre contágio e apresentação de sintomas
    parametros.add('Z',         value = 5.84,   min = 4.83,     max = 6.85,     vary = False )
    #D: dias de transmissão de casos não notificados {assintomáticos}
    parametros.add('D',         value = 3.32,   min = 2.92,     max = 4.04,     vary = False )
    #EQ: dias de transmissão de casos notificados {sintomáticos}
    parametros.add('EQ',        value = 4.82,   min = 3.48,     max = 6.15,     vary = False )
    #NOTA: para a simulação 'mod' o tempo de recuperação da doença é t_UTI + t_rec para os casos fora da UTI
    if tipo_seir == 'mod':
        #t_UTI: dias entre hospitalização e necessidade de UTI
        parametros.add('t_UTI', value = 5.66,   min = 4.98,     max = 6.34,     vary = False )
        #t_rec: dias entre entrada na UTI e recuperação/óbito
        parametros.add('t_rec', value = 16,     min = 1,      max = 20-5.66,  vary = True  ) #lembrar que o t_rec do std é o t_UTI daqui
    elif tipo_seir == 'std':
        #t_rec: dias entre confirmação do caso e recuperação/óbito
        parametros.add('t_rec', value = 20,  min = 0.3, max = 20, vary = False  )
    
    #
    #parâmetros para consideração de medidas de mitigação e distanciamento social
    #
    
    #Consideração de lockdown com porcentagem 'tau_set' para dia 'dia_set' -> (dia, tau) \\ valores da literatura
    #NOTA: o conjunto de variáveis realoca uma porcentagem 'tau_set' da população suscetível para uma população em confinamento a partir do dia 'dia_set'
    for i in range(0, len(lockdown[0]), 2):
        parametros.add('dia_set' + str(int(i / 2 + 1)),  value = lockdown[0][i],     vary = False )
        parametros.add('tau_set' + str(int(i / 2 + 1)),  value = lockdown[0][i + 1], vary = False )
            
    #
    #parâmetros para consideração da variação da mortalidade no tempo
    #
    
    #Consideração da variação da mortalidade a partir do dia 'dob ', gerando uma nova mortalidade 'mort_UTI' na simulação 'std'
    #   ou novos valores de 'xi', 'zeta' e 'mort_UTI' na simulação 'mod' a partir da multiplicação pelo fator ['obt ' ^ (1/3)]
    if control[2]:
        for i in range(control[3]):
            #NOTA: mínimo de 7 para evitar mudanças iniciais bruscas e máximo inicial de pelo menos 10 dias antes do final dos dados para evitar viéses nos dias finais [é mudado remotamente no fitting a partir do parâmetro f]
            parametros.add('dob' + str(i + 1),  value = round(dias_fit/(control[3] + 1) * (i + 1)),    min = 20-int(round(parametros['f'].value)),    max = dias_fit - max(14 - validacao, 0),   vary = True )
            if tipo_seir == 'mod':
                parametros.add('obt' + str(i + 1),      value = 0.0001 + 0.01*control[3]-i,     min = 0.001,   max = 4,    vary = True )
            elif tipo_seir == 'std':
                parametros.add('mort_UTI' + str(i + 1), value = 0.04 + 0.01*control[3]-i,       min = 0.001,   max = 0.4,    vary = True )
    
    parametros.add('ifr', value = ifr, vary = False)
    
    if not control[2]:
        parametros.add('alpha', expr = '1 - (mort_UTI/ifr-1)/(mort_UTI/ifr)', min = 0.0000000001, max = 1)
    else:
        parametros.add('alpha', expr = '1 - (mort_UTI' + str(control[3]) + '/ifr-1)/(mort_UTI' + str(control[3]) + '/ifr)', min = 0.000001, max = 1)
    
    #Consideração de distanciamento social a partir do dia 'dia ', gerando um novo 'beta' [e consequentemente novo 'R'] -> (dia, beta)
    if control[0]:
        for i in range(control[1]):
            #NOTA: mínimo de 7 para evitar mudanças iniciais bruscas e máximo inicial de pelo menos 10 dias antes do final dos dados para evitar viéses nos dias finais [é mudado remotamente no fitting a partir do parâmetro f]
            parametros.add('dia' + str(i + 1),  value = round(dias_fit/(control[1] + 1) * (i + 1)),    min = 1,    max = dias_fit - max(14 - validacao, 0),   vary = True )
            parametros.add('beta' + str(i + 1), value = 1 + 0.01*control[1]-i,  min = 0.0001,   max = 3,    vary = True )
            parametros.add('R' + str(i + 1),    expr = 'beta' + str(i+1) + '*(mu * (1 - alpha) * D + alpha * EQ)')
    
    #
    #parâmetros informativos para o fitting
    #
    
    #R calculado para o primeiro beta
    parametros.add('R', expr = 'beta*(mu * (1 - alpha) * D + alpha * EQ)')
    #mortalidade do modo str a partir dos fatores do modo mod
    if tipo_seir == 'mod':
        parametros.add('m_str_mod', expr = 'xi * zeta * mort_UTI')
    
    return(parametros)
