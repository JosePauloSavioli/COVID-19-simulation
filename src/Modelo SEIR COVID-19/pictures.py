#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 05:19:31 2020

@author: jotape42p
"""

import math
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from datetime import timedelta

font1 = dict(
        family="Courier New, monospace",
        size=14,
        color="#7f7f7f"
        )

font4 = dict(
        family="Courier New, monospace",
        size=18,
        color="#7f7f7f"
        )

#                                                                             #
#####################################-BOOT-####################################
#                                                                             #

def fig_boot(boot, h, info, name_):
    
    fig = px.histogram(boot, x = h, nbins = 50, width = 1500, height = 900)
    fig.write_image(info['path'] + "/Figuras/Boot/0_" + h + name_ + ".png")

#                                                                             #
##################################-POPULAÇÕES-#################################
#                                                                             #

def grafico_populacao(info, data, sim, n = 0):
    
    if n==2:
        iterator = sim['pop_25']
        info['nome'] += '_25'
    elif n:
        iterator = sim['pop_75']
        info['nome'] += '_75'
    else:
        iterator = sim['pop']
    
    fig = make_subplots(rows = 2, cols = 2, 
                        specs = [[{}, {}],
                                 [{'colspan': 2}, None]],
                        horizontal_spacing = 0.05,
                        vertical_spacing = 0.05,
                        subplot_titles = ('Curvas até notificação', 'Curvas pós-notificação', 'Curvas equivalentes SIR'))
    
    x_pp = np.linspace(0, len(sim['pop']['S']) - 1, len(sim['pop']['S']))
    x_p = [data['data_inicial'] + timedelta(days = int(round(sim['result'].params['f'].value)) + x) for x in x_pp]
    
    dc = dict(S= 'Suscetíveis', 
              C= 'Confinados', 
              E= 'Expostos', 
              A= 'Infectados não notificados {assintomáticos}', 
              I= 'Infectados notificados {sintomáticos}', 
              Q= 'Casos leves {quarentenados}', 
              H_tot= 'Casos severos {hospitalizados}', 
              U= 'Casos críticos {UTI}', 
              R= 'Recuperados', 
              D= 'Óbitos', 
              Novos_Casos= 'Casos novos diários', 
              Acum= 'Casos acumulados', 
              Mortes_novas= 'Óbitos novos diários', 
              I_SIR= 'Infectados equivalentes no SIR',
              Cf='Confirmados {não transmitem mais}',
              Ain='Infectados acumulados')
    
    for key, item in iterator.items():
        
        if key in ('E','A','I', 'Novos_Casos'):
            fig.add_trace(go.Scatter(x = x_p,
                                 y = item,
                                 mode = 'lines',
                                 legendgroup = "1",
                                 name = dc[key]),
                        row = 1, col = 1)    
        elif key in ('Cf', 'H_tot', 'U', 'Q', 'D', 'Mortes_novas'):
            fig.add_trace(go.Scatter(x = x_p,
                                 y = item,
                                 mode = 'lines',
                                 legendgroup = "2",
                                 name = dc[key]),
                        row = 1, col = 2)
        elif key in ('S', 'C', 'R', 'D', 'I_SIR', 'Ain'):
            fig.add_trace(go.Scatter(x = x_p,
                                 y = item,
                                 mode = 'lines',
                                 legendgroup = "3",
                                 name = dc[key]),
                        row = 2, col = 1)

    fig.update_layout(title = f"Populações: {info['nome']}", 
                                            font = font4, 
                                            autosize = False, 
                                            width = 2200, 
                                            height = 1500,
                                            legend = dict(font = font1,
                                                          orientation = 'h',
                                                          y = -0.03))
    
    fig.write_image(info['path'] + "/Figuras/pop_" + info['cod'] + '_' + info['nome'] + ".png")

#                                                                             #
###################################-FITTING-###################################
#                                                                             #

def grafico_simulacao(info, data, sim, sl = ()):
    
    if info.get('boot') is not None:
        info['cod'] = 'boot'
        b = "Boot/"
    else:
        b = ''
    
    fig = make_subplots(rows = 2, cols = 2,
                        horizontal_spacing = 0.05,
                        vertical_spacing = 0.05,
                        subplot_titles = ('Casos acumulados','Óbitos'))
    
    if sl:
        if sl[2] != 'boot':
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sl[0]['Acum'],
                                         fill = None,
                                         showlegend=False,
                                         line_color='rgba(255,105,97,0)'),
                 row = 1, col = 1)
            
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sl[1]['Acum'],
                                         fill = 'tonexty',
                                         showlegend=False,
                                         fillcolor = 'rgba(255,105,97,0.3)',
                                         line_color= 'rgba(255,105,97,0)'),
                 row = 1, col = 1)
                 
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sl[0]['D'],
                                         fill = None,
                                         showlegend=False,
                                         line_color='rgba(255,255,255,0)'),
                 row = 1, col = 2)
                 
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sl[1]['D'],
                                         fill = 'tonexty',
                                         showlegend=False,
                                         fillcolor = 'rgba(0,100,80,0.3)',
                                         line_color='rgba(255,255,255,0)'),
                 row = 1, col = 2)
        
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sl[0]['Novos_Casos'],
                                         fill = None,
                                         showlegend=False,
                                         line_color='rgba(255,105,97,0)'),
                 row = 2, col = 1)
            
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sl[1]['Novos_Casos'],
                                         fill = 'tonexty',
                                         showlegend=False,
                                         fillcolor = 'rgba(255,105,97,0.3)',
                                         line_color= 'rgba(255,105,97,0)'),
                 row = 2, col = 1)
                 
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sl[0]['Mortes_novas'],
                                         fill = None,
                                         showlegend=False,
                                         line_color='rgba(255,255,255,0)'),
                 row = 2, col = 2)
                 
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sl[1]['Mortes_novas'],
                                         fill = 'tonexty',
                                         showlegend=False,
                                         fillcolor = 'rgba(0,100,80,0.3)',
                                         line_color='rgba(255,255,255,0)'),
                 row = 2, col = 2)
        
        elif sl[2] == 'boot':
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_max']['Acum'],
                                         fill = None,
                                         showlegend=False,
                                         line_color='rgba(255,105,97,0)'),
                 row = 1, col = 1)
            
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_min']['Acum'],
                                         fill = 'tonexty',
                                         showlegend=False,
                                         fillcolor = 'rgba(255,105,97,0.3)',
                                         line_color= 'rgba(255,105,97,0)'),
                 row = 1, col = 1)
        
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_25']['Acum'],
                                         fill = None,
                                         showlegend=False,
                                         line_color='rgba(255,105,97,0)'),
                 row = 1, col = 1)
            
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_75']['Acum'],
                                         fill = 'tonexty',
                                         showlegend=False,
                                         fillcolor = 'rgba(255,105,97,0.3)',
                                         line_color= 'rgba(255,105,97,0)'),
                 row = 1, col = 1)
        
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_50']['Acum'],
                                         showlegend=False,
                                         line = dict(dash = 'dot'),
                                         line_color='rgba(255,105,97,1)'),
                 row = 1, col = 1)
        
        
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_max']['D'],
                                         fill = None,
                                         showlegend=False,
                                         line_color='rgba(255,255,255,0)'),
                 row = 1, col = 2)
                 
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_min']['D'],
                                         fill = 'tonexty',
                                         showlegend=False,
                                         fillcolor = 'rgba(0,100,80,0.3)',
                                         line_color='rgba(255,255,255,0)'),
                 row = 1, col = 2)
            
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_25']['D'],
                                         fill = None,
                                         showlegend=False,
                                         line_color='rgba(255,105,97,0)'),
                 row = 1, col = 2)
            
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_75']['D'],
                                         fill = 'tonexty',
                                         showlegend=False,
                                         fillcolor = 'rgba(0,100,80,0.3)',
                                         line_color= 'rgba(255,105,97,0)'),
                 row = 1, col = 2)
        
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_50']['D'],
                                         showlegend=False,
                                         line = dict(dash = 'dot'),
                                         line_color='rgba(0,100,80,1)'),
                 row = 1, col = 2)
        
        
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_max']['Novos_Casos'],
                                         fill = None,
                                         showlegend=False,
                                         line_color='rgba(255,105,97,0)'),
                 row = 2, col = 1)
            
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_min']['Novos_Casos'],
                                         fill = 'tonexty',
                                         showlegend=False,
                                         fillcolor = 'rgba(255,105,97,0.3)',
                                         line_color= 'rgba(255,105,97,0)'),
                 row = 2, col = 1)
        
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_25']['Novos_Casos'],
                                         fill = None,
                                         showlegend=False,
                                         line_color='rgba(255,105,97,0)'),
                 row = 2, col = 1)
            
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_75']['Novos_Casos'],
                                         fill = 'tonexty',
                                         showlegend=False,
                                         fillcolor = 'rgba(255,105,97,0.3)',
                                         line_color= 'rgba(255,105,97,0)'),
                 row = 2, col = 1)
        
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_50']['Novos_Casos'],
                                         showlegend=False,
                                         line = dict(dash = 'dot'),
                                         line_color='rgba(255,105,97,1)'),
                 row = 2, col = 1)
        
        
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_max']['Mortes_novas'],
                                         fill = None,
                                         showlegend=False,
                                         line_color='rgba(255,255,255,0)'),
                 row = 2, col = 2)
                 
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_min']['Mortes_novas'],
                                         fill = 'tonexty',
                                         showlegend=False,
                                         fillcolor = 'rgba(0,100,80,0.3)',
                                         line_color='rgba(255,255,255,0)'),
                 row = 2, col = 2)
            
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_25']['Mortes_novas'],
                                         fill = None,
                                         showlegend=False,
                                         line_color='rgba(255,105,97,0)'),
                 row = 2, col = 2)
            
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_75']['Mortes_novas'],
                                         fill = 'tonexty',
                                         showlegend=False,
                                         fillcolor = 'rgba(0,100,80,0.3)',
                                         line_color= 'rgba(255,105,97,0)'),
                 row = 2, col = 2)
        
            fig.add_trace(go.Scatter(x = data['x'],
                                         y = sim['pop_50']['Mortes_novas'],
                                         showlegend=False,
                                         line = dict(dash = 'dot'),
                                         line_color='rgba(0,100,80,1)'),
                 row = 2, col = 2)
    
    a = len(data['casos_acumulados'])
    fig.add_trace(go.Scatter(x = data['x'],
                                 y = sim['pop']['Acum'],
                                 mode = 'lines',
                                 line_color='rgb(255,105,97)',
                                 legendgroup = "1",
                                 name = 'casos acumulados simulados'),
        row = 1, col = 1)
    
    fig.add_trace(go.Bar(x = data['x'][:a],
                                 y = data['casos_ac_no_validation'],
                                 marker = dict(color = 'rgba(211,211,211, 0.9)'),
                                 legendgroup = "1",
                                 name = 'casos acumulados reais'), 
        row = 1, col = 1)
    fig.add_trace(go.Bar(x = data['x'][a:],
                                 y = data['casos_ac_no_validation'][a:],
                                 marker = dict(color = 'rgba(150,150,150, 0.5)'),
                                 showlegend=False), 
        row = 1, col = 1)
    
    fig.add_trace(go.Scatter(x = data['x'],
                                 y = sim['pop']['D'],
                                 mode = 'lines',
                                 line_color='rgb(0,100,80)',
                                 legendgroup = "2",
                                 name = 'óbitos simulados'), 
        row = 1, col = 2)
    
    fig.add_trace(go.Bar(x = data['x'][:a],
                                 y = data['obitos_no_validation'],
                                 marker = dict(color = 'rgba(211,211,211, 0.9)'),
                                 legendgroup = "2",
                                 name = 'óbitos reais'), 
        row = 1, col = 2)
    fig.add_trace(go.Bar(x = data['x'][a:],
                                 y = data['obitos_no_validation'][a:],
                                 marker = dict(color = 'rgba(150,150,150, 0.5)'),
                                 showlegend=False), 
        row = 1, col = 2)
        
    fig.add_trace(go.Scatter(x = data['x'],
                                 y = sim['pop']['Novos_Casos'],
                                 mode = 'lines',
                                 line_color='rgb(255,105,97)',
                                 legendgroup = "1",
                                 name = 'casos acumulados simulados'),
        row = 2, col = 1)
    
    fig.add_trace(go.Bar(x = data['x'][:a],
                                 y = data['casos_novos_no_validation'],
                                 marker = dict(color = 'rgba(211,211,211, 0.9)'),
                                 legendgroup = "1",
                                 name = 'casos acumulados reais'), 
        row = 2, col = 1)
    fig.add_trace(go.Bar(x = data['x'][a:],
                                 y = data['casos_novos_no_validation'][a:],
                                 marker = dict(color = 'rgba(150,150,150, 0.5)'),
                                 showlegend=False), 
        row = 2, col = 1)
    
    fig.add_trace(go.Scatter(x = data['x'],
                                 y = sim['pop']['Mortes_novas'],
                                 mode = 'lines',
                                 line_color='rgb(0,100,80)',
                                 legendgroup = "2",
                                 name = 'óbitos simulados'), 
        row = 2, col = 2)
    
    fig.add_trace(go.Bar(x = data['x'][:a],
                                 y = data['obitos_novos_no_validation'],
                                 marker = dict(color = 'rgba(211,211,211, 0.9)'),
                                 legendgroup = "2",
                                 name = 'óbitos reais'),
        row = 2, col = 2)
    fig.add_trace(go.Bar(x = data['x'][a:],
                                 y = data['obitos_novos_no_validation'][a:],
                                 marker = dict(color = 'rgba(150,150,150, 0.5)'),
                                 showlegend=False),
        row = 2, col = 2)
        
    fig.update_layout(title = f"Ajuste de Curva: {info['nome']} {'Sensibilidade: ' + sl[2] if sl and sl[2] != 'boot' else ''}", 
                                                  font = font4,
                                                  autosize = False, 
                                                  width = 2200, 
                                                  height = 1800,
                                                  legend = dict(font = font1,
                                                                orientation = 'h',
                                                                y = -0.03))
    
    if sl:
        if sl[2] == 'boot':
            fig.write_image(info['path'] + "/Figuras/fitting_boot_" + info['cod'] + '_' + info['nome'] + '_' + sl[2] + ".png")
        else:
            fig.write_image(info['path'] + "/Figuras/Sensibilidade/fitting_" + info['cod'] + '_' + info['nome'] + '_' + sl[2] + ".png")
    else:
        fig.write_image(info['path'] + "/Figuras/" + b + "fitting_" + info['cod'] + '_' + info['nome'] + ".png")

#                                                                             #
#####################################-LOG-#####################################
#                                                                             #

def plot_vert_shape(fig, value, y_min, y_max, dash, color, r = 1, c = 1):
    
    fig.add_shape(dict(type='line', 
                               x0=value,
                               y0=y_min,
                               x1=value,
                               y1=y_max,
                               line = dict(dash = dash,
                                           color = color)),
            row = r,
            col = c)
    
def plot_log(info, data, sim):
    
    fig = make_subplots(rows = 1, cols = 2, subplot_titles = ('Casos acumulados','Óbitos'))
        
    fig.add_trace(go.Scatter(x = data['x'],
                             y = data['casos_ac_no_validation'],
                             mode = 'lines',
                             legendgroup = "1",
                             name = 'casos acumulados reais'),
            row = 1,
            col = 1)
    
    fig.add_trace(go.Scatter(x = data['x'],
                             y = sim['pop']['Acum'],
                             mode = 'lines',
                             legendgroup = "1",
                             name = 'casos acumulados simulados'),
            row = 1,
            col = 1)
            
    fig.add_trace(go.Scatter(x = data['x'],
                             y = data['obitos_no_validation'],
                             mode = 'lines',
                             legendgroup = "2",
                             name = 'óbitos reais'),
            row = 1,
            col = 2)
    
    fig.add_trace(go.Scatter(x = data['x'],
                             y = sim['pop']['D'],
                             mode = 'lines',
                             legendgroup = "2",
                             name = 'óbitos simulados'),
            row = 1,
            col = 2)
    
    r1 = min(min(data['casos_ac_no_validation']) + 0.1, min(sim['pop']['Acum'][:len(data['x'])]) + 0.1)
    s1 = max(max(data['casos_ac_no_validation']), max(sim['pop']['Acum'][:len(data['x'])])) + 50
    r2 = min(min(data['obitos_no_validation']) + 0.1, min(sim['pop']['D'][:len(data['x'])]) + 0.1)
    s2 = max(max(data['obitos_no_validation']), max(sim['pop']['D'][:len(data['x'])])) + 50
    
    for i, value in enumerate(data['lockdown'][0]):
        if not i%2:
            plot_vert_shape(fig, value, r1, s1, 'dot', 'RoyalBlue')
            plot_vert_shape(fig, value, r2, s2, 'dot', 'RoyalBlue', 1, 2)
    for i, value in enumerate(data['lockdown'][2]):
        if not i%2:
            plot_vert_shape(fig, value, r1, s1, 'dash', 'LightSeaGreen')
            plot_vert_shape(fig, value, r2, s2, 'dash', 'LightSeaGreen', 1, 2)
    for i, value in enumerate(data['lockdown'][3]):
        if not i%2:
            plot_vert_shape(fig, value, r1, s1, 'dashdot', 'MediumPurple')
            plot_vert_shape(fig, value, r2, s2, 'dashdot', 'MediumPurple', 1, 2)
    
    fig.update_yaxes(type = 'log')
    fig.update_layout(title = f"Comportamento logarítmico: {info['nome']}", 
                                                            font = font4, 
                                                            autosize = False, 
                                                            width = 2200, 
                                                            height = 900,
                                                            legend_orientation="h")
    
    fig.write_image(info['path'] + "/Figuras/log_" + info['cod'] + '_' + info['nome'] + ".png")

#                                                                             #
#################################-DIFERENÇAS-##################################
#                                                                             #

def plot_dif(df, info, data):
    
    df['Dif_acum_pc'] = df['Dif_acum_pc'].apply(lambda x: x*100)
    df['Dif_obitos_pc'] = df["Dif_obitos_pc"].apply(lambda x: x*100).fillna(0)
    
    fig = make_subplots(rows = 1, cols = 2, subplot_titles = ('Totais','Percentuais'))
    
    fig.add_trace(
            go.Scatter(x = data['x'],
                       y = df['Dif_acum'],
                       mode = 'lines',
                       legendgroup = "1",
                       name = 'Erro total casos acumulados'), 
                       row = 1, col = 1)
    
    fig.add_trace(
            go.Scatter(x = data['x'],
                       y = df['Dif_obitos'],
                       mode = 'lines',
                       legendgroup = "1",
                       name = 'Erro total óbitos'), 
                       row = 1, col = 1)
            
    fig.add_trace(
            go.Scatter(x = data['x'],
                       y = df['Dif_acum_pc'],
                       mode = 'lines',
                       legendgroup = "2",
                       name = 'Erro percentual casos acumulados'), 
                       row = 1, col = 2)
            
    fig.add_trace(
            go.Scatter(x = data['x'],
                       y = df['Dif_obitos_pc'],
                       mode = 'lines',
                       legendgroup = "2",
                       name = 'Erro percentual óbitos'), 
                       row = 1, col = 2)
    
    for i, value in enumerate(data['lockdown'][0]):
        if not i%2:
            plot_vert_shape(fig, value, min(min(df['Dif_acum']), min(df['Dif_obitos'])) - 50, max(max(df['Dif_acum']), max(df['Dif_obitos'])) + 50, 'dot', 'RoyalBlue')
            plot_vert_shape(fig, value, -100, 100, 'dot', 'RoyalBlue', 1, 2)
    for i, value in enumerate(data['lockdown'][2]):
        if not i%2:
            plot_vert_shape(fig, value, min(min(df['Dif_acum']), min(df['Dif_obitos'])) - 50, max(max(df['Dif_acum']), max(df['Dif_obitos'])) + 50, 'dash', 'LightSeaGreen')
            plot_vert_shape(fig, value, -100, 100, 'dash', 'LightSeaGreen', 1, 2)
    for i, value in enumerate(data['lockdown'][3]):
        if not i%2:
            plot_vert_shape(fig, value, min(min(df['Dif_acum']), min(df['Dif_obitos'])) - 50, max(max(df['Dif_acum']), max(df['Dif_obitos'])) + 50, 'dashdot', 'MediumPurple')
            plot_vert_shape(fig, value, -100, 100, 'dashdot', 'MediumPurple', 1, 2)
    
    fig.update_yaxes(range = (min(min(df['Dif_acum']), min(df['Dif_obitos'])) - 50, max(max(df['Dif_acum']), max(df['Dif_obitos'])) + 50))
    fig.update_yaxes(range = (-100, 100), row = 1, col = 2)
    fig.update_layout(title = f"Erros: {info['nome']}", 
                                        font = font4, 
                                        autosize = False, 
                                        width = 2200, 
                                        height = 900,
                                        legend_orientation="h")
    
    fig.write_image(info['path'] + "/Figuras/diferenca_" + info['cod'] + '_' + info['nome'] + ".png")

#                                                                             #
####################################-T_REC-####################################
#                                                                             #

def log_fit(param, ca, ob, title):

    plt.figure(figsize = (15,9))
    
    if 'bestlog' in title:
        plt.title('d = ' + str(param['d']) + ', o = ' + str(param['m']) + ', g = '+ str(param['g']) + ', f = '+ str(param['f']) + ', tr = '+ str(param['tr']))
        plt.plot(ca - param['g'], label = 'ca')
    elif 'bestac' in title:
        plt.title('d = ' + str(param['d']) + ', o = ' + str(1/math.exp(param['m'])) + ', g = '+ str(1/math.exp(param['g'])) + ', f = '+ str(param['f']) + ', tr = '+ str(param['tr']))
        plt.plot(ca / math.exp(param['g']), label = 'ca')
    
    plt.plot(ob, label = 'ob')
    plt.legend()
    plt.savefig(title)
    plt.close()
    
#                                                                             #
##################################-CENÁRIOS-###################################
#                                                                             #
    
def picture(info, data, sim):

    x_pp = np.linspace(0, len(sim['pop']['S']) - 1, len(sim['pop']['S']))
    x_p = [data['data_inicial'] + timedelta(days = int(round(sim['result'].params['f'].value)) + x) for x in x_pp]
    
    fig = make_subplots(rows = 3, cols = 6, 
                        row_heights = [0.4,0.4,0.2],
                        subplot_titles = ('Novos casos', 'Novos óbitos', 'Casos acumulados', 'Óbitos acumulados', 'Ajuste da curva de casos acumulados', 'Ajuste da curva de óbitos', 'Equivalente SIR'),
                        specs = [[{'colspan': 3}, None, None, {'colspan': 3}, None, None],
                                 [{'colspan': 3}, None, None, {'colspan': 3}, None, None],
                                 [{'colspan': 2}, None, {'colspan': 2}, None, {'colspan': 2}, None]],
                        horizontal_spacing = 0.05,
                        vertical_spacing = 0.05)

    fig.add_trace(go.Scatter(x = x_p,
                             y = sim['pop_no_lock']['Novos_Casos'],
                             fill = 'tozeroy',
                             fillcolor = 'rgba(255,105,97, 0.55)',
                             line_color = 'rgb(255,105,97)',
                             legendgroup = '1',
                             name = "Sem quarentena"),
                             row = 1, col = 1)
    fig.add_trace(go.Scatter(x = x_p,
                             y = sim['pop']['Novos_Casos'],
                             fill = 'tozeroy',
                             fillcolor = 'rgba(253, 253, 150, 0.55)',
                             line_color = 'rgb(253, 253, 150)',
                             legendgroup = '1',
                             name = "Atual"),
                             row = 1, col = 1)
    fig.add_trace(go.Scatter(x = x_p,
                             y = sim['pop_wt_lock']['Novos_Casos'],
                             fill = 'tozeroy',
                             fillcolor = 'rgba(0,100,80, 0.55)',
                             line_color = 'rgb(0,100,80)',
                             legendgroup = '1',
                             name = "Com quarentena aumentada"),
                             row = 1, col = 1)
    
    plot_vert_shape(fig, x_p[len(data['x'])], 0, max(sim['pop']['Novos_Casos']), 'dot', 'RoyalBlue')
    
    
    fig.update_yaxes(range = [0, max(sim['pop']['Novos_Casos'])], row=1, col=1)
    
    fig.add_trace(go.Scatter(x = x_p,
                             y = sim['pop_no_lock']['Acum'],
                             fill = 'tozeroy',
                             fillcolor = 'rgba(255,105,97,0.55)',
                             line_color = 'rgb(255,105,97)',
                             showlegend = False),
                             row = 2, col = 1)
    fig.add_trace(go.Scatter(x = x_p,
                             y = sim['pop']['Acum'],
                             fill = 'tozeroy',
                             fillcolor = 'rgba(253, 253, 150,0.55)',
                             line_color = 'rgb(253, 253, 150)',
                             showlegend = False),
                             row = 2, col = 1)
    fig.add_trace(go.Scatter(x = x_p,
                             y = sim['pop_wt_lock']['Acum'],
                             fill = 'tozeroy',
                             fillcolor = 'rgba(0,100,80,0.55)',
                             line_color = 'rgb(0,100,80)',
                             showlegend = False),
                             row = 2, col = 1)

    plot_vert_shape(fig, x_p[len(data['x'])], 0, max(sim['pop']['Acum']), 'dot', 'RoyalBlue', 2, 1)

    fig.update_yaxes(range = [0, max(sim['pop']['Acum'])], row=2, col=1)
    
    fig.add_trace(go.Scatter(x = x_p,
                             y = sim['pop_no_lock']['Mortes_novas'],
                             fill = 'tozeroy',
                             fillcolor = 'rgba(255,105,97,0.55)',
                             line_color = 'rgb(255,105,97)',
                             showlegend = False),
                             row = 1, col = 4)
    fig.add_trace(go.Scatter(x = x_p,
                             y = sim['pop']['Mortes_novas'],
                             fill = 'tozeroy',
                             fillcolor = 'rgba(253, 253, 150,0.55)',
                             line_color = 'rgb(253, 253, 150)',
                             showlegend = False),
                             row = 1, col = 4)
    fig.add_trace(go.Scatter(x = x_p,
                             y = sim['pop_wt_lock']['Mortes_novas'],
                             fill = 'tozeroy',
                             fillcolor = 'rgba(0,100,80,0.55)',
                             line_color = 'rgb(0,100,80)',
                             showlegend = False),
                             row = 1, col = 4)

    plot_vert_shape(fig, x_p[len(data['x'])], 0, max(sim['pop']['Mortes_novas']), 'dot', 'RoyalBlue', 1, 4)

    fig.update_yaxes(range = [0, max(sim['pop']['Mortes_novas'])], row=1, col=4)
    
    fig.add_trace(go.Scatter(x = x_p,
                             y = sim['pop_no_lock']['D'],
                             fill = 'tozeroy',
                             fillcolor = 'rgba(255,105,97,0.55)',
                             line_color = 'rgb(255,105,97)',
                             showlegend = False),
                             row = 2, col = 4)
    fig.add_trace(go.Scatter(x = x_p,
                             y = sim['pop']['D'],
                             fill = 'tozeroy',
                             fillcolor = 'rgba(253, 253, 150,0.55)',
                             line_color = 'rgb(253, 253, 150)',
                             showlegend = False),
                             row = 2, col = 4)
    fig.add_trace(go.Scatter(x = x_p,
                             y = sim['pop_wt_lock']['D'],
                             fill = 'tozeroy',
                             fillcolor = 'rgba(0,100,80,0.55)',
                             line_color = 'rgb(0,100,80)',
                             showlegend = False),
                             row = 2, col = 4)

    plot_vert_shape(fig, x_p[len(data['x'])], 0, max(sim['pop']['D']), 'dot', 'RoyalBlue', 2, 4)

    fig.update_yaxes(range = [0, max(sim['pop']['D'])], row=2, col=4)
    
    a = len(data['casos_confirmados'])
    fig.add_trace(go.Bar(x = data['x'][:a],
                                 y = data['casos_ac_no_validation'],
                                 marker = dict(color = 'rgba(211,211,211,0.9)'),
                                 legendgroup = "2",
                                 name = 'casos acumulados reais'), 
                            row = 3, col = 1)
    fig.add_trace(go.Bar(x = data['x'][a:],
                                 y = data['casos_ac_no_validation'][a:],
                                 marker = dict(color = 'rgba(150,150,150, 0.5)'),
                                 showlegend=False), 
                            row = 3, col = 1)
    fig.add_trace(go.Scatter(x = data['x'],
                                 y = sim['pop']['Acum'],
                                 mode = 'lines',
                                 line_color='rgb(255,105,97)',
                                 legendgroup = "2",
                                 name = 'casos acumulados simulados'),
                            row = 3, col = 1)
    fig.add_trace(go.Bar(x = data['x'][:a],
                                 y = data['obitos_no_validation'],
                                 marker = dict(color = 'rgba(211,211,211,0.9)'),
                                 legendgroup = "3",
                                 name = 'óbitos reais'), 
                            row = 3, col = 3)
    fig.add_trace(go.Bar(x = data['x'][a:],
                                 y = data['obitos_no_validation'][a:],
                                 marker = dict(color = 'rgba(150,150,150, 0.5)'),
                                 showlegend=False), 
                            row = 3, col = 3)
    fig.add_trace(go.Scatter(x = data['x'],
                                 y = sim['pop']['D'],
                                 mode = 'lines',
                                 line_color='rgb(0,100,80)',
                                 legendgroup = "3",
                                 name = 'óbitos simulados'), 
                            row = 3, col = 3)
    
    fig.add_trace(go.Scatter(x = x_p,
                                 y = sim['pop']['S'],
                                 mode = 'lines',
                                 legendgroup = "4",
                                 name = 'Suscetíveis'), 
                            row = 3, col = 5)
    fig.add_trace(go.Scatter(x = x_p,
                                 y = sim['pop']['I_SIR'],
                                 mode = 'lines',
                                 legendgroup = "4",
                                 name = 'Infectados'), 
                            row = 3, col = 5)
    fig.add_trace(go.Scatter(x = x_p,
                                 y = sim['pop']['R'],
                                 mode = 'lines',
                                 legendgroup = "4",
                                 name = 'Recuperados'), 
                            row = 3, col = 5)
    fig.add_trace(go.Scatter(x = x_p,
                                 y = sim['pop']['C'],
                                 mode = 'lines',
                                 legendgroup = "4",
                                 name = 'Confinados'), 
                            row = 3, col = 5)
    
    fig.update_layout(title = f"Cenários futuros de lockdown para {info['nome']}", 
                                                  font = font4,
                                                  autosize = False, 
                                                  width = 2200, 
                                                  height = 1500,
                                                  legend = dict(font = font1,
                                                                orientation = 'h',
                                                                y = -0.03))
    
    fig.write_image(info['path'] + "/Figuras/principal_" + info['cod'] + '_' + info['nome'] + ".png")
    
    