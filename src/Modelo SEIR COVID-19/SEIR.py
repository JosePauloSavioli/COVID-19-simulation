#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 09:49:47 2020

@author: jotape42p
"""

#Função para confinamento
def tau_func(tau, C, S):

    #Resetar o C e S
    pessoas_lock = C
    C = 0
    S = S + pessoas_lock

    #Calcular C e S novos
    C = S * tau
    S = S * (1 - tau)

    return S, C

#Função para inicialização dos lockdowns
def locks(lockdown, control, parametros, ts):
    
    #lockdown do csv
    for i in range(0,len(lockdown[0]),2):
        lockdown[0][i] = max(int(round(parametros['dia_set'+ str(int(i/2) + 1)].value)) - int(round(parametros['f'].value)),1)
        lockdown[0][i+1] = parametros['tau_set'+ str(int(i/2) + 1)].value
    
    #discretização da mortalidade
    if control[2] == 1:
        lockdown[3] = [0] * (2*max(control[3], 0))
        for i in range(0,len(lockdown[3]),2):
            lockdown[3][i] = int(round(parametros['dob' + str(int(i/2 + 1))].value))
            if ts == 'std':
                lockdown[3][i+1] = parametros['mort_UTI'+ str(int(i/2) + 1)].value
            elif ts == 'mod':
                lockdown[3][i+1] = parametros['obt'+ str(int(i/2) + 1)].value
    
    #discretização de beta
    if control[0] == 1:
        lockdown[2] = [0] * (2*max(control[1], 0))
        for i in range(0,len(lockdown[2]),2):
            lockdown[2][i] = int(round(parametros['dia' + str(int(i/2 + 1))].value))
            lockdown[2][i+1] = parametros['beta' + str(int(i/2 + 1))].value
    
    return lockdown



def SEIR(populacao, dias_simulacao, parametros, lockdown, dia_final_fit, control, validacao):
    
    #
    #consideração dos lockdowns
    #
    
    if control[0] == 1:
        for i in range(control[1]):
            parametros['dia' + str(i + 1)].max = dia_final_fit - max(10 - validacao, 0)
    if control[2] == 1:
        for i in range(control[3]):
            parametros['dob' + str(i + 1)].max = dia_final_fit - max(10 - validacao, 0)   
        
    lockdown = locks(lockdown, control, parametros, 'mod')
    
    #
    #inicialização dos parâmetros
    #
    
    I_inicial   = parametros['I_init'].value    #valor inicial de casos confirmados
    
    timestep    = dias_simulacao - 1            #tempo de simulação
    
    #variaveis para casos acumulados
    beta        = parametros['beta'].value      #fator de transmissão
    alpha       = parametros['alpha'].value     #porcentagem de notificação
    Z           = parametros['Z'].value         #dias de incubação
    Da          = parametros['D'].value         #dias de transmissão para não notificados {assintomáticos} [equivalente a sua recuperação]
    EQ          = parametros['EQ'].value        #dias de transmissão para notificados {sintomáticos}
    mu          = parametros['mu'].value        #fator de redução de transmissão de não notificados {assintomáticos}
    
    #variaveis a partir dos casos hospitalizados
    xi          = parametros['xi'].value        #porcentagem de casos severos {hospitalização}
    t_UTI       = parametros['t_UTI'].value     #dias entre hospitalização e necessidade de UTI 
    zeta        = parametros['zeta'].value      #porcentagem de casos críticos {UTI}
    t_r_d       = parametros['t_rec'].value     #dias para recuperação/morte dos casos críticos
    t_rec       = parametros['t_rec'].value + parametros['t_UTI'].value  #dias para recuperação dos casos severos e leves
    mort_UTI    = parametros['mort_UTI'].value  #fator de mortalidade dos casos críticos

    #
    #inicialização das matrizes
    #
    
    #Populações (fazem parte da população, sua soma é a população total)
    I           = [I_inicial * EQ]                  #notificados {Infectados sintomáticos}
    A           = [I[0] * (1 - alpha) / (alpha)]    #não notificados {infectados Assintomáticos}
    E           = [I[0] * Z / alpha]                #Expostos
    H           = [I_inicial]                       #casos severos {Hospitalizados}
    S           = [populacao-I[0]-A[0]-E[0]-H[0]]   #Suscetíveis
    C           = [0]*(timestep+1)                  #Confinados (em lockdown)
    Q           = [0]                               #casos leves {Quarentenados}
    U           = [0]                               #casos críticos {UTI}
    N           = [0]                               #casos severos não críticos {Not in UTI}
    R           = [0]                               #Recuperados
    D           = [0]                               #mortes acumuladas {Deaths}
    
    #Matrizes adicionais
    Acl         = [I_inicial]                       #casos acumulados
    
    #
    #iteração
    #
    
    for t_count in range(timestep):

        #
        #atribuição dos valores de lockdown
        #
        
        #discretização da mortalidade
        if t_count in lockdown[3]:
            obt = parametros['obt' + str(int(lockdown[3].index(t_count)/2 + 1))].value
            xi = min(xi * obt ** (1/3), 1.0) #raiz cúbica por multiplicar a multiplicação dos três (xi, zeta e mort_UTI)
            zeta = max(min(zeta * obt ** (1/3), parametros['zeta'].max), parametros['zeta'].min)
            mort_UTI = max(min(mort_UTI * obt ** (1/3), parametros['mort_UTI'].max), parametros['mort_UTI'].min)
        
        #discretização de beta
        if t_count in lockdown[2]:
            tau_ = parametros['beta'+ str(int(lockdown[2].index(t_count)/2 + 1))].value
            beta = tau_
        
        #tau do csv (lockdown)
        if t_count in lockdown[0]:
            tau = parametros['tau_set'+ str(int(lockdown[0].index(t_count)/2 + 1))].value
            S[t_count], b = tau_func(tau, C[t_count], S[t_count])
            C[t_count:] = [b for _ in C[t_count:]] 

        #tau futuro (não interfere no fitting)
        elif t_count in lockdown[1] and t_count != 0:
            if lockdown[1][lockdown[1].index(t_count) + 1] != 0.0:
                tau = lockdown[1][lockdown[1].index(t_count) + 1]
            else:
                tau = 0.0
            S[t_count], b = tau_func(tau, C[t_count], S[t_count])
            C[t_count:] = [b for _ in C[t_count:]] 
        
        #
        #Runge-Kutta de quarta ordem
        #
        
        ##Runge Kutta 1
        n_exp = beta * (I[t_count] + mu *  A[t_count]) * S[t_count] / populacao
        
        n_as = (1 - alpha) * E[t_count] / Z
        n_inf = alpha * E[t_count] / Z
        
        n_as_r = A[t_count] / Da
        
        n_q = (1 - xi) * I[t_count] / EQ
        n_h = xi * I[t_count] / EQ
        
        n_h_n = (1 - zeta) * H[t_count] / t_UTI
        n_u = zeta * H[t_count] / t_UTI
    
        n_q_r = Q[t_count] / t_rec
        n_u_r = (1 - mort_UTI) * U[t_count] / t_r_d
        n_h_r = N[t_count] / t_r_d
        #n_nu_r = N[t_count] / t_r_d
        
        n_u_d = mort_UTI * U[t_count] / t_r_d
        
        Sk1     = - n_exp
        Ek1     =   n_exp - n_as - n_inf
        Ak1     =           n_as         - n_as_r
        Ik1     =                  n_inf          - n_q - n_h
        Qk1     =                                   n_q       - n_q_r
        Hk1     =                                         n_h         - n_u - n_h_n
        Nk1     =                                                             n_h_n                          - n_h_r
        Uk1     =                                                       n_u                 - n_u_r - n_u_d
        Rk1     =                          n_as_r             + n_q_r                       + n_u_r          + n_h_r
        Dk1     =                                                                                     n_u_d 
        Acl1    =                                   n_q + n_h
        
        ##Runge Kutta 2
        S1 = S[t_count] + Sk1 / 2
        A1 = A[t_count] + Ak1 / 2
        I1 = I[t_count] + Ik1 / 2
        E1 = E[t_count] + Ek1 / 2
        Q1 = Q[t_count] + Qk1 / 2
        H1 = H[t_count] + Hk1 / 2
        N1 = N[t_count] + Nk1 / 2
        U1 = U[t_count] + Uk1 / 2
        
        n_exp = beta * (I1 + mu *  A1) * S1 / populacao
        n_as = (1 - alpha) * E1 / Z
        n_inf = alpha * E1 / Z
        n_as_r = A1 / Da
        n_q = (1 - xi) * I1 / EQ
        n_h = xi * I1 / EQ
        n_h_n = (1 - zeta) * H1 / t_UTI
        n_u = zeta * H1 / t_UTI
        n_q_r = Q1 / t_rec
        n_u_r = (1 - mort_UTI) * U1 / t_r_d
        n_h_r = N1 / t_r_d
        n_u_d = mort_UTI * U1 / t_r_d
        
        Sk2     = - n_exp
        Ek2     =   n_exp - n_as - n_inf
        Ak2     =           n_as         - n_as_r
        Ik2     =                  n_inf          - n_q - n_h
        Qk2     =                                   n_q       - n_q_r
        Hk2     =                                         n_h         - n_u - n_h_n
        Nk2     =                                                             n_h_n                          - n_h_r
        Uk2     =                                                       n_u                 - n_u_r - n_u_d
        Rk2     =                          n_as_r             + n_q_r                       + n_u_r          + n_h_r
        Dk2     =                                                                                     n_u_d 
        Acl2    =                                   n_q + n_h
        
        ##Runge Kutta 3
        S2 = S[t_count] + Sk2 / 2
        A2 = A[t_count] + Ak2 / 2
        I2 = I[t_count] + Ik2 / 2
        E2 = E[t_count] + Ek2 / 2
        Q2 = Q[t_count] + Qk2 / 2
        H2 = H[t_count] + Hk2 / 2
        N2 = N[t_count] + Nk2 / 2
        U2 = U[t_count] + Uk2 / 2
        
        n_exp = beta * (I2 + mu *  A2) * S2 / populacao
        n_as = (1 - alpha) * E2 / Z
        n_inf = alpha * E2 / Z
        n_as_r = A2 / Da
        n_q = (1 - xi) * I2 / EQ
        n_h = xi * I2 / EQ
        n_h_n = (1 - zeta) * H2 / t_UTI
        n_u = zeta * H2 / t_UTI
        n_q_r = Q2 / t_rec
        n_u_r = (1 - mort_UTI) * U2 / t_r_d
        n_h_r = N2 / t_r_d
        n_u_d = mort_UTI * U2 / t_r_d
        
        Sk3     = - n_exp
        Ek3     =   n_exp - n_as - n_inf
        Ak3     =           n_as         - n_as_r
        Ik3     =                  n_inf          - n_q - n_h
        Qk3     =                                   n_q       - n_q_r
        Hk3     =                                         n_h         - n_u - n_h_n
        Nk3     =                                                             n_h_n                          - n_h_r
        Uk3     =                                                       n_u                 - n_u_r - n_u_d
        Rk3     =                          n_as_r             + n_q_r                       + n_u_r          + n_h_r
        Dk3     =                                                                                     n_u_d 
        Acl3    =                                   n_q + n_h
        
        ##Runge Kutta 4
        S3 = S[t_count] + Sk3
        A3 = A[t_count] + Ak3
        I3 = I[t_count] + Ik3
        E3 = E[t_count] + Ek3
        Q3 = Q[t_count] + Qk3
        H3 = H[t_count] + Hk3
        N3 = N[t_count] + Nk3
        U3 = U[t_count] + Uk3
        
        n_exp = beta * (I3 + mu *  A3) * S3 / populacao
        n_as = (1 - alpha) * E3 / Z
        n_inf = alpha * E3 / Z
        n_as_r = A3 / Da
        n_q = (1 - xi) * I3 / EQ
        n_h = xi * I3 / EQ
        n_h_n = (1 - zeta) * H3 / t_UTI
        n_u = zeta * H3 / t_UTI
        n_q_r = Q3 / t_rec
        n_u_r = (1 - mort_UTI) * U3 / t_r_d
        n_h_r = N3 / t_r_d
        n_u_d = mort_UTI * U3 / t_r_d
        
        Sk4     = - n_exp
        Ek4     =   n_exp - n_as - n_inf
        Ak4     =           n_as         - n_as_r
        Ik4     =                  n_inf          - n_q - n_h
        Qk4     =                                   n_q       - n_q_r
        Hk4     =                                         n_h         - n_u - n_h_n
        Nk4     =                                                             n_h_n                          - n_h_r
        Uk4     =                                                       n_u                 - n_u_r - n_u_d
        Rk4     =                          n_as_r             + n_q_r                       + n_u_r          + n_h_r
        Dk4     =                                                                                     n_u_d 
        Acl4    =                                   n_q + n_h
        
        S.append(S[t_count] + (1/6 * (Sk1 + 2*Sk2 + 2*Sk3 + Sk4)))
        E.append(E[t_count] + (1/6 * (Ek1 + 2*Ek2 + 2*Ek3 + Ek4)))
        A.append(A[t_count] + (1/6 * (Ak1 + 2*Ak2 + 2*Ak3 + Ak4)))
        I.append(I[t_count] + (1/6 * (Ik1 + 2*Ik2 + 2*Ik3 + Ik4)))
        Q.append(Q[t_count] + (1/6 * (Qk1 + 2*Qk2 + 2*Qk3 + Qk4)))
        H.append(H[t_count] + (1/6 * (Hk1 + 2*Hk2 + 2*Hk3 + Hk4)))
        U.append(U[t_count] + (1/6 * (Uk1 + 2*Uk2 + 2*Uk3 + Uk4)))
        N.append(N[t_count] + (1/6 * (Nk1 + 2*Nk2 + 2*Nk3 + Nk4)))
        R.append(R[t_count] + (1/6 * (Rk1 + 2*Rk2 + 2*Rk3 + Rk4)))
        D.append(D[t_count] + (1/6 * (Dk1 + 2*Dk2 + 2*Dk3 + Dk4)))
        
        Acl.append(Acl[t_count] + (1/6 * (Acl1 + 2*Acl2 + 2*Acl3 + Acl4)))

    #
    #retorna as populações [populações adicionais como casos novos e mortes novas são calculadas fora da função do modelo]
    #

    return dict(tuple(zip(['S', 'C', 'E', 'A', 'I', 'Q', 'H', 'U', 'N', 'R', 'D', 'Acum'], [list(map(round, x)) for x in [S,C,E,A,I,Q,H,U,N,R,D,Acl]])))

def Warmup(populacao, parametros, I_inicial):
    
    #
    #inicialização dos parâmetros
    #
    
    beta        = parametros['beta'].value      #fator de transmissão
    alpha       = parametros['alpha'].value     #porcentagem de notificação
    Z           = parametros['Z'].value         #dias de incubação
    Da          = parametros['D'].value         #dias de transmissão para não notificados {assintomáticos} [equivalente a sua recuperação]
    EQ          = parametros['EQ'].value        #dias de transmissão para notificados {sintomáticos}
    mu          = parametros['mu'].value        #fator de redução de transmissão de não notificados {assintomáticos}
    
    t_dea       = parametros['t_rec'].value     #dias para recuperação/morte dos casos confirmados
    mort        = parametros['mort_UTI'].value  #fator de mortalidade

    #
    #inicialização das matrizes
    #
    
    I           = [1]                      #notificados {Infectados sintomáticos}
    A           = [0]        #não notificados {infectados Assintomáticos}
    E           = [0]                    #Expostos
    Cf          = [0]                           #Confirmados, é necessário para separar quem passa a doença de quem não passa a doença
    R           = [0]                                   #Recuperados
    D           = [0]                                   #Mortes acumuladas {Deaths}
    S           = [populacao-I[0]-A[0]-E[0]-Cf[0]-R[0]-D[0]]      #Suscetíveis
    
    #matrizes adicionais
    Acl         = [Cf[0]]                           #casos acumulados
    Ain         = [I[0]+A[0]]                           #infectados totais acumulados

    #
    #iteração
    #
    
    t_count = 0
    
    while Acl[t_count] < I_inicial:

        #
        #Runge-Kutta de quarta ordem
        #
        
        ##Runge Kutta 1
        n_exp  = beta * (I[t_count] + mu *  A[t_count]) * S[t_count] / populacao
        
        n_as   = (1 - alpha) * E[t_count] / Z
        n_inf  = alpha * E[t_count] / Z
        
        n_as_r = A[t_count] / Da
        
        n_in_c = I[t_count] / EQ
        n_cf_r = (1 - mort) * Cf[t_count] / t_dea
        
        n_cf_d = mort * Cf[t_count] / t_dea

        Sk1     = - n_exp
        Ek1     =   n_exp - n_as - n_inf
        Ak1     =           n_as         - n_as_r
        Ik1     =                  n_inf          - n_in_c
        Cfk1    =                                   n_in_c - n_cf_r - n_cf_d
        Rk1     =                          n_as_r          + n_cf_r
        Dk1     =                                                   + n_cf_d
        Acl1    =                                   n_in_c
        Ain1    =           n_as + n_inf
        
        ##Runge Kutta 2
        S1 = S[t_count] + Sk1 / 2
        A1 = A[t_count] + Ak1 / 2
        I1 = I[t_count] + Ik1 / 2
        Cf1= Cf[t_count]+ Cfk1/ 2
        E1 = E[t_count] + Ek1 / 2
        
        n_exp = beta * (I1 + mu *  A1) * S1 / populacao
        n_as = (1 - alpha) * E1 / Z
        n_inf = alpha * E1 / Z
        n_as_r = A1 / Da
        n_in_c = I1 / EQ
        n_cf_r = (1 - mort) * Cf1 / t_dea
        n_cf_d = mort * Cf1 / t_dea

        Sk2     = - n_exp
        Ek2     =   n_exp - n_as - n_inf
        Ak2     =           n_as         - n_as_r
        Ik2     =                  n_inf          - n_in_c
        Cfk2    =                                   n_in_c - n_cf_r - n_cf_d
        Rk2     =                          n_as_r          + n_cf_r
        Dk2     =                                                   + n_cf_d
        Acl2    =                                   n_in_c
        Ain2    =           n_as + n_inf
        
        ##Runge Kutta 3
        S2 = S[t_count] + Sk2 / 2
        A2 = A[t_count] + Ak2 / 2
        I2 = I[t_count] + Ik2 / 2
        Cf2= Cf[t_count]+ Cfk2/ 2
        E2 = E[t_count] + Ek2 / 2
        
        n_exp = beta * (I2 + mu *  A2) * S2 / populacao
        n_as = (1 - alpha) * E2 / Z
        n_inf = alpha * E2 / Z
        n_as_r = A2 / Da
        n_in_c = I2 / EQ
        n_cf_r = (1 - mort) * Cf2 / t_dea
        n_cf_d = mort * Cf2 / t_dea

        Sk3     = - n_exp
        Ek3     =   n_exp - n_as - n_inf
        Ak3     =           n_as         - n_as_r
        Ik3     =                  n_inf          - n_in_c
        Cfk3    =                                   n_in_c - n_cf_r - n_cf_d
        Rk3     =                          n_as_r          + n_cf_r
        Dk3     =                                                   + n_cf_d
        Acl3    =                                   n_in_c
        Ain3    =           n_as + n_inf
        
        ##Runge Kutta 4
        S3 = S[t_count] + Sk3
        A3 = A[t_count] + Ak3
        I3 = I[t_count] + Ik3
        Cf3= Cf[t_count]+ Cfk3
        E3 = E[t_count] + Ek3
        
        n_exp = beta * (I3 + mu *  A3) * S3 / populacao
        n_as = (1 - alpha) * E3 / Z
        n_inf = alpha * E3 / Z
        n_as_r = A3 / Da
        n_in_c = I3 / EQ
        n_cf_r = (1 - mort) * Cf3 / t_dea
        n_cf_d = mort * Cf3 / t_dea

        Sk4     = - n_exp
        Ek4     =   n_exp - n_as - n_inf
        Ak4     =           n_as         - n_as_r
        Ik4     =                  n_inf          - n_in_c
        Cfk4    =                                   n_in_c - n_cf_r - n_cf_d
        Rk4     =                          n_as_r          + n_cf_r
        Dk4     =                                                   + n_cf_d
        Acl4    =                                   n_in_c
        Ain4    =           n_as + n_inf
        
        S.append(S[t_count] + (1/6 * (Sk1 + 2*Sk2 + 2*Sk3 + Sk4)))
        E.append(E[t_count] + (1/6 * (Ek1 + 2*Ek2 + 2*Ek3 + Ek4)))
        A.append(A[t_count] + (1/6 * (Ak1 + 2*Ak2 + 2*Ak3 + Ak4)))
        I.append(I[t_count] + (1/6 * (Ik1 + 2*Ik2 + 2*Ik3 + Ik4)))
        Cf.append(Cf[t_count] + (1/6 * (Cfk1 + 2*Cfk2 + 2*Cfk3 + Cfk4)))
        R.append(R[t_count] + (1/6 * (Rk1 + 2*Rk2 + 2*Rk3 + Rk4)))
        D.append(D[t_count] + (1/6 * (Dk1 + 2*Dk2 + 2*Dk3 + Dk4)))
        Acl.append(Acl[t_count] + (1/6 * (Acl1 + 2*Acl2 + 2*Acl3 + Acl4)))
        Ain.append(Ain[t_count] + (1/6 * (Ain1 + 2*Ain2 + 2*Ain3 + Ain4)))
        
        t_count += 1
        
        if t_count > 1000:
            break
        
        #print(f"{S[t_count+1]:.2f}\n{E[t_count+1]:.2f}\n{I[t_count+1]:.2f}\n{A[t_count+1]:.2f}\n{Cf[t_count+1]:.2f}\n{Acl[t_count+1]:.2f}\n{beta} {alpha} {mort}")
        #input()
        
    #
    #retorna as populações [populações adicionais como casos novos e mortes novas são calculadas fora da função do modelo]
    #
    
    #print(S[t_count], E[t_count], I[t_count], A[t_count], Cf[t_count], R[t_count], D[t_count], Acl[t_count], Ain[t_count])
    
    return S[t_count], E[t_count], I[t_count], A[t_count], Cf[t_count], R[t_count], D[t_count], Acl[t_count], Ain[t_count]


def SEIR_std(populacao, dias_simulacao, parametros, lockdown, dia_final_fit, control, validacao):

    #
    #consideração dos lockdowns
    #
    
    if control[0]:
        for i in range(control[1]):
            parametros['dia' + str(i + 1)].max = dia_final_fit - max(10 - validacao, 0)
    if control[2]:
        for i in range(control[3]):
            parametros['dob' + str(i + 1)].max = dia_final_fit - max(10 - validacao, 0)
    
    lockdown = locks(lockdown, control, parametros, 'std')
    
    #
    #inicialização dos parâmetros
    #
    
    I_inicial   = parametros['I_init'].value    #valor inicial de casos confirmados
    #print(I_inicial, '\n')
    timestep    = dias_simulacao - 1            #tempo de simulação
    
    beta        = parametros['beta'].value      #fator de transmissão
    alpha       = parametros['alpha'].value     #porcentagem de notificação
    Z           = parametros['Z'].value         #dias de incubação
    Da          = parametros['D'].value         #dias de transmissão para não notificados {assintomáticos} [equivalente a sua recuperação]
    EQ          = parametros['EQ'].value        #dias de transmissão para notificados {sintomáticos}
    mu          = parametros['mu'].value        #fator de redução de transmissão de não notificados {assintomáticos}
    
    t_dea       = parametros['t_rec'].value     #dias para recuperação/morte dos casos confirmados
    mort        = parametros['mort_UTI'].value  #fator de mortalidade

    #
    #inicialização das matrizes
    #
    
    I           = [0]                      #notificados {Infectados sintomáticos}
    A           = [0]        #não notificados {infectados Assintomáticos}
    E           = [0]                    #Expostos
    Cf          = [0]                           #Confirmados, é necessário para separar quem passa a doença de quem não passa a doença                        
    C           = [0]*(timestep+1)                      #Confinados (em lockdown)
    R           = [0]                                   #Recuperados
    D           = [0]                                   #Mortes acumuladas {Deaths}
    S           = [0]      #Suscetíveis
    
    #matrizes adicionais
    Acl         = [0]                           #casos acumulados
    Ain         = [0]                           #infectados totais acumulados

    S[0], E[0], I[0], A[0], Cf[0], R[0], D[0], Acl[0], Ain[0] = Warmup(populacao, parametros, I_inicial)
    
    #
    #iteração
    #
    
    for t_count in range(timestep):
        
        #
        #atribuição dos valores de lockdown
        #
        
        #discretização da mortalidade
        if t_count in lockdown[3]:
            obt = parametros['mort_UTI' + str(int(lockdown[3].index(t_count)/2 + 1))].value
            mort = obt
        
        #discretização de beta
        if t_count in lockdown[2]:
            tau_ = parametros['beta'+ str(int(lockdown[2].index(t_count)/2 + 1))].value
            beta = tau_
        
        #tau do csv (lockdown)
        if t_count in lockdown[0]:
            tau = parametros['tau_set'+ str(int(lockdown[0].index(t_count)/2 + 1))].value
            S[t_count], b = tau_func(tau, C[t_count], S[t_count])
            C[t_count:] = [b for _ in C[t_count:]] 
        #tau futuro (não interfere no fitting)
        elif t_count in lockdown[1] and t_count != 0:
            if lockdown[1][lockdown[1].index(t_count) + 1] != 0.0:
                tau = lockdown[1][lockdown[1].index(t_count) + 1]
            else:
                tau = 0.0
            S[t_count], b = tau_func(tau, C[t_count], S[t_count])
            C[t_count:] = [b for _ in C[t_count:]] 
        
        #
        #Runge-Kutta de quarta ordem
        #
        
        ##Runge Kutta 1
        n_exp  = beta * (I[t_count] + mu *  A[t_count]) * S[t_count] / populacao
        
        n_as   = (1 - alpha) * E[t_count] / Z
        n_inf  = alpha * E[t_count] / Z
        
        n_as_r = A[t_count] / Da
        
        n_in_c = I[t_count] / EQ
        n_cf_r = (1 - mort) * Cf[t_count] / t_dea
        
        n_cf_d = mort * Cf[t_count] / t_dea

        Sk1     = - n_exp
        Ek1     =   n_exp - n_as - n_inf
        Ak1     =           n_as         - n_as_r
        Ik1     =                  n_inf          - n_in_c
        Cfk1    =                                   n_in_c - n_cf_r - n_cf_d
        Rk1     =                          n_as_r          + n_cf_r
        Dk1     =                                                   + n_cf_d
        Acl1    =                                   n_in_c
        Ain1    =           n_as + n_inf
        
        ##Runge Kutta 2
        S1 = S[t_count] + Sk1 / 2
        A1 = A[t_count] + Ak1 / 2
        I1 = I[t_count] + Ik1 / 2
        Cf1= Cf[t_count]+ Cfk1/ 2
        E1 = E[t_count] + Ek1 / 2
        
        n_exp = beta * (I1 + mu *  A1) * S1 / populacao
        n_as = (1 - alpha) * E1 / Z
        n_inf = alpha * E1 / Z
        n_as_r = A1 / Da
        n_in_c = I1 / EQ
        n_cf_r = (1 - mort) * Cf1 / t_dea
        n_cf_d = mort * Cf1 / t_dea

        Sk2     = - n_exp
        Ek2     =   n_exp - n_as - n_inf
        Ak2     =           n_as         - n_as_r
        Ik2     =                  n_inf          - n_in_c
        Cfk2    =                                   n_in_c - n_cf_r - n_cf_d
        Rk2     =                          n_as_r          + n_cf_r
        Dk2     =                                                   + n_cf_d
        Acl2    =                                   n_in_c
        Ain2    =           n_as + n_inf
        
        ##Runge Kutta 3
        S2 = S[t_count] + Sk2 / 2
        A2 = A[t_count] + Ak2 / 2
        I2 = I[t_count] + Ik2 / 2
        Cf2= Cf[t_count]+ Cfk2/ 2
        E2 = E[t_count] + Ek2 / 2
        
        n_exp = beta * (I2 + mu *  A2) * S2 / populacao
        n_as = (1 - alpha) * E2 / Z
        n_inf = alpha * E2 / Z
        n_as_r = A2 / Da
        n_in_c = I2 / EQ
        n_cf_r = (1 - mort) * Cf2 / t_dea
        n_cf_d = mort * Cf2 / t_dea

        Sk3     = - n_exp
        Ek3     =   n_exp - n_as - n_inf
        Ak3     =           n_as         - n_as_r
        Ik3     =                  n_inf          - n_in_c
        Cfk3    =                                   n_in_c - n_cf_r - n_cf_d
        Rk3     =                          n_as_r          + n_cf_r
        Dk3     =                                                   + n_cf_d
        Acl3    =                                   n_in_c
        Ain3    =           n_as + n_inf
        
        ##Runge Kutta 4
        S3 = S[t_count] + Sk3
        A3 = A[t_count] + Ak3
        I3 = I[t_count] + Ik3
        Cf3= Cf[t_count]+ Cfk3
        E3 = E[t_count] + Ek3
        
        n_exp = beta * (I3 + mu *  A3) * S3 / populacao
        n_as = (1 - alpha) * E3 / Z
        n_inf = alpha * E3 / Z
        n_as_r = A3 / Da
        n_in_c = I3 / EQ
        n_cf_r = (1 - mort) * Cf3 / t_dea
        n_cf_d = mort * Cf3 / t_dea

        Sk4     = - n_exp
        Ek4     =   n_exp - n_as - n_inf
        Ak4     =           n_as         - n_as_r
        Ik4     =                  n_inf          - n_in_c
        Cfk4    =                                   n_in_c - n_cf_r - n_cf_d
        Rk4     =                          n_as_r          + n_cf_r
        Dk4     =                                                   + n_cf_d
        Acl4    =                                   n_in_c
        Ain4    =           n_as + n_inf
        
        S.append(S[t_count] + (1/6 * (Sk1 + 2*Sk2 + 2*Sk3 + Sk4)))
        E.append(E[t_count] + (1/6 * (Ek1 + 2*Ek2 + 2*Ek3 + Ek4)))
        A.append(A[t_count] + (1/6 * (Ak1 + 2*Ak2 + 2*Ak3 + Ak4)))
        I.append(I[t_count] + (1/6 * (Ik1 + 2*Ik2 + 2*Ik3 + Ik4)))
        Cf.append(Cf[t_count] + (1/6 * (Cfk1 + 2*Cfk2 + 2*Cfk3 + Cfk4)))
        R.append(R[t_count] + (1/6 * (Rk1 + 2*Rk2 + 2*Rk3 + Rk4)))
        D.append(D[t_count] + (1/6 * (Dk1 + 2*Dk2 + 2*Dk3 + Dk4)))
        Acl.append(Acl[t_count] + (1/6 * (Acl1 + 2*Acl2 + 2*Acl3 + Acl4)))
        Ain.append(Ain[t_count] + (1/6 * (Ain1 + 2*Ain2 + 2*Ain3 + Ain4)))
        
        #print(f"{S[t_count+1]:.2f}\n{E[t_count+1]:.2f}\n{I[t_count+1]:.2f}\n{A[t_count+1]:.2f}\n{Cf[t_count+1]:.2f}\n{Acl[t_count+1]:.2f}\n{beta} {alpha} {mort}")
        #input()
        
    #
    #retorna as populações [populações adicionais como casos novos e mortes novas são calculadas fora da função do modelo]
    #

    return dict(tuple(zip(['S', 'C', 'E', 'A', 'I', 'R', 'D', 'Acum', 'Cf', 'Ain'], [list(map(round, x)) for x in [S,C,E,A,I,R,D,Acl,Cf,Ain]])))

