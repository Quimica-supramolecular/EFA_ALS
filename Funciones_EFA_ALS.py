# -*- coding: utf-8 -*-
"""
Author: Karen L. Ochoa Lara (karen.ochoa@unison.mx)
Author: José O. Juaréz Sanchez (octavio.juarez@unison.mx)
Author: Pedro J. Gomez Vega (pedro.gomez@unison.mx)
Copyright (c) 2021 Universidad de Sonora 

**Licencia MIT**

Por la presente se concede permiso, libre de cargos, a cualquier persona 
que obtenga una copia de este software y de los archivos de documentación 
asociados (el "Software"), a utilizar el Software sin restricción, incluyendo 
sin limitación los derechos a usar, copiar, modificar, fusionar, publicar, 
distribuir, sublicenciar, y/o vender copias del Software, y a permitir a 
las personas a las que se les proporcione el Software a hacer lo mismo, sujeto 
a las siguientes condiciones: 
    
El aviso de copyright anterior y este aviso de permiso se incluirán en todas 
las copias o partes sustanciales del Software. 

EL SOFTWARE SE PROPORCIONA "COMO ESTÁ", SIN GARANTÍA DE NINGÚN TIPO, EXPRESA 
O IMPLÍCITA, INCLUYENDO, PERO NO LIMITADO A GARANTÍAS DE COMERCIALIZACIÓN, 
IDONEIDAD PARA UN PROPÓSITO PARTICULAR E INCUMPLIMIENTO. EN NINGÚN CASO LOS 
AUTORES O PROPIETARIOS DE LOS DERECHOS DE AUTOR SERÁN RESPONSABLES DE NINGUNA 
RECLAMACIÓN, DAÑOS U OTRAS RESPONSABILIDADES, YA SEA EN UNA ACCIÓN DE CONTRATO, 
AGRAVIO O CUALQUIER OTRO MOTIVO, DERIVADAS DE, FUERA DE O EN CONEXIÓN CON EL 
SOFTWARE O SU USO U OTRO TIPO DE ACCIONES EN EL SOFTWARE. 
"""

# =============================================================================
# En este módulo no debería cambiarse nada. Sin embargo, si se presentan dificul-
# -tades con algún ajuste podría variar solo las variables que tienen una nota 
# descriptiva de la variable enmarcada de manera similar a la presente.
# =============================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
from scipy import optimize
from scipy.misc import derivative
import warnings
warnings.filterwarnings("ignore")

class modelo_ajuste:
    def __init__(self, C, A, Y, H, G, K, lb, ub, nw, nc, ns, graph):
        self.C = C
        self.A = A
        self.Y = Y
        self.H = H
        self.G = G
        self.K = K
        self.lb = lb
        self.ub = ub
        self.nw = nw
        self.nc = nc
        self.ns = ns
        self.graph = graph

    def model_1_1(self):
        C = self.C
        A = self.A
        Y = self.Y
        H = self.H
        G = self.G
        K = self.K
        lb = self.lb
        ub = self.ub
        nw = self.nw
        nc = self.nc
        ns = self.ns
        graph = self.graph
        
        ssq  = 1e2
        it = 0
        r_mu = []
        r_ssq = list()
# =============================================================================
#       Si desea proporcionar más ciclos iterativos incremente el número de 
#       esta "it < 150" parte del código.
# =============================================================================
        while it < 150 or ssq < 1e-5:
            it += 1
            G_f = G - (C.T[:,1])
            # minimización de las componentes lineales
            A_0 = []
            for l in range(0, nw):
                
                def eps(e): 
                    f = Y.T[:,l] - (C.T @ e)
                    return np.sum(f * f)
                
                bnds = sp.optimize.Bounds(lb, np.inf, keep_feasible=True)
                aj_0 = optimize.minimize(eps, A[:,l], method = "slsqp", bounds= bnds, tol=1e-30)
                A_0.append(aj_0.x)                              
                
                A_1 = np.array(A_0)
                
            # minimización de las componentes no lineales
            C_0 = []
            for i in range(0, len(C.T)):
                
                """
                Las siguientes funciones son funciones de penalización que se
                aplicarán a la función objetivo. 
                """
                def p1(x): return  G[i] - (x[1] + G_f[i])
            
                def p2(x): return H[i] - (x[0] + x[1])
                        
                def p3(x):
                    return x[1] - (((1/K) + H[i] + G[i]) - ((-(1/K) - H[i] - G[i])**2 - 4*H[i]*G[i])**0.5)/2
        
                def p4(x): return x[1] - (K * x[0] * G_f[i])
                
                def p5(G_f):
                    return G_f[i]**2 - (G_f[i]*(G[i] - H[i] - (1/K))) - (G[i]/K)

                def bounds(C):
                    return [(lb,ub)]*len(C)
                
                bounds = bounds(C)

                pg = 1e30 / (1 * 10**(it-1))

                if pg < 1e10:
                    pg = 1e10
                    
                #Función objetivo
                def C_cal(c): 
                    pen = 1e14 
                    fun = sum(Y[:,i] - (c @ A_1.T))**2 + pg*(p1(c)**2 + p2(c)**2 
                                                        + p3(c)**2 + p4(c)**2
                                                        + p5(G_f)**2)
                    return fun

                C_aj = optimize.minimize(C_cal, C[:,i], method = "bfgs", jac = derivative(C_cal, C[:,i]))
          
                C_0.append(C_aj.x)
    
            C_1 = np.array(C_0)
# =============================================================================
#           En algunas ocasiones será necesario limitar el valor de C_1 a ser
#           mayor que 0 en todo momento, para ello quite el simbolo # de la 
#           siguiente línea.
# =============================================================================
            #C_1[C_1 < 0] = 0

            G_free = G - (C_1[:,1])
            
            k = sum(C_1[:,1]) / sum(C_1[:,0] * G_free)
            
            y_cal = C_1 @ A_1.T
            d = y_cal - Y.T
            ssq_1 = np.sum(d*d)
            mu = abs((ssq - ssq_1)/ ssq)
            lof = ((sum(sum((d**2))) / sum(sum((Y**2)))))**0.5 * 100
            MAE = abs(sum(sum(d)) / nw)
            dif_en_ct = round(max(100 - (np.sum(C_1, 1) * 100 / max(H))), 2)
            

            if graph == True:
                
                plt.plot(G, C_1, ":*")
                plt.xlabel("[G], M", size = "xx-large")
                plt.ylabel("[H_libre], M", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
                
                plt.plot(range(0, nw), A_1, "-")
                plt.xlabel("# canales", size = "xx-large")
                plt.ylabel("Matriz A de propiedad obervada", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
                
                plt.plot(range(0, nw), y_cal.T, "k:")
                plt.plot(range(0, nw), Y, "k-", alpha = 0.5)
                plt.xlabel("# canales", size = "xx-large")
                plt.ylabel("Propiedad observada (Y)", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
            
            def conteo(lista, u_lista):
                cont = 0
                for ele in lista: 
                    if (ele == u_lista): 
                        cont = cont + 1
                return cont
            
            c_ssq = conteo(r_ssq, ssq_1)

            print("="*50)
            print("Suma de cuadrados: ",ssq_1)
            print("ssq se ha repetido: ",c_ssq)
            print("\u03BC: ",mu)
            print("Falta de ajuste (%): ",lof)
            print("Error absoluto medio: ",MAE)
            print("Constante de asociación :",round(k, 2))
            print("diferencia en C total (%): ", dif_en_ct)
            print("="*50)
            
            stats = np.array([ssq_1, mu, lof, MAE, dif_en_ct])

            if mu < 1e-4 and lof < 1 and ssq_1 < 5e-3 and dif_en_ct <= 5:
                print("#"*50)
                print("Se ha logrado la convergencia")
                print("#"*50)
                break
            
            if it == 100:
                print("#"*50)
                print("No se ha logrado la convergencia. Considere cambiar K.")
                print("#"*50)
                break
            
            if ssq > 1e-3 and c_ssq == 10:
                print("#"*50)
                print("No se ha logrado la convergencia, la funcion se ha estancado en un mínimo local. Considere cambiar K, el valor pen de penalización o el parametro tol_k")
                print("#"*50)
                break
            
            K = np.array(k)
            C = C_1.T
            A = np.linalg.pinv(C.T) @ Y.T 
            ssq = ssq_1
            r_ssq.append(ssq)
            r_mu.append(mu)
            
        return C, A, K, stats
        

    def model_1_2(self):
        C = self.C
        A = self.A
        Y = self.Y
        H = self.H
        G = self.G
        K = self.K
        lb = self.lb
        ub = self.ub
        nw = self.nw
        nc = self.nc
        ns = self.ns
        graph = self.graph
    
# =============================================================================
#     La variable tol_k se utiliza para limitar la variación en el valor de las 
#     constantes de asociación, por default se mantiene a un 25%. Antes de 
#     considerar variar este parametro varie el valor de K. 
# =============================================================================
        tol_k = 0.25
        K1_sup = K[0] + (K[0] * tol_k)
        K1_inf = K[0] - (K[0] * tol_k)
        K2_sup = K[1] + (K[1] * tol_k)
        K2_inf = K[1] - (K[1] * tol_k)
        
        ssq  = 1e2
        it = 0
        r_mu = []
        r_ssq = []
# =============================================================================
#       Si desea proporcionar más ciclos iterativos incremente el número de 
#       esta parte del código "it < 150".
# =============================================================================
        while it < 150:
            G_f = G - (C.T[:,1] + (2*C.T[:,2]))
            it += 1
            
            # minimización de las componentes lineales
            A_0 = []
            for l in range(0, nw):
                
                def eps(e): 
                    f = Y.T[:,l] - (C.T @ e)
                    return np.sum(f * f)
                
                bnds = sp.optimize.Bounds(lb, np.inf, keep_feasible=True)
                aj_0 = optimize.minimize(eps, A[:,l], method = "slsqp", tol=1e-30)
                A_0.append(aj_0.x)
                
            A_1 = np.array(A_0) 
            
            #minimización de las componentes no lineales.
            C_0 = []
            for i in range(0, len(C.T)):
                """
                Las siguientes funciones son funciones de penalización que se
                aplicarán a la función objetivo. 
                """
        
                def p1(x): return H[i] - (x[0] + x[1] + x[2]) 
                    
                def p2(x): return G[i] - (G_f[i] + x[1] + (2*x[2])) 
    
                def p3(x): return x[1] - (K[0] * x[0] * G_f[i]) 
                
                def p4(x): return x[2] - (K[1] * x[1] * G_f[i]) 
            
                def p5(x): return lb - x[0] 
    
                def p6(x): return lb - x[1] 
    
                def p7(x): return lb - x[2]
                
                def p11(K): return K[0] - K[1] #no viable
        
                def p8(Gf):
                    p = 1/K[1] + 2 * H[i] - G[i]
                    q = 1/(K[0] * K[1]) + H[i]/K[1] - G[i]/K[1]
                    r = -G[i] / (K[0] * K[1])
                    return 0 - (Gf**3 + p*Gf + q*Gf + r)
                
                def p9(K): return  - K[0]
                
                def p10(K): return - K[1] 
                def p12(x): return x[0] - ub
                def p13(x): return x[1] - ub
                def p14(x): return x[2] - ub
                     
                def bounds(C):
                    return [(lb,ub)]*len(C)
                
                bounds = bounds(C) 

                pg = 1e30 / (1 * 10**(it-1))
 
                if pg < 1e12:
                    pg = 1e12
                
                #Función ojetivo
                def C_cal(c): 
                    pen = 1e14
                    fun = sum(Y[:,i] - (c @ A_1.T))**2 + pg*(p1(c)**2 + p2(c)**2 + p3(c)**2 + p4(c)**2
                              + max(0, p5(c))**2 + max(0, p6(c))**2 + max(0, p7(c))**2
                              + p8(G_f[i])**2 + max(0, p9(K))**2 + max(0, p10(K))**2 
                              + max(0, p12(c))**2 + max(0, p13(c))**2 + max(0, p14(c))**2)
                    return fun     

                C_aj = optimize.minimize(C_cal, C[:,i], method = "bfgs", jac = derivative(C_cal, C[:,i]))
    
                C_0.append(C_aj.x)
                
                
            C_1 = np.array(C_0)
# =============================================================================
#           En algunas ocasiones será necesario limitar el valor de C_1 a ser
#           mayor que 0 en todo momento, para ello quite el simbolo # de la 
#           siguiente línea.
# =============================================================================
            #C_1[C_1 < 0] = 0

            H_free = C_1[:,0]
            HG_free = C_1[:,1]
            HG2_free = C_1[:,2]
            G_free = G - (HG_free + (2*HG2_free))
        
            k1_as = np.array(sum(HG_free) / sum(H_free * G_free))
            k2_as = np.array((sum(HG2_free) / sum(H_free * G_free**2)) / k1_as)

            k1_as[k1_as <= K1_inf] = K1_inf
            k2_as[k2_as <= K2_inf] = K2_inf
            k1_as[k1_as >= K1_sup] = K1_sup
            k2_as[k2_as >= K2_sup] = K2_sup
            
            k = [k1_as, k2_as]
            
            y_cal  = C_1 @ A_1.T
            d = y_cal - Y.T
            ssq_1 = np.sum(d*d)
            mu = abs((ssq - ssq_1)/ ssq)
            lof = ((sum(sum((d**2))) / sum(sum((Y**2)))))**0.5 * 100
            MAE = abs(sum(sum(d)) / nw)
            dif_en_ct = round(max(100 - (np.sum(C_1, 1) * 100 / max(H))), 2)
            
            if graph == True:
                
                plt.plot(G, C_1, ":*")
                plt.xlabel("[G], M", size = "xx-large")
                plt.ylabel("[H_libre], M", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
                
                plt.plot(range(0, nw), A_1, "-")
                plt.xlabel("# canales", size = "xx-large")
                plt.ylabel("Matriz A de propiedad obervada", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
                
                plt.plot(range(0, nw), y_cal.T, "k:")
                plt.plot(range(0, nw), Y, "k-", alpha = 0.5)
                plt.xlabel("# canales", size = "xx-large")
                plt.ylabel("Propiedad observada (Y)", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
                

            def conteo(lista, u_lista):
                cont = 0
                for ele in lista: 
                    if (ele == u_lista): 
                        cont = cont + 1
                return cont
            
            c_ssq = conteo(r_ssq, ssq_1)

            print("="*50)
            print("Suma de cuadrados: ",ssq_1)
            print("ssq se ha repetido: ",c_ssq)
            print("\u03BC: ",mu)
            print("Falta de ajuste (%): ",lof)
            print("Error absoluto medio: ",MAE)
            print("Constante de asociación :",k)
            print("diferencia en C total (%): ",dif_en_ct)
            print("="*50)
            
            stats = np.array([ssq_1, mu, lof, MAE, dif_en_ct])
        
            if mu < 1e-4 and lof < 1 and ssq < 5e-3 and dif_en_ct <= 5:
                print("#"*50)
                print("Se ha logrado la convergencia")
                print("#"*50)
                break

            if it == 100:
                print("#"*50)
                print("No se ha logrado la convergencia, considere cambiar K.")
                print("#"*50)
                break

            if ssq > 1e-3 and c_ssq == 10:
                print("#"*50)
                print("No se ha logrado la convergencia, la funcion se ha estancado en un mínimo local. Considere cambiar K.")
                print("#"*50)
                break
            
            K = k
            C = C_1.T 
            A = np.linalg.pinv(C.T) @ Y.T
            ssq = ssq_1
            r_ssq.append(ssq)
            r_mu.append(mu)
                
        return C, A, K, stats

    def model_1_3(self):
        C = self.C
        A = self.A
        Y = self.Y
        H = self.H
        G = self.G
        K = self.K
        lb = self.lb
        ub = self.ub
        nw = self.nw
        nc = self.nc
        ns = self.ns
        graph = self.graph  

# =============================================================================
#     La variable tol_k se utiliza para limitar la variación en el valor de las 
#     constantes de asociación, por default se mantiene a un 25%. Antes de 
#     considerar variar este parametro varie el valor de K.
# =============================================================================
        tol_k = 0.25
        K1_sup = K[0] + (K[0] * tol_k)
        K1_inf = K[0] - (K[0] * tol_k)
        K2_sup = K[1] + (K[1] * tol_k)
        K2_inf = K[1] - (K[1] * tol_k)
        K3_sup = K[2] + (K[2] * tol_k)
        K3_inf = K[2] - (K[2] * tol_k)

        ssq  = 1e2
        it = 0
        r_mu = []
        r_ssq = []
# =============================================================================
#       Si desea proporcionar más ciclos iterativos incremente el número de 
#       esta parte del código "it < 150".
# =============================================================================
        while it < 150 or ssq < 1e-5:
            G_f = G - (C.T[:,1] + (2*C.T[:,2]) + (3*C.T[:,3]))
            it += 1
            
            # minimización de las componentes lineales
            A_0 = []
            for l in range(0, nw):
                
                def eps(e): 
                    f = Y.T[:,l] - (C.T @ e)
                    return np.sum(f * f)
                
                bnds = sp.optimize.Bounds(lb, np.inf, keep_feasible=True)
                aj_0 = optimize.minimize(eps, A[:,l], method = "slsqp", tol=1e-30)
                A_0.append(aj_0.x)
                
            A_1 = np.array(A_0) 
            
            # minimización de las componentes no lineales
            C_0 = []
            for i in range(0, len(C.T)):
                """
                Las siguientes funciones son funciones de penalización que se
                aplicarán a la función objetivo. 
                """
        
                def p1(x): return H[i] - (x[0] + x[1] + x[2] + x[3]) 
                    
                def p2(x): return G[i] - (G_f[i] + x[1] + (2*x[2]) + (3*x[3])) 

                def p3(x): return x[1] - (K[0] * x[0] * G_f[i]) 
                
                def p4(x): return x[2] - (K[1] * x[1] * G_f[i]) 
                
                def p5(x): return x[3] - (K[2] * x[2] * G_f[i]) 
            
                def p6(x): return lb - x[0]
                def p7(x): return lb - x[1] 
                def p8(x): return lb - x[2]
                def p9(x): return lb - x[3]
                def p10(K): return lb - K[0]
                def p11(K): return lb - K[1] 
                def p12(K): return lb - K[2] 
                def p13(x): return x[0] - ub
                def p14(x): return x[1] - ub
                def p15(x): return x[2] - ub
                def p16(x): return x[3] - ub
                
                def bounds(C): return [(lb,ub)]*len(C)
                bounds = bounds(C)
                
                pg = 1e50 / (1 * 10**(it-1))

                if pg < 1e12:
                    pg = 1e12
                    
                def C_cal(c): 
                    fun = sum(Y[:,i] - (c @ A_1.T))**2 + sum(c)**2 + pg*(p1(c)**2 + p2(c)**2 + p3(c)**2 + p4(c)**2 + p5(c)**2
                                                    + max(0, p6(c))**2 + max(0, p7(c))**2 + max(0, p8(c))**2
                                                    + max(0, p9(c))**2 + max(0, p10(K))**2 + max(0, p11(K))**2
                                                    + max(0, p12(K))**2 + max(0, p13(c))**2 + max(0, p14(c))**2
                                                    + max(0, p15(c))**2 + max(0, p16(c))**2)
                    return fun
                
                C_aj = optimize.minimize(C_cal, C[:,i], method = "bfgs",  jac = derivative(C_cal,C[:,i]))
                C_0.append(C_aj.x)
            
            
            C_1 = np.array(C_0)
# =============================================================================
#           En algunas ocasiones será necesario limitar el valor de C_1 a ser
#           mayor que 0 en todo momento, para ello quite el simbolo # de la 
#           siguiente línea.
# =============================================================================
            #C_1[C_1 < 0] = 0

            H_free = C_1[:,0]
            HG_free = C_1[:,1]
            HG2_free = C_1[:,2]
            HG3_free = C_1[:,3]
            G_free = G - (HG_free + (2*HG2_free) + (3*HG3_free))
        
            k1_as = np.array(sum(HG_free) / sum(H_free * G_free))
            k2_as = np.array((sum(HG2_free) / sum(HG_free * G_free))) 
            k3_as = np.array((sum(HG3_free) / sum(HG2_free * G_free)))
            
            k1_as[k1_as <= K1_inf] = K1_inf
            k2_as[k2_as <= K2_inf] = K2_inf
            k3_as[k3_as <= K3_inf] = K3_inf
            k1_as[k1_as >= K1_sup] = K1_sup
            k2_as[k2_as >= K2_sup] = K2_sup
            k3_as[k3_as >= K3_sup] = K3_sup
            
            k = np.array([k1_as, k2_as, k3_as])
            k[k < 0] = 0

            y_cal = C_1 @ A_1.T
            d = y_cal - Y.T
            ssq_1 = np.sum(d*d)
            mu = abs((ssq - ssq_1)/ ssq)
            lof = ((sum(sum((d**2))) / sum(sum((Y**2)))))**0.5 * 100
            MAE = abs(sum(sum(d)) / nw)
            dif_en_ct = round(max(100 - (np.sum(C_1, 1) * 100 / max(H))), 2)
            
            if graph == True:

                plt.plot(G, C_1, ":*")
                plt.xlabel("[G], M", size = "xx-large")
                plt.ylabel("[H_libre], M", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
                
                plt.plot(range(0, nw), A_1, "-")
                plt.xlabel("# canales", size = "xx-large")
                plt.ylabel("Matriz A de propiedad obervada", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
                
                plt.plot(range(0, nw), y_cal.T, "k:")
                plt.plot(range(0, nw), Y, "k-", alpha = 0.5)
                plt.xlabel("# canales", size = "xx-large")
                plt.ylabel("Propiedad observada (Y)", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
                
            def conteo(lista, u_lista):
                cont = 0
                for ele in lista: 
                    if (ele == u_lista): 
                        cont = cont + 1
                return cont
            
            c_ssq = conteo(r_ssq, ssq_1)

            print("="*50)
            print("Suma de cuadrados: ",ssq_1)
            print("ssq se ha repetido: ",c_ssq)
            print("\u03BC: ",mu)
            print("Falta de ajuste (%): ",lof)
            print("Error absoluto medio: ",MAE)
            print("Constante de asociación :",k)
            print("diferencia en C total (%): ", dif_en_ct)
            print("="*50)
            
            stats = np.array([ssq_1, mu, lof, MAE, dif_en_ct])
            
            if mu < 1e-4 and lof < 1 and ssq_1 < 5e-3 and dif_en_ct <= 5:
                print("#"*50)
                print("Se ha logrado la convergencia")
                print("#"*50)
                break

            if it == 100:
                print("#"*50)
                print("No se ha logrado la convergencia. Considere cambiar K.")
                print("#"*50)
                break

            if ssq > 1e-3 and c_ssq == 10:
                print("#"*50)
                print("No se ha logrado la convergencia, la funcion se ha estancado en un mínimo local, considere cambiar K.")
                print("#"*50)
                break

            K = k
            C = C_1.T 
            A = np.linalg.pinv(C.T) @ Y.T
            ssq = ssq_1
            r_ssq.append(ssq)
            r_mu.append(mu)
        return C, A, K, stats
    
    def model_1_4(self):
        C = self.C
        A = self.A
        Y = self.Y
        H = self.H
        G = self.G
        K = self.K
        lb = self.lb
        ub = self.ub
        nw = self.nw
        nc = self.nc
        ns = self.ns
        graph = self.graph  

# =============================================================================
#     La variable tol_k se utiliza para limitar la variación en el valor de las 
#     constantes de asociación, por default se mantiene a un 25%. Antes de 
#     considerar variar este parametro varie el valor de K. 
# =============================================================================
        tol_k = 0.25
        K1_sup = K[0] + (K[0] * tol_k)
        K1_inf = K[0] - (K[0] * tol_k)
        K2_sup = K[1] + (K[1] * tol_k)
        K2_inf = K[1] - (K[1] * tol_k)
        K3_sup = K[2] + (K[2] * tol_k)
        K3_inf = K[2] - (K[2] * tol_k)
        K4_sup = K[3] + (K[3] * tol_k)
        K4_inf = K[3] - (K[3] * tol_k)

        ssq  = 1e2
        it = 0
        r_mu = []
        r_ssq = []
# =============================================================================
#       Si desea proporcionar más ciclos iterativos incremente el número de 
#       esta parte del código "it < 150".
# =============================================================================
        while it < 150 or ssq < 1e-5:
            G_f = G - (C.T[:,1] + (2*C.T[:,2]) + (3*C.T[:,3]) + (4*C.T[:,4]))
            it += 1
            
            # minimización de las componentes lineales
            A_0 = []
            for l in range(0, nw):
                
                def eps(e): 
                    f = Y.T[:,l] - (C.T @ e)
                    return np.sum(f * f)
                
                bnds = sp.optimize.Bounds(lb, np.inf, keep_feasible=True)
                aj_0 = optimize.minimize(eps, A[:,l], method = "slsqp", tol=1e-30)
                A_0.append(aj_0.x)
                
            A_1 = np.array(A_0) 
            
            # minimización de las componentes no lineales
            C_0 = []
            for i in range(0, len(C.T)):
                """
                Las siguientes funciones son funciones de penalización que se
                aplicarán a la función objetivo. 
                """
        
                def p1(x): return H[i] - (x[0] + x[1] + x[2] + x[3] + x[4]) 
                    
                def p2(x): return G[i] - (G_f[i] + x[1] + (2*x[2]) + (3*x[3]) + (4*x[4])) 

                def p3(x): return x[1] - (K[0] * x[0] * G_f[i]) 
                
                def p4(x): return x[2] - (K[1] * x[1] * G_f[i]) 
                
                def p5(x): return x[3] - (K[2] * x[2] * G_f[i])

                def p6(x): return x[4] - (K[3] * x[3] * G_f[i]) 
            
                def p7(x): return lb - x[0]
                def p8(x): return lb - x[1] 
                def p9(x): return lb - x[2]
                def p10(x): return lb - x[3]
                def p11(x): return lb - x[4]
                def p12(K): return lb - K[0]
                def p13(K): return lb - K[1] 
                def p14(K): return lb - K[2]
                def p15(K): return lb - K[3]  
                def p16(x): return x[0] - ub
                def p17(x): return x[1] - ub
                def p18(x): return x[2] - ub
                def p19(x): return x[3] - ub
                def p20(x): return x[4] - ub    
                
                def bounds(C): return [(lb,ub)]*len(C)
                bounds = bounds(C)
                
                pg = 1e50 / (1 * 10**(it-1))

                if pg < 1e12:
                    pg = 1e12
                    
                def C_cal(c): 
                    fun = sum(Y[:,i] - (c @ A_1.T))**2 + sum(c)**2 + pg*(p1(c)**2 + p2(c)**2 + p3(c)**2 + p4(c)**2 + p5(c)**2
                                                    + p6(c)**2 + max(0, p7(c))**2 + max(0, p8(c))**2 + max(0, p9(c))**2
                                                    + max(0, p10(c))**2 + max(0, p11(c))**2 + max(0, p12(K))**2 + max(0, p13(K))**2
                                                    + max(0, p14(K))**2  + max(0, p15(K))**2 + max(0, p16(c))**2 + max(0, p17(c))**2
                                                    + max(0, p18(c))**2 + max(0, p19(c))**2 + max(0, p20(c))**2)
                    return fun
                
                C_aj = optimize.minimize(C_cal, C[:,i], method = "bfgs",  jac = derivative(C_cal,C[:,i]))
                C_0.append(C_aj.x)
                
            C_1 = np.array(C_0)
# =============================================================================
#           En algunas ocasiones será necesario limitar el valor de C_1 a ser
#           mayor que 0 en todo momento, para ello quite el simbolo # de la 
#           siguiente línea.
# =============================================================================
            C_1[C_1 < 0] = 0
            
            H_free = C_1[:,0]
            HG_free = C_1[:,1]
            HG2_free = C_1[:,2]
            HG3_free = C_1[:,3]
            HG4_free = C_1[:,4]
            G_free = G - (HG_free + (2*HG2_free) + (3*HG3_free) + (4*HG4_free))
        
            k1_as = np.array(sum(HG_free) / sum(H_free * G_free))
            k2_as = np.array((sum(HG2_free) / sum(HG_free * G_free))) 
            k3_as = np.array((sum(HG3_free) / sum(HG2_free * G_free)))
            k4_as = np.array((sum(HG4_free) / sum(HG3_free * G_free)))

            k1_as[k1_as <= K1_inf] = K1_inf
            k2_as[k2_as <= K2_inf] = K2_inf
            k3_as[k3_as <= K3_inf] = K3_inf
            k4_as[k4_as <= K4_inf] = K4_inf
            k1_as[k1_as >= K1_sup] = K1_sup
            k2_as[k2_as >= K2_sup] = K2_sup
            k3_as[k3_as >= K3_sup] = K3_sup
            k4_as[k4_as >= K4_sup] = K4_sup

            k = np.array([k1_as, k2_as, k3_as, k4_as])
            k[k < 0] = 0

            y_cal = C_1 @ A_1.T
            d = y_cal - Y.T
            ssq_1 = np.sum(d*d)
            mu = abs((ssq - ssq_1)/ ssq)
            lof = ((sum(sum((d**2))) / sum(sum((Y**2)))))**0.5 * 100
            MAE = abs(sum(sum(d)) / nw)
            dif_en_ct = round(max(100 - (np.sum(C_1, 1) * 100 / max(H))), 2)
            
            if graph == True:

                plt.plot(G, C_1, ":*")
                plt.xlabel("[G], M", size = "xx-large")
                plt.ylabel("[H_libre], M", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
                
                plt.plot(range(0, nw), A_1, "-")
                plt.xlabel("# canales", size = "xx-large")
                plt.ylabel("Matriz A de propiedad obervada", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
                
                plt.plot(range(0, nw), y_cal.T, "k:")
                plt.plot(range(0, nw), Y, "k-", alpha = 0.5)
                plt.xlabel("# canales", size = "xx-large")
                plt.ylabel("Propiedad observada (Y)", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
                
            def conteo(lista, u_lista):
                cont = 0
                for ele in lista: 
                    if (ele == u_lista): 
                        cont = cont + 1
                return cont
            
            c_ssq = conteo(r_ssq, ssq_1)

            print("="*50)
            print("Suma de cuadrados: ",ssq_1)
            print("ssq se ha repetido: ",c_ssq)
            print("\u03BC: ",mu)
            print("Falta de ajuste (%): ",lof)
            print("Error absoluto medio: ",MAE)
            print("Constante de asociación :",k)
            print("diferencia en C total (%): ",dif_en_ct)
            print("="*50)
            
            stats = np.array([ssq_1, mu, lof, MAE, dif_en_ct])

            if mu < 1e-4 and lof < 1 and ssq_1 < 5e-3 and dif_en_ct <= 5:
                print("#"*50)
                print("Se ha logrado la convergencia")
                print("#"*50)
                break

            if it == 100:
                print("#"*50)
                print("No se ha logrado la convergencia. Considere cambiar K.")
                print("#"*50)
                break

            if ssq > 1e-3 and c_ssq == 10:
                print("#"*50)
                print("No se ha logrado la convergencia, la funcion se ha estancado en un mínimo local, considere cambiar K.")
                print("#"*50)
                break


            K = k
            C = C_1.T 
            A = np.linalg.pinv(C.T) @ Y.T
            ssq = ssq_1
            r_ssq.append(ssq)
            r_mu.append(mu)
        return C, A, K, stats

    
    def model_1_5(self):
        C = self.C
        A = self.A
        Y = self.Y
        H = self.H
        G = self.G
        K = self.K
        lb = self.lb
        ub = self.ub
        nw = self.nw
        nc = self.nc
        ns = self.ns
        graph = self.graph  

# =============================================================================
#     La variable tol_k se utiliza para limitar la variación en el valor de las 
#     constantes de asociación, por default se mantiene a un 25%. Antes de 
#     considerar variar este parametro varie el valor de K. 
# =============================================================================
        tol_k = 0.25
        K1_sup = K[0] + (K[0] * tol_k)
        K1_inf = K[0] - (K[0] * tol_k)
        K2_sup = K[1] + (K[1] * tol_k)
        K2_inf = K[1] - (K[1] * tol_k)
        K3_sup = K[2] + (K[2] * tol_k)
        K3_inf = K[2] - (K[2] * tol_k)
        K4_sup = K[3] + (K[3] * tol_k)
        K4_inf = K[3] - (K[3] * tol_k)
        K5_sup = K[4] + (K[4] * tol_k)
        K5_inf = K[4] - (K[4] * tol_k)

        ssq  = 1e2
        it = 0
        r_mu = []
        r_ssq = []
# =============================================================================
#       Si desea proporcionar más ciclos iterativos incremente el número de 
#       esta parte del código "it < 150".
# =============================================================================
        while it < 150 or ssq < 1e-5:
            G_f = G - (C.T[:,1] + (2*C.T[:,2]) + (3*C.T[:,3]) + (4*C.T[:,4]) + (5*C.T[:,5]))
            it += 1
            
            # minimización de las componentes lineales
            A_0 = []
            for l in range(0, nw):
                
                def eps(e): 
                    f = Y.T[:,l] - (C.T @ e)
                    return np.sum(f * f)
                
                bnds = sp.optimize.Bounds(lb, np.inf, keep_feasible=True)
                aj_0 = optimize.minimize(eps, A[:,l], method = "slsqp", tol=1e-30)
                A_0.append(aj_0.x)
                
            A_1 = np.array(A_0)
            
            # minimización de las componentes no lineales
            C_0 = []
            for i in range(0, len(C.T)):
                """
                Las siguientes funciones son funciones de penalización que se
                aplicarán a la función objetivo. 
                """
        
                def p1(x): return H[i] - (x[0] + x[1] + x[2] + x[3] + x[4] + x[5]) 
                    
                def p2(x): return G[i] - (G_f[i] + x[1] + (2*x[2]) + (3*x[3]) + (4*x[4]) + (5*x[5])) 

                def p3(x): return x[1] - (K[0] * x[0] * G_f[i]) 
                
                def p4(x): return x[2] - (K[1] * x[1] * G_f[i]) 
                
                def p5(x): return x[3] - (K[2] * x[2] * G_f[i])

                def p6(x): return x[4] - (K[3] * x[3] * G_f[i]) 
            
                def p7(x): return x[5] - (K[4] * x[4] * G_f[i])

                def p8(x): return lb - x[0]
                def p9(x): return lb - x[1] 
                def p10(x): return lb - x[2]
                def p11(x): return lb - x[3]
                def p12(x): return lb - x[4]
                def p13(x): return lb - x[5]
                def p14(K): return lb - K[0]
                def p15(K): return lb - K[1] 
                def p16(K): return lb - K[2]
                def p17(K): return lb - K[3]  
                def p18(K): return lb - K[4]
                def p19(x): return x[0] - ub
                def p20(x): return x[1] - ub
                def p21(x): return x[2] - ub
                def p22(x): return x[3] - ub
                def p23(x): return x[4] - ub
                def p24(x): return x[5] - ub    
                
                def bounds(C): return [(lb,ub)]*len(C)
                bounds = bounds(C)
                
                pg = 1e50 / (1 * 10**(it-1))

                if pg < 1e12:
                    pg = 1e12
                    
                #Función objetivo
                def C_cal(c): 
                    fun = sum(Y[:,i] - (c @ A_1.T))**2 + pg*(p1(c)**2 + p2(c)**2 + p3(c)**2 + p4(c)**2 + p5(c)**2
                                                    + p6(c)**2 + p7(c)**2 + max(0, p8(c))**2 + max(0, p9(c))**2 + max(0, p10(c))**2
                                                    + max(0, p11(c))**2 + max(0, p12(c))**2 + max(0, p13(c))**2 + max(0, p14(K))**2 
                                                    + max(0, p15(K))**2 + max(0, p16(K))**2 + max(0, p17(K))**2 + max(0, p18(K))**2 
                                                    + max(0, p19(c))**2 + max(0, p20(c))**2 + max(0, p21(c))**2 + max(0, p22(c))**2 
                                                    + max(0, p23(c))**2 + max(0, p24(c))**2)
                    return fun
                
                C_aj = optimize.minimize(C_cal, C[:,i], method = "bfgs",  jac = derivative(C_cal,C[:,i]))
                C_0.append(C_aj.x)
                
            C_1 = np.array(C_0)
# =============================================================================
#           En algunas ocasiones será necesario limitar el valor de C_1 a ser
#           mayor que 0 en todo momento, para ello quite el simbolo # de la 
#           siguiente línea.
# =============================================================================
            C_1[C_1 < 0] = 0
        
            H_free = C_1[:,0]
            HG_free = C_1[:,1]
            HG2_free = C_1[:,2]
            HG3_free = C_1[:,3]
            HG4_free = C_1[:,4]
            HG5_free = C_1[:,5]
            G_free = G - (HG_free + (2*HG2_free) + (3*HG3_free) + (4*HG4_free) + (5*HG5_free))
        
            k1_as = np.array(sum(HG_free) / sum(H_free * G_free))
            k2_as = np.array((sum(HG2_free) / sum(HG_free * G_free)))
            k3_as = np.array((sum(HG3_free) / sum(HG2_free * G_free)))
            k4_as = np.array((sum(HG4_free) / sum(HG3_free * G_free)))
            k5_as = np.array((sum(HG5_free) / sum(HG4_free * G_free)))

            k1_as[k1_as <= K1_inf] = K1_inf
            k2_as[k2_as <= K2_inf] = K2_inf
            k3_as[k3_as <= K3_inf] = K3_inf
            k4_as[k4_as <= K4_inf] = K4_inf
            k5_as[k5_as <= K5_inf] = K5_inf
            k1_as[k1_as >= K1_sup] = K1_sup
            k2_as[k2_as >= K2_sup] = K2_sup
            k3_as[k3_as >= K3_sup] = K3_sup
            k4_as[k4_as >= K4_sup] = K4_sup
            k5_as[k5_as >= K5_sup] = K5_sup

            k = np.array([k1_as, k2_as, k3_as, k4_as, k5_as])
            k[k < 0] = 0
            
            y_cal = C_1 @ A_1.T
            d = y_cal - Y.T
            ssq_1 = np.sum(d*d)
            mu = abs((ssq - ssq_1)/ ssq)
            lof = ((sum(sum((d**2))) / sum(sum((Y**2)))))**0.5 * 100
            MAE = abs(sum(sum(d)) / nw)
            dif_en_ct = round(max(100 - (np.sum(C_1, 1) * 100 / max(H))), 2)
            
            if graph == True:
                
                plt.plot(G, C_1, ":*")
                plt.xlabel("[G], M", size = "xx-large")
                plt.ylabel("[H_libre], M", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
                
                plt.plot(range(0, nw), A_1, "-")
                plt.xlabel("# canales", size = "xx-large")
                plt.ylabel("Matriz A de propiedad obervada", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
                
                plt.plot(range(0, nw), y_cal.T, "k:")
                plt.plot(range(0, nw), Y, "k-", alpha = 0.5)
                plt.xlabel("# canales", size = "xx-large")
                plt.ylabel("Propiedad observada (Y)", size = "xx-large")
                plt.xticks(size = "large")
                plt.yticks(size = "large")
                plt.show()
                
            def conteo(lista, u_lista):
                cont = 0
                for ele in lista: 
                    if (ele == u_lista): 
                        cont = cont + 1
                return cont
            
            c_ssq = conteo(r_ssq, ssq_1)

            print("="*50)
            print("Suma de cuadrados: ",ssq_1)
            print("ssq se ha repetido: ",c_ssq)
            print("\u03BC: ",mu)
            print("Falta de ajuste (%): ",lof)
            print("Error absoluto medio: ",MAE)
            print("Constante de asociación :",k)
            print("diferencia en C total (%): ",dif_en_ct)
            print("="*50)
            
            stats = np.array([ssq_1, mu, lof, MAE, dif_en_ct])
            
            if mu < 1e-4 and lof < 1 and ssq_1 < 5e-3 and dif_en_ct <= 5:
                print("#"*50)
                print("Se ha logrado la convergencia")
                print("#"*50)
                break

            if it == 100:
                print("#"*50)
                print("No se ha logrado la convergencia. Considere cambiar K.")
                print("#"*50)
                break

            if ssq > 1e-3 and c_ssq == 10:
                print("#"*50)
                print("No se ha logrado la convergencia, la funcion se ha estancado en un mínimo local, considere cambiar K.")
                print("#"*50)
                break

            K = k
            C = C_1.T 
            A = np.linalg.pinv(C.T) @ Y.T
            ssq = ssq_1
            r_ssq.append(ssq)
            r_mu.append(mu)
        return C, A, K, stats