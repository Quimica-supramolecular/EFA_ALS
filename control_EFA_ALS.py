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

"""
======================
¿Qué hace EFA-ALS?
======================

El siguiente programa sirve para el ajuste de datos obtenidos por las técnicas 
de UV-Vis, fluorescencia, IR o cualquier otra que siga la relación que se obser-
-va en la ecuación 1. 

                         Y = CA^T + E            Ec.1

Donde Y es la matriz de la propiedad observada (absorbancia, fluorescencia,etc.), 
C es la matriz de concentraciones de las especies, A^T es la matriz transpuesta
de las absortividades molares y E es el error asociado a la medición. Para poder
utilizar la ecuación 1 es indispensable tener un estimado inicial del perfil de
concentraciones o de la absortividad molar de las especies en solución. En este
caso, para este programa se ha eligido obtener un estimado inicial del perfil 
de concentraciones utilizando la metodología de evolving factor análysis (EFA)
de manera similar a lo descrito por Gampp et al., 1985[1]. Una vez obtenido el 
perfil de concentraciones iniciales se estima un modelo químico coherente con el
número de especies en solución y se utiliza el método de minimización de 
alternating least squares(ALS) aplicando funciones de penalización de manera 
similar a lo reportado por Gemperline y Cash en 2003 [2].

Cabe señalar que el cáculo de equilibrios químicos en este programa esta 
limitado a modelos receptor-huesped desde 1:1 a 1:5.

Referencias 

1.- Harald Gampp, Marcel Maeder, Charles J. Meyer and Andreas D. Zuberbühler.
    Talanta, Vol. 32, No. 12, pp. 1133-1139, 1985.
2.- Paul J. Gemperline and Eric Cash. Anal. Chem. 2003, 75, 4236-4243.

======================
Como usar el programa.
======================

Este es el módulo a ejecutar. Primeramente, se le solicitarán algunos 
datos tales como:
    
    - Nombre de archivo sin usar la extensión. xlsx, el programa ya lo añade 
      automáticamente.
      
    - El nombre de las hojas de donde se encuentran los datos del experimento.
      Primero se solicitará el nombre de la hoja que contiene los espectros 
      de las titulaciones de UV-Vis, Fluorescencia, IR, DC, etc. Los datos de
      los espectros deberán estar acomodados en el orden de obtención 
      experimental (libre, 1era adición, 2da adición,…) y en columnas con las 
      longitudes de onda como filas.
      
    - Se solicitará el nombre de la hoja donde se encuentran las concentracio-
      -nes con las que se realizó la titulación. Estas concentraciones deben 
      estar acomodadas en columnas en el siguiente orden. 
          -- Primera columna: concentración total del hospedero o receptor
          -- Segunda columna: concentración del huésped.
          -- En las filas se tendrán las mezclas receptor-huésped utilizadas 
             en cada paso de la titulación. 
          
    - Se le preguntará cuantos autovalores desea incluir en el cálculo. Este 
      paso sirve para reducir el tamaño de los datos a analizar, eliminar ruido 
      de los espectros, define el número de especies en solución para el cálcu-
      -lo y al mismo tiempo define el modelo que se utilizará. Si no tiene una
      estimación del modelo de equilibrio químico a utilizar se recomienda hacer
      dos cálculos, el primero utilizando los autovalores claramente separados
      del resto. El segundo cálculo realícelo sumando un autovalor, elija el que
      sea más coherente con su sistema químico y obtenga el mejor ajuste a los
      datos observados. NOTA: Debe ser un número mayor a 1 y menor de 6.
      
    - Deberá elegir un nivel de tolerancia para separar el ruido del análisis 
      de EFA y los autovalores importantes. Por default el valor es 0.25 y 
      deberá oprimir el numero cero para trabajar con ese valor, el cual será 
      útil la mayoría de las veces. Sin embargo, habrá ocasiones en el que ese
      valor no sea el más adecuado y deberá analizar la gráfica de EFA para 
      establecer ese valor de manera más exacta.
     
    - Por último, se le preguntará si desea ver las gráficas obtenidas de los 
      cálculos en cada paso iterativo. Gran parte del presente software fue 
      realizado para aprovechar las utilidades del entorno de desarrollo Spyder,
      el cual permite ver las gráficas en un panel especial sin interrumpir la
      ejecución del código. Por ende, si está utilizando Spyder puede elegir 
      ver las gráficas de manera iterativa sin problema alguno. Sin embargo, si
      está utilizando otro entorno de desarrollo sin un panel similar al de 
      Spyder, considere que la ejecución del programa se detendrá cuando cada 
      gráfica aparezca y deberá cerrarla para que la ejecución continúe, para 
      estos casos cuando se le pregunte si desea ver las gráficas oprima Enter 
      y que de este modo la ejecución de los cálculos sea más fluida.
     
    - El archivo de salida será un archivo xlsx y tendrá el mismo nombre que el 
      archivo de entrada con el añadido de "_salida". En este archivo se 
      incluirán datos como los espectros calculados, las concentraciones de las
      especies en el equilibrio, la o las constantes de asociación, las absor-
      -tividades molares de las especies en solución y los estadísticos del 
      ajuste realizado. 
     
    - Por otra parte, considere que para ejecutar el programa deberá tener ins-
      -talado los módulos de numpy, scipy, pandas y matplotlib. 
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# Datos para el analisis 

n_archivo = input(str("Ingrese el nombre del archivo .xlsx sin incluir la extensión: ")) #ingrese el nombre del archivo de freeoffice o excel
datos = n_archivo + ".xlsx"
espectros = input(str("Ingrese el nombre de la hoja con espectros: ")) # ingrese el nombre de la hoja con los espectros
concentraciones = input(str("Ingrese el nombre de la hoja con las concentraciones utilizadas: ")) # ingrese el nombre de la hoja con las concentraciones utilizadas
spec = np.array(pd.read_excel(datos, espectros, header=0, index_col=0))
concentracion = np.array(pd.read_excel(datos,concentraciones, header=0))

C_T = concentracion[:,0:2] 
G = C_T[:,1]
H = C_T[:,0]
nc = len(C_T)
nw = len(spec)

u, s, v = np.linalg.svd(spec, full_matrices=False)

plt.plot(range(0, nc), np.log10(s), "o")
plt.ylabel("log(EV)", size = "xx-large")
plt.xlabel("# de autovalores", size = "xx-large")
plt.xticks(size = "large")
plt.yticks(size = "large")
plt.show()

EV = int(input("¿Cuantos autovalores incluirá en el cálculo?: ", ))

Y = u[:,0:EV] @ np.diag(s[0:EV:]) @ v[0:EV:]

#EFA fijo
 
L = range(1,(nc + 1), 1)
L2 = range(0, nc, 1)

X = []
for i in L:
    uj, sj, vj = np.linalg.svd(spec.T[:i,:], full_matrices=False)
    X.append(sj)

ev_s = pd.DataFrame(X)
ev_s0 = np.array(ev_s)

X2 = []
for i in L2:
    ui, si, vi = np.linalg.svd(spec.T[i:,:], full_matrices=False)
    X2.append(si)

ev_s1 = pd.DataFrame(X2)
ev_s10 = np.array(ev_s1) 

plt.figure()
plt.plot(G, np.log10(ev_s0), "k-o")
plt.plot(G, np.log10(ev_s10), "b:o")
plt.ylabel("log(EV)", size = "xx-large")
plt.xlabel("[G], M", size = "xx-large")
plt.xticks(size = "large")
plt.yticks(size = "large")
plt.show()
    
C1 = (ev_s.iloc[:,0:EV]).fillna(0)
C2 = (ev_s1.iloc[:,0:EV]).fillna(0) 
 
EFA0 = []
for i in range(0, EV):
    EFA1 = np.array([C1.iloc[:,i], C2.iloc[:,-1-i]]) 
    EFA2 = np.min(EFA1, 0)
    EFA0.append(EFA2)
    
EFA = np.array(EFA0) 

plt.plot(G, EFA.T, ":o")
plt.xlabel("[G], M", size = "xx-large")
plt.ylabel("u.a.", size = "xx-large")
plt.xticks(size = "large")
plt.yticks(size = "large")
plt.show()

print("\n Escriba la tolerancia para separar los autovalores del ruido. Para obtener el valor por default escriba 0", end="")
tolerancia = float(input("¿Cuál es el nivel de tolerancia deseada?: ", ))
if tolerancia == 0:
    tolerancia = 0.25

EFA = np.log10(abs(EFA / EFA[EFA != 0].min())) - tolerancia

c_e = EFA / EFA.max()
c_e[c_e < 0] = 0

C = c_e * max(H)

plt.plot(G, C.T, ":o")
plt.ylabel("[H_libre], M", size = "xx-large")
plt.xlabel("[G], M", size = "xx-large")
plt.xticks(size = "large")
plt.yticks(size = "large")
plt.show()

n_K = EV - 1

if n_K == 1:
    k_e = float(input("Indique un valor estimado para la constante de asociación: ",))
else:
    k_e = []
    for i in range(n_K):
        print("K" + str(i+1) + ":", end="")
        i = float(input("Indique un valor estimado para esta constante de asociación: ",))
        k_e.append(i)

K = np.array(k_e)
A = (np.linalg.pinv(C.T) @ Y.T) 

graph = str(input("¿Quiere ver las gráficas en cada paso iterativo? Oprima s para si o enter para no: ", ))
if graph == "s":
    graph = True

ns = len(C.T)
nc = len(C)
nw = len(A.T)
lb = 0
ub = max(H)

# Gráficas

def graficas(G, nw, c_sp, A_cal):
    
    plt.plot(G, c_sp.T, ":o")
    plt.xlabel("[G], M", size = "xx-large")
    plt.ylabel("[H_libre], M", size = "xx-large")
    plt.xticks(size = "large")
    plt.yticks(size = "large")
    plt.show()
    
    plt.plot(range(0, nw), A_cal.T, "-")
    plt.xlabel("# canales", size = "xx-large")
    plt.ylabel("Matriz A de propiedad obervada", size = "xx-large")
    plt.xticks(size = "large")
    plt.yticks(size = "large")
    plt.show()
    
    y_cal = c_sp.T @ A_cal

    plt.plot(range(0, nw), y_cal.T, "k:")
    plt.plot(range(0, nw), Y, "k-", alpha = 0.5)
    plt.xlabel("# canales", size = "xx-large")
    plt.ylabel("Propiedad observada (Y)", size = "xx-large")
    plt.xticks(size = "large")
    plt.yticks(size = "large")
    plt.show()

# Implemetación de ALS. 

from Funciones_EFA_ALS import *

if n_K == 1:
    ajuste_datos = modelo_ajuste(C, A, Y, H, G, K, lb, ub, nw, nc, ns, graph)
    c_sp, A_cal, K_ev, stats = ajuste_datos.model_1_1()

if n_K == 2:
    ajuste_datos = modelo_ajuste(C, A, Y, H, G, K, lb, ub, nw, nc, ns, graph)
    c_sp, A_cal, K_ev, stats = ajuste_datos.model_1_2()
    
if n_K == 3:
    ajuste_datos = modelo_ajuste(C, A, Y, H, G, K, lb, ub, nw, nc, ns, graph)
    c_sp, A_cal, K_ev, stats = ajuste_datos.model_1_3()

if n_K == 4:
    ajuste_datos = modelo_ajuste(C, A, Y, H, G, K, lb, ub, nw, nc, ns, graph)
    c_sp, A_cal, K_ev, stats = ajuste_datos.model_1_4()

if n_K == 5:
    ajuste_datos = modelo_ajuste(C, A, Y, H, G, K, lb, ub, nw, nc, ns, graph)
    c_sp, A_cal, K_ev, stats = ajuste_datos.model_1_5()
    

# Elaboracion de gráficas de salida 
graficas(G, nw, c_sp, A_cal)

# Archivo de salida

Y_cal = c_sp.T @ A_cal
Y_cal = pd.DataFrame(Y_cal)
stats = pd.DataFrame(stats, index= ["Suma de cuadrados", "\u03BC", "Falta de ajuste (%)",\
                                    "Error absoluto medio", "Diferencia en C total (%)"])

if n_K == 1:
    K_ev = pd.DataFrame([K_ev], index = ["K"])
    c_sp = pd.DataFrame(c_sp, index = ["[H_libre]", "[HG]"])
    A_cal = pd.DataFrame(A_cal, index = ["H_libre", "HG"])

if n_K == 2:
    K_ev = pd.DataFrame(K_ev, index = ["K1", "K2"])
    c_sp = pd.DataFrame(c_sp, index = ["[H_libre]", "[HG]", "[HG\u2082]"])
    A_cal = pd.DataFrame(A_cal, index = ["H_libre", "HG", "HG\u2082"])

if n_K == 3:
    K_ev = pd.DataFrame(K_ev, index = ["K1", "K2", "K3"])
    c_sp = pd.DataFrame(c_sp, index = ["[H_libre]", "[HG]", "[HG\u2082]",\
                                       "[HG\u2083]"])
    A_cal = pd.DataFrame(A_cal, index = ["H_libre", "HG", "HG\u2082",\
                                        "HG\u2083"])

if n_K == 4:
    K_ev = pd.DataFrame(K_ev, index = ["K1", "K2", "K3", "K4"])
    c_sp = pd.DataFrame(c_sp, index = ["[H_libre]", "[HG]", "[HG\u2082]",\
                                       "[HG\u2083]", "[HG\u2084]"])
    A_cal = pd.DataFrame(A_cal, index = ["H_libre", "HG", "HG\u2082",\
                                         "HG\u2083", "HG\u2084"])
    
if n_K == 5:
    K_ev = pd.DataFrame(K_ev, index = ["K1", "K2", "K3", "K4", "K5"])
    c_sp = pd.DataFrame(c_sp, index = ["[H_libre]", "[HG]", "[HG\u2082]",\
                                       "[HG\u2083]", "[HG\u2084]", "[HG\u2085]"])
    A_cal = pd.DataFrame(A_cal, index = ["H_libre", "HG", "HG\u2082",\
                                       "HG\u2083", "HG\u2084", "HG\u2085"])
    
with pd.ExcelWriter(n_archivo + "_salida" + ".xlsx") as writer:
    c_sp.to_excel(writer, sheet_name="c_especies")
    A_cal.to_excel(writer, sheet_name="A_calculada")
    K_ev.to_excel(writer, sheet_name="Constantes de asociación")
    Y_cal.to_excel(writer, sheet_name="Y_cal")
    stats.to_excel(writer, sheet_name="Estadísticos")