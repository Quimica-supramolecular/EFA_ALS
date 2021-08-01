# EFA-ALS

¿Que hace EFA-ALS?
==================

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
      sea mas coherente con su sistema químico y obtenga el mejor ajuste a los
      datos observados. NOTA: Debe ser un numero mayor a 1 y menor de 6.
      
    - Deberá elegir un nivel de tolerancia para separar el ruido del análisis 
      de EFA y los autovalores importantes. Por default el valor es 0.25 y 
      deberá oprimir el numero cero para trabajar con ese valor, el cual será 
      útil la mayoría de las veces. Sin embargo, habrá ocasiones en el que ese
      valor no sea el mas adecuado y deberá analizar la gráfica de EFA para 
      establecer ese valor de manera mas exacta.
     
    - Por último, se le preguntará si desea ver las gráficas obtenidas de los 
      cálculos en cada paso iterativo. Gran parte del presente software fue 
      realizado para aprovechar las utilidades del entorno de desarrollo Spyder,
      el cual permite ver las gráficas en un panel especial sin interrumpir la
      ejecución del código. Por ende, si está utilizando Spyder puede elegir 
      ver las gráficas de manera iterativa sin problema alguno. Sin embargo, si
      esta utilizando otro entorno de desarrollo sin un panel similar al de 
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
