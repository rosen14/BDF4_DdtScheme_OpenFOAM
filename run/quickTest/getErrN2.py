import re
import os
import numpy as np
from getCaseVectorT import extraer_T_vector
from getCaseNx import extraer_divisiones_x
from getExactGaussian import gaussiana_discretizada
import matplotlib.pyplot as plt
# Este script recorre los directorios que comienzan con 'case_bdf4_quick_dt_',
# extrae el vector T del campo escalar T y la cantidad de divisiones en el eje x
# de cada caso, y calcula el error en norma 2 respecto a la solución exacta discretizada.

for caseDir in os.listdir('.'):
    if os.path.isdir(caseDir) and caseDir.startswith('case_bdf4_quick_dt_'):
        ruta_blockMeshDict = os.path.join(caseDir, 'system', 'blockMeshDict')
        nx = extraer_divisiones_x(ruta_blockMeshDict)
        T_s = extraer_T_vector(caseDir)
        # Calcula la solución exacta discretizada
        x, T_e = gaussiana_discretizada(nx, centro=0.6, ancho=0.05)
        # Calcula el error en norma 2
        error_n2 = np.linalg.norm(T_s - T_e) / np.sqrt(nx)
        print(f"Case: {caseDir}, Nx: {nx}, Error N2: {error_n2:.16f}")
        # Graficar los resultados
        '''
        plt.figure(figsize=(10, 6))
        plt.title(f"{caseDir}")
        plt.xlabel('X')
        plt.ylabel('f')
        plt.plot(x, T_e, label='Solución Exacta', color='blue')
        plt.plot(x, T_s, label='Solución Numérica', color='red', linestyle='--')
        plt.legend()
        plt.grid()
        plt.show()
        '''

