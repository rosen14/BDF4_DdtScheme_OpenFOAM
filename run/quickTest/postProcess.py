from utility import (
    get_case_dirs,
    get_time_dirs,
    get_internal_field,
    get_deltaT,
    extraer_divisiones_x,
    gaussiana_discretizada,
    is_float
)
import os
import numpy as np
import pandas as pd
pd.set_option('display.precision', 16)

# Recorrer directorios de casos
data = []
for caseDirName in get_case_dirs('.'):
    # Para cada caso, buscar el directorio del último paso de tiempo
    sortedTimes = sorted([float(d) for d in get_time_dirs(caseDirName) if is_float(d)])
    lastTime = sortedTimes[-1]
    T_s = get_internal_field(os.path.join(caseDirName, str(lastTime), 'T'))
    
    deltaT = get_deltaT(os.path.join(caseDirName, '0/uniform','time'))

    ruta_blockMeshDict = os.path.join(caseDirName, 'system', 'blockMeshDict')
    nx = extraer_divisiones_x(ruta_blockMeshDict)
    x, T_e = gaussiana_discretizada(nx, centro=0.6, ancho=0.05)
    error_n2 = np.linalg.norm(T_s - T_e) / np.sqrt(nx)
    data.append({
        'caseDirName': caseDirName,
        'deltaT': deltaT,
        'deltaX': 1.0 / nx,  # Asumiendo que el dominio es [0,1]
        'error': error_n2
    })

df = pd.DataFrame(data)
df = df.sort_values(by="deltaX").reset_index(drop=True)
orderP = np.log2(df["error"].iloc[1:].values / df["error"].iloc[:-1].values)

# Crear la columna 'Order P', con NaN en el último valor
df["Order P"] = np.append(orderP, np.nan)

print(df)