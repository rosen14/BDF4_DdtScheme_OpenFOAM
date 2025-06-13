import os
import re
import numpy as np
import pandas as pd
pd.set_option('display.precision', 16)
# Solución exacta
def exact_solution(t, lamb):
    return np.cos(lamb * t)

def is_case_dir(path):
    if 'case_' in os.path.basename(path):
        return True
    return False

def is_float(cadena):
    try:
        float(cadena)
        return True
    except ValueError:
        return False
    

def get_time_dirs(path):
    return [
        d for d in os.listdir(path)
        if os.path.isdir(os.path.join(path, d)) and is_float(d)]
            

def get_case_dirs(path):
    return [
        d for d in os.listdir(path)
        if os.path.isdir(os.path.join(path, d)) and is_case_dir(d)]

def get_internal_field_uniform(ruta_archivo):
    with open(ruta_archivo, "r") as archivo:
        contenido = archivo.read()

    # Buscar la línea con 'internalField   uniform <valor>;'
    match = re.search(r"internalField\s+uniform\s+([-\d\.Ee+]+);", contenido)
    if match:
        return float(match.group(1))
    else:
        raise ValueError("No se encontró un valor 'internalField uniform' válido.")

# Recorrer directorios de casos
data = []
for caseDirName in get_case_dirs('.'):
    # Para cada caso, buscar el directorio del último paso de tiempo
    sortedTimes = sorted([float(d) for d in get_time_dirs(caseDirName) if is_float(d)])
    lastTime = sortedTimes[-1]
    lastValue = get_internal_field_uniform(os.path.join(caseDirName, str(lastTime), 'phi'))
    
    deltaT = sortedTimes[1] - sortedTimes[0]
    e_s = exact_solution(lastTime, 1000.0)
    data.append({
        'caseDirName': caseDirName,
        'deltaT': deltaT,
        'lastTime': lastTime,
        'exact_solution': e_s,  # Usando lambda = 1000.0
        'simulated_solution': lastValue,
        'error': abs(lastValue - e_s)
    })

df = pd.DataFrame(data)
df = df.sort_values(by="deltaT").reset_index(drop=True)
orderP = np.log2(df["error"].iloc[1:].values / df["error"].iloc[:-1].values)

# Crear la columna 'Order P', con NaN en el último valor
df["Order P"] = np.append(orderP, np.nan)

print(df)


