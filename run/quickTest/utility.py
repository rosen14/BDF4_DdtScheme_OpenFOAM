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

def extraer_divisiones_x(ruta_blockMeshDict):
    with open(ruta_blockMeshDict, 'r') as f:
        contenido = f.read()

    # Busca una línea con un bloque hex (...) (Nx Ny Nz)
    patron = re.search(r'hex\s*\([^\)]+\)\s*\(\s*(\d+)\s+(\d+)\s+(\d+)\s*\)', contenido)

    if patron:
        divisiones_x = int(patron.group(1))
        return divisiones_x
    else:
        print("No se encontró una entrada 'hex' con divisiones en el archivo.")
        return None

def get_deltaT(ruta_archivo):
    with open(ruta_archivo, 'r') as archivo:
        contenido = archivo.read()
    
    # Buscar la línea que contiene: deltaT <número>;
    match = re.search(r'\bdeltaT\s+([-\d\.Ee+]+);', contenido)
    if match:
        return float(match.group(1))
    else:
        raise ValueError("No se encontró un valor 'deltaT' válido.")

def gaussiana_discretizada(N, centro=0.6, ancho=0.05):
    """
    Retorna un vector de tamaño N con los valores de una gaussiana
    centrada en `centro` y con desviación estándar `ancho`, en [0,1].
    
    Parámetros:
    - N: número de puntos (int)
    - centro: posición del centro de la gaussiana (float)
    - ancho: desviación estándar (float)

    Retorna:
    - x: posiciones en [0,1] (np.ndarray)
    - g: valores de la gaussiana (np.ndarray)
    """
    dx = 1.0 / N
    x = np.linspace(dx/2, 1 - dx/2, N)
    g = np.exp(-((x - centro) ** 2) / (2 * ancho ** 2))
    return x, g
    
def get_internal_field(ruta_archivo):

    with open(ruta_archivo, 'r') as f:
        lineas = f.readlines()

    datos = []
    dentro_de_lista = False
    for linea in lineas:
        linea = linea.strip()
        if linea.startswith("internalField") and "nonuniform" in linea:
            dentro_de_lista = True
            continue
        if dentro_de_lista:
            if linea.isdigit():
                continue  # tamaño del campo
            elif linea == "(":
                continue
            elif linea == ");":
                break
            else:
                try:
                    datos.append(float(linea))
                except ValueError:
                    pass
    return datos