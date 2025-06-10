import os
import re

def es_tiempo(d):
    """Determina si un nombre de carpeta representa un tiempo numérico"""
    return re.match(r"^\d+(\.\d+)?$", d)

def encontrar_ultimo_tiempo(directorio):
    """Encuentra la carpeta de último tiempo dentro del directorio del caso"""
    tiempos = [d for d in os.listdir(directorio) if os.path.isdir(os.path.join(directorio, d)) and es_tiempo(d)]
    tiempos_ordenados = sorted(tiempos, key=lambda x: float(x))
    return tiempos_ordenados[-1] if tiempos_ordenados else None

def extraer_T_vector(ruta_caso):
    """Extrae el vector interno del campo escalar T"""
    carpeta_tiempo = encontrar_ultimo_tiempo(ruta_caso)
    if carpeta_tiempo is None:
        print("No se encontraron carpetas de tiempo.")
        return []

    ruta_T = os.path.join(ruta_caso, carpeta_tiempo, "T")
    if not os.path.isfile(ruta_T):
        print(f"No se encontró el archivo T en {carpeta_tiempo}.")
        return []

    with open(ruta_T, 'r') as f:
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

if __name__ == "__main__":
    ruta_caso = "."  # Asume que estás en el directorio del caso
    T_valores = extraer_T_vector(ruta_caso)
    print("Vector T:", T_valores)
