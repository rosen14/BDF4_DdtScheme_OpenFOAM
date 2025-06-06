import os
import re
import math

def es_tiempo(nombre):
    return re.match(r"^\d+(\.\d+)?$", nombre)


def solucion_exacta(t, c=1.0):
    return math.exp(-c*t)


def extraer_phi_uniform(carpeta_caso):
    tiempos = sorted(
        [d for d in os.listdir(carpeta_caso) if os.path.isdir(os.path.join(carpeta_caso, d)) and es_tiempo(d)],
        key=lambda x: float(x)
    )

    datos = []

    for tiempo in tiempos:
        archivo_phi = os.path.join(carpeta_caso, tiempo, "phi")
        if not os.path.isfile(archivo_phi):
            continue

        with open(archivo_phi) as f:
            for linea in f:
                if "internalField" in linea and "uniform" in linea:
                    valor = float(linea.split("uniform")[1].replace(";", "").strip())
                    t_val = float(tiempo)
                    phi_ex = solucion_exacta(t_val, 5.0)
                    datos.append((t_val, valor, phi_ex))
                    break
    return datos

def main():
    base = os.getcwd()
    carpetas = [d for d in os.listdir(base) if os.path.isdir(d)]

    for caso in carpetas:
        datos = extraer_phi_uniform(caso)
        if datos:
            with open(f"phi_vs_time_{caso}.csv", "w") as f:
                f.write("time,phi,phi_exact\n")
                for t, phi, phi_ex in datos:
                    f.write(f"{t},{phi},{phi_ex}\n")
            print(f"Guardado: phi_vs_time_{caso}.csv")
        else:
            print(f"No se encontró 'phi' en el caso '{caso}'")

if __name__ == "__main__":
    main()

