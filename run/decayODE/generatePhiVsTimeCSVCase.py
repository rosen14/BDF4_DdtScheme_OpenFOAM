import os
import re

def es_tiempo(nombre):
    return re.match(r"^\d+(\.\d+)?$", nombre)

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
                    datos.append((float(tiempo), valor))
                    break
    return datos

def main():
    base = os.getcwd()
    carpetas = [d for d in os.listdir(base) if os.path.isdir(d)]

    for caso in carpetas:
        datos = extraer_phi_uniform(caso)
        if datos:
            with open(f"phi_vs_time_{caso}.csv", "w") as f:
                f.write("time,phi\n")
                for t, phi in datos:
                    f.write(f"{t},{phi}\n")
            print(f"Guardado: phi_vs_time_{caso}.csv")
        else:
            print(f"No se encontró 'phi' en el caso '{caso}'")

if __name__ == "__main__":
    main()

