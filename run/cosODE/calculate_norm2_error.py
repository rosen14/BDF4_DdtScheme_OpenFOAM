import csv
import math
import sys

def calcular_norma_error(nombre_archivo, skip=300):
    tiempos = []
    phi_vals = []
    phi_exact_vals = []

    with open(nombre_archivo, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Saltar encabezado
        for i, row in enumerate(reader):
            if i < skip:
                continue  # Saltar primeros pasos de tiempo

            try:
                tiempo = float(row[0])
                phi = float(row[1])
                phi_exact = float(row[2])

                tiempos.append(tiempo)
                phi_vals.append(phi)
                phi_exact_vals.append(phi_exact)
            except ValueError:
                continue  # Saltar líneas con errores

    if not phi_vals:
        print("No hay suficientes datos luego de descartar los primeros pasos.")
        return None

    errores = [phi - exact for phi, exact in zip(phi_vals, phi_exact_vals)]
    norma = math.sqrt(sum(e**2 for e in errores))
    #print(1/(len(errores)-1))
    #print((tiempos[1] - tiempos[0])/(tiempos[-1] - tiempos[0]))
    error = norma*math.sqrt(1/(len(errores)-1))  # Norma L2
    return error

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python calcular_error.py phi_vs_time_casoX.csv")
    else:
        archivo = sys.argv[1]
        norma_error = calcular_norma_error(archivo)
        if norma_error is not None:
            print(f"Norma del error (después de descartar n pasos): {norma_error:.6e}")
