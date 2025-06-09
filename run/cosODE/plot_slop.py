import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

def extraer_dt(nombre_archivo):
    """Extrae el valor de dt desde el nombre del archivo, asume formato *_dt_00003125.csv"""
    match = re.search(r"_dt_(\d+)\.csv$", nombre_archivo)
    if match:
        dt_str = match.group(1)
        # Insertamos el punto decimal delante del primer dígito no nulo
        dt_float = float("0." + dt_str.lstrip("0")) if dt_str.lstrip("0") else 0.0
        return dt_float
    return None

def calcular_error_csv(nombre_archivo, saltar=0):
    """Calcula el error L2 entre phi y phi_exact, descartando los primeros `saltar` pasos"""
    tiempos, phi, phi_exact = [], [], []

    with open(nombre_archivo, 'r') as f:
        next(f)  # saltar encabezado
        for linea in f:
            columnas = linea.strip().split(",")
            if len(columnas) < 3:
                continue
            tiempos.append(float(columnas[0]))
            phi.append(float(columnas[1]))
            phi_exact.append(float(columnas[2]))

    # Convertir a numpy y descartar los primeros pasos
    phi = np.array(phi[saltar:])
    phi_exact = np.array(phi_exact[saltar:])

    # Norma 2 del error
    error = np.linalg.norm(phi - phi_exact) / np.sqrt(len(phi))
    return error

def main():
    archivos = [f for f in os.listdir() if f.startswith("phi_vs_time_case_bdf4_dt_") and f.endswith(".csv")]
    dt_list = []
    error_list = []

    for archivo in archivos:
        dt = extraer_dt(archivo)
        if dt is not None:
            error = calcular_error_csv(archivo)
            dt_list.append(dt)
            error_list.append(error)
            print(f"{archivo}: dt = {dt:.8f}, error = {error:.2e}")

    # Ordenar por dt
    dt_list, error_list = zip(*sorted(zip(dt_list, error_list)))

    # Log-log para estimar orden
    log_dt = np.log10(dt_list)
    log_error = np.log10(error_list)
    slope, intercept, r_value, _, _ = linregress(log_dt, log_error)

    print(f"\nOrden estimado del método: {abs(slope):.4f} (R² = {r_value**2:.6f})")

    # Graficar
    plt.figure(figsize=(6,4))
    plt.plot(log_dt, log_error, 'o-', label=f'Orden ≈ {abs(slope):.2f}')
    plt.xlabel('log_10(Delta t)')
    plt.ylabel('log_{10}(error)')
    plt.grid(True)
    plt.legend()
    plt.title("Orden de convergencia temporal")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

