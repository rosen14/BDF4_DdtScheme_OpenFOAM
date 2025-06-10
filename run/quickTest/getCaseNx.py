import re
import os

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

if __name__ == "__main__":
    ruta = os.path.join("case_bdf4_quick_dt_0001_dx_001","system", "blockMeshDict")
    if os.path.exists(ruta):
        nx = extraer_divisiones_x(ruta)
        if nx is not None:
            print(f"Cantidad de divisiones en el eje x: {nx}")
    else:
        print("No se encontró el archivo system/blockMeshDict")
