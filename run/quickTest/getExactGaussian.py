import numpy as np

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

# Ejemplo de uso
if __name__ == "__main__":
    N = 10
    x, g = gaussiana_discretizada(N, centro=0.6, ancho=0.05)
    print("x =", x)
    print("gaussiana =", g)
