import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from pregunta_2 import *  # Importar funcion de pregunta #2

def L_k(x, k, x_points):
    """
    L_k: Calcula el polinomio base de Lagrange para un índice.

    Parámetros:
        x : Variable simbólica para evaluar el polinomio.
        k : int, Índice del polinomio.
        x_points : Lista o arreglo con los puntos x donde se interpola.

    Retorna:
            Polinomio simbólico L_k(x) que es 1 en x_points[k] y 0 en x_points[j].
    """
    result = 1
    x_k = x_points[k]
    for j, x_j in enumerate(x_points):
        if j != k:
            result *= (x - x_j) / (x_k - x_j)
    return sp.simplify(result)

def lagrange_interpolante(x_points, y_points):
    """
    lagrange_interpolante: Construye el polinomio interpolante de Lagrange simbólico para los puntos dados.

    Parámetros:
        x_points : list o array
        y_points : list o array

    Retorna:
        Polinomio interpolante simbólico p_n(x) y la variable simbólica x
    """
    x = sp.symbols('x')
    p_n = 0
    n = len(x_points)
    for k in range(n):
        p_n += y_points[k] * L_k(x, k, x_points)
    return sp.simplify(p_n), x

# Generar resultado con edo2
h = 1
p = lambda x: -1/x
q = lambda x: (1/(4*x**2)) - 1
r = lambda x: 0
x_num, y_num = edo2(p, q, r, 1, 6, 1, 0, h)

# Interpolación de Lagrange
p_lagrange, x_symbol = lagrange_interpolante(x_num, y_num)
print("Polinomio de Lagrange simbólico:")
print(p_lagrange)

# Convertir polinomio a función numérica
p_lagrange_func = sp.lambdify(x_symbol, p_lagrange, 'numpy')

# Puntos para graficar el polinomio
x_dense = np.linspace(1, 6, 300)
y_lagrange_dense = p_lagrange_func(x_dense)

# Graficar ambas curvas
plt.plot(x_num, y_num, 'o-', label='Solución edo2 (h=1)', color='blue')
plt.plot(x_dense, y_lagrange_dense, '--', label='Polinomio de Lagrange', color='red')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Comparación de edo2 vs. interpolante de Lagrange')
plt.legend()
plt.grid(True)
plt.show()
