import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from pregunta_2 import *  # para edo2 de la pregunta 2

def calcular_coeficientes_trazador_cubico(x, y):
    """
    Calcula los coeficientes de los polinomios cúbicos naturales.
    Parámetros:
        x : array
        y : array

    Retorna:
        a, b, c, d : Coeficientes de los polinomios cúbicos en cada intervalo
        x :          Los nodos originales (sin modificar).
    """
    n = len(x) - 1  # número de intervalos
    h = np.diff(x)

    A = np.zeros(n-1)
    B = np.zeros(n-1)
    C = np.zeros(n-1)
    D = np.zeros(n-1)

    for j in range(1, n):
        i = j - 1  # índice del sistema reducido
        B[i] = 2 * (h[j-1] + h[j])
        D[i] = 6 * ((y[j+1] - y[j]) / h[j] - (y[j] - y[j-1]) / h[j-1])
        if i > 0:
            A[i-1] = h[j-1]
        if i < (n-1) -1:
            C[i] = h[j]

    # Resolver sistema para M internos (funcion originalmente en el archivo de Pregunta2)
    M_internal = thomas_algorithm(A, B, C, D)
    M = np.zeros(n+1)
    M[1:n] = M_internal  # condiciones

    # Calcular coeficientes del trazador
    a = (M[1:] - M[:-1]) / (6 * h)
    b = M[:-1] / 2
    c = (y[1:] - y[:-1]) / h - (2*h*M[:-1] + h*M[1:]) / 6
    d = y[:-1]

    return a, b, c, d, x

def evaluar_trazador_cubico(a, b, c, d, x_nodes, x_eval):
    """
    Evalúa el trazador cúbico definido por coeficientes (a,b,c,d) en los puntos dados.

    Parámetros:
        a, b, c, d :    Coeficientes de los polinomios spline para cada intervalo.
        x_nodes :       Nodos que definen los intervalos de los trazadores.
        x_eval :        Puntos donde se desea evaluar el trazador.

    Retorna:
        y_eval : Valores del trazador evaluado en x_eval.
    """
    n = len(a)
    y_eval = np.zeros_like(x_eval)
    for i, xv in enumerate(x_eval):
        j = np.searchsorted(x_nodes, xv) - 1
        if j < 0:
            j = 0
        elif j >= n:
            j = n - 1
        dx = xv - x_nodes[j]
        y_eval[i] = a[j]*dx**3 + b[j]*dx**2 + c[j]*dx + d[j]
    return y_eval

def polinomios_trazador_cubico_simbolico(a, b, c, d, x_nodes):
    """
    Construye el trazaador cúbico simbólicamente para cada intervalo.

    Parámetros:
        a, b, c, d :    Coeficientes para cada intervalo.
        x_nodes :       Nodos usados para definir los intervalos.

    Retorna:
        Lista con los trazador cubico para cada intervalo.
    """
    x = sp.symbols('x')
    splines = []
    for j in range(len(a)):
        dx = x - x_nodes[j]
        S = a[j]*dx**3 + b[j]*dx**2 + c[j]*dx + d[j]
        splines.append(sp.simplify(S))
    return splines


# Generar solución edo2
h = 1
p = lambda x: -1/x
q = lambda x: (1/(4*x**2)) - 1
r = lambda x: 0
x_num, y_num = edo2(p, q, r, 1, 6, 1, 0, h)

# Calcular coeficientes trazador cúbico
a, b, c, d, x_nodes = calcular_coeficientes_trazador_cubico(x_num, y_num)

# Obtener trazador simbólicos
cubic_sym = polinomios_trazador_cubico_simbolico(a, b, c, d, x_nodes)

print("Polinomios cúbicos spline por intervalo:")
for i, S in enumerate(cubic_sym):
    print(f"Intervalo [{x_nodes[i]}, {x_nodes[i+1]}]:")
    sp.pprint(S)
    print()

# Evaluar el trazador para graficar
x_densidad = np.linspace(1, 6, 300)
y_cubico = evaluar_trazador_cubico(a, b, c, d, x_nodes, x_densidad)

# Graficar
plt.plot(x_num, y_num, 'o-', label='Solución edo2 (h=1)', color='blue')
plt.plot(x_densidad, y_cubico, '--', label='Trazador cúbico', color='green')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Comparación edo2 y trazador cúbico')
plt.legend()
plt.grid(True)
plt.show()
