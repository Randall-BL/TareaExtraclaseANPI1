#Ejercicio 2: Diferencias finitas y método de Thomas
import numpy as np
import matplotlib.pyplot as plt

def thomas_algorithm(a, b, c, d):
    """
    thomas_algorithm: Resuelve un sistema de ecuaciones lineales tridiagonal.

    Este algoritmo es una aplicación de la eliminación gaussiana optimizada
    para sistemas tridiagonales.

    Parámetros:
        a :  Diagonal inferior de la matriz tridiagonal (tamaño n-1).
        b :  Diagonal principal de la matriz tridiagonal (tamaño n).
        c :  Diagonal superior de la matriz tridiagonal (tamaño n-1).
        d :  Vector del lado derecho del sistema de ecuaciones (tamaño n).

    Retorna:
        x :  Vector solución del sistema de ecuaciones (tamaño n).
    """
    n = len(d)
    c_prime = np.zeros(n-1)
    d_prime = np.zeros(n)
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]

    for i in range(1, n-1):
        temp = b[i] - a[i-1]*c_prime[i-1]
        c_prime[i] = c[i]/temp
        d_prime[i] = (d[i] - a[i-1]*d_prime[i-1]) / temp

    d_prime[n-1] = (d[n-1] - a[n-2]*d_prime[n-2]) / (b[n-1] - a[n-2]*c_prime[n-2])

    x = np.zeros(n)
    x[-1] = d_prime[-1]
    for i in reversed(range(n-1)):
        x[i] = d_prime[i] - c_prime[i]*x[i+1]
    return x

def edo2(p, q, r, a, b, y0, yn, h):
    """
    edo2: Resuelve una Ecuación Diferencial Ordinaria (EDO) lineal de segundo orden
          con condiciones de frontera utilizando el método de diferencias finitas.

    La EDO tiene la forma: y'' + p(x)y' + q(x)y = r(x)
    Las condiciones de frontera son y(a_interval) = y0, y(b_interval) = yn.

    Parámetros:
        p : function, Coeficiente de y' en la EDO.
        q : function, Coeficiente de y en la EDO.
        r : function, Término independiente de la EDO.
        a_interval : float, Límite inferior del intervalo de la solución.
        b_interval : float, Límite superior del intervalo de la solución.
        y0 : float, Condición de frontera en x = a_interval.
        yn : float, Condición de frontera en x = b_interval.
        h : float, Tamaño del paso para la discretización.

    Retorna:
        x_points : numpy.ndarray, Puntos x donde se calcula la solución.
        y_solution : numpy.ndarray, Solución numérica y(x) en los puntos x_points.
    """
    n = int((b - a) / h)
    x = np.linspace(a, b, n+1)
    A = np.zeros(n-1)
    B = np.zeros(n-1)
    C = np.zeros(n-1)
    D = np.zeros(n-1)

    for i in range(1, n):
        xi = x[i]
        A[i-1] = -1 - h*p(xi)/2
        B[i-1] = 2 + h**2*q(xi)
        C[i-1] =- 1 + h*p(xi)/2
        D[i-1] = -h**2*r(xi)

    D[0] -= A[0]*y0
    D[-1] -= C[-1]*yn

    y_inner = thomas_algorithm(A[1:], B, C[:-1], D)
    y = np.concatenate([[y0], y_inner, [yn]])
    return x, y

# Validación con problema especificado
if __name__ == "__main__":
    # Validación con problema especificado
    p = lambda x: -1/x
    q = lambda x: (1/(4*x**2)) - 1
    r = lambda x: 0

    h_vals = [1, 0.5, 0.2, 0.1, 0.01]
    x_exact = np.linspace(1, 6, 1000)
    y_exact = np.sin(6 - x_exact) / (np.sin(5) * np.sqrt(x_exact))

    for h in h_vals:
        x_num, y_num = edo2(p, q, r, 1, 6, 1, 0, h)
        plt.plot(x_num, y_num, label=f'h={h}')

    plt.plot(x_exact, y_exact, 'k', label='Exacta', linewidth=2)
    plt.legend()
    plt.title('Aproximación por diferencias finitas')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()


