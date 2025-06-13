import numpy as np

def cuarta_derivada_aproximada(f, x, h=1e-3):
    return (f(x - 2*h) - 4*f(x - h) + 6*f(x) - 4*f(x + h) + f(x + 2*h)) / h**4

def max_cuarta_derivada(f, a, b, muestras=1000):
    print("Calculando la cuarta derivada...")
    print(f"Buscando el máximo de |f^(4)(x)| en [{a}, {b}]...")
    x_vals = np.linspace(a, b, muestras)
    derivadas = [abs(cuarta_derivada_aproximada(f, x)) for x in x_vals]
    max_val = max(derivadas)
    print(f"max|f^(4)(x)| = {max_val:.6f}")
    return max_val

def regla_simpson_compuesta(f, a, b, n):
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = [f(xi) for xi in x]
    suma = y[0] + y[-1]
    suma += 4 * sum(y[i] for i in range(1, n, 2))
    suma += 2 * sum(y[i] for i in range(2, n-1, 2))
    return (h / 3) * suma

def cota_simpson_puntos(f, a, b, tol):
    alpha_max = max_cuarta_derivada(f, a, b)
    # Calcular n real (aún no entero)
    h_real = ((180 * tol) / ((b - a) * alpha_max))**(1/4)
    n_real = (b - a) / h_real
    if n_real % 2 != 0:
        n = int(np.ceil(n_real / 2) * 2)  # Redondear al siguiente par
    else:
        n = int(np.ceil(n_real))
    
    # Calcular la cota final de error
    h = (b - a) / n
    error_cota = ((b - a) * h**4 / 180) * alpha_max

    print("\nRESULTADOS:")
    print(f"n calculado (real): {n_real:.2f}")
    print(f"n ajustado (entero par): {n}")
    print(f"Cota de error con n={n}: {error_cota:.2e}")
    return n

# -------------------- PRUEBA --------------------

f = lambda x: np.exp(x**2 - 10*x + 26)
a, b = 5, 5.5
tol = 1e-8

n = cota_simpson_puntos(f, a, b, tol)
integral_aprox = regla_simpson_compuesta(f, a, b, n)
print(f"\nIntegral aproximada con Simpson y n={n}: {integral_aprox:.10f}")
