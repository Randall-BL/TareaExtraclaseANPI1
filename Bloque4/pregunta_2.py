import numpy as np
import math
from sympy import symbols, diff, lambdify, exp



# Función para calcular el valor de n necesario para cumplir la tolerancia
def cota_simpson_puntos(f_simbolica, a, b, tol):
    """
    Calcula el número de subintervalos 'n' (par) necesario para que la cota de error
    de la regla compuesta de Simpson sea menor que una tolerancia dada 'tol'.

    Args:
        f_simbolica: La función f(x) definida simbólicamente (usando sympy.symbols).
        a: Límite inferior del intervalo de integración.
        b: Límite superior del intervalo de integración.
        tol: La tolerancia deseada para la cota de error.

    Returns:
        El número de subintervalos 'n' (entero par) requerido.
    """
    x = symbols('x')

    # Derivada cuarta de f (simbólica) usando sympy
    f_4_sym = diff(f_simbolica(x), x, 4)

    # Convertimos la derivada simbólica a una función numérica evaluable
    f_4_num = lambdify(x, f_4_sym, 'numpy')
    num_puntos_eval = 1000 # Número de puntos para la estimación de alpha_max
    xs = np.linspace(a, b, num_puntos_eval)
    alpha_max = max(abs(f_4_num(xi)) for xi in xs)

    factor = ((b - a)**5 * alpha_max) / (180 * tol)
    n_float = factor ** 0.25 # n mínimo requerido, puede ser no entero

    n = math.ceil(n_float)

    # Asegurar que n sea par, ya que la regla de Simpson compuesta requiere un número par de subintervalos.
    if n % 2 == 1:
        n += 1

    return int(n) # Convertimos a entero por si acaso

# Regla compuesta de Simpson para n+1 puntos (es decir, n subintervalos)
def simpson_compuesto(f_numerica, a, b, n):
    """
    Calcula la aproximación de una integral definida usando la regla compuesta de Simpson.

    Args:
        f_numerica: La función numérica f(x) a integrar.
        a: Límite inferior del intervalo de integración.
        b: Límite superior del intervalo de integración.
        n: Número de subintervalos (debe ser par).

    Returns:
        La aproximación de la integral.
    """
    if n % 2 == 1:
        raise ValueError("El número de subintervalos 'n' para la regla de Simpson debe ser par.")

    h = (b - a) / n
    x_puntos = np.linspace(a, b, n + 1) # n+1 puntos para n subintervalos
    y_valores = f_numerica(x_puntos)

    S = y_valores[0] + y_valores[-1] # f(a) + f(b)
    S += 4 * np.sum(y_valores[1:n:2]) # 4 * (f(x1) + f(x3) + ... + f(xn-1))
    S += 2 * np.sum(y_valores[2:n-1:2]) # 2 * (f(x2) + f(x4) + ... + f(xn-2))

    return (h / 3) * S

def f_simbolica(x):
    return exp(x) * (26 - 10*x + x**2)

# --- Parámetros del Problema ---
a_val = 5
b_val = 5.5
tolerancia = 1e-8

print("--- Calculando n para la Cota de Error de Simpson ---")
# Obtener el n necesario para cumplir la tolerancia
n_requerido = cota_simpson_puntos(f_simbolica, a_val, b_val, tolerancia)
print(f"\nValor de n necesario para que la cota de error sea menor que {tolerancia}: {n_requerido}")


x_sym = symbols('x') # Necesitamos definir 'x' como un símbolo aquí también para lambdify
f_numerica = lambdify(x_sym, f_simbolica(x_sym), 'numpy')

print("\n--- Verificación de la Aproximación de la Integral ---")
# Aproximación de la integral con Simpson usando el 'n' calculado
try:
    aprox_integral = simpson_compuesto(f_numerica, a_val, b_val, n_requerido)
    print(f"Aproximación de la integral de f(x) en [{a_val}, {b_val}] usando n = {n_requerido}: {aprox_integral:.10f}")
except ValueError as e:
    print(f"Error al calcular la integral con Simpson: {e}")