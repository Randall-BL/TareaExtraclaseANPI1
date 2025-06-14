import numpy as np

def nodos_pesos(k):
    """
    Retorna los nodos (x) y pesos (w) de la cuadratura de Gauss-Legendre de la pagina brindada.

    Los valores están tabulados (no calculados) para órdenes entre 1 y 10 inclusive,
    y corresponden al intervalo [-1, 1].

    Parámetros:
        k : orden de la cuadratura

    Retorna:
        dupla: (w, x), pesos y nodos
    """
    nodos_pesos_dict = {
        1: (
            np.array([2.0]),    #Vacio no se usa
            np.array([0.0])
        ),
        2: (
            np.array([1.0000000000000000, 1.0000000000000000]),
            np.array([-0.5773502691896257, 0.5773502691896257])
        ),
        3: (
            np.array([0.8888888888888888, 0.5555555555555556, 0.5555555555555556]),
            np.array([0.0000000000000000, -0.7745966692414834, 0.7745966692414834])
        ),
        4: (
            np.array([0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538]),
            np.array([-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526])
        ),
        5: (
            np.array([0.5688888888888889, 0.4786286704993665, 0.4786286704993665,
                      0.2369268850561891, 0.2369268850561891]),
            np.array([0.0000000000000000, -0.5384693101056831, 0.5384693101056831,
                      -0.9061798459386640, 0.9061798459386640])
        ),
        6: (
            np.array([0.3607615730481386, 0.3607615730481386, 0.4679139345726910,
                      0.4679139345726910, 0.1713244923791704, 0.1713244923791704]),
            np.array([-0.6612093864662645, 0.6612093864662645, -0.2386191860831969,
                      0.2386191860831969, -0.9324695142031521, 0.9324695142031521])
        ),
        7: (
            np.array([0.4179591836734694, 0.3818300505051189, 0.3818300505051189,
                      0.2797053914892766, 0.2797053914892766, 0.1294849661688697,
                      0.1294849661688697]),
            np.array([0.0000000000000000, -0.4058451513773972, 0.4058451513773972,
                      -0.7415311855993945, 0.7415311855993945, -0.9491079123427585,
                      0.9491079123427585])
        ),
        8: (
            np.array([0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                      0.3137066458778873, 0.2223810344533745, 0.2223810344533745,
                      0.1012285362903763, 0.1012285362903763]),
            np.array([-0.1834346424956498, 0.1834346424956498, -0.5255324099163290,
                      0.5255324099163290, -0.7966664774136267, 0.7966664774136267,
                      -0.9602898564975363, 0.9602898564975363])
        ),
        9: (
            np.array([0.3302393550012598, 0.1806481606948574, 0.1806481606948574,
                      0.0812743883615744, 0.0812743883615744, 0.3123470770400029,
                      0.3123470770400029, 0.2606106964029354, 0.2606106964029354]),
            np.array([0.0000000000000000, -0.8360311073266358, 0.8360311073266358,
                      -0.9681602395076261, 0.9681602395076261, -0.3242534234038089,
                      0.3242534234038089, -0.6133714327005904, 0.6133714327005904])
        ),
        10: (
            np.array([0.2955242247147529, 0.2955242247147529, 0.2692667193099963,
                      0.2692667193099963, 0.2190863625159820, 0.2190863625159820,
                      0.1494513491505806, 0.1494513491505806, 0.0666713443086881,
                      0.0666713443086881]),
            np.array([-0.1488743389816312, 0.1488743389816312, -0.4333953941292472,
                      0.4333953941292472, -0.6794095682990244, 0.6794095682990244,
                      -0.8650633666889845, 0.8650633666889845, -0.9739065285171717,
                      0.9739065285171717])
        )
    }

    if k < 1 or k > 10:
        raise ValueError("k debe estar entre 1 y 10")

    return nodos_pesos_dict[k]


def cuad_gauss(f, a, b, k):
    """
    Aplica la cuadratura Gaussiana simple para aproximar la integrales definidas 
    utilizando el intervalo [-1, 1].

    Parámetros:
        f : función a integrar
        a : límite inferior del intervalo
        b : límite superior del intervalo
        k : orden de la cuadratura (nodos/pesos)

    Retorna:
        Aproximación numérica de la integral definida
    """
    w, x = nodos_pesos(k)  # pesos y nodos en [-1,1]
    # Transformar nodos
    x_trans = (1/2) * ((b - a) * x + (b + a))
    #Constante + sumatoria
    return (1/2) * (b - a) * np.sum(w * f(x_trans))


def cuad_gauss_comp(f, a, b, k, n):
    """
    Aplica la cuadratura Gaussiana de forma compuesta sobre el intervalo [a, b].

    Parámetros:
        f : función a integrar
        a : límite inferior del intervalo
        b : límite superior del intervalo
        k : orden de la cuadratura
        n : número de subintervalos

    Retorna:
        aproximación numérica de la integral definida
    """
    h = (b - a) / n
    aprox = 0.0 #Sumatoria empieza en 0
    for i in range(n):
        xi = a + i * h  # Expresion necesaria
        bi = xi + h     # Expresion necesaria

        #Sumatoria de cuad_gauss simple
        aprox += cuad_gauss(f, xi, bi, k)
    return aprox


def cuad_gauss_iter(f, a, b, k, tol=1e-8, max_n=1000000):
    """
    Implementa una versión iterativa de la cuadratura Gaussiana compuesta.

    Parámetros:
        f : función a integrar
        a : límite inferior del intervalo
        b : límite superior del intervalo
        k : orden de la cuadratura
        tol : tolerancia de convergencia 
        max_n : máximo número de subintervalos 

    Retorna:
        aproximación numérica de la integral cumpliendo tolerancia
    """
    n = 1
    Sn = cuad_gauss_comp(f, a, b, k, n) # Usa la compuesta para hacer calculos iterativos
    while n <= max_n:
        n += 1
        print(n)
        Sn1 = cuad_gauss_comp(f, a, b, k, n)
        if abs(Sn1 - Sn) < tol:
            return Sn1
        Sn = Sn1
    raise RuntimeError("No se alcanzó la tolerancia con n <= 1e-8")

##########
#Aplicación
##########
f = lambda x: np.cos(x) * np.exp(x)
a = 2
b = 5
k = 7

print("Simple:", cuad_gauss(f, a, b, k))
print("Compuesta:", cuad_gauss_comp(f, a, b, k, 20))
iter = cuad_gauss_iter(f, a, b, k)
print(f"Iterativa: {iter}")
