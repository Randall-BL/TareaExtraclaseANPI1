#Ejercicio 2: Diferencias finitas y método de Thomas
import numpy as np
import matplotlib.pyplot as plt

def thomas_algorithm(a, b, c, d):
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
    n = int((b - a) / h)
    x = np.linspace(a, b, n+1)
    A = np.zeros(n-1)
    B = np.zeros(n-1)
    C = np.zeros(n-1)
    D = np.zeros(n-1)

    for i in range(1, n):
        xi = x[i]
        A[i-1] = 1 - h*p(xi)/2
        B[i-1] = -2 + h**2*q(xi)
        C[i-1] = 1 + h*p(xi)/2
        D[i-1] = -h**2*r(xi)

    D[0] -= A[0]*y0
    D[-1] -= C[-1]*yn

    y_inner = thomas_algorithm(A[1:], B, C[:-1], D)
    y = np.concatenate([[y0], y_inner, [yn]])
    return x, y

# Validación con problema especificado
p = lambda x: 1/x
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

