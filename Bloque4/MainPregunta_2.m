f = @(x) exp(x.^2 - 10*x + 26);
a = 5;
b = 5.5;
tol = 1e-8;

n = cota_simpson_puntos(f, a, b, tol);

