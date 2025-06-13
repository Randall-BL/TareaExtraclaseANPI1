function n = cota_simpson_puntos(f, a, b, tol)
  printf("Calculando la cuarta derivada...\n");
  printf("Buscando el máximo de |f⁽⁴⁾(x)| en [%.1f, %.1f]...\n", a, b);

  muestras = 1000;
  h_der = 1e-3;
  x_vals = linspace(a, b, muestras);
  derivadas = zeros(1, muestras);

  for i = 1:muestras
    x = x_vals(i);
    derivadas(i) = abs((f(x - 2*h_der) - 4*f(x - h_der) + 6*f(x) - 4*f(x + h_der) + f(x + 2*h_der)) / h_der^4);
  end

  alpha_max = max(derivadas);
  printf("max|f⁽⁴⁾(x)| = %.6f\n", alpha_max);

  h_teorico = ((180 * tol) / ((b - a) * alpha_max))^(1/4);
  n_real = (b - a) / h_teorico;
  n = 2 * ceil(n_real / 2);  % Par más cercano mayor o igual

  h = (b - a) / n;
  error_estimado = ((b - a) * h^4 / 180) * alpha_max;

  printf("\nRESULTADOS:\n");
  printf("n calculado (real): %.2f\n", n_real);
  printf("n ajustado (entero par): %d\n", n);
  printf("Cota de error con n=%d: %.2e\n", n, error_estimado);

  integral_aprox = regla_simpson_compuesta(f, a, b, n);
  printf("\nIntegral aproximada con Simpson y n=%d: %.10f\n", n, integral_aprox);
endfunction

function I = regla_simpson_compuesta(f, a, b, n)
  h = (b - a) / n;
  x = linspace(a, b, n + 1);
  y = arrayfun(f, x);

  suma = y(1) + y(end);
  suma += 4 * sum(y(2:2:end-1));
  suma += 2 * sum(y(3:2:end-2));

  I = (h / 3) * suma;
endfunction


