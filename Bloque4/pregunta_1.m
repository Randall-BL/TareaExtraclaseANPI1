function Ej1
  pkg load symbolic
  syms x

  f = log(asin(x)) / log(x);  % Función simbólica
  a = 0.2;
  b = 0.8;
  n = 2;  % Número de subintervalos (n > 0)

  err_trap = cota_error_trapecio(f, a, b, n);
  err_simp = cota_error_simpson(f, a, b, n);

  fprintf('Cota de error (Trapecio) con n=%d: %.10e\n', n, err_trap);
  fprintf('Cota de error (Simpson) con n=%d: %.10e\n', n, err_simp);
endfunction

function dmax = max_abs_deriv(fsym, a, b, orden)
  syms x
  f_deriv = diff(fsym, x, orden);
  f_abs = matlabFunction(abs(f_deriv), 'vars', x);

  % fminbnd minimiza, así que usamos el negativo para hallar el máximo
  [~, minval] = fminbnd(@(x) -f_abs(x), a, b);
  dmax = abs(minval);
endfunction

function err = cota_error_trapecio(fsym, a, b, n)
  dmax = max_abs_deriv(fsym, a, b, 2);  % Segunda derivada
  err = ((b - a)^3 / (12)) * dmax;
endfunction

function err = cota_error_simpson(fsym, a, b, n)
  if mod(n, 2) ~= 0
    error('Simpson requiere un número par de subintervalos');
  end
  dmax = max_abs_deriv(fsym, a, b, 4);  % Cuarta derivada
  err = ((b - a)^5 / (180 * n^4)) * dmax;
endfunction

