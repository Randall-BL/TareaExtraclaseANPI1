function pregunta_1
  clc; clear; close all;
  pkg load symbolic
  syms x

  f = log(asin(x)) / log(x);  % Función simbólica

  % Límite inferior y superior de la integral
  a = 0.2;
  b = 0.8;

  n = 2;  % Número de subintervalos (n > 0)

  err_trap = cota_error_trapecio(f, a, b, n);
  err_simp = cota_error_simpson(f, a, b, n);

  % Imprime las cotas de error
  fprintf('Cota de error (Trapecio): %.10e\n', err_trap);
  fprintf('Cota de error (Simpson): %.10e\n', err_simp);
endfunction

function dmax = max_abs_deriv(fsym, a, b, orden)
  % max_abs_deriv: Calcula el valor máximo de la derivada absoluta de una función simbólica
  % en el intervalo [a, b] para el orden de derivada especificado.
  %
  % Parámetros:
  % fsym  : La función simbólica que representa la ecuación a derivar.
  % a     : Límite inferior del intervalo.
  % b     : Límite superior del intervalo.
  % orden : El orden de la derivada .
  %
  % Retorna:
  % dmax : El valor máximo absoluto de la derivada en el intervalo [a, b].
  %
  % Este valor se usa para calcular la cota de error de las reglas del trapecio y Simpson.

  syms x
  f_deriv = diff(fsym, x, orden); %Derivada
  f_abs = matlabFunction(abs(f_deriv), 'vars', x);

  % fminbnd minimiza, así que usamos el negativo para hallar el máximo
  [~, minval] = fminbnd(@(x) -f_abs(x), a, b);
  dmax = abs(minval);
endfunction

function err = cota_error_trapecio(fsym, a, b, n)
  % cota_error_trapecio: Calcula la cota de error para la regla del trapecio.
  %
  % Parámetros:
  % fsym : La función simbólica que se va a integrar.
  % a    : Límite inferior del intervalo.
  % b    : Límite superior del intervalo.
  % n    : Número de subintervalos.
  %
  % Retorna:
  % err  : La cota de error para el método del trapecio.

  dmax = max_abs_deriv(fsym, a, b, 2);  % Segunda derivada
  err = ((b - a)^3 / (12)) * dmax;
endfunction

function err = cota_error_simpson(fsym, a, b, n)
  % cota_error_simpson: Calcula la cota de error para la regla de Simpson.
  %
  % Parámetros:
  % fsym  : La función simbólica que se va a integrar.
  % a     : Límite inferior del intervalo.
  % b     : Límite superior del intervalo.
  % n     : Número de subintervalos (debe ser par).
  %
  % Retorna:
  % err (float): La cota de error para el método de Simpson.
  if mod(n, 2) ~= 0
    error('Simpson requiere un número par de subintervalos');
  end
  dmax = max_abs_deriv(fsym, a, b, 4);  % Cuarta derivada
  err = ((b - a)^5 / (2880)) * dmax;
endfunction

