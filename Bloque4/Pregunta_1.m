pkg load symbolic

% Función para calcular la cota de error de la regla del trapecio
function error_bound = trapezoid_error_bound(f_expr, a, b, n)
    % f_expr: expresión simbólica de la función
    % a, b: límites de integración
    % n: número de subintervalos

    % Declarar variable simbólica
    syms x

    % Calcular la segunda derivada de f
    f_second_deriv = diff(f_expr, x, 2);

    % Convertir a función numérica para evaluación
    f_second_numeric = matlabFunction(abs(f_second_deriv));

    % Encontrar el máximo de |f''(x)| en [a,b]
    % Usamos fminbnd con función negativa para encontrar el máximo
    try
        [~, max_f_second] = fminbnd(@(x) -f_second_numeric(x), a, b);
        max_f_second = -max_f_second;
    catch
        % Si hay problemas con fminbnd, evaluar en puntos discretos
        x_vals = linspace(a, b, 1000);
        f_vals = arrayfun(f_second_numeric, x_vals);
        max_f_second = max(f_vals);
    end

    % Calcular la cota de error del trapecio
    h = (b - a) / n;
    error_bound = (b - a) * h^2 * max_f_second / 12;

    fprintf('Regla del Trapecio:\n');
    fprintf('h = %.6f\n', h);
    fprintf('max|f''''(x)| en [%.2f, %.2f] = %.6f\n', a, b, max_f_second);
    fprintf('Cota de error = %.8e\n\n', error_bound);
end

% Función para calcular la cota de error de la regla de Simpson
function error_bound = simpson_error_bound(f_expr, a, b, n)
    % f_expr: expresión simbólica de la función
    % a, b: límites de integración
    % n: número de subintervalos (debe ser par)

    if mod(n, 2) ~= 0
        error('n debe ser par para la regla de Simpson');
    end

    % Declarar variable simbólica
    syms x

    % Calcular la cuarta derivada de f
    f_fourth_deriv = diff(f_expr, x, 4);

    % Convertir a función numérica para evaluación
    f_fourth_numeric = matlabFunction(abs(f_fourth_deriv));

    % Encontrar el máximo de |f^(4)(x)| en [a,b]
    try
        [~, max_f_fourth] = fminbnd(@(x) -f_fourth_numeric(x), a, b);
        max_f_fourth = -max_f_fourth;
    catch
        % Si hay problemas con fminbnd, evaluar en puntos discretos
        x_vals = linspace(a, b, 1000);
        f_vals = arrayfun(f_fourth_numeric, x_vals);
        max_f_fourth = max(f_vals);
    end

    % Calcular la cota de error de Simpson
    h = (b - a) / n;
    error_bound = (b - a) * h^4 * max_f_fourth / 180;

    fprintf('Regla de Simpson:\n');
    fprintf('h = %.6f\n', h);
    fprintf('max|f^(4)(x)| en [%.2f, %.2f] = %.6f\n', a, b, max_f_fourth);
    fprintf('Cota de error = %.8e\n\n', error_bound);
end

% Validación con la integral dada: ∫[0.2 to 0.8] ln(sin^(-1)(x))/ln(x) dx
fprintf('=== VALIDACIÓN DE LAS IMPLEMENTACIONES ===\n\n');

% Definir la función simbólica
syms x
f_sym = log(asin(x)) / log(x);

% Parámetros de integración
a = 0.2;
b = 0.8;
n_trap = 10;  % número de subintervalos para trapecio
n_simp = 10;  % número de subintervalos para Simpson (debe ser par)

fprintf('Función: f(x) = ln(sin^(-1)(x))/ln(x)\n');
fprintf('Intervalo: [%.1f, %.1f]\n', a, b);
fprintf('Número de subintervalos: %d\n\n', n_trap);

% Calcular cotas de error
fprintf('COTAS DE ERROR:\n');
error_trap = trapezoid_error_bound(f_sym, a, b, n_trap);
error_simp = simpson_error_bound(f_sym, a, b, n_simp);

% Implementación adicional de los métodos numéricos para comparación
function I = trapezoid_rule(f_func, a, b, n)
    h = (b - a) / n;
    x = a:h:b;
    y = arrayfun(f_func, x);
    I = h * (y(1)/2 + sum(y(2:end-1)) + y(end)/2);
end

function I = simpson_rule(f_func, a, b, n)
    if mod(n, 2) ~= 0
        error('n debe ser par para Simpson');
    end
    h = (b - a) / n;
    x = a:h:b;
    y = arrayfun(f_func, x);

    odd_sum = sum(y(2:2:end-1));   % índices impares
    even_sum = sum(y(3:2:end-1));  % índices pares

    I = h/3 * (y(1) + 4*odd_sum + 2*even_sum + y(end));
end

% Convertir función simbólica a numérica para integración
f_numeric = @(x) log(asin(x)) ./ log(x);

% Calcular aproximaciones numéricas
I_trap = trapezoid_rule(f_numeric, a, b, n_trap);
I_simp = simpson_rule(f_numeric, a, b, n_simp);

fprintf('RESULTADOS DE INTEGRACIÓN NUMÉRICA:\n');
fprintf('Aproximación por trapecio: %.8f\n', I_trap);
fprintf('Aproximación por Simpson:  %.8f\n', I_simp);

% Intentar calcular el valor "exacto" usando quad para comparación
try
    I_exact = quad(f_numeric, a, b, 1e-10);
    fprintf('Valor de referencia (quad): %.8f\n\n', I_exact);

    fprintf('ERRORES REALES vs COTAS:\n');
    fprintf('Error real trapecio:  %.8e (cota: %.8e)\n', abs(I_trap - I_exact), error_trap);
    fprintf('Error real Simpson:   %.8e (cota: %.8e)\n', abs(I_simp - I_exact), error_simp);

    if abs(I_trap - I_exact) <= error_trap
        fprintf('✓ Cota del trapecio es válida\n');
    else
        fprintf('✗ Cota del trapecio podría ser insuficiente\n');
    end

    if abs(I_simp - I_exact) <= error_simp
        fprintf('✓ Cota de Simpson es válida\n');
    else
        fprintf('✗ Cota de Simpson podría ser insuficiente\n');
    end
catch
    fprintf('No se pudo calcular valor de referencia con quad\n');
end
