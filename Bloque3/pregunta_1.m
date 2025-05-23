function script_principal()
    intervalos = [2 10];
    m_vals = [10 20 50 100 250];
    y0 = 4;

    figure;
    hold on;
    colors = lines(length(m_vals));

    for i = 1:length(m_vals)
        m = m_vals(i);
        [x, y] = runge_kutta_6(intervalos(1), intervalos(2), y0, m);
        plot(x, y, 'Color', colors(i,:), 'DisplayName', ["Aprox m = " num2str(m)]);
    end

    x_exact = linspace(intervalos(1), intervalos(2), 1000);
    y_exact = x_exact .* log(x_exact / 2) + 2 * x_exact;
    plot(x_exact, y_exact, 'k--', 'DisplayName', 'Solución exacta');

    xlabel('x');
    ylabel('y');
    title('Comparación de soluciones: Runge-Kutta orden 6 vs Exacta');
    legend show;
    grid on;
    hold off;
end

function [x, y] = runge_kutta_6(a, b, y0, m)
    h = (b - a) / (m - 1);
    x = linspace(a, b, m);
    y = zeros(1, m);
    y(1) = y0;

    f = @(x, y) (x + y) / x;

    for i = 1:m-1
        k1 = h * f(x(i), y(i));
        k2 = h * f(x(i) + h/3, y(i) + k1/3);
        k3 = h * f(x(i) + 2*h/5, y(i) + (4*k1 + 6*k2)/25);
        k4 = h * f(x(i) + h, y(i) + (k1 - 12*k2 + 15*k3)/4);
        k5 = h * f(x(i) + 2*h/3, y(i) + (6*k1 + 90*k2 - 50*k3 + 8*k4)/81);
        k6 = h * f(x(i) + 4*h/5, y(i) + (6*k1 + 36*k2 + 10*k3 + 8*k4 + 7*k5)/75);

        y(i+1) = y(i) + (23*k1 + 125*k2 - 81*k5 + 125*k6) / 192;
    end
end

% Ejecutar automáticamente si se compila
script_principal();

