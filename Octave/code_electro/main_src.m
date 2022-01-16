clear all; clc;

%% Parameters
%

% Cylinders
nb_cylinders = 10; c_radius = 1e-1;
% Incident wave parameters
alpha = 0; lambda = 10; k = 2 * pi / lambda;

% Computation parameters

x_lim = [-10 10 300]; y_lim = [-10 10 300]; % intervals + precision

% Generating cylinders
cylinders = random_cylinders(nb_cylinders, c_radius);

%% Util inline/lambda/anonymous (whatever you call it) functions
%

u_i = @(k, x, y, a) exp(-1i * k * (x .* cos(a) .+ y .* sin(a)));

coeff_scat = coeff_scatter(cylinders, k, c_radius, alpha);

x = linspace(x_lim(1), x_lim(2), x_lim(3));
y = linspace(y_lim(1), y_lim(2), y_lim(3));
[X Y] = meshgrid(x, y);

E = u_i(k, X, Y, alpha);

dist_to_cylinders = zeros(nb_cylinders);

for ii = 1:(x_lim(3) * y_lim(3));

    if any(vecnorm([X(ii) Y(ii)] - cylinders(:, :), 2, 2) < c_radius);
        E(ii) = 0;
    else
        dist_to_cylinders = vecnorm([X(ii) Y(ii)] .- cylinders, 2, 2);

        for jj = 1:nb_cylinders
            E(ii) = E(ii) + ...
                coeff_scat(jj) .* besselh(0, 1, k * abs(dist_to_cylinders(jj)));
        end

    end

end

figure(1, 'Position', [100 100 1400 700]);
scat_plot = pcolor(X, Y, real(E));
set(scat_plot, 'EdgeColor', 'none');
colorbar;

figure(2, 'Position', [100 100 1400 700]);
scat_plot = pcolor(X, Y, imag(E));
set(scat_plot, 'EdgeColor', 'none');
colorbar;
