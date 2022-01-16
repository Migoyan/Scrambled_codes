clear all; clc;

%% Parameters
%

%% Cylinders
%nb_cylinders = 10; c_radius = 1e-1; size_box = [3 5];
%cylinders = random_cylinders(nb_cylinders, c_radius, size_box(1), size_box(2));

nb_cylinders_xy = [5 10]; c_radius = 5e-2; size_box = [3 5];
cylinders = ordered_cylinders(
	nb_cylinders_xy(1), nb_cylinders_xy(2), c_radius, size_box(1), size_box(2)
);
nb_cylinders = prod(nb_cylinders_xy);
% Incident wave parameters
alpha = 0; lambda = 10; k = 2 * pi / lambda;

% Computation parameters

x_lim = [-10 10 500]; y_lim = [-10 10 500]; % intervals + precision



%% Util inline/lambda/anonymous (whatever you call it) functions
%

u_i = @(k, x, y, a) exp(-1i * k * (x .* cos(a) .+ y .* sin(a)));

coeff_scat = coeff_scatter(cylinders, k, c_radius, alpha);

x = linspace(x_lim(1), x_lim(2), x_lim(3));
y = linspace(y_lim(1), y_lim(2), y_lim(3));
[X Y] = meshgrid(x, y);

E = u_i(k, X, Y, alpha);

dist_to_cylinders = zeros(x_lim(3), y_lim(3), nb_cylinders);

for ii = 1:x_lim(3);

	for jj = 1:y_lim(3);
		dist_to_cylinders(ii, jj, :) = vecnorm([X(ii, jj) Y(ii, jj)] .- cylinders, 2, 2);
	end

end

for ii = 1:x_lim(3);

	for jj = 1:y_lim(3);
		E(ii, jj) = E(ii, jj) + sum(coeff_scat .* ...
			squeeze(besselh(0, 1, k * abs(dist_to_cylinders(ii, jj, :))))
		);
	end

end

E = E .* prod(dist_to_cylinders > c_radius, 3);
figure(1, 'Position', [100 100 1400 700]);
scat_plot = pcolor(X, Y, real(E));
set(scat_plot, 'EdgeColor', 'none');
colorbar;

figure(2, 'Position', [100 100 1400 700]);
scat_plot = pcolor(X, Y, imag(E));
set(scat_plot, 'EdgeColor', 'none');
colorbar;
