clear all; clc;

%% Parameters
%

% Cylinders

nb_cylinders_xy = [3 8]; c_radius = 1e-3; size_box = [2 4];
cylinders = ordered_cylinders(
	nb_cylinders_xy(1), nb_cylinders_xy(2), c_radius, size_box(1), size_box(2)
);
nb_cylinders = prod(nb_cylinders_xy);

% Incident wave parameters
alpha = 0; list_lambda = linspace(1, 40, 100);

% Segment
S_x = 3; S_y = [-2 2 50];

Segment = linspace(S_y(1), S_y(2), S_y(3))';

%% Util inline/lambda/anonymous (whatever you call it) functions
%

u_i = @(k, x, y, a) exp(-1i * k * (x .* cos(a) .+ y .* sin(a)));

c = 299792458;
mu0 = 12.566370614e-7;
for jj = 1:length(list_lambda);
	lambda = list_lambda(jj);
	fprintf("Lambda = %d\n", lambda);
	k = 2 * pi / lambda;
	coeff_scat = coeff_scatter(cylinders, k, c_radius, alpha);

	P = 0;

	omega = k * c;

	U = u_i(k, S_x, Segment, alpha);
	du = -1i .* k .* U;

	dist_to_cylinders = zeros(S_y(3), nb_cylinders);
	for ii = 1:S_y(3);

		dist_to_cylinders(ii, :) = vecnorm([S_x Segment(ii)] .- cylinders, 2, 2);

	end

	for ii = 1:S_y(3);
		U(ii) = U(ii) + sum(coeff_scat .*
			besselh(0, 1, k * abs(dist_to_cylinders(ii, :)))'
		);
		du(ii) = du(ii) .+ sum(coeff_scat .*
			dbesselh(0, 1, k * abs(dist_to_cylinders(ii, :)))'
		);

		P = P + 1/(2*omega*mu0) * imag(conj(U(ii)) * du(ii));
	end
	Ptg(jj) = P * (S_y(2) - S_y(1)) / S_y(3);
end

plot(list_lambda, Ptg);
