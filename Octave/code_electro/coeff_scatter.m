function s0_i = coeff_scatter (cylinders, k, c_radius, a)
	%% Util inline/lambda/anonymous (whatever you call it) functions
	%
	S0_i = @(k, radius) -1 / (1 + bessely(0, k * radius) / besselj(0, k * radius));
	u_i = @(k, x, y, a) exp(-1i * k * (x .* cos(a) + y .* sin(a)));

	scatter_matrix_i = inverse(scattering_matrix(cylinders, k, c_radius));
	S0 = -S0_i(k, c_radius);
	Q_matrix = S0 .* u_i(k, cylinders(:, 1), cylinders(:, 2), a);
	s0_i = scatter_matrix_i * Q_matrix;
end
