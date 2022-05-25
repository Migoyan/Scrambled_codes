function scatter_matrix = scattering_matrix (cylinders, k, c_radius)

	%% Util inline/lambda/anonymous (whatever you call it) functions
	%
	S0_i = @(k, radius) -1 / (1 + bessely(0, k * radius) / besselj(0, k * radius));
	H_ij = @(k, c_1, c_2) besselh(0, 1, k * abs(vecnorm(c_1 - c_2,2)));
	nb_cylinders = size(cylinders)(1);
	scatter_matrix = ones(nb_cylinders, nb_cylinders);

	S0 = -S0_i(k, c_radius);

	for ii = 1:nb_cylinders;

		for jj = 1:nb_cylinders;
			if ii == jj;
				scatter_matrix(ii, jj) = 1;
			else
				scatter_matrix(ii, jj) = S0 * ...
				H_ij(k, cylinders(ii, :), cylinders(jj, :));
			end
			
		end

	end

end
