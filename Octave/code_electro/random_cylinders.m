function cylinders_pos = random_cylinders (nb, radius, size_x, size_y)
	%% This function dispose nb cylinders at random position in 1x1 box and check
	%% that there is no overlap between cylinders

	cylinders_pos = zeros(nb, 2); % 1 = pos_x, 2 = pos_y

	for ii = 1:nb;

		while 1; % do while equivalent
			cylinders_pos(ii, :) = rand(1, 2); overlap = false;

			% Check that the new cylinder doesn't overlap with an old one
			for jj = 1:ii - 1;

				if vecnorm(cylinders_pos(ii, :) - cylinders_pos(jj, :), 2, 2) <= 2 * radius;
					overlap = true;
					break
				end

			end

			if !overlap; % While condition
				break
			end

		end

	end

	% Centering the 1x1 box around Zero
	cylinders_pos = cylinders_pos .- 0.5;
end
