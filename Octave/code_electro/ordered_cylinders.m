function cylinders_pos = ordered_cylinders (nb_x, nb_y, radius, size_x, size_y)
	dx = size_x/(nb_x-1); dy = size_y/(nb_y-1);
	if dx < 2*radius || dy < 2*radius;
		error('radius too large');
	end
	x = 0; y = 0;
	for i = 1:nb_x*nb_y;
		cylinders_pos(i,:) = [x y];
		if mod(i, nb_y) == 0;
			y = 0;
			x = x + dx;
			
		else
			y = y + dy;
		end
	end
	cylinders_pos(:,1) = cylinders_pos(:,1) -  size_x/2;
	cylinders_pos(:,2) = cylinders_pos(:,2) -  size_y/2;
end
