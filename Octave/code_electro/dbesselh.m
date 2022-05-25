function [retval] = dbesselh (n, k, x)
  retval = n .* besselh(n, k, x) / x .- besselh(n+1, k, x);
endfunction
