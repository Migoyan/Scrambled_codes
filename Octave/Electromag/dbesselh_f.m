## Copyright (C) 2021 temporaire invite24
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {Function File} {@var{retval} =} dbesselh_f (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: temporaire invite24 <invite24@sa01.intern>
## Created: 2021-10-01

function [retval] = dbesselh_f (n, x)
  retval = 0.5 * (besselh(n-1, x) - besselh(n+1, x));
endfunction
