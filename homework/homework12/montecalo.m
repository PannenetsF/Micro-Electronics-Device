## Copyright (C) 2020 pannenetsf
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{velocity} =} montecalo (@var{E}, @var{t_total})
##
## @seealso{}
## @end deftypefn

## Author: pannenetsf <pannenets.f@foxmail.com>
## Created: 2020-11-30

function velocity = montecalo (E, tau, t_total)

Tau = 1/tau;

m0 = 9.10938356e-31;
m_eff = 1.09*m0; % 300K Si
q = 1.602176634e-19;
a = E * q / m_eff;

% The integration of probality is 1.
mean_delta_t = 1/Tau;
num = 10 * round(t_total/mean_delta_t) + 1;

E_size = size(E,2);
r = rand(E_size,num);
delta_t = -log(r) / Tau;
l = 0.5 * a' .* (delta_t.^2);
L = sum(l,2);

velocity = L ./ sum(delta_t,2);

endfunction
