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
## @deftypefn {} {@var{retval} =} plot_montecalo (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: pannenetsf <pannenetsf@manjaro>
## Created: 2020-11-30

E = 0:10000:1e5;
N = 5;
D = 3;
t_total = 1e-6;
tau = 10.^-(1:D:(1+D*N));
k = tau * 0;
subplot(2,1,1)
for i = 1:N
  velocity = montecalo(E, tau(i), t_total);
  plot(E, velocity);
  k(i) = max(velocity ./ E')
  hold on;
endfor

subplot(2,1,2)
plot(tau, k)
