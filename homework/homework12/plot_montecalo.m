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

E = 0:10:1000;
N = 7;
D = 1;
t_total = 1e-6;
tau = 1e-10 * (1:N);
k = tau * 0;
eps = 1e-10;
m0 = 9.10938356e-31;
m_eff = 1.09*m0; % 300K Si
q = 1.602176634e-19;

subplot(3,1,1)
for i = 1:N
  velocity = montecalo(E, tau(i), t_total);
  plot(E, velocity,'o-o');
  k(i) = mean(velocity) / mean(E+eps);
  hold on;
endfor

subplot(3,1,2)
plot(tau, k,'o-o')
subplot(3,1,3)
plot(tau, q/m_eff*tau,'x-x')
