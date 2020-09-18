T = 300;
k = 8.617e-5;
e0 = 8.85e-14;
q = 1.602e-19;
KS = 11.8;
ni = 1e10;
EG = 1.12;

NB = logspace(14,17);
VA = [0.5, 0, -10];

Vbi = EG/2 + k*T .* log(NB./ni);
W = zeros(3, size(NB,2));
W = 1.0e4 * sqrt(2 * KS * e0 / q .* (Vbi - VA').*(1./NB));

close
loglog(NB, W, '-'); grid
axis([1.0e14, 1.0e17, 1.0e-1, 1.0e1])
xlabel('NA or ND (cm^-3)');
ylabel('W (um)');
set(gca, 'DefaultTextUnits', 'normalized')
text(.38, .26, 'VA=0.5V');
text(.38, .5, 'VA=0V');
text(.38, .75, 'VA=-10V');
text(.77, .82, 'Si 300K');
text(.77, .79, 'p+/n and n+/p');
set(gca, 'DefaultTextUnits', 'data') % gca Return a handle to the current axes object.
