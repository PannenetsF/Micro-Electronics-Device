T = 300;
k = 8.617e-5;
e0 = 8.85e-14;
q = 1.602e-19;
KS = 11.8;
ni = 1e10;
EG = 1.12;

xleft = -3.5e-4;
xright = -xleft;
NA = input('Please enter p-side doping (cm^-3), NA = ');
ND = input('Please enter n-side doping (cm^-3), ND = ');

Vbi = k*T*log((NA*ND)/ni^2);
xN = sqrt(2 * KS * e0 / q * NA * Vbi / (ND * (NA+ND)));
xP = sqrt(2 * KS * e0 / q * ND * Vbi / (NA * (NA+ND)));
x = linspace(xleft, xright, 200);

Vx1 = (Vbi - q * ND .* (xN-x).^2/(2*KS*e0).*(x<=xN).*(x>=0));
Vx2 = 0.5 * q * NA .*(xP + x).^2/(KS * e0).*(x >= -xP & x < 0);
Vx = Vx1 + Vx2;
VMAX = 3;
EF = Vx(1) + VMAX/2-k*T*log(NA/ni);


close
plot(x, -Vx+EG/2+VMAX/2); grid
axis([xleft, xright, 0, VMAX])
axis('off');hold on 

plot(x, -Vx+VMAX/2, 'w:');
plot([xleft, xright], [EF, EF], 'w');
plot([0 0], [0.15, VMAX-0.5], 'w--');

text(xleft*1.08, (-Vx(1)+EG/2+VMAX/2-0.05),'Ec');
text(xright*1.02, (-Vx(200)+EG/2+VMAX/2-0.5), 'Ec');
text(xleft*1.08, (-Vx(1)-EG/2+VMAX/2-0.05),'Ev');
text(xright*1.02, (-Vx(200)-EG/2+VMAX/2-0.5), 'Ev');
text(xleft*1.08, (-Vx(1)+VMAX/2-0.05),'Ei');
text(xright*1.02, (EF-0.5), 'EF');
set(gca, 'DefaultTextUnits', 'normalized')
text(.18, 0, 'pside');
text(.47, .0, 'x=0');
text(.75, .0, 'nside');
set(gca, 'DefaultTextUnits', 'data') % gca Return a handle to the current axes object.

