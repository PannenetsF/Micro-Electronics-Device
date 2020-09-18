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
VA = input('Please enter VA (V), VA = ');
%NA = 1e18;
%ND = 1e16;

Vbi=k*T*log((NA*ND)/ni^2);
Vbi=Vbi-VA;
xN=sqrt(2*KS*e0/q*NA*Vbi/(ND*(NA+ND)));    % Depletion width n-side
xP=sqrt(2*KS*e0/q*ND*Vbi/(NA*(NA+ND)));    % Depletion width p-side
x = linspace(xleft, xright, 200);
% Vx1=(Vbi-q*ND.* (xN-x).^2/(2*KS*e0).*(x<=xN).*(x>=0));
Vx1=(Vbi-q*ND.*(xN-x).^2/(2*KS*e0).*(x<=xN)).*(x>=0); 
Vx2=0.5*q*NA.*(xP+x).^2/(KS*e0).*( x>=-xP & x<0 );
Vx = Vx1 + Vx2;
VMAX = 3;
EF = Vx(1) + VMAX/2-k*T*log(NA/ni);


close

subplot(5,1,1);

str_title = sprintf('ND = %e, NA = %e Enegry Band', ND, NA);
title(str_title);

plot(x, -Vx+EG/2+VMAX/2); grid
axis([xleft, xright, 0, VMAX])
axis('off');hold on 
plot ( x, -Vx-EG/2+VMAX/2);
plot(x, -Vx+VMAX/2, 'w:');
plot([xleft, xright], [EF, EF], 'w');
plot([0 0], [0.15, VMAX-0.5], 'w--');

text(xleft*1.08, (-Vx(1)+EG/2+VMAX/2-0.05),'Ec');
text(xright*1.02, (-Vx(200)+EG/2+VMAX/2-0.05), 'Ec');
text(xleft*1.08, (-Vx(1)-EG/2+VMAX/2-0.05),'Ev');
text(xright*1.02, (-Vx(200)-EG/2+VMAX/2-0.5), 'Ev');
text(xleft*1.08, (-Vx(1)+VMAX/2-0.05),'Ei');
text(xright*1.02, (EF-0.05), 'EF');
set(gca, 'DefaultTextUnits', 'normalized')
text(.18, 0, 'pside');
text(.47, .0, 'x=0');
text(.75, .0, 'nside');
set(gca, 'DefaultTextUnits', 'data') % gca Return a handle to the current axes object.
title(str_title);


subplot(5,1,2);

str_title = sprintf('ND = %e, NA = %e Distro of Impurities', ND, NA);
title(str_title);

hold on;

axis([xleft, xright, -20, 20]);
plot(x, -log10(NA*(x<0)));
plot(x, log10(ND*(x>=0)));
xlabel('x axis');
ylabel('ND-NA in log10');


subplot(5,1,3);
eps = 1e-16;
str_title = sprintf('ND = %e, NA = %e desity of charge', ND, NA);
title(str_title);
hold on;
axis(2*[10*(-xP), xN]);
mask_p = (x < 0) & (x >= -xP-eps);
mask_n = (x >= 0) & (x <= xN);
plot(x, -log10(q*NA*mask_p));
plot(x, log10(q*ND*mask_n+eps));

xlabel('x axis');
ylabel('\rho in log10');




subplot(5,1,4);

str_title = sprintf('ND = %e, NA = %e electric field', ND, NA);
title(str_title);
hold on;
axis(2*[(-xP), xN]);

plot(x, (-q*NA*mask_p/KS/e0.*(xP+x)));
plot(x, (-q*ND*mask_n/KS/e0.*(xN-x)));

xlabel('x axis');
ylabel('E');



subplot(5,1,5);

str_title = sprintf('ND = %e, NA = %e potential', ND, NA);
title(str_title);
hold on;
axis([(-xP), xN]);

plot(x, (q*NA*mask_p/2/KS/e0.*(xP+x).^2));
plot(x, (Vbi-q*ND*mask_n/KS/e0.*(xN).^2));

xlabel('x axis');
ylabel('E');