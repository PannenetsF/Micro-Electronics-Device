% Equilibrium Energy Band Diagram Generator
%(Si, 300K, nondegenerately doped step junction)

%Constants
T=300;        % Temperature in Kelvin
k=8.617e-5;   % Boltzmann constant (eV/K)
e0=8.85e-14;  % permittivity of free space (F/cm)
q=1.602e-19;  % charge on an electron (coul)
KS=11.8;      % Dielectric constant of Si at 300K
ni=1.0e10;    % intrinsic carrier conc. in Silicon at 300K (cm^-3)
EG=1.12;      % Silicon band gap (eV)

%Control constants
xleft = -3.5e-4;   % Leftmost x position
xright = -xleft;   % Rightmost x position
% NA=input ('Please enter p-side doping (cm^-3), NA = ');
% ND=input ('Please enter n-side doping (cm^-3), ND = ');
NA = 1e17;
ND = 1e14;
%Computations
Vbi=k*T*log((NA*ND)/ni^2); 
xN=sqrt(2*KS*e0/q*NA*Vbi/(ND*(NA+ND)));    % Depletion width n-side
xP=sqrt(2*KS*e0/q*ND*Vbi/(NA*(NA+ND)));    % Depletion width p-side
x = linspace(xleft, xright, 200);
Vx1=(Vbi-q*ND.*(xN-x).^2/(2*KS*e0).*(x<=xN)).*(x>=0); 
Vx2=0.5*q*NA.*(xP+x).^2/(KS*e0).*( x>=-xP & x<0 );
Vx=Vx1+Vx2;     % V as a function of x
VMAX = 3;                        % Maximium Plot Voltage
EF=Vx(1)+VMAX/2-k*T*log(NA/ni);  % Fermi level

%Plot Diagram
close
plot ( x, -Vx+EG/2+VMAX/2);
axis ([xleft xright 0 VMAX]);
axis ('off');  hold on
plot ( x, -Vx-EG/2+VMAX/2);
plot ( x, -Vx+VMAX/2,'w:');
plot ( [xleft xright], [ EF EF ], 'w' );
plot ( [ 0 0 ], [ 0.15 VMAX-0.5 ], 'w--' );
text(xleft*1.08,(-Vx(1)+EG/2+VMAX/2-.05),'Ec');
text(xright*1.02,(-Vx(200)+EG/2+VMAX/2-.05),'Ec');
text(xleft*1.08,(-Vx(1)-EG/2+VMAX/2-.05),'Ev');
text(xright*1.02,(-Vx(200)-EG/2+VMAX/2-.05),'Ev');
text(xleft*1.08,(-Vx(1)+VMAX/2-.05),'Ei');
text(xright*1.02, EF-.05,'EF');
set(gca,'DefaultTextUnits','normalized')
text(.18, 0,'p-side');
text ( .47, 0, 'x=0');
text(.75, 0,'n-side');
set(gca,'DefaultTextUnits','data')
hold off



