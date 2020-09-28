T=300;        % Temperature in Kelvin
k=8.617e-5;   % Boltzmann constant (eV/K)
e0=8.85e-14;  % permittivity of free space (F/cm)
q=1.602e-19;  % charge on an electron (coul)
KS=11.8;      % Dielectric constant of Si at 300K
ni=1.0e10;    % intrinsic carrier conc. in Silicon at 300K (cm^-3)
EG=1.12;      % Silicon band gap (eV) 

NA = 1e16;
ND = 1e19;

NDref = 1.3e17;
NAref = 2.35e17;
unmin = 92;
upmin = 54.3;
un0 = 1268;
up0 = 406.9;
an = 0.91;
ap = 0.88;

un = unmin + un0 / (1 + (ND/NDref)^an);
up = upmin + up0 / (1 + (NA/NAref)^ap);

Dn = k * T * un;
Dp = k * T * up;

tp = 1e-6;
tn = tp;

Ln = sqrt(tn * Dn);
Lp = sqrt(tp * Dp);

J0 = q * (Dn*ni^2/Ln/NA + Dp*ni^2/Lp/ND);


Vbi=k*T*log((NA*ND)/ni^2); 
xN=sqrt(2*KS*e0/q*NA*Vbi/(ND*(NA+ND)));    % Depletion width n-side
xP=sqrt(2*KS*e0/q*ND*Vbi/(NA*(NA+ND)));    % Depletion width p-side

W = xN + xP;

J = J0 * e^(VA / k / T) - J0