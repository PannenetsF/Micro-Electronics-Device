EG = 1.12;
kT = 0.0259;
ni = 1.0e10;

ND=logspace(14,17);
Vbi = EG/2 + kT.*log(ND./ni);

close 
semilogx(ND, Vbi);grid
axis([1.0e14,1.0e17,0.75,1]);
xlabel('NA or ND (cm^-3)');
ylabel('Vbi (V)');
text(1e16,0.8,'Si 300k');
text(1e16, 0.78, 'p+/n and + /p diodes');