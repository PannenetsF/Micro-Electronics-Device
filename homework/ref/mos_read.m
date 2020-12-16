 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%                                                   %%  
         %%     1D Drift Diffusion Model for pn Diodes        %%  
         %%     Equilibrium and Non Equilibrium Solver        %%
         %%                                                   %%  
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
clear all;
close all;


% Defining the Fundamental and Material Constants %
% 定义了相关的常量

q     = 1.602E-19;        % C or [J/eV] 元电荷量 单位为 C
kb    = 1.38E-23;         % [J/K] 波尔兹曼常量 单位为 J/K
eps   = 1.05E-12;         % This includes the eps  = 11.7 for Si [F/cm] Si的绝对节点常数
T     = 300;              % [K] 本次仿真的环境温度
ni    = 1.5E10;           % Intrinsic carrier concentration [1/cm^3] 本征载流子浓度
Vt    = kb*T/q;           % [eV] 此温度下的热电压
RNc   = 2.8E19;           % This is 2.8e20 in the FORTRAN file 2
TAUN0 = 0.1E-6;           % Electron SRH life time 
TAUP0 = 0.1E-6;           % Hole SRH life time
mun0   = 1500;            % Electron Mobility in cm2/V-s 电子迁移率
mup0   = 1000;             % Hole Mobility in cm2/V-s 空穴迁移率
                        
                        
dEc = Vt*log(RNc/ni);

% Define Doping Values %

Na = 1E16;             % [1/cm^3] 受主掺杂浓度 单位为 1/cm^3
Nd = 1E17;             % [1/cm^3] 施主掺杂浓度 单位为 1/cm^3

% Calculate relevant parameters for the simulation %
% 计算相关参数

Vbi = Vt*log(Na*Nd/(ni*ni)); % 计算内建电势
W   = sqrt(2*eps*(Na+Nd)*Vbi/(q*Na*Nd))     % [cm] 计算耗尽层的宽度 单位为 cm
Wn  = W*sqrt(Na/(Na+Nd))                    % [cm] 计算 n 型一侧的宽度 单位为 cm
Wp  = W*sqrt(Nd/(Na+Nd))                    % [cm] 计算 p 型一侧的宽度 单位为 cm
Wone = sqrt(2*eps*Vbi/(q*Na))               % [cm] 
E_p = q*Nd*Wn/eps                           % [V/cm] pn结交界处电场强度
Ldn = sqrt(eps*Vt/(q*Nd));                  
Ldp = sqrt(eps*Vt/(q*Na));
Ldi = sqrt(eps*Vt/(q*ni))

% Calculate relevant parameters in an input file %

% Write to a file
save input_params.txt Na Nd Vbi W Wn Wp E_p Ldn Ldp

%Material_Constants    %Define some material constants
      
% Setting the size of the simulation domain based 
% on the analytical results for the width of the depletion regions
% for a simple pn-diode %

x_max = 0;
if(x_max < Wn)
    x_max = Wn;
end
if(x_max < Wp)
    x_max = Wp;
end
x_max = 20*x_max

% Setting the grid size based on the extrinsic Debye lengths %

dx = Ldn;
if(dx > Ldp)
    dx=Ldp;
end
dx = dx/20;

% Calculate the required number of grid points and renormalize dx %

n_max = x_max/dx
n_max = round(n_max);

dx = dx/Ldi;    % Renormalize lengths with Ldi

% Set up the doping C(x)=Nd(x)-Na(x) that is normalized with ni %

for i = 1:n_max
    if(i <= n_max/2)
        dop(i) = - Na/ni;
    elseif(i > n_max/2)
        dop(i) = Nd/ni;
    end
end

% Initialize the potential based on the requirement of charge
% neutrality throughout the whole structure

for i = 1: n_max
    zz = 0.5*dop(i);
    if(zz > 0)
        xx = zz*(1 + sqrt(1+1/(zz*zz)));
    elseif(zz <  0)
        xx = zz*(1 - sqrt(1+1/(zz*zz)));
    end
    fi(i) = log(xx);
    n(i) = xx;
    p(i) = 1/xx;
end

delta_acc = 1E-5;               % Preset the Tolerance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                      %%
%%               EQUILIBRIUM  SOLUTION PART BEGINS                      %%
%%                                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%(A) Define the elements of the coefficient matrix for the internal nodes and
%    initialize the forcing function

dx2 = dx*dx;
for i = 1: n_max
    a(i) = 1/dx2;
    c(i) = 1/dx2;
    b(i) = -(2/dx2+exp(fi(i))+exp(-fi(i)));
    f(i) = exp(fi(i)) - exp(-fi(i)) - dop(i) - fi(i)*(exp(fi(i))+exp(-fi(i)));
 end


%(B) Define the elements of the coefficient matrix and initialize the forcing
%    function at the ohmic contacts 

a(1) = 0;
c(1) = 0;
b(1) = 1;
f(1) = fi(1);
a(n_max) = 0;
c(n_max) = 0;
b(n_max) = 1;
f(n_max) = fi(n_max);

%(C)  Start the iterative procedure for the solution of the linearized Poisson
%     equation using LU decomposition method:

flag_conv = 0;		           % convergence of the Poisson loop
k_iter= 0;
while(~flag_conv)            
    k_iter = k_iter + 1; 
    
    alpha(1) = b(1);
    for i=2:n_max
        beta(i)=a(i)/alpha(i-1);
        alpha(i)=b(i)-beta(i)*c(i-1);
    end
    
% Solution of Lv = f %    

    v(1) = f(1);
    for i = 2:n_max
        v(i) = f(i) - beta(i)*v(i-1);
    end
     
% Solution of U*fi = v %    

    temp = v(n_max)/alpha(n_max);
    delta(n_max) = temp - fi(n_max);
    fi(n_max)=temp;
    for i = (n_max-1):-1:1       %delta%
        temp = (v(i)-c(i)*fi(i+1))/alpha(i);
        delta(i) = temp - fi(i);
        fi(i) = temp;
    end
    
    delta_max = 0;
    
    for i = 1: n_max
        xx = abs(delta(i));
        if(xx > delta_max)
            delta_max=xx;
        end
        %sprintf('delta_max = %d',delta_max)      %'k_iter = %d',k_iter,'
        
    end

     %delta_max=max(abs(delta));
   
% Test convergence and recalculate forcing function and 
% central coefficient b if necessary
    
    if(delta_max < delta_acc)
        flag_conv = 1;
    else
        for i = 2: n_max-1
            b(i) = -(2/dx2 + exp(fi(i)) + exp(-fi(i)));
            f(i) = exp(fi(i)) - exp(-fi(i)) - dop(i) - fi(i)*(exp(fi(i)) + exp(-fi(i)));
        end
    end
end


xx1(1) = dx*1e4;
for i = 2:n_max-1 
    Ec(i) = dEc - Vt*fi(i);     %Values from the second Node%
    ro(i) = -ni*(exp(fi(i)) - exp(-fi(i)) - dop(i));
    el_field1(i) = -(fi(i+1) - fi(i))*Vt/(dx*Ldi);
    el_field2(i) = -(fi(i+1) - fi(i-1))*Vt/(2*dx*Ldi);
    n(i) = exp(fi(i));
    p(i) = exp(-fi(i));
    xx1(i) = xx1(i-1) + dx*Ldi*1e4;
end


Ec(1) = Ec(2);
Ec(n_max) = Ec(n_max-1);
xx1(n_max) = xx1(n_max-1) + dx*Ldi*1e4;
el_field1(1) = el_field1(2);
el_field2(1) = el_field2(2);
el_field1(n_max) = el_field1(n_max-1);
el_field2(n_max) = el_field2(n_max-1);
nf = n*ni;
pf = p*ni;
ro(1) = ro(2);
ro(n_max) = ro(n_max-1);

figure(1)
plot(xx1, Vt*fi,'r','LineWidth',2)
xlabel('x [um]');
ylabel('Potential [eV]');
title('Potential vs Position - at Equilibrium');

figure(2)
plot(xx1, el_field1,'r','LineWidth',2)
hold on;
plot(xx1, el_field2,'r','LineWidth',2)
xlabel('x [um]');
ylabel('Electric Field [V/cm]');
title('Field Profile vs Position - at Equilibrium');

figure(3)
%plot(xx1, nf,'g','LineWidth',2)
semilogy(xx1, nf,'g','LineWidth',2)
hold on;
%plot(xx1, pf,'r','LineWidth',2)
semilogy(xx1, pf,'r','LineWidth',2)
xlabel('x [um]');
ylabel('Electron & Hole Densities [1/cm^3]');
title('Electron & Hole Densities vs Position - at Equilibrium');
legend('n','p');
%axis([0 6.75 0 10.2e17])

figure(4)
%plot(xx1, ro,'r','LineWidth',2)
plot(xx1, q*ro,'r','LineWidth',2)
xlabel('x [um]');
%ylabel('Total Charge Density [1/cm^3]');
ylabel('Total Charge Density [C/cm^3]');
title('Total Charge Density vs Position - at Equilibrium');
%axis([0.5 5 -3e17 8e17])

figure(5)
plot(xx1, Ec,'r','LineWidth',2)
xlabel('x [um]');
%ylabel('Total Charge Density [1/cm^3]');
ylabel('Conduction Band Energy (eV)');
title('Conduction Band vs Position - at Equilibrium');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 END OF EQUILIBRIUM  SOLUTION PART                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                      %%
%%               NON-EQUILIBRIUM  SOLUTION PART BEGINS                  %%
%%                                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              1. Calculate Low filed mobility                         %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Prameters for Low field mobility calculation %%

    TL = 300;                    % Temp in Kelvin
    N  = Na + Nd;                % Local (total) impurity concentration

    MU1N.CAUG   = 55.24;         % cm2/(V.s)
    MU2N.CAUG   = 1429.23;       % cm2/(V.s)
    ALPHAN.CAUG = 0.0;           % unitless
    BETAN.CAUG  = -2.3;          % unitless
    GAMMAN.CAUG = -3.8;          % unitless
    DELTAN.CAUG = 0.73;          % unitless
    NCRITN.CAUG = 1.072*10^17;   % cm-3

    MU1P.CAUG   = 49.7;          % cm2/(V.s)
    MU2P.CAUG   = 479.37;        % cm2/(V.s)
    ALPHAP.CAUG = 0.0;           % unitless
    BETAP.CAUG  = -2.2;          % unitless
    GAMMAP.CAUG = 13.7;          % unitless
    DELTAP.CAUG = 0.70;          % unitless
    NCRITP.CAUG = 1.606*10^17;   % cm-3
    BETAN = 2.0;
    BETAP = 1.0;

    
% %     mun0 = ( MU1N.CAUG*((TL/300)^ALPHAN.CAUG) ) ...
% %       + (( (MU2N.CAUG*((TL/300)^BETAN.CAUG)) - (MU1N.CAUG*((TL/300)^ALPHAN.CAUG)) ) ... 
% %            / ( 1 + ((TL/300)^GAMMAN.CAUG) * ((N/NCRITN.CAUG)^DELTAN.CAUG) ))
% %     
% %     
% %     mup0 = ( MU1P.CAUG*((TL/300)^ALPHAP.CAUG) ) ... 
% %       + (( (MU2P.CAUG*((TL/300)^BETAP.CAUG)) - (MU1P.CAUG*((TL/300)^ALPHAP.CAUG)) ) ... 
% %            / ( 1 + ((TL/300)^GAMMAP.CAUG) * ((N/NCRITP.CAUG)^DELTAP.CAUG) ))
    
    VSATN = (2.4*10^7) / (1 + 0.8*exp(TL/600));  % Saturation Velocity of Electrons
    VSATP = VSATN;                               % Saturation Velocity of Holes

%%%%%%%%%%%%%%%%%%% END of Low Field Mobility Calculation %%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   2. Start the main Loop to increment the Anode voltage by Vt=KbT/q  %% 
%%      till it reaches 0.625V.                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vindex=0;

for VA = 0:0.33*Vt:0.625                % Start VA increment loop
    VA 
    
    Each_Step   = 0.33*Vt 
    Total_Steps = 0.625/(0.33*Vt)
    vindex = vindex +1
    Vplot(vindex) = VA;
    
    fi(1) = fi(1) + VA;            % Apply potential to Anode (1st node)  
    %fi(1)
    
    flag_conv2 = 0;		           % Convergence of the Poisson loop
    k_itern= 0;
    
    %% Initialize the First and Last Node for Poisson's eqn
    a(1) = 0;
    c(1) = 0;
    b(1) = 1;
    f(1) = fi(1);
    a(n_max) = 0;
    c(n_max) = 0;
    b(n_max) = 1;
    f(n_max) = fi(n_max);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3. Start the Poisson equation solver loop to calculate the        %%
    %%    potential for each Anode voltage increase                      %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    while(~flag_conv2)             % Start Poisson's eqn
        
        k_itern = k_itern + 1

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% 3.1 . Calculate Field Dependant Mobility for each value of 'fi'   %% 
        %%       at each node point of the PN diode.                         %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%%% To test with Constant Mobility without field dependancy.        
% %         for i = 1:n_max           % Start Loop for Field Dep Mobility 
% %             mup(i) = mup0;
% %             mun(i) = mun0;
% %         end   

        %% Calculate the Electric Field at each Node
        for i = 2:n_max-1
            Ef(i) = abs(fi(i) - fi(i+1))*Vt/(dx*Ldi);
        end
            
            Ef(1)     = Ef(2); 
            Ef(n_max) = Ef(n_max-1);

        %% Calculate the Field Dependant Mobility at each Node
        for i = 1:n_max           

              pdeno  = (mup0 * Ef(i) / VSATP) ^ BETAP;
              mup(i)   = mup0 * ( (1/(1 + pdeno)) ^(1/BETAP)); 
                            
              ndeno  = (mun0 * Ef(i) / VSATN) ^ BETAN;
              mun(i) = mun0 * ( (1/(1 + ndeno)) ^(1/BETAN)); 
                            
        end   
        
        mup(1)     = mup(2);
        mup(n_max) = mup(n_max-1);
        
        mun(1)     = mun(2);
        mun(n_max) = mun(n_max-1);
        
        %%%%%%%%%%% END of FIELD Dependant Mobility Calculation %%%%%%%%%%% 
      
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %% 3.2 Solve Continuity Equation for Electron and Holes using LU Decomposition %%                                    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %(A) Define the elements of the coefficient matrix and initialize the forcing
        %    function at the ohmic contacts for ELECTRON and HOLE Continuity Eqns 


        an(1) = 0;              %Co-ef for electron at Anode
        bn(1) = 1;              %Co-ef for electron at Anode
        cn(1) = 0;              %Co-ef for electron at Anode
        ap(1) = 0;              %Co-ef for hole     at Anode
        bp(1) = 1;              %Co-ef for hole     at Anode
        cp(1) = 0;              %Co-ef for hole     at Anode
        %fnp(1) = (Ldi*Ldi*dx2/Vt) * ( p(1)*n(1) - 1 ) / ( TAUP0*(n(1) + 1 ) + TAUN0*(p(1) + 1 ) );
        fn(1) = n(1);
        fp(1) = p(1);
        
        an(n_max) = 0;          %Co-ef for electron at Cathode
        bn(n_max) = 1;          %Co-ef for electron at Cathode
        cn(n_max) = 0;          %Co-ef for electron at Cathode
        ap(n_max) = 0;          %Co-ef for hole     at Cathode
        bp(n_max) = 1;          %Co-ef for hole     at Cathode
        cp(n_max) = 0;          %Co-ef for hole     at Cathode
        %fnp(n_max) = (Ldi*Ldi*dx2/Vt) * ( p(n_max)*n(n_max) - 1 ) / ( TAUP0*(n(n_max) + 1) + TAUN0*(p(n_max) + 1) );
        fn(n_max) = n(n_max);
        fp(n_max) = p(n_max);


        %(B) Define the elements of the coefficient matrix for the internal nodes and
        %    initialize the forcing function


        for i = 2: n_max-1
            munim1by2 = (mun(i-1)+mun(i))/2;         
            munip1by2 = (mun(i)+mun(i+1))/2;         
            mupim1by2 = (mup(i-1)+mup(i))/2;         
            mupip1by2 = (mup(i)+mup(i+1))/2;
            
            %% Co-efficients for HOLE Continuity eqn
                cp(i) = mupip1by2 * BER(fi(i) - fi(i+1));
                ap(i) = mupim1by2 * BER(fi(i) - fi(i-1));
           bp(i) = -( mupim1by2 * BER(fi(i-1) - fi(i)) + mupip1by2 * BER(fi(i+1) - fi(i)));
            %% Co-efficients for ELECTRON Continuity eqn
                cn(i) = munip1by2 * BER(fi(i+1) - fi(i));
                an(i) = munim1by2 * BER(fi(i-1) - fi(i));
           n(i) = -( munim1by2 * BER(fi(i) - fi(i-1)) + munip1by2 * BER(fi(i) - fi(i+1)));

            %% Forcing Function for ELECTRON and HOLE Continuity eqns
       fn(i) = (Ldi*Ldi*dx2/Vt) * ( p(i)*n(i) - 1 ) / ( TAUP0*(n(i) + 1) + TAUN0*(p(i)+1));
       fp(i) = (Ldi*Ldi*dx2/Vt) * ( p(i)*n(i) - 1 ) / ( TAUP0*(n(i) + 1) + TAUN0*(p(i)+1));
        end

        
        %(C)  Start the iterative procedure for the solution of the linearized Continuity
        %     equation for "ELECTRONS" using LU decomposition method:


            alphan(1) = bn(1);
            for i=2:n_max
                betan(i)=an(i)/alphan(i-1);
                alphan(i)=bn(i)-betan(i)*cn(i-1);
            end

        % Solution of Lv = f %    

            vn(1) = fn(1);
            for i = 2:n_max
                vn(i) = fn(i) - betan(i)*vn(i-1);
            end

        % Solution of U*fi = v %    

            tempn = vn(n_max)/alphan(n_max);
            %deltan(n_max) = tempn - n(n_max);    
            n(n_max)=tempn;
            for i = (n_max-1):-1:1       %delta%
                tempn = (vn(i)-cn(i)*n(i+1))/alphan(i);
              %  deltan(i) = tempn - n(i);   
                n(i) = tempn;
            end
       
       %%%%%%%%%%%%%%%%%%%%%%% END of ELECTRON Continuty Solver %%%%%%%%%%%  

            
        %(D)  Start the iterative procedure for the solution of the linearized Continuity
        %     equation for "HOLES" using LU decomposition method:


            alphap(1) = bp(1);
            for i=2:n_max
                betap(i)=ap(i)/alphap(i-1);
                alphap(i)=bp(i)-betap(i)*cp(i-1);
            end

        % Solution of Lv = f %    

            vp(1) = fp(1);
            for i = 2:n_max
                vp(i) = fp(i) - betap(i)*vp(i-1);
            end

        % Solution of U*fi = v %    

            tempp = vp(n_max)/alphap(n_max);
            %deltap(n_max) = tempp - p(n_max);    
            p(n_max)=tempp;
            for i = (n_max-1):-1:1       %delta%
                tempp = (vp(i)-cp(i)*p(i+1))/alphap(i);
             %   deltap(i) = tempp - p(i);   
                p(i) = tempp;
            end
    
       %%%%%%%%%%%%%%%%%%%%%%% END of HOLE Continuty Solver %%%%%%%%%%%  

       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       %% 3.3 Calculate potential fi again with new values of "n" and "p"%%
       %%     and check for convergence                                  %%     
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       % Recalculate forcing function and central coefficient b for fi 

        for i = 2: n_max-1
            b(i) = -(2/dx2 + n(i) + p(i));
            f(i) = n(i) - p(i) - dop(i) - (fi(i)*(n(i) + p(i)));  
%% here values of n(i) and p(i) are used in place of exp(fi(i))
        end

        
        % Solve for Updated potential given the new value of Forcing 
        % Function using LU decomposition 
        
        alpha(1) = b(1);
        for i=2:n_max
            beta(i)=a(i)/alpha(i-1);
            alpha(i)=b(i)-beta(i)*c(i-1);
        end

        % Solution of Lv = f %

        v(1) = f(1);
        for i = 2:n_max
            v(i) = f(i) - beta(i)*v(i-1);
        end

        % Solution of U*fi = v %

        temp = v(n_max)/alpha(n_max);
        delta(n_max) = temp - fi(n_max);
        fi(n_max)=temp;
        for i = (n_max-1):-1:1       %delta%
            temp = (v(i)-c(i)*fi(i+1))/alpha(i);
            delta(i) = temp - fi(i);
            fi(i) = temp;
        end

        delta_max = 0;

        for i = 1: n_max
            xx = abs(delta(i));
            if(xx > delta_max)
                delta_max=xx;        %% Calculate the max error
            end
        end

        %delta_max=max(abs(delta));

        % Test convergence and start the loop if necessary else increase
        % the applied potential
        %% delta_max
        if(delta_max < delta_acc)
            flag_conv2 = 1;
        end

    end   % End of WHILE Loop for Poisson's eqn solver

  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                        CALCULATE CURRENT                             %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          % Electron Current       
    for i=2:n_max-1
        Jnim1by2(vindex,i) = (q*mun(i)*Vt/(dx*Ldi)) * ni*( n(i)*BER(fi(i)-fi(i-1)) - n(i-1)*BER(fi(i-1)-fi(i)) );
        Jnip1by2(vindex,i) = (q*mun(i)*Vt/(dx*Ldi)) * ni*( n(i+1)*BER(fi(i+1)-fi(i))  - n(i)*BER(fi(i)-fi(i+1)) );
        Jelec(vindex,i) = (Jnip1by2(vindex,i) + Jnim1by2(vindex,i))/2;
        
% %         % Electron Current with only one node
% %         Jnim1by2 = (q*mun(2)*Vt/dx*Ldi) * ( (n(2)*ni)*BER((fi(2)-fi(2-1))) 
% %                        - (n(2-1)*ni)*BER((fi(2-1)-fi(2))) );
% %         Jnip1by2 = (q*mun(2)*Vt/dx*Ldi) * ( (n(2+1)*ni)*BER((fi(2+1)-fi(2))) 
% %                        - (n(2)*ni)*BER((fi(2)-fi(2+1))) );
% %         Jelec(vindex) = (Jnip1by2 + Jnim1by2)/2;

          % Hole Current   
% %         Jpim1by2 = (q*mup(2)*Vt/dx*Ldi) * ( (p(2)*ni)*BER((fi(2-1)-fi(2))) 
% %                        - (p(2-1)*ni)*BER((fi(2)-fi(2-1))) );
% %         Jpip1by2 = (q*mup(2)*Vt/dx*Ldi) * ( (p(2+1)*ni)*BER((fi(2)-fi(2+1))) 
% %                        - (p(2)*ni)*BER((fi(2+1)-fi(2))) );
% %         Jhole(vindex) = (Jpip1by2 + Jpim1by2)/2;

        Jpim1by2(vindex,i) = (q*mup(i)*Vt/(dx*Ldi)) * ni*( p(i)*BER((fi(i-1)-fi(i))) - p(i-1)*BER((fi(i)-fi(i-1))) );
        Jpip1by2(vindex,i) = (q*mup(i)*Vt/(dx*Ldi)) * ni*( p(i+1)*BER((fi(i)-fi(i+1))) - p(i)*BER((fi(i+1)-fi(i))) );
        Jhole(vindex,i) = (Jpip1by2(vindex,i) + Jpim1by2(vindex,i))/2;

    end

        Jtotal = Jelec + Jhole;
%%         Jtotal(vindex) = Jelec;
        
        
end   % End of main FOR loop for VA increment.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 END OF NON-EQUILIBRIUM  SOLUTION PART                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Write the results of the simulation in files %

xx1(1) = dx*1e4;
for i = 2:n_max-1 
    Ec(i) = dEc - Vt*fi(i);     %Values from the second Node%
    ro(i) = -ni*(n(i) - p(i) - dop(i));
    el_field1(i) = -(fi(i+1) - fi(i))*Vt/(dx*Ldi);
    el_field2(i) = -(fi(i+1) - fi(i-1))*Vt/(2*dx*Ldi);
    xx1(i) = xx1(i-1) + dx*Ldi*1e4;
end

Jtotal(:,1) = Jtotal(:,2);
Jelec(:,1) = Jelec(:,2);
Jhole(:,1) = Jhole(:,2);
Jtotal(:,n_max) = Jtotal(:,(n_max-1));
Jelec(:,n_max) = Jelec(:,(n_max-1));
Jhole(:,n_max) = Jhole(:,(n_max-1));

Ec(1) = Ec(2);
Ec(n_max) = Ec(n_max-1);
xx1(n_max) = xx1(n_max-1) + dx*Ldi*1e4;
el_field1(1) = el_field1(2);
el_field2(1) = el_field2(2);
el_field1(n_max) = el_field1(n_max-1);
el_field2(n_max) = el_field2(n_max-1);
nf = n*ni;
pf = p*ni;
ro(1) = ro(2);
ro(n_max) = ro(n_max-1);

%% Calculate Quasi Fermi Level - Efn Efp
for i = 1:n_max
    Ei(i)   = Ec(i) - 0.56;
    Efn(i)  = Ei(i) + Vt*log(nf(i)/ni);
    Efp(i)  = Ei(i) - Vt*log(pf(i)/ni);
end
    Ev = Ec - 1.12;

figure(14)
plot(xx1, Ec,'black','LineWidth',2.5);
hold on;
plot(xx1, Ev,'black','LineWidth',2.5);
hold on;
plot(xx1, Ei,'--black','LineWidth',2.5);
hold on;
plot(xx1, Efn,'r','LineWidth',2.5);
hold on;
plot(xx1, Efp,'b','LineWidth',2.5);
xlabel('x [um]');
ylabel('Energy [eV]');
title('Quasi Fermi Levels (Efn & Efp) vs Position - at Applied Bias(0.625V)');
legend('Ec','Ev','Ei','Efn','Efp');
axis([0 7 -1 1]);


figure(6)
plot(xx1, Ec,'b','LineWidth',2)
xlabel('x [um]');
ylabel('Conduction Band Energy (eV)');
title('Conduction Band vs Position - at Applied Bias (0.625)');

figure(7)
plot(xx1, Vt*fi,'b','LineWidth',2)
xlabel('x [um]');
ylabel('Potential [eV]');
title('Potential vs Position - at Applied Bias(0.625V)');

figure(8)
plot(xx1, el_field1,'b','LineWidth',2)
hold on;
plot(xx1, el_field2,'b','LineWidth',2)
xlabel('x [um]');
ylabel('Electric Field [V/cm]');
title('Field Profile vs Position - at Applied Bias(0.625V)');

figure(9)
%plot(xx1, nf,'g','LineWidth',2)
semilogy(xx1, nf,'g','LineWidth',2)
hold on;
%plot(xx1, pf,'b','LineWidth',2)
semilogy(xx1, pf,'b','LineWidth',2)
xlabel('x [um]');
ylabel('Electron & Hole Densities [1/cm^3]');
title('Electron & Hole Densities vs Position - at Applied Bias(0.625V)');
legend('n','p');
%axis([0 6.75 0 10.2e17])

figure(10)
%plot(xx1, ro,'b','LineWidth',2)
plot(xx1, q*ro,'b','LineWidth',2)
xlabel('x [um]');
%ylabel('Total Charge Density [1/cm^3]');
ylabel('Total Charge Density [C/cm^3]');
title('Total Charge Density vs Position - at Applied Bias(0.625V)');
%axis([0.5 5 -3e17 8e17])


figure(11)
plot(Vplot, Jtotal(:,2),'r','LineWidth',2)
hold on
plot(Vplot, Jhole(:,2),'g','LineWidth',2)
hold on
plot(Vplot, Jelec(:,2),'b','LineWidth',2)
xlabel('VA [V]');
ylabel('Total Current Density [Amp/cm^2]');
title('I vs V Plot');
legend('Jtotal','Jhole','Jelec','2');


figure(12)
plot(Vplot, Jtotal(:,2),'r','LineWidth',2)
xlabel('VA [V]');
ylabel('Total Current Density [Amp/cm^2]');
title('I vs V Plot');
%legend('Jtotal','Jhole','Jelec','2');

figure(13)
plot(xx1,Jtotal((round((Total_Steps)-1)),:),'b','LineWidth',2)
xlabel('x [um]');
ylabel('Total Current Density [A/cm^2]');
title('Total Current Density vs Position - at Applied Bias(0.625V)');
axis([0 7 0 6]);


%figure(5)
%plot(xx1, n)
%hold all
%plot(xx1, p)

save cond_band xx1 Ec;
save tot_charge xx1 ro;
save el_field xx1 el_field1 el_field2;
save np_data xx1 nf pf;
save pot_1 fi;
