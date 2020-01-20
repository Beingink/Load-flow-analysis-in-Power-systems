% Calling the excel sheets and getting the line and bus datas 
clc;
linedata = xlsread('linedata');
busdata = xlsread('busdata');

j=sqrt(-1);  i = sqrt(-1);

fb = linedata(:,1);      % From bus number
tb = linedata(:,2);      % To bus number
R = linedata(:,3);       % Resistance, R
X = linedata(:,4);       % Reactance, X...
b = j*linedata(:,5);     % Shunt Admittance, B/2 

Z = R + j*X;             % Z matrix
y= 1./Z;   %branch admittance

nbranch=length(linedata(:,1)); % no. of branches
nbus = max(max(fb), max(tb));  % no. of buses

%  Forming the Y Bus Matrix

for n = 1:nbranch
Ybus=zeros(nbus,nbus);     % initialize Ybus to zero
               
% Formation of the off diagonal elements
for k=1:nbranch;
       Ybus(fb(k),tb(k))=Ybus(fb(k),tb(k))-y(k);
       Ybus(tb(k),fb(k))=Ybus(fb(k),tb(k));
    end
end

% Formation of the diagonal elements
for  n=1:nbus
     for k=1:nbranch
         if fb(k)==n | tb(k)==n
         Ybus(n,n) = Ybus(n,n)+y(k) + b(k);
            else, end
     end
end
fprintf('\n\n\t\t\t --------[POWER FLOW SOLUTIONS]---------');
fprintf('\n\n\t\t\t             <--METHODS-->');
fprintf('\n\n\t *1) GAUSS SIDEL METHOD  \t   *2) NEWTON RAPHSON METHOD \t *3) FAST DECOUPLED METHOD');
fprintf('\n\n\t          [Tolerance: 0.0001]     [Acceleration Factor:1.6]')


casi=input('\n\n   ENTER THE METHOD YOU PREFER :');

switch(casi)
    
    case 1 
% Gauss-Seidel method

basemva = 100;        %Base MVA 
tolerance = 0.0001;   %Tolerance 
mi = 80;              %Maximum Iterations
af = 1.6;             %Acceleration factor

% Keys for check purposes
Vm=0; delta=0; yload=0; deltad =0;

nbus = length(busdata(:,1));
for k=1:nbus
n=busdata(k,1);
bt(n)=busdata(k,2); 
Vm(n)=busdata(k,3); 
delta(n)=busdata(k, 4);
Pd(n)=busdata(k,5); 
Qd(n)=busdata(k,6); 
Pg(n)=busdata(k,7); 
Qg(n) = busdata(k,8);

    if Vm(n) <= 0  
            Vm(n) = 1.0; V(n) = 1 + j*0;
    else delta(n) = pi/180*delta(n);
         V(n) = Vm(n)*(cos(delta(n)) + j*sin(delta(n)));
         P(n)=(Pg(n)-Pd(n))/basemva;
         Q(n)=(Qg(n)-Qd(n))/basemva;
         S(n) = P(n) + j*Q(n);
    end
DV(n)=0;
end
num = 0; AcurBus = 0; converge = 1;
Vc = zeros(nbus,1)+j*zeros(nbus,1); Sc = zeros(nbus,1)+j*zeros(nbus,1);
iter=0;
maxerror=10;
sumc1=0;
sumc2=0;
sumc3=0;
sumc4=0;
while maxerror >= tolerance & iter <= mi
iter=iter+1;
  for n = 1:nbus;
  YV = 0+j*0;
    for L = 1:nbranch;
            if fb(L) == n, k=tb(L);
            YV = YV + Ybus(n,k)*V(k);
            elseif tb(L) == n, k=fb(L);
            YV = YV + Ybus(n,k)*V(k);
            end
    end
       Sc = conj(V(n))*(Ybus(n,n)*V(n) + YV) ;
       Sc = conj(Sc);
       DP(n) = P(n) - real(Sc);
       DQ(n) = Q(n) - imag(Sc);
         if bt(n) == 1
         S(n) =Sc; P(n) = real(Sc); Q(n) = imag(Sc); DP(n) =0; DQ(n)=0;
         Vc(n) = V(n);
         elseif bt(n) == 2
         Q(n) = imag(Sc); S(n) = P(n) + j*Q(n);

           
             Qgc = Q(n)*basemva + Qd(n) ;
             
         end
       if bt(n) ~= 1
       Vc(n) = (conj(S(n))/conj(V(n)) - YV )/ Ybus(n,n);
       else, end
          if bt(n) == 0
          V(n) = V(n) + af*(Vc(n)-V(n));
          elseif bt(n) == 2
          VcI = imag(Vc(n));
          VcR = sqrt(Vm(n)^2 - VcI^2);
          Vc(n) = VcR + j*VcI;
           V(n) = V(n) + af*(Vc(n) -V(n));
          end
   end
  maxerror=max( max(abs(real(DP))), max(abs(imag(DQ))) );
   if iter == mi & maxerror > tolerance
   fprintf('\nWARNING: Iterative solution did not converged after ')
   fprintf('%g', iter), fprintf(' iterations.\n\n')
   fprintf('Press Enter to terminate the iterations and print the results \n')
   converge = 0; pause, else, end
   
end
if converge ~= 1
   tech= ('                      ITERATIVE SOLUTION DID NOT CONVERGE'); else, 
   tech=('                   Power Flow Solution - Gauss-Seidel Method');
end   
k=0;
for n = 1:nbus
  Vm(n) = abs(V(n)); deltad(n) = angle(V(n))*180/pi;
     if bt(n) == 1
     S(n)=P(n)+j*Q(n);
     Pg(n) = P(n)*basemva + Pd(n);
     Qg(n) = Q(n)*basemva + Qd(n) ;
     k=k+1;
     Pgg(k)=Pg(n);
     elseif  bt(n) ==2
     k=k+1;
     Pgg(k)=Pg(n);
     S(n)=P(n)+j*Q(n);
     Qg(n) = Q(n)*basemva + Qd(n) ;
     end
        sumc1 = sumc1 + Pg(n);  
       sumc2 = sumc2+ Qg(n); 
       sumc3 = sumc3 + Pd(n);  
      sumc4 = sumc4 + Qd(n);
     yload(n) = (Pd(n)- j*Qd(n))/(basemva*Vm(n)^2);
end

busdata(:,3)=Vm'; busdata(:,4)=deltad';



fprintf('\n\n\t\t\t\t\t GAUSS SIDEL SOLUTION');

fprintf('\n\n\t 1) Y BUS');
fprintf('\n\t 2) LINE FLOW SOLUTION');
fprintf('\n\t 3) LINE LOSSES SOLUTION');
fprintf('\n\t 4) EXIT');

opt=input('\n\n Choose your option : ');

if(opt==1)
    
%  DISPLAYING Y BUS
fprintf('                               Y BUS  \n\n')

display(Ybus);

end
if (opt==2)
    
%  DISPLAYING POWER FLOW SOLUTIONS UPTO A VALUE OF 3 DECIMAL PLACES

disp(tech)
fprintf('                               %g Iterations  \n\n', iter)
head =['    Bus   Voltage   Angle     ------Load------     ---Generation--- '
       '    No.   Mag.      Degree      MW       Mvar        MW       Mvar  '
       '                                                                    '];
disp(head)
for n=1:nbus
     fprintf(' |%5g', n), fprintf(' |%7.3f', Vm(n)),
     fprintf(' |%8.3f', deltad(n)), fprintf(' |%9.3f', Pd(n)),
     fprintf(' |%9.3f', Qd(n)),  fprintf(' |%9.3f', Pg(n)),
     fprintf(' |%9.3f \n', Qg(n)), 
end
    fprintf('      \n'), fprintf('    Total              ')
    fprintf('    %9.3f', sumc3), fprintf('    %9.3f', sumc4),
    fprintf('   %9.3f', sumc1), fprintf('  %9.3f\n\n',sumc2),
end


if(opt==3)

% CALCULATING LINE FLOW LOSSES 

SLT = 0;
fprintf('\n')
fprintf('                           Line Flow and Losses \n\n')
fprintf('     --Line--    Power at bus & line flow    --Line loss-- \n')
fprintf('     from  to    MW      Mvar     MVA         MW      Mvar   \n')

for n = 1:nbus
busprt = 0;
   for L = 1:nbranch;
       if busprt == 0
       fprintf('   \n'), fprintf('%6g', n), fprintf('      %9.3f', P(n)*basemva)
       fprintf('%9.3f', Q(n)*basemva), fprintf('%9.3f\n', abs(S(n)*basemva))

       busprt = 1;
       else, end
       if fb(L)==n      k = tb(L);
       In = (V(n) - V(k))*y(L) + b(L);
       Ik = (V(k) - V(n))*y(L) + b(L)*V(k);
       Snk = V(n)*conj(In)*basemva;
       Skn = V(k)*conj(Ik)*basemva;
       SL  = Snk + Skn;
       SLT = SLT + SL;
       elseif tb(L)==n  k = fb(L);
       In = (V(n) - V(k))*y(L) + b(L)*V(n);
       Ik = (V(k) - V(n))*y(L) + b(L);
       Snk = V(n)*conj(In)*basemva;
       Skn = V(k)*conj(Ik)*basemva;
       SL  = Snk + Skn;
       SLT = SLT + SL;
       else, end
         if fb(L)==n | tb(L)==n
         fprintf('%12g', k),
         fprintf('%9.3f', real(Snk)), fprintf('%9.3f', imag(Snk))
         fprintf('%9.3f', abs(Snk)),
         fprintf('%9.3f', real(SL)),
             if fb(L) ==n 
             fprintf('%9.3f \n', imag(SL))
             else, fprintf('%9.3f\n', imag(SL))
             end
         else, end
  end
end
SLT = SLT/2;
fprintf('   \n'), fprintf('    Total loss                         ')
fprintf('%9.3f', real(SLT)), fprintf('%9.3f\n', imag(SLT))
clear Ik In SL SLT Skn Snk

end

if(opt==4)
fprintf('\n\n\t    Have a good day! ');
end



break

    case 2

%  Newton-Raphson method

basemva = 100;        %Base MVA 
tolerance = 0.0001;   %Tolerance 
mi = 80;              %Maximum Iterations

% Keys for check purposes
ns=0; ng=0; Vm=0; delta=0; yload=0; deltad=0;

nbus = length(busdata(:,1));
for k=1:nbus
n=busdata(k,1);             % Bus Number
bt(n)=busdata(k,2);         % Bus type
Vm(n)=busdata(k,3);         % Magnitude of bus voltage
delta(n)=busdata(k, 4);     % Bus voltage Angle
Pd(n)=busdata(k,5);         % Power Demand
Qd(n)=busdata(k,6);         % Reactive Power Demand  
Pg(n)=busdata(k,7);         % Power generated
Qg(n) = busdata(k,8);       % Reactive Power Generated

% Making Flat Bus Voltage Profile
    
     if Vm(n) <= 0  
        Vm(n) = 1.0; 
         V(n) = 1 + j*0;  %Rectangular Form
     else
         delta(n) = pi/180*delta(n);
         
         V(n) = Vm(n)*(cos(delta(n)) + j*sin(delta(n))); %Polar form
         
% Converting Powers into Per Unit        
         
         P(n)=(Pg(n)-Pd(n))/basemva;
         Q(n)=(Qg(n)-Qd(n))/basemva;
         
         S(n) = P(n) + j*Q(n);
    end
end

% Identifying different Bus types

for k=1:nbus
if bt(k) == 1, ns = ns+1; else, end
if bt(k) == 2  ng = ng+1; else, end

ngs(k) = ng;
nss(k) = ns;
end

% Collecting datas from the previously obtained Y Bus

Ym=abs(Ybus); t = angle(Ybus);
m=2*nbus-ng-2*ns; % Segmenting the iteration factor to distinct the bus types

maxerror = 1; 
converge=1;
iter = 0;

% Start of iterations

clear A  DC   J  DX

while maxerror >= tolerance & iter <= mi         % Test for maximum power mismatch
for i=1:m
for k=1:m
   A(i,k)=0;      %Initializing Jacobian matrix
end, end
iter = iter+1;
for n=1:nbus
nn=n-nss(n);
lm=nbus+n-ngs(n)-nss(n)-ns;

J11=0; J22=0; J33=0; J44=0; % The jacobian matrix elements

for i=1:nbranch
     if fb(i) == n | tb(i) == n
        if fb(i) == n,  l = tb(i); end
        if tb(i) == n,  l = fb(i); end
        J11=J11+ Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
        J33=J33+ Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
        if bt(n)~=1
        J22=J22+ Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
        J44=J44+ Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
        else, end
        if bt(n) ~= 1  & bt(l) ~=1
        lk = nbus+l-ngs(l)-nss(l)-ns;
        ll = l -nss(l);
      % off diagonalelements of J1
        A(nn, ll) =-Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
              if bt(l) == 0  % off diagonal elements of J2
              A(nn, lk) =Vm(n)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));end
              if bt(n) == 0  % off diagonal elements of J3
              A(lm, ll) =-Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n)+delta(l)); end
              if bt(n) == 0 & bt(l) == 0  % off diagonal elements of  J4
              A(lm, lk) =-Vm(n)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));end
        else end
     else , end
   end
   Pk = Vm(n)^2*Ym(n,n)*cos(t(n,n))+J33;
   Qk = -Vm(n)^2*Ym(n,n)*sin(t(n,n))-J11;
   if bt(n) == 1 P(n)=Pk; Q(n) = Qk; end   % Reference bus P
     if bt(n) == 2  Q(n)=Qk;
        
           Qgc = Q(n)*basemva + Qd(n) -  (n);     end
   if bt(n) ~= 1
     A(nn,nn) = J11;  %diagonal elements of J1
     
     DC(nn) = P(n)-Pk; %Final power
   end
   if bt(n) == 0
     A(nn,lm) = 2*Vm(n)*Ym(n,n)*cos(t(n,n))+J22;  %diagonal elements of J2
     A(lm,nn)= J33;        %diagonal elements of J3
     A(lm,lm) =-2*Vm(n)*Ym(n,n)*sin(t(n,n))-J44;  %diagonal of elements of J4
     
     DC(lm) = Q(n)-Qk; % Final Reactive Power
   end
end
%Taking inverse

DX=A\DC';

for n=1:nbus
  nn=n-nss(n);
  lm=nbus+n-ngs(n)-nss(n)-ns;
    if bt(n) ~= 1
    delta(n) = delta(n)+DX(nn); end
    if bt(n) == 0
    Vm(n)=Vm(n)+DX(lm); end
 end
  maxerror=max(abs(DC));
     if iter == mi & maxerror > tolerance 
   fprintf('\nWARNING: Iterative solution did not converged after ')
   fprintf('%g', iter), fprintf(' iterations.\n\n')
   fprintf('Press Enter to terminate the iterations and print the results \n')
   converge = 0; pause, else, end
   
end

if converge ~= 1
   tech= ('                      ITERATIVE SOLUTION DID NOT CONVERGE'); else, 
   tech=('                   Power Flow Solution - Newton-Raphson Method');
end   
V = Vm.*cos(delta)+j*Vm.*sin(delta);
deltad=180/pi*delta;
i=sqrt(-1);
k=0;
sumc1=0;
sumc2=0;
sumc3=0;
sumc4=0;

for n = 1:nbus
     if bt(n) == 1
     k=k+1;
     S(n)= P(n)+j*Q(n);
     Pg(n) = P(n)*basemva + Pd(n);
     Qg(n) = Q(n)*basemva + Qd(n);
     Pgg(k)=Pg(n);
     Qgg(k)=Qg(n); 
  
     elseif  bt(n) ==2
     k=k+1;
     S(n)=P(n)+j*Q(n);
     Qg(n) = Q(n)*basemva + Qd(n);
     Pgg(k)=Pg(n);
     Qgg(k)=Qg(n); 
    
     end
   sumc1 = sumc1 + Pg(n);  
    sumc2 = sumc2+ Qg(n); 
       sumc3 = sumc3 + Pd(n);  
     sumc4 = sumc4 + Qd(n);
yload(n) = (Pd(n)- j*Qd(n))/(basemva*Vm(n)^2);
end
busdata(:,3)=Vm'; busdata(:,4)=deltad';




fprintf('\n\n\t\t\t\t\t NEWTON RAPHSON SOLUTION');

fprintf('\n\n\t 1) Y BUS');
fprintf('\n\t 2) LINE FLOW SOLUTION');
fprintf('\n\t 3) LINE LOSSES SOLUTION');
fprintf('\n\t 4) EXIT');

opt=input('\n\n Choose your option : ');

if(opt==1)

%  DISPLAYING Y BUS
fprintf('                               Y BUS  \n\n')

display(Ybus);

end

if(opt==2)
%  DISPLAYING POWER FLOW SOLUTIONS UPTO A VALUE OF 3 DECIMAL PLACES

disp(tech)
fprintf('                               %g Iterations  \n\n', iter)
head =['    Bus   Voltage   Angle     ------Load------     ---Generation--- '
       '    No.   Mag.      Degree      MW       Mvar        MW       Mvar  '
       '                                                                    '];
disp(head)
for n=1:nbus
     fprintf(' |%5g', n), fprintf(' |%7.3f', Vm(n)),
     fprintf(' |%8.3f', deltad(n)), fprintf(' |%9.3f', Pd(n)),
     fprintf(' |%9.3f', Qd(n)),  fprintf(' |%9.3f', Pg(n)),
     fprintf(' |%9.3f \n', Qg(n)), 
end
    fprintf('      \n'), fprintf('    Total              ')
    fprintf('    %9.3f', sumc3), fprintf('    %9.3f', sumc4),
    fprintf('   %9.3f', sumc1), fprintf('  %9.3f\n\n',sumc2),
end


if(opt==3)

% CALCULATING LINE FLOW LOSSES 

SLT = 0;
fprintf('\n')
fprintf('                           Line Flow and Losses \n\n')
fprintf('     --Line--    Power at bus & line flow    --Line loss-- \n')
fprintf('     from  to    MW      Mvar     MVA         MW      Mvar   \n')

for n = 1:nbus
busprt = 0;
   for L = 1:nbranch;
       if busprt == 0
       fprintf('   \n'), fprintf('%6g', n), fprintf('      %9.3f', P(n)*basemva)
       fprintf('%9.3f', Q(n)*basemva), fprintf('%9.3f\n', abs(S(n)*basemva))

       busprt = 1;
       else, end
       if fb(L)==n      k = tb(L);
       In = (V(n) - V(k))*y(L) + b(L);
       Ik = (V(k) - V(n))*y(L) + b(L)*V(k);
       Snk = V(n)*conj(In)*basemva;
       Skn = V(k)*conj(Ik)*basemva;
       SL  = Snk + Skn;
       SLT = SLT + SL;
       elseif tb(L)==n  k = fb(L);
       In = (V(n) - V(k))*y(L) + b(L)*V(n);
       Ik = (V(k) - V(n))*y(L) + b(L);
       Snk = V(n)*conj(In)*basemva;
       Skn = V(k)*conj(Ik)*basemva;
       SL  = Snk + Skn;
       SLT = SLT + SL;
       else, end
         if fb(L)==n | tb(L)==n
         fprintf('%12g', k),
         fprintf('%9.3f', real(Snk)), fprintf('%9.3f', imag(Snk))
         fprintf('%9.3f', abs(Snk)),
         fprintf('%9.3f', real(SL)),
             if fb(L) ==n 
             fprintf('%9.3f \n', imag(SL))
             else, fprintf('%9.3f\n', imag(SL))
             end
         else, end
  end
end
SLT = SLT/2;
fprintf('   \n'), fprintf('    Total loss                         ')
fprintf('%9.3f', real(SLT)), fprintf('%9.3f\n', imag(SLT))
clear Ik In SL SLT Skn Snk
end

if(opt==4)
fprintf('\n\n\t    Have a good day! ');
end

   ;
    break
    
    
    case 3
        %   Fast Decoupled method

ns=0; Vm=0; delta=0; yload=0; deltad=0; % Keys for check purposes

basemva = 100;        %Base MVA 
tolerance = 0.0001;   %Tolerance 
mi = 80;              %Maximum Iterations
af = 1.6;             %Acceleration factor


nbus = length(busdata(:,1));
for k=1:nbus
n=busdata(k,1);
bt(n)=busdata(k,2); 
Vm(n)=busdata(k,3); 
delta(n)=busdata(k,4);
Pd(n)=busdata(k,5); 
Qd(n)=busdata(k,6); 
Pg(n)=busdata(k,7); 
Qg(n)= busdata(k,8);

    if Vm(n) <= 0  Vm(n) = 1.0; V(n) = 1 + j*0;
    else delta(n) = pi/180*delta(n);
         V(n) = Vm(n)*(cos(delta(n)) + j*sin(delta(n)));
         P(n)=(Pg(n)-Pd(n))/basemva;
         Q(n)=(Qg(n)-Qd(n))/basemva;
         S(n) = P(n) + j*Q(n);
    end
if bt(n) == 1, ns = ns+1; else, end
nss(n) = ns;
end
Ym = abs(Ybus); t = angle(Ybus);
ii=0;
for ib=1:nbus
     if bt(ib) == 0 | bt(ib) == 2
     ii = ii+1;
      jj=0;
         for jb=1:nbus
             if bt(jb) == 0 | bt(jb) == 2
             jj = jj+1;
             B1(ii,jj)=imag(Ybus(ib,jb));
             else,end
         end
     else, end
end

ii=0;
for ib=1:nbus
     if bt(ib) == 0
     ii = ii+1;
      jj=0;
         for jb=1:nbus
             if bt(jb) == 0
             jj = jj+1;
             B2(ii,jj)=imag(Ybus(ib,jb));
             else,end
         end
     else, end
end

B1inv=inv(B1); B2inv = inv(B2);
maxerror = 1; converge = 1; 
iter = 0;
sumc1=0;
sumc2=0;
sumc3=0;
sumc4=0;

% Start of iterations
while maxerror >= tolerance & iter <= mi % Test for max. power mismatch
iter = iter+1;
id=0; iv=0;
for n=1:nbus
nn=n-nss(n);
J11=0;  J33=0;
   for i=1:nbranch
     if fb(i) == n | tb(i) == n
        if fb(i) == n,  l = tb(i); end
        if tb(i) == n,  l = fb(i); end
        J11=J11+ Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
        J33=J33+ Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
     else , end
   end
   Pk = Vm(n)^2*Ym(n,n)*cos(t(n,n))+J33;
   Qk = -Vm(n)^2*Ym(n,n)*sin(t(n,n))-J11;
   if bt(n) == 1 P(n)=Pk; Q(n) = Qk; end   % Swing bus P
     if bt(n) == 2  Q(n)=Qk;
         Qgc = Q(n)*basemva + Qd(n);
         
     end
   if bt(n) ~= 1
   id = id+1;
     DP(id) = P(n)-Pk;
     DPV(id) = (P(n)-Pk)/Vm(n);
   end
   if bt(n) == 0
   iv=iv+1;
     DQ(iv) = Q(n)-Qk;
     DQV(iv) = (Q(n)-Qk)/Vm(n);
   end
end
Dd=-B1\DPV';
DV=-B2\DQV';
id=0;iv=0;
  for n=1:nbus
    if bt(n) ~= 1
    id = id+1;
    delta(n) = delta(n)+Dd(id); end
    if bt(n) == 0
    iv = iv+1;
    Vm(n)=Vm(n)+DV(iv); end
  end
    maxerror=max(max(abs(DP)),max(abs(DQ)));
   if iter ==mi & maxerror > tolerance
   fprintf('\nWARNING: Iterative solution did not converged after ')
   fprintf('%g', iter), fprintf(' iterations.\n\n')
   fprintf('Press Enter to terminate the iterations and print the results \n')
   converge = 0; pause, else, end
   
end
if converge ~= 1
   tech= ('                      ITERATIVE SOLUTION DID NOT CONVERGE'); else, 
   tech=('                   Power Flow Solution - Fast Decoupled Method');
end   
k=0;
V = Vm.*cos(delta)+j*Vm.*sin(delta);
deltad=180/pi*delta;
clear A; clear DC; clear DX
i=sqrt(-1);
for n = 1:nbus
     if bt(n) == 1
     S(n)=P(n)+j*Q(n);
     Pg(n) = P(n)*basemva + Pd(n);
     Qg(n) = Q(n)*basemva + Qd(n);
     k=k+1;
     Pgg(k)=Pg(n);
     elseif  bt(n) ==2
     S(n)=P(n)+j*Q(n);
     Qg(n) = Q(n)*basemva + Qd(n); 
     k=k+1;
     Pgg(k)=Pg(n);
     end
     sumc1 = sumc1 + Pg(n);  
    sumc2 = sumc2+ Qg(n); 
       sumc3 = sumc3 + Pd(n);  
     sumc4 = sumc4 + Qd(n);
yload(n) = (Pd(n)- j*Qd(n))/(basemva*Vm(n)^2);
end
busdata(:,3)=Vm'; busdata(:,4)=deltad';

clear Pk Qk  DP DQ J11 J33 B1 B1inv B2 B2inv DPV  DQV Dd delta ib id ii iv jb jj

fprintf('\n\n\t\t\t\t\t FAST DECOUPLED SOLUTION');

fprintf('\n\n\t 1) Y BUS');
fprintf('\n\t 2) LINE FLOW SOLUTION');
fprintf('\n\t 3) LINE LOSSES SOLUTION');
fprintf('\n\t 4) EXIT');

opt=input('\n\n Choose your option : ');

if(opt==1)



%  DISPLAYING Y BUS
fprintf('                               Y BUS  \n\n')

display(Ybus);
end

if(opt==2)

%  DISPLAYING POWER FLOW SOLUTIONS UPTO A VALUE OF 3 DECIMAL PLACES

disp(tech)
fprintf('                               %g Iterations  \n\n', iter)
head =['    Bus   Voltage   Angle     ------Load------     ---Generation--- '
       '    No.   Mag.      Degree      MW       Mvar        MW       Mvar  '
       '                                                                    '];
disp(head)
for n=1:nbus
     fprintf(' |%5g', n), fprintf(' |%7.3f', Vm(n)),
     fprintf(' |%8.3f', deltad(n)), fprintf(' |%9.3f', Pd(n)),
     fprintf(' |%9.3f', Qd(n)),  fprintf(' |%9.3f', Pg(n)),
     fprintf(' |%9.3f \n', Qg(n)), 
end
    fprintf('      \n'), fprintf('    Total              ')
    fprintf('    %9.3f', sumc3), fprintf('    %9.3f', sumc4),
    fprintf('   %9.3f', sumc1), fprintf('  %9.3f\n\n',sumc2),


end

if(opt==3)

% CALCULATING LINE FLOW LOSSES 

SLT = 0;
fprintf('\n')
fprintf('                           Line Flow and Losses \n\n')
fprintf('     --Line--    Power at bus & line flow    --Line loss-- \n')
fprintf('     from  to    MW      Mvar     MVA         MW      Mvar   \n')

for n = 1:nbus
busprt = 0;
   for L = 1:nbranch;
       if busprt == 0
       fprintf('   \n'), fprintf('%6g', n), fprintf('      %9.3f', P(n)*basemva)
       fprintf('%9.3f', Q(n)*basemva), fprintf('%9.3f\n', abs(S(n)*basemva))

       busprt = 1;
       else, end
       if fb(L)==n      k = tb(L);
       In = (V(n) - V(k))*y(L) + b(L);
       Ik = (V(k) - V(n))*y(L) + b(L)*V(k);
       Snk = V(n)*conj(In)*basemva;
       Skn = V(k)*conj(Ik)*basemva;
       SL  = Snk + Skn;
       SLT = SLT + SL;
       elseif tb(L)==n  k = fb(L);
       In = (V(n) - V(k))*y(L) + b(L)*V(n);
       Ik = (V(k) - V(n))*y(L) + b(L);
       Snk = V(n)*conj(In)*basemva;
       Skn = V(k)*conj(Ik)*basemva;
       SL  = Snk + Skn;
       SLT = SLT + SL;
       else, end
         if fb(L)==n | tb(L)==n
         fprintf('%12g', k),
         fprintf('%9.3f', real(Snk)), fprintf('%9.3f', imag(Snk))
         fprintf('%9.3f', abs(Snk)),
         fprintf('%9.3f', real(SL)),
             if fb(L) ==n 
             fprintf('%9.3f \n', imag(SL))
             else, fprintf('%9.3f\n', imag(SL))
             end
         else, end
  end
end
SLT = SLT/2;
fprintf('   \n'), fprintf('    Total loss                         ')
fprintf('%9.3f', real(SLT)), fprintf('%9.3f\n', imag(SLT))
clear Ik In SL SLT Skn Snk
end

if(opt==4)
fprintf('\n\n\t    Have a good day! ');
end

break
end
 