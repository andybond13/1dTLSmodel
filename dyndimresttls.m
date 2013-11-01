%Same as dyndimrest but we solve it with the TLS.
% it is a guitar string in tension 
% for which suddeny damage is added at the center.
% and a brittle damage model Y = Yc is used
clear all;
%faire converger le cas mineaire
% this function gives d in terms of phi
function res = dval(phi,lc)
x = phi/lc;
%cas lin d = phi/lc
if   (x < 0)  res = 0;
elseif (x > 1) res = 1;
else 
%res = phi/lc; %lin
res = 2 * x - x*x ;  %quad
%res = 3 * x - 3*x*x + x*x*x ;  %cubic
%res = x*x;   %sqrt
%if (x<=0.5) res = 2*x^2; else res = -2*x^2 + 4 * x - 1; endif% s shape
endif
endfunction

% this function is the derivative of d with respect to phi
function res = dp(phi,lc)
res = 0;
x = phi/lc;
if   (x >= 0 && x <= 1)  
%res = 1/lc;           %lin
res = 2/lc * (1-x) ;  %quad
%res = 1/lc * (3-6*x+3*x*x); %cubic
%res = 1/lc * 2*x; %sqrt
%if (x<=0.5) res = 4*x/lc; else res = (-4 * x + 4)/lc; endif% s shape
endif
endfunction

%second derivative of d with respect to phi
function res = dpp(phi,lc)
res = 0;
x = phi/lc;
%cas lin d = phi/lc
%cas quad d = 2 * (phi/lc) -(phi/lc)**2 
if   (x >= 0 && x <= 1)  
res = -2/(lc*lc) ; % quad
%res = (1/lc)*(1/lc) * 6 * (-1+x); %cubic
%res = (1/lc)*(1/lc) * 2; %sqrt
%if (x<=0.5) res = 4/(lc*lc); else res = -4/(lc*lc); endif% s shape
endif
endfunction

Ntim = 200;
Nelt = 200;
Nnod = Nelt+1;
E = 1; %(beton)
rho = 1; %(beton)
A = 1; % barre de 10cm sur 10cm
c = sqrt(E/rho);
L = 1;
h = 1/Nelt; % 
cfl = 1.;
dt = cfl * h/c;% 
Yc = E/10;
ec = sqrt(2 * Yc / E);
sigc = E * ec;
lc = 0.2;
% two gauss point on the element
pg(1) = (1-sqrt(3)/3)/2;
pg(2) = (1+sqrt(3)/3)/2;

d = zeros(Ntim, Nelt);
s = zeros(Ntim, Nelt);
e = zeros(Ntim, Nelt);
energy = zeros(Ntim, Nelt);
Y = zeros(Ntim, Nelt);
strain_energy = zeros(Ntim, 1);
kinetic_energy = zeros(Ntim, 1);
dissip_energy = zeros(Ntim, 1);



bpos = zeros(Ntim, Nelt);
grad = zeros(Ntim, Nelt);

%s(1, 1:Nelt) = ones(1,Nelt);
%e(1, 1:Nelt) = ones(1,Nelt);

u = zeros(Ntim, Nnod);
ustat = zeros(Ntim, Nnod);
Ystat = zeros(Ntim, Nelt);
v = zeros(Ntim, Nnod);
a = zeros(Ntim, Nnod);
phi = zeros(Ntim, Nnod);
YmYc= zeros(Ntim, Nelt);
phidot = zeros(Ntim,1);
ddotbar = zeros(Ntim,1);
dissip = zeros(Ntim,1);
nbiter = zeros(Ntim,1);

    
m = rho * h * A * ones(Nnod);
m(1) = m(1)/2;
m(Nnod) = m(Nnod)/2;

%initialization
for j=1:Nnod;
  x(j) = (j-1) * h;
endfor;
for j=1:Nelt;
  xe(j) = 0.5 * (x(j) + x(j+1));
endfor;
for i=1:Ntim;
  t(i) = (i-1) * dt;
endfor;

% Initially the bar is loaded and all elements are at Yc
% A tls is placed on the first element, obviously it satifies
%   Ybar = Yc
% This extra damage will create new stress that are less than in the next element
% and ill thus give an unloading wave.

for j=1:Nnod;
  u(1,j) = x(j) * ec * L * 0.999;
  ustat(1,j) = x(j) * ec * L * 0.999;
  phi(1,j) = h-x(j);
endfor

phidot(1) = 0;

for j=1:Nelt;
  e(1,j) = (u(1,j+1)-u(1,j))/h;
  if (phi(1,j) > 0  && phi(1,j+1) > 0)
    for k=1:2
      philoc = pg(k)*phi(1,j)+ (1-pg(k))*phi(1,j+1);
      dloc = dval(philoc,lc);
      s(1,j) += 0.5 * (1-dloc) * E * e(1,j);
    endfor
  elseif  (phi(1,j) <= 0 && phi(1,j+1) <= 0)
    s(1,j) = E * e(1,j);
  elseif  (phi(1,j) > 0 && phi(1,j+1) <= 0)
    delta = abs(phi(1,j))/(abs(phi(1,j))+abs(phi(1,j+1)));
    sloc = 0;
    for k=1:2
      philoc = pg(k)*phi(1,j);
      dloc = dval(philoc,lc);
      sloc += 0.5 * (1-dloc) * E * e(1,j);
    endfor    
    s(1,j) = delta *  sloc +  (1-delta) * E * e(1,j);
  endif
endfor

%acceleration
a(1,1) = 0; 
a(1,Nnod) = 0;
for j=2:Nnod-1;
  a(1,j) = A*(s(1,j) - s(1,j-1)) /m(j); 
endfor

nbiter(1) = 0;

for i=2:Ntim;

% prediction
v(i,:)= v(i-1,:) + 0.5*dt*a(i-1,:);  
u(i,:)= u(i-1,:) + dt*v(i-1,:) + 0.5*dt*dt*a(i-1,:);  

%def computation and Y update.

for j=1:Nelt;
  e(i,j) = (u(i,j+1)-u(i,j))/(h);
  %b=0.5*E*e(i,j)*e(i,j)-Yc;
endfor

% moving the localization front
% we compute a = integral (Yn+1 - Yc) d' in the current non-local zone
% then we compute b = (Yn+1-Yc) d' on the front
% the shift in level set if the ratio of the two.

for j=1:Nnod;
  phi(i,j) = phi(i-1,j);
endfor

err_crit = 1e15;
nbiter(i) = 0;
while (err_crit > 1.e-6)
  nbiter(i) += 1;
  residu_Y = 0; tangent_Y = 0;
  loop_residu = 0;
  loop_tangent = 0;
  for j=1:Nelt;
    if (phi(i,j) > 0 && phi(i,j+1) > 0)
      for k=1:2
	philoc = pg(k)*phi(i,j) + (1-pg(k))*phi(i,j+1);
	residu_Y += h * 0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dp(philoc,lc);
	tangent_Y += h * 0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dpp(philoc,lc);
      endfor
      loop_residu += 1;
    elseif  (phi(i,j) > 0 && phi(i,j+1) <= 0)
      delta = h * abs(phi(i,j))/(abs(phi(i,j))+abs(phi(i,j+1)));
      for k=1:2
        philoc = pg(k)*phi(i,j);
	residu_Y += delta *  0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dp(philoc,lc);
	tangent_Y += delta *  0.5 * (0.5 * E * e(i,j) * e(i,j) - Yc) * dpp(philoc,lc);
      endfor
      loop_residu += 1;
      if (delta < h) tangent_Y += (0.5 * E * e(i,j) * e(i,j) - Yc)* dp(0,lc);
      else tangent_Y += (0.5 * E * e(i,j+1) * e(i,j+1) - Yc)  * dp(0,lc);
      endif
      loop_tangent +=1;
    endif
  endfor

% law Ybar = Yc
if 1
  %i
  %loop_residu
  %residu_Y
  %loop_tangent
  %tangent_Y
  %i
  %dval(phi(i-1,1),lc)
  %err_crit = abs(residu_Y)/(Yc*dval(phi(i-1,1),lc));
  %if (abs(tangent_Y) <= 1.e-10) tangent_Y = 0.02; endif 
  %tangent_Y
  %dphi = -residu_Y/tangent_Y;
  %dphi
  %i
  %dval(phi(i,1),lc)
  YbarmYc = residu_Y/(dval(phi(i,1),lc));
  residu = YbarmYc/Yc;
  err_crit = abs(residu);
  tangent = tangent_Y/(Yc*dval(phi(i,1),lc)) - (dp(phi(i,1),lc)/dval(phi(i,1),lc)^2) * (YbarmYc/Yc);
  %if (abs(tangent) <= 1.e-10) err_crit = 0.; dphi = 0; 
  %else
  dphi = - residu/tangent;
  %endif

endif

%delay effect law.
if 0
tauc=0.4;
coef_a=1;
l = max(phi(i,1),lc);
dpbar = dval(phi(i,1),lc)/l;
dphi_iter = phi(i,1) - phi(i-1,1);
YbarmYc = residu_Y/(dval(phi(i,1),lc));
if (YbarmYc < 0) YbarmYc = 0.; endif
residu_delay = - dpbar * dphi_iter * tauc/dt + (1-exp(-coef_a*YbarmYc/Yc));
err_crit = abs(residu_delay);

if (YbarmYc < 0) tangent = 0;
else
tangent = tangent_Y/(Yc*dval(phi(i,1),lc)) - (dp(phi(i,1),lc)/dval(phi(i,1),lc)^2) * (YbarmYc/Yc);
endif

tangent_delay = ((-dp(phi(i,1),lc)/l+dval(phi(i,1),lc)/l^2) * dphi_iter - dval(phi(i,1),lc)/l ) * tauc/dt + coef_a*exp(-coef_a*YbarmYc/Yc)*tangent;

dphi = -residu_delay/tangent_delay;
endif

for j=1:Nnod;
  phi(i,j) = phi(i,j) + dphi;
endfor

%err_crit = 0.

endwhile
  
phidot(i) = (phi(i,1) - phi(i-1,1))/dt;
    


%updating the stress
for j=1:Nelt;
  if (phi(i,j) > 0 && phi(i,j+1) > 0)
    for k=1:2
      philoc = pg(k)*phi(i,j)+ (1-pg(k))*phi(i,j+1);
      dlocg(k) = dval(philoc,lc);
      s(i,j) += 0.5 * (1-dlocg(k)) * E * e(i,j);
    endfor
  elseif  (phi(i,j) <= 0 && phi(i,j+1) <= 0)
    s(i,j) = E * e(i,j);
    dlocg(1) = dlocg(2) = 0;
  elseif  (phi(i,j) > 0 && phi(i,j+1) <= 0)
    delta = abs(phi(i,j))/(abs(phi(i,j))+abs(phi(i,j+1)));
    sloc = 0;
    for k=1:2
      philoc = pg(k)*phi(i,j);
      dlocg(k) = dval(philoc,lc);
      sloc += 0.5 * (1-dlocg(k)) * E * e(i,j);
    endfor    
    s(i,j) = delta * sloc +  (1-delta) * E * e(i,j);
  endif
  d(i,j) = 0.5*(dlocg(1)+dlocg(2));
endfor

%acceleration
%if (phi(i,1)/lc < 1) 
if (phi(i,1) <= lc) a(i,1) = 0; 
else a(i,1) =  A*s(i,1) /m(1);
endif 
a(i,Nnod) = 0;
for j=2:Nnod-1;
  a(i,j) = A*(s(i,j) - s(i,j-1)) /m(j); 
endfor


%correction
v(i,:)= v(i,:) + 0.5*dt*a(i,:);  


endfor

for i = 1:Ntim;
  for j=1:Nelt;
    Y(i,j) = 0.5*E*e(i,j)^2;
    YmYc(i,j) = Y(i,j)/Yc - 1;
    if (i > 1) 
      dissip(i) += h * Y(i,j) * (d(i,j)-d(i-1,j)); 
      kinetic_energy(i) += 0.5 * h * rho * 0.5 * (((u(i,j) - u(i-1,j))^2)+((u(i,j+1) - u(i-1,j+1))^2))/dt^2; 
      %ustat(i,j+1) = ustat(i,j) + h*s(1,Nelt)/(E*(1-d(i,j)));
    endif;
    energy(i,j)=h*Y(i,j)*(1-d(i,j));
  endfor
  %ustat(i,:) *= u(1,Nnod)/ustat(i,Nnod);
  for j=1:Nelt;
    Ystat(i,j) = 0.5*E*((ustat(i,j+1)-ustat(i,j))/h)^2;
  endfor;
  if (i > 1) dissip_energy(i) = dissip_energy(i-1) + dissip(i); endif
  strain_energy(i) = sum (energy(i,:));
  tot_energy(i) = strain_energy(i) + kinetic_energy(i) + dissip_energy(i);
endfor




