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
    %if (x<=0.5) res = 2*x^2; else res = -2*x^2 + 4 * x - 1; end% s shape
end
end