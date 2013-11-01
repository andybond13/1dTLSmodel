% this function is the derivative of d with respect to phi
function res = dp(phi,lc)
res = 0;
x = phi/lc;
if   (x >= 0 && x <= 1)
    %res = 1/lc;           %lin
    res = 2/lc * (1-x) ;  %quad
    %res = 1/lc * (3-6*x+3*x*x); %cubic
    %res = 1/lc * 2*x; %sqrt
    %if (x<=0.5) res = 4*x/lc; else res = (-4 * x + 4)/lc; end% s shape
end
end