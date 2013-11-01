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
    %if (x<=0.5) res = 4/(lc*lc); else res = -4/(lc*lc); end% s shape
end
end