function[phi] = nucleate(x,phi,xnuc,phinuc)

h = 0;
loc = 0;
for i=1:length(x)-1
    if ((xnuc >= x(i)) && (xnuc < x(i+1)))
        loc = i;
        h = x(i+1)-x(i);
        delta = (xnuc - x(i))/h;
        break;
    end
end

assert(loc ~= 0);

phi(loc) = phinuc - delta*h;
phi(loc+1) = phinuc - (1-delta)*h;