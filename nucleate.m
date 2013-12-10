function[phi] = nucleate(t,x,phi,xnuc,phinuc)
%x      -mesh
%phi    -level-set calculated at this time-step
%xnuc   -location(s) of localizations to be nucleated
%phinuc -amount of the level-set to be set at nucleated localizations

assert(length(xnuc) == length(phinuc));
assert(length(xnuc) > 0);

for j = 1:length(xnuc)
    
    h = 0;
    loc = 0;
    for i=1:length(x)-1
        if ((xnuc(j) >= x(i)) && (xnuc(j) < x(i+1)))
            loc = i;
            h = x(i+1)-x(i);
            delta = (xnuc(j) - x(i))/h;
            break;
        end
    end
    
    assert(loc ~= 0);
    
    phi(loc) = phinuc(j) - delta*h;
    phi(loc+1) = phinuc(j) - (1-delta)*h;
    sprintf('crack nucleated, t = %f, x = %f',t,xnuc(j))
end

