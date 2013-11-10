function[phinew,segment]=analyzeDamage(x,phi,h)

%produce:
%new phi based on distances - maxima

list_max = [];
value_max = [];
%direction = []; %1 or -1: slope
for i = 1:length(phi)-1
    if ((abs(phi(i)-phi(i+1)) > h*1e-8) && (i > 1) && (i < length(phi)-1)) %local maxima/minima
        if ((phi(i) > phi(i-1)) && (phi(i+1)>phi(i+2)))
            delta = (phi(i+1)-phi(i) + h)/2;
            list_max(end+1) = x(i) + delta;
            value_max(end+1) = phi(i) + delta;
        end
    end
    if ((i == 1) && phi(1)>phi(2))
                    list_max(end+1) = x(i);
                    value_max(end+1) = phi(i);
    end
    if ((i == length(phi)-1) && phi(i+1)>phi(i))
                    list_max(end+1) = x(i+1);
                    value_max(end+1) = phi(i+1);
    end

end

segment = {};
for j=1:length(list_max)
    segment{j} = [];
end

assert(length(x) >= 1);

for i = 1:length(phi)
   phinew(i) = -min(-value_max+abs(x(i) -list_max));
   dir = 0;
   
   for j=1:length(list_max)
       if (phinew(i) == -min(-value_max(j)+abs(x(i) -list_max(j))))
           %phinew(i) = phinew(i) + value_max(j);
           %correct sign
           %phinew(i) = (x(i)-list_zero(j)) * direction(j);
           %add to segment
           segment{j}(end+1) = i;
           break;
       end
   end
   phinew(i) = max(phinew(i),phi(i));
end

%join segments together if same peak
for j=1:length(list_max)-1
    if (abs(list_max(j)-list_max(j+1)) < 1e-8)
        segment{j} = [segment{j} segment{j+1}];
        segment{j+1} = [];
    end
end