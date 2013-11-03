function[phi,segment]=analyzeDamage(x,phi,h)

%produce:
%new phi based on distances

list_zero = [];
direction = []; %1 or -1: slope
for i = 1:length(phi)-1
    if (phi(i)*phi(i+1) <= 0)
        phi1 = phi(i);
        phi2 = phi(i+1);
        delta = phi1/(phi1+phi2)*h;
        list_zero(end+1) = x(i) + delta;
        dir = 1;
        if (phi2 < phi1)
            dir = -1;
        end
        direction(end+1) = dir;
    end
end

segment = {};
for j=1:length(list_zero)
    segment{j} = [];
end

for i = 1:length(phi)
   phinew(i) = min(abs(x(i) -list_zero));
   dir = 0;
   
   for j=1:length(list_zero)
       if (phinew(i) == abs(x(i)-list_zero(j)))
           %correct sign
           phinew(i) = (x(i)-list_zero(j)) * direction(j);
           %add to segment
           segment{j}(end+1) = i;
           break;
       end
   end
end

%join segments together if same peak
for j=1:length(list_zero)-1
    if ((direction(j) == 1)  && (direction(j+1) == -1))
        segment{j} = [segment{j} segment{j+1}];
        segment{j+1} = [];
    end
    
end
