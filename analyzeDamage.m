function[phinew,newsegment]=analyzeDamage(x,phi,h)

%produce:
%new phi based on distances - maxima

if (sum(phi) == -1*length(phi))
    phinew = phi;
    newsegment={};
    return;
end

list_max = [];
value_max = [];
%direction = []; %1 or -1: slope
for i = 1:length(phi)-1
    if ((abs(phi(i)-phi(i+1)) > h*1e-8) && (i > 1) && (i < length(phi)-1)) %local maxima/minima
        if ((phi(i) > phi(i-1)) && (phi(i+1)>phi(i+2)))
            delta = (phi(i+1)-phi(i) + h)/2;
            if (phi(i)+delta < 0)
                continue;
            end
            list_max(end+1) = x(i) + delta;
            value_max(end+1) = phi(i) + delta;
        end
    end
    if ((i == 1) && (phi(1)-phi(2)> -eps))

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

if (length(value_max) == length(list_max) && length(value_max) == 0)
    phinew = phi;
    return;
end

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
    if (abs(list_max(j)-list_max(j+1)) < h)
        segment{j} = [segment{j} segment{j+1}];
        segment{j+1} = [];
        segment{j} = sort(segment{j});
    end
end


%new
%split hat segments into two
newsegment = {};
list_maxnew = [];
value_maxnew = [];
for i=1:length(segment)
    indices = segment{i};
    if (length(indices) == 0)
        continue; %don't copy empties
    end
    assert(length(indices) > 1);
    if ((phinew(indices(1)) < phinew(indices(2))) && (phinew(indices(end)) < phinew(indices(end-1))))
        %hat - duplicate
        list_maxnew(end+1) = list_max(i);
        list_maxnew(end+1) = list_max(i);
        value_maxnew(end+1) = value_max(i);
        value_maxnew(end+1) = value_max(i);
        iphimax = intersect([indices],find(phi==max(phi(indices))));
        newsegment{end+1} = indices(find(indices<=iphimax));
        newsegment{end+1} = indices(find(indices>iphimax));
    elseif ((phinew(indices(1)) > phinew(indices(2))))% && (phinew(indices(end)) < phinew(indices(end-1))))
        %not hat - copy as is
        list_maxnew(end+1) = list_max(i);
        value_maxnew(end+1) = value_max(i);
        newsegment{end+1} = segment{i};
    elseif ((phinew(indices(1)) < phinew(indices(2))) && (phinew(indices(end)) > phinew(indices(end-1))))
        %hot hat - copy as is
        list_maxnew(end+1) = list_max(i);
        value_maxnew(end+1) = value_max(i);
        newsegment{end+1} = segment{i};        
    else
        [list_max;value_max]
        [indices;x(indices);phinew(indices)]
        plot(x(indices),phinew(indices))
        assert(1==0);
    end
end

for i=1:length(newsegment)
    indices = newsegment{i};
    for j = indices
        phinew(j) = value_maxnew(i)-abs(x(j) -list_maxnew(i));
    end
end


1+1;