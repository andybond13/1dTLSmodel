function [nfrags,fragmentList] = findFragments(x,phi,d)

fragmentList={};

sbegin = 0;
send = 0;
for i=1:length(x)-1
    if (sbegin == 0)
        if (i == 1)
            sbegin = i;
        end
    end
    
    if (send == 0)
        if (i == length(x)-1)
            send = i+1;
        end
        if ((i > 1) && (i < length(x)))
            if (phi(i)>phi(i-1) && phi(i)>phi(i+1) && d(i) == 1)
                send = i;
            end
        end
        
    end
    
    if (sbegin*send ~= 0)
        list = [sbegin send];
        fragmentList{end+1} = list;
        sbegin = send + 1;
        send = 0;
    end
end

%calculate total number of fragments (removing symmetry simplification)
nfrags = length(fragmentList)*2;
if (d(1) < 1)
    nfrags = nfrags - 1;
end

1+1;