function phi = checkFailureCriteria(t,x,phi,criterion,elemOrNodal,qty,absOrAsIs,phiPos,failvalue)
%x              -mesh
%phi            -level-set calculated at this time-step
%criterion      -criterion to compare against for failure
%elemOrNodal    -either 'elem' or 'nodal' - is criterion elemental or nodal
%qty            -the quantity to be compared to the criterion -e.g. s,Y
%absOrAsIs      -whether to compare (0) qty or (1) abs(qty) to criterion
%phiPos         -whether(1) or not(0) failure cannot occur depending on if phi>0
%failvalue      -what to call phi at localization zone if created - e.g. h
%failure if qty > criterion

assert(strcmp(elemOrNodal,'elem') || strcmp(elemOrNodal,'nodal'));
assert(absOrAsIs == 0 || absOrAsIs == 1);
if (strcmp(elemOrNodal,'nodal'))
    assert(length(x) == length(qty));
else
    assert(length(x) == length(qty)+1);
end

assert((length(qty) == length(criterion)) || (length(criterion) == 1))
if (length(criterion) == 1)
    criterion = criterion*ones(size(qty));
end
    
xlist = [];

for i=1:length(qty)
    if (phiPos == 0) %can't fail if phi>0
        if (strcmp(elemOrNodal,'nodal'))
            if (phi(i) > 0)
                continue;
            end
        else
            if (phi(i)>0 || phi(i+1)>0)
                continue;
            end
        end
    end
    
    
    qtyc = qty(i);
    if (absOrAsIs == 1)
        qtyc = abs(qtyc);
    end
    if (qtyc > criterion(i))
        if (strcmp(elemOrNodal,'nodal'))
            xlist(end+1) = x(i);
        else
            %assume middle of element
            xlist(end+1) = 0.5*(x(i)+x(i+1));
        end
    end
end
%nucleate list
for i=1:length(xlist)
    phi = nucleate(t,x,phi,xlist(i),failvalue);
end