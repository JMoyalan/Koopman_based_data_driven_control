function [value,isterminal,direction] = eventfun(t,x,xlim)
value = 1;
for i1 = 1 : length(x)
    if abs(x(i1)) >= xlim(i1)
        value = 0;
    end
end
isterminal = 1;
direction = 0;
end