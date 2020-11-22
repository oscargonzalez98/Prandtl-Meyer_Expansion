function [value] = h(x, E, theta, H)

if(x<=E)
    value = H;
end

if(x>E)
    value = H + (x - E)*tand(theta);
end

end