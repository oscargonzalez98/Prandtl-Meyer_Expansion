function [value] = y_s(x, E, theta)

if(x<=E)
    value = 0;
end

if(x<E)
    value = -1*(x-E)*tand(theta);
end


end