function [result] = parcialNU_parcial_x(NU, theta, h, E)

if(x<= E)
    
    result = 0;
end

if(x > E)
    
    result = (1 - NU)*(tand(theta) / h);
    
end

end