function [x,y] = positionsetzero(x,y)
% this function is to set the starting point of the xy coordinate to
% zeros 
%

x0 = x(x>0);
y0 = y(y>0);

x = x - min(x0);
y = y - min(y0);

end