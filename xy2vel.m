function [vel, dist] = xy2vel(x,y,pixsiz)
% input:
% x: x position value in pixel
% y: y position value in pixel
% pixsiz: pixel size cm/pixel

% output:
% vel: velocity in cm/s
% dist: distance in cm

Fs = 10; % assume frame rate is 10 Hz

dist = pixsiz*sqrt((x(2:end)-x(1:(end-1))).^2 + (y(2:end)-y(1:(end-1))).^2);
vel = dist*Fs;

dist = [0;dist(:)];
vel  = [0;vel(:)];

end