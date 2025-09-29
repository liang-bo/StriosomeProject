function dataset = add_smooth_vel_dist(dataset)

mnum = numel(dataset.mouse);
pixsiz = 38/380; % cm/pix

for m = 1:mnum
    
    for s = 1:numel(dataset.mouse(m).file)
        
        x = dataset.mouse(m).file(s).behavior.x;
        y = dataset.mouse(m).file(s).behavior.y;
        
        [vel, ~] = xy2vel(x,y,pixsiz);
        velfilt = smooth(1:length(vel),vel,10,'rloess');
        velfilt(velfilt<0) = 0;
        
        dataset.mouse(m).file(s).behavior.velocity = velfilt;
        
        fprintf('mouse %d, session %d is done. \n', m, s);
    end
    
end





end