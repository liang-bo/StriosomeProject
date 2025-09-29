function [x,y] = bin2patch(bin)
% convert the binary data to patch vectors
% debug data
% bin = mousedata.behavior(1).obj1{1, 1};
bin = double(bin);
L = length(bin);%tt = 1:L;
binedge = gradient(bin);

binedgeup = find(binedge>0)+1;binedgeup(binedgeup>L) = L;binedgeup(binedgeup<1) = 1;
binedgedown = find(binedge<0);binedgedown(binedgedown>L) = L;binedgedown(binedgedown<1) = 1;

binedgeup = binedgeup(1:2:end);
binedgedown = binedgedown(1:2:end);

egnumup = length(binedgeup);
egnumdown = length(binedgedown);

if ~(isempty(binedgeup)||isempty(binedgedown))
    binedgeup1 = binedgeup(1);%
    % binedgeupend = binedgeup(end);
    binedgedown1 = binedgedown(1);%
    % binedgedownend = binedgedown(end);
    
    switch egnumup-egnumdown
        
        case 1
            binedgedown = [binedgedown;L];
            egnum = egnumup;
        case -1
            binedgeup = [1;binedgeup];
            egnum = egnumdown;
        case 0
            
            if binedgeup1 > binedgedown1
                binedgeup = [1;binedgeup];
                binedgedown = [binedgedown;L];
                egnum = egnumup+1;
            else
                egnum = egnumup;
            end
            
        otherwise
            printf('something is wrong')
    end
    
    x = zeros(4,egnum);
    y = zeros(4,egnum);
    for ii = 1:egnum
        x(:,ii) = [binedgeup(ii);binedgeup(ii);binedgedown(ii);binedgedown(ii);];
        y(:,ii) = [0;1;1;0;];
    end
else
    x = [];
    y = [];
end
% For debugging purpose
% figure
% plot(tt,bin,tt(binedgeup),ones(egnumnew,1),'^',tt(binedgedown),ones(egnumnew,1),'v')
% hold on;
% patch('XData',x,'YData',y)

end