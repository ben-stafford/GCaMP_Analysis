function [CV, DirCV, Ori_pref, L_ori, L_dir] = get_CV_DirCV(Rxi, theta)
%theta - radian
%Rxi = responses at each theta
%all values shoule be above 0, so may need to normalize 

if max(theta)>(2*pi)
    theta_r = deg2rad(theta); 
else
    theta_r = theta;
end

if isstruct(Rxi)
    Rxi = cell2mat(struct2cell(Rxi));
end

%combined Fourier method: 
Z = sum ( Rxi .* exp((-2i*theta))); % -2i is different from r in circ_r calculation
Ori_pref = abs(rad2deg(angle(Z))); 

Znorm = Z/sum(Rxi); % normalize
L_ori = abs((Znorm)); 
%CV = L_ori ; %1-CV = L_ori, no need to change
%CV = 1 - L_ori; 

%broken down by components: 
CV = ( sqrt( sum(Rxi.*sin(2*theta_r))^2 + (sum(Rxi.*cos(2*theta_r)))^2)) / sum(Rxi);




Zd = sum ( Rxi .* exp((-1i*theta))); % (same as r in circ_r calculation)
Zdnorm = Zd/sum(Rxi); % normalize
L_dir = abs((Zdnorm)); % (circ_r function does not use 'real' function, just takes abs)
%DirCV = L_dir; %1-dirCV = L_dir, no need to change
%DirCV = 1 - L_dir; 

DirCV = ( sqrt( sum(Rxi.*sin(1*theta_r))^2 + (sum(Rxi.*cos(1*theta_r)))^2)) / sum(Rxi);


end

