function [x_data_interp, y_data_interp, x_fit, y_fit, pref_dir, fwhm, osi_fit, dsi_fit, resnorm] = fit_sum_von_mises(x_data, y_data, interp_points, plot_it, visible)


%x_data = deg2rad(theta); %theta = [0:30:330];
%y_data = double(ori_peak_r_x); 
%interp_points = 359;% keep as deg

% % For testing
% x_data = tmp_rads;
% y_data = tmp_all_dirs_pref_sf_tf_mat;
% interp_points = 360;

%Flag to turn on or off plotting
if ~exist('plot_it','var')
    plot_it = 0;
end

%Flag to make plots visible or not; default is to display the plots
if ~exist('visible','var')
    visible = 'on';
end

% Force all dfof values to positive to prevent large DSIs/OSIs
if min(y_data) < 0
    y_data = y_data-min(y_data);
end

% Normalize the data being fit
y_data = y_data./max(y_data);

sum_of_vonMises = @(v, x_data)v(1) + (v(2)*exp(v(3)*((cos(x_data - v(4)^2) - 1)))) + (v(5)*exp( v(6)*((cos(x_data - v(7)^2) - 1))));

x0(1) = min(y_data); %v(1) = baseline
[x_amp, max_theta_ind] = max(y_data); %v(2) = Amp1
x0(2) = x_amp*0.8;
x0(3) = 1.0; %v(3) = k1, width 1, if k=2, FWHM= 50 deg
x0(4) = x_data(max_theta_ind); %v(4) = peak Ori 1
x0(5) = 0.8*(x0(2)); %v(5) = Amp2
x0(6) = 1.5; %v(6) = k2, width 2
x0(7) = (x_data(max_theta_ind) + pi); %v(7) = Peak Ori 2

% if x0(7)>2*pi % peak was towards the end of theta range
%     x0(5) = x0(2); 
%     x0(7) = x0(4);
%     x0(2) = 0.8* (x0(2)); 
%     x0(4) = (x_data(max_theta_ind) - pi); 
% end

x0= double(x0); 

options = optimoptions('lsqcurvefit', 'MaxIter', 500, 'Display', 'off');

lb = [ -1 0 0.1 0 0 0.1 0]; %KS added to try to remove neg going peak amps
ub = [ 0.5 15 10 359 15 10 359]; %KS added to try to remove neg going peak amps

%interpolate data: 
x_data_interp = linspace(0, 2*pi, interp_points);
y_data(end+1) = y_data(1); %wrap it
x_data(end+1) = 2*pi; 
y_data_interp = interp1(x_data, y_data, x_data_interp, 'linear', 'extrap');

%fit: %KS commented out interpolating before fitting, should interp after
%fitting. Need to make sure the y_fit will interp to the function itself in
%line 47
[x_fit, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(sum_of_vonMises, x0, x_data, y_data, lb, ub, options);
y_fit = sum_of_vonMises(x_fit, x_data_interp);

% Now try the fitting again if not well fit
if resnorm < 10
    x0(1) = min(y_data); %v(1) = baseline
    [x0(2), max_theta_ind] = max(y_data); %v(2) = Amp1
    x0(3) = 1.5; %v(3) = k1, width 1, if k=2, theta(0.5)= 50 deg
    x0(4) = x_data(max_theta_ind); %v(4) = peak Ori 1
    x0(5) = 0.5*(x0(2)); %v(5) = Amp2
    x0(6) = 1.5; %v(6) = k2, width 2
    x0(7) = (x_data(max_theta_ind) + pi); %v(7) = Peak Ori 2

    if x0(7)>2*pi % peak was towards the end of theta range
        x0(5) = x0(2); 
        x0(7) = x0(4);
        x0(2) = 0.5* (x0(2)); 
        x0(4) = (x_data(max_theta_ind) - pi); 
    end

    x0 = double(x0);
    %fit:
    [x_fit, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(sum_of_vonMises, x0, x_data_interp, y_data_interp, [], [], options);
    y_fit = sum_of_vonMises(x_fit, x_data_interp);
end

%Get fit params: 
k1 = x_fit(3);
k2 = x_fit(6);

% Find preferred and orthogonal directions from the fitted data
[pref_amp, pref_idx] = max(y_fit);
pref_dir = rad2deg(x_data_interp(pref_idx));
orth1_idx = mod(floor(pref_dir) + 90, interp_points)+1;
orth1_dir = rad2deg(x_data_interp(orth1_idx));
orth1_amp = y_fit(orth1_idx);
orth2_idx = mod(floor(pref_dir) - 90, interp_points)+1;
orth2_dir = rad2deg(x_data_interp(orth2_idx));
orth2_amp = y_fit(orth2_idx);

% Look for second peak in the fitted data
mins = find(islocalmin(y_fit)); % to look for other peak
if isempty(mins)
    no_min = 1; 
    mins(1) = nan;
    mins(2) = nan;
    peak = 1;
elseif size(mins) == 1 % Likely only one peak, because only one minimum
    peak = 1;
    no_min = 0; 
else
    % Preferred direction is between the minimums, so look outside them for
    % the second peak
    if pref_idx >= mins(1) && pref_idx <= mins(2)
        peak = 2;
        shift_val = interp_points - mins(2); % value to shift the fitted data to be able to continuously search for second peak
        no_min = 0;
    % Preferred direction is outside the minimums
    else
        peak = 2;
        no_min = 0;
    end
end

% Find FWHM of preferred direction peak
half_max = pref_amp/2;
range = (ceil(interp_points/2))/2;
% If the peak is too close to either side of the unit circle, shift 180
% degrees to make values continuous
% Shift the fitted data so that the pref response is centered at 180
y_fit_shift = circshift(y_fit,(ceil(interp_points/2)-pref_idx));
% if pref_idx <= range || pref_idx >= interp_points - range
%     y_fit_shift = circshift(y_fit,ceil(interp_points/2));
% else
%     y_fit_shift = y_fit;
% end
[max_shift, shift_idx] = max(y_fit_shift);
min_range = mod(shift_idx-range,interp_points);
max_range = mod(shift_idx+range,interp_points);
% Check below and above the max value separately in case fit isn't
% perfectly symmetric. Might be unnecessary.
y_fit_shift_below = y_fit_shift(min_range:shift_idx);
[~,closest_idx_below] = min(abs(y_fit_shift_below-half_max));
y_fit_shift_above = y_fit_shift(shift_idx:max_range);
[~,closest_idx_above] = min(abs(y_fit_shift_above-half_max));
y_fit_shift_vals = y_fit_shift(min_range:max_range);
[~,closest_idx_all] = min(abs(y_fit_shift_vals-half_max));
fwhm = closest_idx_above + (length(y_fit_shift_below)-closest_idx_below);

% I think this is calculating the half width at half max
% FWHM = rad2deg(acos((log((1+exp(-2*k))/2)/k) + 1));

% Now determine the null direction. This changes depending in whether there
% are one or two peaks.
if peak == 1
    % If FWHM is really broad, using the fit could give spurious results,
    % so just use the value 180 degress from the max value
    if fwhm < 180
        if (pref_idx + 180)+round(fwhm/2) > interp_points
            max_shift = (pref_idx + 180)+round(fwhm/2);
            shift_diff = max_shift - interp_points;
            y_fit_null_shift = circshift(y_fit,-(shift_diff+2));
            pref_idx_shift = pref_idx - (shift_diff+2);
            min_val = mod((pref_idx_shift + 180)-round(fwhm/2),interp_points);
            max_val = mod((pref_idx_shift + 180)+round(fwhm/2),interp_points);
            [null_amp, null_idx] = max(y_fit_null_shift(min_val:max_val));
            null_idx = mod(min_val + shift_diff,interp_points)+1;
        else
            min_val = mod((pref_idx + 180)-round(fwhm/2),interp_points);
            max_val = mod((pref_idx + 180)+round(fwhm/2),interp_points);
            [null_amp, null_idx] = max(y_fit(min_val:max_val));
            null_idx = min_val + null_idx;
        end
    else
        null_amp = y_fit(mod(pref_idx + 180,interp_points)+1);
        null_idx = mod(pref_idx + 180,interp_points)+1;
    end
elseif peak == 2
    if no_min == 0
        % If there are two minimum values, and the preferred direction did not fall between them check between them
        if length(mins) > 1 && ~exist('shift_val','var')
            [null_amp,null_idx] = max(y_fit(mins(1):mins(2)));
            null_idx = mod(mins(1) + null_idx,interp_points)+1;
        % Two minimum values, but preferred direction is between them
        elseif length(mins) > 1 && exist('shift_val','var')
            y_fit_peak2_shift = circshift(y_fit,shift_val+2); % add a buffer for visualization and so you don't get zero values later
            min_val = mod(mins(2)+shift_val+2,interp_points);
            max_val = mod(mins(1)+shift_val+2,interp_points);
            [null_amp,null_idx] = max(y_fit_peak2_shift(min_val:max_val));
            null_idx = mod(null_idx - shift_val,interp_points)+1;
        % If there's only one min, search around it using the fwhm/2
        else
            min_val = mod((mins(2) + 180)-round(fwhm/2),interp_points);
            max_val = mod((mins(2) + 180)+round(fwhm/2),interp_points);
            [null_amp, null_idx] = max(y_fit(min_val:max_val));
            null_idx = mod(min_val + null_idx,interp_points)+1;
        end
    elseif no_min == 1 % If there are no minimums, but a second peak is detected, just look 180 degrees +/- fwhm/2 from pref_dir
        min_val = mod((pref_idx + 180)-round(fwhm/2),interp_points);
        max_val = mod((pref_idx + 180)+round(fwhm/2),interp_points);
        [null_amp, null_idx] = max(y_fit(min_val:max_val));
        null_idx = mod(min_val + null_idx,interp_points)+1;
    end
end

null_dir = rad2deg(x_data_interp(null_idx));

% If the null direction is 'too close' to the preferred direction, it's
% either a fitting error (or a bad roi) so force null direction to be 180
% degrees from preferred
if abs(null_dir-pref_dir) < 90
    null_amp = y_fit(mod(pref_idx + 180,interp_points));
    null_idx = mod(pref_idx + 180,interp_points)+1;
end

% Calculate OSI and DSI from fitted data
osi_fit = abs((pref_amp - mean([orth1_amp orth2_amp]))/(pref_amp + mean([orth1_amp orth2_amp])));
dsi_fit = (pref_amp - null_amp)/(pref_amp + null_amp);

if plot_it == 1
    % Plot the data wtih the fit overlaid:
%     figure(1); clf;
    f = figure('visible',visible); 
    raw_data = plot(x_data_interp, y_data_interp, 'k-');
    hold on
    fitted_data = plot(x_data_interp, y_fit, 'r-', 'LineWidth', 2.5);
    if pref_amp < 0 % should never happen
        ylim([min(y_data_interp) max(y_data_interp) + abs(max(y_data_interp)*0.4)]);
    else
        ylim([0 1.4]);
    end
    % xl(0,360);
    pref = plot(x_data_interp(pref_idx),pref_amp,'cx','Markersize',10, 'Linewidth', 5);
    null = plot(x_data_interp(null_idx),null_amp,'cx','Markersize',10, 'Linewidth', 5);
    orth1 = plot(x_data_interp(orth1_idx),orth1_amp,'bx','Markersize',10, 'Linewidth', 5);
    orth2 = plot(x_data_interp(orth2_idx),orth2_amp,'bx','Markersize',10, 'Linewidth', 5);
    legend('raw', 'vonMises fit')
    legend box off
    xl(deg2rad(0),deg2rad(360));
    xticks(deg2rad(0:30:330))
    xticklabels({0:30:300})
    xlabel('Degrees')
    ylabel('dF/F')
    str = ['Dir pref = ', num2str(pref_dir,3),  char(176), ', FWHM = ', num2str(fwhm,3), char(176), newline, 'Dir 2 peak = ', num2str(null_dir,3), char(176), newline, 'OSI = ', num2str(osi_fit, 3), newline,  'DSI = ', num2str(dsi_fit, 3)];
    dim = [.125 .825 .5 .1];  
    annotation('textbox',dim,'String',str,'FitBoxToText', 'on', 'EdgeColor','none')
    set(gcf, 'Position', [680 639 454 339]);
end

clear shift_val

end