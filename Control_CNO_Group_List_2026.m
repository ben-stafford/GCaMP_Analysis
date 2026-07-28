% Script for generating mat file to feed into Tuning_Analysis code so that
% several control or CNO datasets can be fed into
% Control_CNO_Speed_Group_Compare.m to make comparisons between control and
% CNO conditions. Currently works only for analysis of speed tuning.
clear

%% Nr5a1 Datasets
%% Speed Stimulus
data = struct();

data.animal_ids = {
    '69';
    '71';
    '66';
    '7823'
    '129'
    };
data.file_names = {
    {'day1_000_001_speed_control.mat','day1_004_005_speed_cno.mat'};
    {'day1_000_001_speed_control.mat','day1_004_005_speed_cno.mat'};
    {'day1_000_001_speed_control.mat','day1_004_005_speed_cno.mat'};
    {'day1_000_001_speed_control.mat','day1_004_005_speed_cno.mat'};
    {'day1_000_001_speed_control.mat','day1_004_005_speed_cno.mat'};
    };
data.stim_type = {
    'speed';
    'speed';
    'speed';
    'speed';
    'speed';
    };

if ismac
    cd ~/Dropbox/Salk/GCaMP/Data/
elseif ispc
    cd C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data
end

save('Control_CNO_Speed_Nr5a1_List', 'data');

%% Direction Stimulus
data = struct();

data.animal_ids = {
    '69';
    '71';
    '66';
    '7823'
    '129'
    };
data.file_names = {
    {'day1_002_direction_control.mat','day1_003_direction_cno.mat'};
    {'day1_002_direction_control.mat','day1_003_direction_cno.mat'};
    {'day1_002_direction_control.mat','day1_003_direction_cno.mat'};
    {'day1_002_direction_control.mat','day1_006_direction_cno.mat'};
    {'day1_002_direction_control.mat','day1_006_direction_cno.mat'};
    };
data.stim_type = {
    'direction';
    'direction';
    'direction';
    'direction';
    'direction';
    };

if ismac
    cd ~/Dropbox/Salk/GCaMP/Data/
elseif ispc
    cd C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data
end

save('Control_CNO_Direction_Nr5a1_List', 'data');

%% SF and TF Direction Stimulus
data = struct();

data.animal_ids = {
    '69';
    '71';
    '66';
    '129';
    '7823';
    };
data.file_names = {
    {'day2_000_001_002_sf_tf_direction_control.mat','day2_004_005_006_sf_tf_direction_cno.mat'};
    {'day2_000_001_002_sf_tf_direction_control.mat','day2_004_005_006_sf_tf_direction_cno.mat'};
    {'day2_000_001_002_sf_tf_direction_control.mat','day2_004_005_006_sf_tf_direction_cno.mat'};
    {'day2_000_001_002_sf_tf_direction_control.mat','day2_003_004_005_sf_tf_direction_cno.mat'};
    {'day3_000_001_002_sf_tf_direction_control.mat','day3_003_004_005_sf_tf_direction_cno.mat'};
    };
data.stim_type = {
    'sf_tf_direction';
    'sf_tf_direction';
    'sf_tf_direction';
    'sf_tf_direction';
    'sf_tf_direction';
    };

if ismac
    cd ~/Dropbox/Salk/GCaMP/Data/
elseif ispc
    cd C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data
end

save('Control_CNO_SF_TF_Direction_Nr5a1_List', 'data');

%% Nr5a1 No DREADDs Datasets
%% Speed Stimulus
data = struct();

data.animal_ids = {
    % '281'
    '277'
    '557'
    '563'
    '559'
    };
data.file_names = {
    % {'day1_000_001_speed_control.mat','day1_003_004_speed_cno.mat'};
    {'day2_000_001_speed_control.mat','day2_004_005_speed_cno.mat'};
    {'day1_000_001_speed_control.mat','day1_003_004_speed_cno.mat'};
    {'day1_000_001_speed_control.mat','day1_003_004_speed_cno.mat'};
    {'day2_000_001_speed_control.mat','day2_003_004_speed_cno.mat'};
    };
data.stim_type = {
    % 'speed';
    'speed';
    'speed';
    'speed';
    'speed';
    };

if ismac
    cd ~/Dropbox/Salk/GCaMP/Data/
elseif ispc
    cd C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data
end

save('Control_No_DREADDs_Speed_Nr5a1_List', 'data');

%% Direction Stimulus
data = struct();

data.animal_ids = {
    '277'
    '281'
    };
data.file_names = {
    {'day1_002_direction_control.mat','day1_005_direction_cno.mat'};
    {'day1_002_direction_control.mat','day1_005_direction_cno.mat'};
    };
data.stim_type = {
    'direction';
    'direction';
    };

if ismac
    cd ~/Dropbox/Salk/GCaMP/Data/
elseif ispc
    cd C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data
end

save('Control_No_DREADDs_Direction_Nr5a1_List', 'data');

%% SF and TF Direction Stimulus
data = struct();

data.animal_ids = {
    '277';
    '281';
    '559';
    };
data.file_names = {
    {'day3_000_001_002_sf_tf_direction_control.mat','day3_004_005_006_sf_tf_direction_cno.mat'};
    {'day2_000_001_002_sf_tf_direction_control.mat','day2_003_004_005_sf_tf_direction_cno.mat'};
    {'day1_000_001_002_sf_tf_direction_control.mat','day1_003_004_005_sf_tf_direction_cno.mat'};
    };
data.stim_type = {
    'sf_tf_direction';
    'sf_tf_direction';
    'sf_tf_direction';
    };

if ismac
    cd ~/Dropbox/Salk/GCaMP/Data/
elseif ispc
    cd C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data
end

save('Control_No_DREADDs_SF_TF_Direction_Nr5a1_List', 'data');


%% Nr5a1 No CNO Datasets
%% Speed Stimulus
data = struct();

data.animal_ids = {
    '129'
    };
data.file_names = {
    {'day1_000_001_speed_control.mat','day1_004_005_speed_cno.mat'};
    };
data.stim_type = {
    'speed';
    };

if ismac
    cd ~/Dropbox/Salk/GCaMP/Data/
elseif ispc
    cd C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data
end

save('Control_No_CNO_Speed_Nr5a1', 'data');

%% Tlx3 Datasets
%% Speed Stimulus
data = struct();

data.animal_ids = {
    '7411';
    '7412';
    '33'; 
    '110';
    '109';
    };
data.file_names = {
    {'day2_000_001_speed_control.mat','day2_003_004_speed_cno.mat'};
    {'day3_000_001_speed_control.mat','day3_003_004_speed_cno.mat'};
    {'day1_000_001_speed_control.mat','day1_004_005_speed_cno.mat'};
    {'day1_000_001_speed_control.mat','day1_004_005_speed_cno.mat'};
    {'day1_000_001_speed_control.mat','day1_004_005_speed_cno.mat'};
    };
data.stim_type = {
    'speed';
    'speed';
    'speed';
    'speed';
    'speed';
    };

if ismac
    cd ~/Dropbox/Salk/GCaMP/Data/
elseif ispc
    cd C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data
end

save('Control_CNO_Speed_Tlx3_List', 'data');

%% SF x TF Direction Stimulus
data = struct();

data.animal_ids = {
    '7411';
    '7412';
    '109';
    '110';
    };
data.file_names = {
    {'day3_000_001_002_sf_tf_direction_control.mat','day3_003_004_005_sf_tf_direction_cno.mat'};
    {'day2_000_001_002_sf_tf_direction_control.mat','day2_003_004_005_sf_tf_direction_cno.mat'};
    {'day2_000_001_002_sf_tf_direction_control.mat','day2_004_005_006_sf_tf_direction_cno.mat'};
    {'day2_000_001_002_sf_tf_direction_control.mat','day2_004_005_006_sf_tf_direction_cno.mat'};    
    };
data.stim_type = {
    'sf_tf_direction';
    'sf_tf_direction';
    'sf_tf_direction';
    'sf_tf_direction';
    };

if ismac
    cd ~/Dropbox/Salk/GCaMP/Data/
elseif ispc
    cd C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data
end

save('Control_CNO_SF_TF_Direction_Tlx3_List', 'data');

%% Tlx3 no DREADDs
% Speed Stimulus
data = struct();

data.animal_ids = {
    '255';
    '267';
    '272';
    };
data.file_names = {
    {'day1_000_001_speed_control.mat','day1_003_004_speed_cno.mat'};
    {'day3_000_001_speed_control.mat','day3_003_004_speed_cno.mat'};
    {'day1_000_001_speed_control.mat','day1_003_004_speed_cno.mat'};
    
    };
data.stim_type = {
    'speed';
    'speed';
    'speed';
    };

if ismac
    cd ~/Dropbox/Salk/GCaMP/Data/
elseif ispc
    cd C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data
end

save('Control_No_DREADDs_Speed_Tlx3_List', 'data');

%% Crh Datasets
%% Speed Stimulus
data = struct();

data.animal_ids = {
    '541';
    '447';
    };
data.file_names = {
    {'day1_000_001_speed_control.mat','day1_003_004_speed_cno.mat'};
    {'day1_000_001_speed_control.mat','day1_003_004_speed_cno.mat'};
    };
data.stim_type = {
    'speed';
    'speed';
    };

if ismac
    cd ~/Dropbox/Salk/GCaMP/Data/
elseif ispc
    cd C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data
end

save('Control_CNO_Speed_Crh_List', 'data');

%% SF TF Direction Stimulus
data = struct();

data.animal_ids = {
    '541';
    '447';
    };
data.file_names = {
    {'day2_000_001_002_sf_tf_direction_control.mat','day2_003_004_005_sf_tf_direction_cno.mat'};
    {'day2_000_001_002_sf_tf_direction_control.mat','day2_003_004_005_sf_tf_direction_cno.mat'};
    };
data.stim_type = {
    'sf_tf_direction';
    'sf_tf_direction';
    };

if ismac
    cd ~/Dropbox/Salk/GCaMP/Data/
elseif ispc
    cd C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data
end

save('Control_CNO_SF_TF_Direction_Crh_List', 'data');