clc
clear all
close all

%% ~~~~~~~~~~ specify the path
dataset_root = fullfile(pwd); % /Volumes/mumi-ext-001/mumi-data/Project-Glioma-FrequencyFluctuations/Codes
addpath(genpath(fullfile(pwd)))
addpath(genpath(fullfile(dataset_root,'Tools/DPABI_V7.0_230110/DPABI_V7.0_230110/Subfunctions')));

% path for the folders
dataset_root = dataset_root(1:end-6);
code_path   = fullfile(dataset_root,'Codes');
data_path   = fullfile(dataset_root,'Data');
info_path   = fullfile(dataset_root,'Infos');
atlas_path  = fullfile(dataset_root,'Atlas');
result_path  = fullfile(dataset_root,'Results');

% extract all the files in the data folder
data_list   = dir(data_path);
data_list   = data_list(4:end); % first three rows are unnecessary
 
% Load packages
addpath(genpath(fullfile('/mnt/storage/tier2/MUMI-EXT-001/mumi-data','USERS','Manuela_Moretto','TOOLS','MATLAB_UTILS','NIfTI_20140122')))
addpath(genpath(fullfile('/mnt/storage/tier2/MUMI-EXT-001/mumi-data','USERS','Manuela_Moretto','TOOLS','MATLAB_UTILS','github_repo')))
addpath(genpath(pwd));

%% Compute static ALFF
% loop for the subjects
for jj = 1:length(data_list) %1:57

    % Get the path for anatomy and function for one subject
    SUBJ = data_list(jj).name;
    subj_folder = fullfile(data_path,SUBJ); 
    anat_folder = dir(fullfile(subj_folder,'*.anat'));
    anat_path = fullfile(subj_folder,anat_folder.name);
    func_folder = dir(fullfile(subj_folder,'*.fMRI'));
    func_path = fullfile(subj_folder,func_folder.name);
   
    % Define the path of functional data
    preproc_func_name = dir(fullfile(func_path,'*nonBPF-MNI-smooth.*'));
    preproc_func_name = preproc_func_name(end).name; % only last row is necessary
    
    % read functional data
    nifti_file = fullfile(func_path,preproc_func_name); % specify the data path
    % preproc_func       = niftiread(nifti_file);
    %preproc_func_header = spm_vol(nifti_file(1:end-3)); % get the header from nii.file, not nii.gz.

    new_result_path = fullfile(result_path,SUBJ);
    mkdir(new_result_path);
  
    % ~~~~~~~~~ Compute static ALFF
    ASamplePeriod = 2.6; % TR(s)
    HighCutoff = 0.1; % low pass
    LowCutoff = 0.01; % high pass
    AMaskFilename = fullfile(anat_path,'T1_gmedge_to_MNI_nonlin.nii.gz'); % GM mask
    AResultFilename = fullfile(new_result_path,'ALFF'); 
    TemporalMask = [];
    ScrubbingMethod = [];
    CUTNUMBER = 10;
    [ALFFBrain, fALFFBrain, Header] = y_alff_falff_KM(nifti_file,ASamplePeriod, HighCutoff, LowCutoff, AMaskFilename, AResultFilename, TemporalMask, ScrubbingMethod, CUTNUMBER);
end


%% ~~~~~~~~~ Compute dynamic ALFF 

% ~~~~ get the HGG index, read info.excel file
excel_info_path = fullfile(info_path,"lnif_prep_status.xlsx");
excel_info_data = readtable(excel_info_path);
excel_info_data(58:end,:) = [];
str_HGG = 'LGG';
str_exclusion = 'DeTe_93';
HGG_index = strcmp(excel_info_data.Tumor_grade, str_HGG);
exclusion_index = ~strcmp(excel_info_data.ID, str_exclusion);
HGG_exclusion_index = HGG_index & exclusion_index;
data_list = data_list(HGG_exclusion_index,:);

%%
for SUB_LOOP = 11:length(data_list) % loop for subject 1:41

    % Get the path for anatomy and function for one subject
    SUBJ = data_list(SUB_LOOP).name;
    %
    subj_folder = fullfile(data_path,SUBJ); 
    anat_folder = dir(fullfile(subj_folder,'*.anat'));
    anat_path = fullfile(subj_folder,anat_folder.name);
    func_path = fullfile(subj_folder,'func','sdc_unc');
   
    % ~~~~~~Define the path of functional data
    preproc_func_name = dir(fullfile(func_path,'*nonBPF-MNI-smooth.nii'));
    preproc_func_name = preproc_func_name(end).name; % only last row is necessary
    
    % read functional data
    nifti_file = fullfile(func_path,preproc_func_name); % specify the data path
    preproc_dfunc = niftiread(nifti_file); % data for dynamic
    preproc_dfunc_header = spm_vol(nifti_file);

    % specify the parameters
    window_size = 23-1; % time window=23, -1 is necessary give for loop
    lastTR = 270-23; % 247, necessary for loop

    ASamplePeriod = 2.6; % TR(s)
    HighCutoff = 0.1; % low pass
    LowCutoff = 0.01; % high pass
    AMaskFilename = fullfile(anat_path,'T1_gmedge_to_MNI_nonlin.nii.gz'); % GM mask
    TemporalMask = [];
    ScrubbingMethod = [];
    CUTNUMBER = 10;

    % create file for dALFF
    dALFF_path = fullfile(result_path,SUBJ,'dALFF');
    mkdir(dALFF_path); % create a folder (dALFF) in the results folder


    for step_size = 2:3 % loop for step size
        % ~~~~~~~~~ 2 step size 
        if step_size == 2
            stepsize_path = fullfile(dALFF_path,sprintf('Stepsize_%d',step_size));
            mkdir(stepsize_path); % create a folder (stepsize_2) in the dALFF folder
            disp(['Folder "', step_size, '" created at: ', stepsize_path]);
            for ii = 1:step_size:lastTR
                preproc_dfunc_window = preproc_dfunc(:,:,:,ii:ii+window_size); % get the data for each window
                AResultFilename = fullfile(stepsize_path,sprintf('dALFF_%d',ii)); % create new file for each window
                Header = preproc_dfunc_header(ii);
                [ALFFBrain, fALFFBrain, Header] = y_alff_falff_KM(preproc_dfunc_window,ASamplePeriod, HighCutoff, LowCutoff, AMaskFilename, AResultFilename, TemporalMask, ScrubbingMethod, Header, CUTNUMBER);
                disp(['Iteration',num2str(ii)])
            end
        end
    
        % ~~~~~~~~~ 3 step size
        if step_size == 3
            stepsize_path = fullfile(dALFF_path,sprintf('Stepsize_%d',step_size));
            mkdir(stepsize_path); % create a folder (stepsize_3) in the dALFF folder
            disp(['Folder "', step_size, '" created at: ', stepsize_path]);
            for ii = 1:step_size:lastTR
                preproc_dfunc_window = preproc_dfunc(:,:,:,ii:ii+window_size); % get the data for each window
                AResultFilename = fullfile(stepsize_path,sprintf('dALFF_%d',ii)); % create new file for each window
                Header = preproc_dfunc_header(ii);
                [ALFFBrain, fALFFBrain, Header] = y_alff_falff_KM(preproc_dfunc_window,ASamplePeriod, HighCutoff, LowCutoff, AMaskFilename, AResultFilename, TemporalMask, ScrubbingMethod, Header, CUTNUMBER);
                disp(['Iteration',num2str(ii)])
            end
        end

%     % ~~~~~~~~~ 4 step size
%     if step_size == 4
%         stepsize_path = fullfile(dALFF_path,sprintf('Stepsize_%d',step_size));
%         mkdir(stepsize_path); % create a folder (stepsize_4) in the dALFF folder
%         disp(['Folder "', step_size, '" created at: ', stepsize_path]);
%         for ii = 1:step_size:lastTR
%             preproc_dfunc_window = preproc_dfunc(:,:,:,ii:ii+window_size); % get the data for each window
%             AResultFilename = fullfile(stepsize_path,sprintf('dALFF_%d',ii)); % create new file for each window
%             Header = preproc_func_header(ii);
%             [ALFFBrain, fALFFBrain, Header] = y_alff_falff(preproc_dfunc_window,ASamplePeriod, HighCutoff, LowCutoff, AMaskFilename, AResultFilename, TemporalMask, ScrubbingMethod, Header, CUTNUMBER);
%             disp(['Iteration',num2str(ii)])
%         end
%     end
%     % ~~~~~~~~~ 5 step size
%     if step_size == 5
%         stepsize_path = fullfile(dALFF_path,sprintf('Stepsize_%d',step_size));
%         mkdir(stepsize_path); % create a folder (stepsize_5) in the dALFF folder
%         disp(['Folder "', step_size, '" created at: ', stepsize_path]);
%         for ii = 1:step_size:lastTR
%             preproc_dfunc_window = preproc_dfunc(:,:,:,ii:ii+window_size); % get the data for each window
%             AResultFilename = fullfile(stepsize_path,sprintf('dALFF_%d',ii)); % create new file for each window
%             Header = preproc_func_header(ii);
%             [ALFFBrain, fALFFBrain, Header] = y_alff_falff(preproc_dfunc_window,ASamplePeriod, HighCutoff, LowCutoff, AMaskFilename, AResultFilename, TemporalMask, ScrubbingMethod, Header, CUTNUMBER);
%             disp(['Iteration',num2str(ii)])

    end % end of step size loop

    disp(['Done ',SUBJ])
end % Subject loop


%% ~~~~~~~ Normalization

% ~~~~ ALFF&fALFF  % each voxel / global mean (GM mask)
for ii = 1:length(data_list) % 1:57 = subject loop

    % Get each subject name
    SUB_filename = data_list(ii).name;
    SUB_result_path = fullfile(result_path,SUB_filename);

    % read nii.file
    % ALFF
    ALFF_path = fullfile(SUB_result_path,'ALFF.nii');
    staticALFF_header = niftiinfo(ALFF_path);
    df_staticALFF = niftiread(ALFF_path);
    %fALFF
    fALFF_path = fullfile(SUB_result_path,'fALFF.nii');
    staticfALFF_header = niftiinfo(fALFF_path);
    df_staticfALFF = niftiread(fALFF_path);
    
    % compute mean and devide by it
    % ALFF
    Nonzeroindex_ALFF = (df_staticALFF ~= 0);
    GloMean_ALFF = mean(df_staticALFF(Nonzeroindex_ALFF));
    normal_df_staticALFF = df_staticALFF./GloMean_ALFF;
    % fALFF
    Nonzeroindex_fALFF = (df_staticfALFF ~= 0);
    GloMean_fALFF = mean(df_staticfALFF(Nonzeroindex_fALFF));
    normal_df_staticfALFF = df_staticfALFF./GloMean_fALFF;
    
    % nii.file output
    % ALFF
    filename = 'normalized_ALFF.nii';
    fullpath = fullfile(SUB_result_path,filename);
    niftiwrite(normal_df_staticALFF,fullpath,staticALFF_header);
    % fALFF
    filename = 'normalized_fALFF.nii';
    fullpath = fullfile(SUB_result_path,filename);
    niftiwrite(normal_df_staticfALFF,fullpath,staticfALFF_header);
    
    %in command, show the working subject
    disp(sprintf('%s is done', data_list(ii).name))
end % Subject loop

%% ~~~~ Normalization for dALFF&dfALFF  % each voxel / global mean (GM mask)
% Subject loop
% ~~~~ get the HGG index, read info.excel file
excel_info_path = fullfile(info_path,"lnif_prep_status.xlsx");
excel_info_data = readtable(excel_info_path);
excel_info_data(58:end,:) = [];
str_HGG = 'LGG';
str_exclusion = 'DeTe_93';
HGG_index = strcmp(excel_info_data.Tumor_grade, str_HGG);
exclusion_index = ~strcmp(excel_info_data.ID, str_exclusion);
HGG_exclusion_index = HGG_index & exclusion_index;
data_list = data_list(HGG_exclusion_index,:);

for SUB_LOOP = 1:length(data_list)

    SUBJ = data_list(SUB_LOOP).name;
    subj_folder = fullfile(data_path,SUBJ)
    anat_folder = dir(fullfile(subj_folder,'*.anat'));
    anat_path = fullfile(subj_folder,anat_folder.name);
    func_path = fullfile(subj_folder,'func','sdc_unc');
    result_folder = fullfile(result_path,SUBJ); 
    
% Stepsize loop
    for STEP_LOOP = 2:3
        % ~~~~~ read nii.gz file
        % dALFF
        step_path_dynamic = sprintf('/dALFF/Stepsize_%d/4D_dALFF.nii.gz',STEP_LOOP);
        dALFF_step_path = fullfile(result_folder,step_path_dynamic);
        dynamicALFF_header = niftiinfo(dALFF_step_path);
        df_dynamicALFF = niftiread(dALFF_step_path);

        % fdALFF
        step_path_fdynamic = sprintf('/dALFF/Stepsize_%d/4D_fdALFF.nii.gz',STEP_LOOP);
        fdALFF_step_path = fullfile(result_folder,step_path_fdynamic);
        fdynamicALFF_header = niftiinfo(fdALFF_step_path);
        df_fdynamicALFF = niftiread(fdALFF_step_path);
            
    
        % ~~~~~ compute mean and devide by it
        % ALFF
        normal_df_dALFF = [];
        for ii = 1:size(df_dynamicALFF,4) %1:the number of window
            current_df_dynamicALFF = df_dynamicALFF(:,:,:,ii);
            Nonzeroindex_dALFF = (current_df_dynamicALFF ~= 0);
            GloMean_dALFF = mean(current_df_dynamicALFF(Nonzeroindex_dALFF));
            normal_df_dALFF(:,:,:,ii) = current_df_dynamicALFF./GloMean_dALFF;
        end
         normal_df_dALFF = single(normal_df_dALFF);

        % fALFF
        normal_df_fdALFF = [];
        for ii = 1:size(df_fdynamicALFF,4) %1:the number of window
            current_df_fdynamicALFF = df_fdynamicALFF(:,:,:,ii);
            Nonzeroindex_fdALFF = (current_df_fdynamicALFF ~= 0);
            GloMean_fdALFF = mean(current_df_fdynamicALFF(Nonzeroindex_fdALFF));
            normal_df_fdALFF(:,:,:,ii) = current_df_fdynamicALFF./GloMean_fdALFF;
        end
        normal_df_fdALFF = single(normal_df_fdALFF);
    
        % ~~~~~ nii.file output
        step_path = sprintf('/dALFF/Stepsize_%d',STEP_LOOP);
        step_path = fullfile(result_path,SUBJ,step_path);
        % ALFF
        filename_dALFF = 'normalized_4D_dALFF.nii';
        fullpath_dALFF = fullfile(step_path,filename_dALFF);
        niftiwrite(normal_df_dALFF,fullpath_dALFF,dynamicALFF_header);
        % fALFF
        filename_fdALFF = 'normalized_4D_fdALFF.nii';
        fullpath_fdALFF = fullfile(step_path,filename_fdALFF);
        niftiwrite(normal_df_fdALFF,fullpath_fdALFF,fdynamicALFF_header);

    end % stepsize end
end % subject loop end



%% ~~~~~~~~ Each time course of dALFF and fdALFF

% ~~~~ get the HGG index, read info.excel file
excel_info_path = fullfile(info_path,"lnif_prep_status.xlsx");
excel_info_data = readtable(excel_info_path);
excel_info_data(58:end,:) = [];
% str_HGG = 'HGG';
% HGG_index = strcmp(excel_info_data.Tumor_grade, str_HGG);
% HGG_exclusion_index = HGG_index & exclusion_index;
str_exclusion = {'ToF_02_54','DeTe_93'};
exclusion_index = ~ismember(excel_info_data.ID, str_exclusion);
data_list = data_list(exclusion_index,:);
excel_info_data = excel_info_data(exclusion_index,:);

% Functional networks template
parcellation_path = fullfile(atlas_path, 'Schaefer2018_200Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz');
parcellation_data = niftiread(parcellation_path);
% Functional networks csv file
FNcsv_path = fullfile(atlas_path,'Schaefer2018_200Parcels_Kong2022_17Networks_order.csv');
FNcsv_data = readmatrix(FNcsv_path);
FNcsv_data_table = readtable(FNcsv_path);
FNcsv_data(101:200,7) = FNcsv_data(101:200,7)+17;
FNcsv_data_table.NETWORKID(101:200) = FNcsv_data_table.NETWORKID(101:200)+17;


df_timecourse = excel_info_data;


% read normalized_4D_dALFF and fdALFF.nii
for SUB_LOOP = 1:size(df_timecourse,1)

    SUBJ = data_list(SUB_LOOP).name;
    subj_folder = fullfile(data_path,SUBJ);
    anat_folder = dir(fullfile(subj_folder,'*.anat'));
    anat_path = fullfile(subj_folder,anat_folder.name);
    func_path = fullfile(subj_folder,'func','sdc_unc');
    result_folder = fullfile(result_path,SUBJ); 
    
    % ~~~~~~~ import dALFF map, FN template,network csv file    
    % import dALFF (4D map)
    step_path = fullfile(result_folder,'dALFF','Stepsize_2','normalized_4D_fdALFF.nii.gz');
    img_data = niftiread(step_path);



% ~~~~~~ extract voxles corresponding 200 parcels

    for WINDOW_LOOP = 1:size(img_data,4)
        % get one window
        df = img_data(:,:,:,WINDOW_LOOP);
    
        garage_parcel = {};
        tmp = unique(parcellation_data(:));
        tmp = tmp(tmp~=0);
        for current_parcel_index = 1:max(tmp) % 1:200
            parcel_mask = (parcellation_data == current_parcel_index); % Extract the mask for the current parcel   
            parcel_data = df(parcel_mask);
            parcel_data = parcel_data(parcel_data ~= 0); % remove 0
            garage_parcel{current_parcel_index,1} = parcel_data; % Store parcel_data
            garage_parcel{current_parcel_index,2} = mean(parcel_data); % store mean of each parcel
            parcel_names = FNcsv_data_table.ROIName(current_parcel_index);
            garage_parcel{current_parcel_index,3} = parcel_names;
        end
        ALLWINDOW_parcel{WINDOW_LOOP,1} = garage_parcel;
    
    
        % ~~~~~~ classify 200 parcels into 17*2(Each hem) networks
        garage_network = {};
        netmp = unique(FNcsv_data(:,7));
    
        for network_index = 1:max(netmp) % 1:34
            network_mask = (FNcsv_data(:,7) == network_index); % logical array
            network_data = garage_parcel(network_mask,1); % voxels for a network
            network_data = cell2mat(network_data);
            % store into variables
            garage_network{network_index,1} = data_list(SUB_LOOP).name; % column 1, values for a network 
            garage_network{network_index,2} = mean(network_data); % Column 2, mean of a network
            networknames = FNcsv_data_table.NETWORKNAME(FNcsv_data_table.NETWORKID == network_index); % extract network name
            garage_network{network_index,3} = networknames(1); % column 3, each network name
            if network_index < 18
               garage_network{network_index,4} = 'LEFT';
                if strcmp(excel_info_data.Tumor_lateralization(SUB_LOOP), 'LEFT')
                   garage_network{network_index,5} = 'Ipsilesional';
                else
                   garage_network{network_index,5} = 'Contralesional';
                end
            else
                garage_network{network_index,4} = 'RIGHT';
                if strcmp(excel_info_data.Tumor_lateralization(SUB_LOOP), 'RIGHT')
                   garage_network{network_index,5} = 'Ipsilesional';
                else
                   garage_network{network_index,5} = 'Contralesional';
                end
            end
        end
        ALLWINDOW_networks{WINDOW_LOOP,1} = garage_network;
    
    end % Window_LOOP end

df_timecourse.dALFF(SUB_LOOP) = {ALLWINDOW_networks};
disp(sprintf('%s is done', data_list(SUB_LOOP).name))
end %SUB_LOOP end
 

filename = 'ALLSUB_fdALFF_34FN';
save(filename, 'df_timecourse');


%% ~~~~~~~~ 17 Functional networks_ALFF
% ~~~~ read info.excel file
excel_info_path = fullfile(info_path,"lnif_prep_status.xlsx");
excel_info_data = readtable(excel_info_path);
excel_info_data(58:end,:) = [];
% str_HGG = 'HGG';
% HGG_index = strcmp(excel_info_data.Tumor_grade, str_HGG);
% HGG_exclusion_index = HGG_index & exclusion_index;
str_exclusion = {'ToF_02_54','DeTe_93'};
exclusion_index = ~ismember(excel_info_data.ID, str_exclusion);
data_list = data_list(exclusion_index,:);


for SUB_LOOP = 1:length(data_list)

    SUBJ = data_list(SUB_LOOP).name;
    subj_folder = fullfile(data_path,SUBJ)
    anat_folder = dir(fullfile(subj_folder,'*.anat'));
    anat_path = fullfile(subj_folder,anat_folder.name);
    func_path = fullfile(subj_folder,'func','sdc_unc');
    result_folder = fullfile(result_path,SUBJ); 

    
    % ~~~~~~~ import ALFF map, FN template,network csv file    
    % import ALFF (3D map)
    step_path = fullfile(result_folder,'normalized_fALFF.nii.gz');
    img_data = niftiread(step_path);
    % Functional networks template
    parcellation_path = fullfile(atlas_path, 'Schaefer2018_200Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz');
    parcellation_data = niftiread(parcellation_path);
    % Functional networks csv file
    FNcsv_path = fullfile(atlas_path,'Schaefer2018_200Parcels_Kong2022_17Networks_order.csv');
    FNcsv_data = readmatrix(FNcsv_path);
    FNcsv_data_table = readtable(FNcsv_path);
    
    
    % ~~~~~~ extract voxles corresponding 200 parcels
    garage_parcel = {};
    tmp = unique(parcellation_data(:));
    tmp = tmp(tmp~=0);
    for current_parcel_index = 1:max(tmp) % 1:200
        parcel_mask = (parcellation_data == current_parcel_index); % Extract the mask for the current parcel   
        parcel_data = img_data(parcel_mask);
        parcel_data = parcel_data(parcel_data ~= 0); % remove 0
        garage_parcel{current_parcel_index,1} = parcel_data; % Store parcel_data
        garage_parcel{current_parcel_index,2} = mean(parcel_data); % store mean of each parcel
        parcel_names = FNcsv_data_table.ROIName(current_parcel_index);
        garage_parcel{current_parcel_index,3} = parcel_names;
    end
    ALLSUB_parcel{SUB_LOOP,1} = garage_parcel;

    % ~~~~~~ classify 200 parcels into 17 netowrks
    garage_network = {};
    netmp = unique(FNcsv_data(:,7));
    for network_index = 1:max(netmp); % 1:17
        network_mask = (FNcsv_data(:,7) == network_index); % logical array
        network_data = garage_parcel(network_mask,1); % voxels for a network
        network_data = cell2mat(network_data);
        % store into variables
        garage_network{network_index,1} = network_data; % column 1, values for a network 
        garage_network{network_index,2} = mean(network_data); % Column 2, mean of a network
        networknames = FNcsv_data_table.NETWORKNAME(FNcsv_data_table.NETWORKID == network_index); % extract network name
        garage_network{network_index,3} = networknames(1); % column 3, each network name
    end
    ALLSUB_networks{SUB_LOOP,1} = garage_network;
 

    disp(sprintf('%s is done',data_list(SUB_LOOP).name))
end % SUB_LOOP end

filename = 'ALLSUB_fALFF_200parcels';
save(filename,'ALLSUB_parcel')
% save(filename,'ALLSUB_parcel','ALLSUB_networks') % List all variables to save

%% ~~~~~~~~ 17 Functional networks_dALFF
% ~~~~ get the HGG index, read info.excel file
excel_info_path = fullfile(info_path,"lnif_prep_status.xlsx");
excel_info_data = readtable(excel_info_path);
excel_info_data(58:end,:) = [];
% str_HGG = 'HGG';
% HGG_index = strcmp(excel_info_data.Tumor_grade, str_HGG);
% HGG_exclusion_index = HGG_index & exclusion_index;
str_exclusion = {'ToF_02_54','DeTe_93'};
exclusion_index = ~ismember(excel_info_data.ID, str_exclusion);
data_list = data_list(exclusion_index,:);


for SUB_LOOP = 1:length(data_list)

    SUBJ = data_list(SUB_LOOP).name;
    subj_folder = fullfile(data_path,SUBJ)
    anat_folder = dir(fullfile(subj_folder,'*.anat'));
    anat_path = fullfile(subj_folder,anat_folder.name);
    func_path = fullfile(subj_folder,'func','sdc_unc');
    result_folder = fullfile(result_path,SUBJ); 

    % stepsize loop
    for STEP_LOOP = 2:3
        
        % ~~~~~~~ import dALFF avg map, FN template,network csv file    
        % import std of dALFF (3D map)
        step_path = sprintf('/dALFF/Stepsize_%d/nmz_std_fdALFF.nii.gz',STEP_LOOP);
        step_path = fullfile(result_folder,step_path);
        img_data = niftiread(step_path);
        % Functional networks template
        parcellation_path = fullfile(atlas_path, 'Schaefer2018_200Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz');
        parcellation_data = niftiread(parcellation_path);
        % Functional networks csv file
        FNcsv_path = fullfile(atlas_path,'Schaefer2018_200Parcels_Kong2022_17Networks_order.csv');
        FNcsv_data = readmatrix(FNcsv_path);
        FNcsv_data_table = readtable(FNcsv_path);
        
        
        % ~~~~~~ extract voxles corresponding 200 parcels
        garage_parcel = {};
        tmp = unique(parcellation_data(:));
        tmp = tmp(tmp~=0);
        for current_parcel_index = 1:max(tmp) % 1:200
            parcel_mask = (parcellation_data == current_parcel_index); % Extract the mask for the current parcel   
            parcel_data = img_data(parcel_mask);
            parcel_data = parcel_data(parcel_data ~= 0); % remove 0
            garage_parcel{current_parcel_index,1} = parcel_data; % Store parcel_data
            garage_parcel{current_parcel_index,2} = mean(parcel_data); % store mean of each parcel
            parcel_names = FNcsv_data_table.ROIName(current_parcel_index);
            garage_parcel{current_parcel_index,3} = parcel_names;
        end
        ALLSUB_parcel{SUB_LOOP,[STEP_LOOP-1]} = garage_parcel;
    
        % ~~~~~~ classify 200 parcels into 17 netowrks
        garage_network = {};
        netmp = unique(FNcsv_data(:,7));
        for network_index = 1:max(netmp); % 1:17
            network_mask = (FNcsv_data(:,7) == network_index); % logical array
            network_data = garage_parcel(network_mask,1); % voxels for a network
            network_data = cell2mat(network_data);
            % store into variables
            garage_network{network_index,1} = network_data; % column 1, values for a network 
            garage_network{network_index,2} = mean(network_data); % Column 2, mean of a network
            networknames = FNcsv_data_table.NETWORKNAME(FNcsv_data_table.NETWORKID == network_index); % extract network name
            garage_network{network_index,3} = networknames(1); % column 3, each network name
        end
        ALLSUB_networks{SUB_LOOP,[STEP_LOOP-1]} = garage_network;
     
    end % STEP_LOOP end

    disp(sprintf('%s is done',data_list(SUB_LOOP).name))
end % SUB_LOOP end


filename = 'ALLSUB_fdALFF_200parcels';
save(filename,'ALLSUB_parcel') % List all variables to save






% Below, I don't use




%% ~~~~~ make a table for plot
varNames = {'dALFF','Networks'};
Y = [];
X = [];
for ss = 1:17
    ALFF_value = squeeze(garage_network{ss,1});
    Network_name = repmat({garage_network{ss,3}},length(ALFF_value), 1);
    Y = [Y; ALFF_value];
    X = [X; Network_name];
end
TableforPlot = table(Y,X, 'VariableNames', varNames);


% ~~~~~ plot distribution
functionalnetworkname = string(garage_network(:,3))';
x = categorical(string(TableforPlot.Networks),functionalnetworkname);
y = TableforPlot.dALFF;
ax1 = nexttile;
swarmchart(ax1,x,y,20,x,'.','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('dALFF');
title(sprintf('dALFF (stepsize %d)',STEP_LOOP),'FontSize',20);

%% ~~~~ ALLSUB 17 FN Average of static ALFF&fALFF 

% 17 functional networks
% Functionan networks template
parcellation_path = fullfile(atlas_path, 'Schaefer2018_200Parcels_17Networks_order_FSLMNI152_2mm.nii.gz');
parcellation_data = niftiread(parcellation_path);
% Functional networks csv file
FNcsv_path = fullfile(atlas_path,'Schaefer2018_200Parcels_17Networks_order_FSLMNI152_2mm.Centroid_RAS.csv');
FNcsv_data = readmatrix(FNcsv_path);
FNcsv_data_table = readtable(FNcsv_path);

% ~~~~ Extract the 17 Fcuntioan netowork value for each sub
Garage_ALLSUB_parcel = [];
Garage_ALLSUB_network = [];
for ii = 1:length(data_list) % 1:57 = subject loop

    % Get each subject name
    SUB_filename = data_list(ii).name;
    SUB_result_path = fullfile(result_path,SUB_filename);

    % read nii.file
    % ALFF
    ALFF_path = fullfile(SUB_result_path,'normalized_ALFF.nii');
    staticALFF_header = niftiinfo(ALFF_path);
    df_staticALFF = niftiread(ALFF_path);

    %fALFF
    fALFF_path = fullfile(SUB_result_path,'normalized_fALFF.nii');
    staticfALFF_header = niftiinfo(fALFF_path);
    df_staticfALFF = niftiread(fALFF_path);

    % ~~~~~~ extract voxles corresponding 200 parcels
    garage_parcel = {};
    tmp = unique(parcellation_data(:));
    tmp = tmp(tmp~=0);
    for current_parcel_index = 1:max(tmp) % 1:200
        parcel_mask = (parcellation_data == current_parcel_index); % Extract the mask for the current parcel   
        parcel_data = df_staticALFF(parcel_mask);
        parcel_data = parcel_data(parcel_data ~= 0); % remove 0
        garage_parcel{current_parcel_index,1} = parcel_data; % Store parcel_data
        garage_parcel{current_parcel_index,2} = mean(parcel_data); % store mean of each parcel  
    end
    % store each subject data into ALLSUB variable
    Garage_ALLSUB_parcel{ii,1} = garage_parcel;


    % ~~~~~~ classify 200 parcels into 17 netowrks
    garage_network = {};
    netmp = unique(FNcsv_data(:,7));
    for network_index = 1:max(netmp); % 1:17
        network_mask = (FNcsv_data(:,7) == network_index); % logical array
        network_data = garage_parcel(network_mask,1); % voxels for a network
        network_data = cell2mat(network_data);
        % store into variables
        garage_network{network_index,1} = network_data; % column 1, values for a network 
        garage_network{network_index,2} = mean(network_data); % Column 2, mean of a network
        networknames = FNcsv_data_table.NETWORKNAME(FNcsv_data_table.NETWORKID == network_index); % extract network name
        garage_network{network_index,3} = networknames(1); % column 3, each network name
    end
    % store each subject data into ALLSUB variable
    Garage_ALLSUB_network{ii,1} = garage_network;
    
    % outputing the process
    disp(sprintf('%s is done',data_list(ii).name))

  
end % Subject loop


%% ~~~~~~~~ Plot 17 FN Average of static ALFF&fALFF
% Prepare the data and variable
HGG_string = {'HGG'};
LGG_string = {'LGG'};
excel_info_path = fullfile(info_path,"lnif_prep_status.xlsx");
excel_info_data = readtable(excel_info_path);
excel_info_data(58:end,:) = [];

% Get HGG index
HGG_patient = strcmp(excel_info_data.Tumor_grade, 'HGG');
HGG_patient = find(HGG_patient);

% Make the vector containing (1) mean for each FN (2) FN name (3) Glioma grade
df_ALL = [];
for ii = 1:length(Garage_ALLSUB_network) % 1:57 subject loop
    df = Garage_ALLSUB_network{ii, 1}(:,2:3); % Get the mean and corresponding FN
    if any(HGG_patient == ii) % True if HGG patient
        glioma_name = repmat(HGG_string,17,1); % Write HGG for 17 times
    else % if LGG
        glioma_name = repmat(LGG_string,17,1);% Write LGG for 17 times
    end
    df(:,3) = glioma_name % Create the 3rd column = (3) glioma grade

    % Store each subject into ALL subject variable
    df_ALL = [df_ALL; df];

end % end subject loop (1:57)

% ~~~~~ Prepare the dataset for plot
x = categorical(string(df_ALL(:,2)));
y = cell2mat(df_ALL(:,1));
Glioma_Index = categorical(string(df_ALL(:,3)));

% Plot
ax1 = nexttile;
swarmchart(ax1,x,y,50,Glioma_Index,'filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor',[0 .5 .5]);
ylabel('static ALFF mean per subjects');
xlabel('Functional networks');
title('57 Glioma Patients (HGG vs. LGG)', 'FontSize',15);

%% ~~~~~ dALFF&fdALFF  % each volxel for each window / global mean (GM mask) for each window

% Stepsize loop
for ss = 2:5

    % ~~~~~ read nii.file
    % dALFF
    step_path_dynamic = sprintf('/dALFF/Stepsize_%d/4D_dALFF.nii.gz',ss);
    dALFF_path = fullfile(result_path,step_path_dynamic);
    dynamicALFF_header = niftiinfo(dALFF_path);
    df_dynamicALFF = niftiread(dALFF_path);
    %fdALFF
    step_path_fdynamic = sprintf('/dALFF/Stepsize_%d/4D_fdALFF.nii.gz',ss);
    fdALFF_path = fullfile(result_path,step_path_fdynamic);
    fdynamicALFF_header = niftiinfo(fdALFF_path);
    df_fdynamicALFF = niftiread(fdALFF_path);
    

    % ~~~~~ compute mean and devide by it
    % ALFF
    normal_df_dALFF = [];
    for ii = 1:size(df_dynamicALFF,4) %1:the number of window
        current_df_dynamicALFF = df_dynamicALFF(:,:,:,ii);
        Nonzeroindex_dALFF = (current_df_dynamicALFF ~= 0);
        GloMean_dALFF = mean(current_df_dynamicALFF(Nonzeroindex_dALFF));
        normal_df_dALFF(:,:,:,ii) = current_df_dynamicALFF./GloMean_dALFF;
    end
     normal_df_dALFF = single( normal_df_dALFF);
    % fALFF
     normal_df_fdALFF = [];
    for ii = 1:size(df_fdynamicALFF,4) %1:the number of window
        current_df_fdynamicALFF = df_fdynamicALFF(:,:,:,ii);
        Nonzeroindex_fdALFF = (current_df_fdynamicALFF ~= 0);
        GloMean_fdALFF = mean(current_df_fdynamicALFF(Nonzeroindex_fdALFF));
        normal_df_fdALFF(:,:,:,ii) = current_df_fdynamicALFF./GloMean_fdALFF;
    end
     normal_df_fdALFF = single(normal_df_fdALFF);

    % ~~~~~ nii.file output
    step_path = sprintf('/dALFF/Stepsize_%d',ss);
    step_path = fullfile(result_path,step_path)
    % ALFF
    filename_dALFF = 'normalized_4D_dALFF.nii';
    fullpath_dALFF = fullfile(step_path,filename_dALFF);
    niftiwrite(normal_df_dALFF,fullpath_dALFF,dynamicALFF_header);
    % fALFF
    filename_fdALFF = 'normalized_4D_fdALFF.nii';
    fullpath_fdALFF = fullfile(step_path,filename_fdALFF);
    niftiwrite(normal_df_fdALFF,fullpath_fdALFF,fdynamicALFF_header);

end

%% ~~~~~~~ Compute the standard deviation of dynamic ALFF and FALFF
% read nii.file
% step2_path = fullfile(result_path,'dALFF','Stepsize_2','4D_dALFF.nii.gz');
% step2_dALFF = niftiread(step2_path);
% step2_dALFF_header = niftiinfo(step2_path);
% % extract two voxels
% Onevoxel = step2_dALFF(45, 48, 50,:);
% Onevoxel = Onevoxel(:);
% Othervoxel = step2_dALFF(47, 49, 50,:);
% Othervoxel = Othervoxel(:);
% % calculate SD
% std_Onevoxel = std(Onevoxel);
% std_Othervoxel = std(Othervoxel);
% % plot time series
% time = 1:step2_dALFF_header.ImageSize(4); % get the number of slices
% figure;
% plot(time,[Onevoxel,Othervoxel])
% xlabel('TR');
% ylabel('ALFF');
% xlim([-1,127]);
% title('dALFF (2 step-size)','Fontsize', 16);
% legend('Voxel_n1','Voxel_n2');

