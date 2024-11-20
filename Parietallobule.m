% Parietal lobule analysis

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
result_path  = fullfile(dataset_root,'Results_ALLSUB');
csv_path   = fullfile(dataset_root,'CSV_file');

% extract all the files in the data folder
data_list   = dir(data_path);
data_list   = data_list(4:end); % first three rows are unnecessary



%% ~~~~~~~~~~ fALFF
% clearvars -except ALLSUB_dALFF_EachHemNetworks ALLSUB_ALFF_EachHemNetworks;
% Download the data
load('ALLSUB_ALFF_200parcels.mat')
% clearvars -except ALLSUB_parcel
% Merge the parcels into Functional Networks per hemisphere

% Functional networks csv file
FNcsv_path = fullfile(atlas_path,'Schaefer2018_200Parcels_Kong2022_parietallobule.csv');
FNcsv_data = readmatrix(FNcsv_path);
FNcsv_data_table = readtable(FNcsv_path);

%% ~~~~~~ Assign the index for functional networks for each hemisphere, from 1:34
% ~~~~~ Create new data structure from 55(SUB)*1 cell inside 34(Each network)*3(Parcels, NETWORK_Names, Network_ID) 
% Only stepsize2
% Transfrom from cell to table
NetworkNames = FNcsv_data_table.NETWORKNAME;
ALLSUB_parcels_table = {};
for SUB_LOOP = 1:length(ALLSUB_parcel)
    Sub_data = ALLSUB_parcel{SUB_LOOP,1};
    Sub_table = array2table(Sub_data, 'VariableNames', {'Parcels', 'Mean','NETWORK_Names'});
    Sub_table.Parietal_ID = FNcsv_data_table.ParietalID;
    ALLSUB_parcels_table{SUB_LOOP,1} = Sub_table;
end

%% ~~~~~ Merge percels for each hemisphere Network
% ~~~~~ ALL SUBJECTS
EachHemNetworks = [];
for SUB_LOOP = 1:length(ALLSUB_parcels_table)

    df = ALLSUB_parcels_table{SUB_LOOP, 1};
    x_vectors = [];
        % ~~~~~~ A SUBJECT 
        for NETWORK_LOOP = 1:12
            % OneNetwork_parcels_index = strcmp(string(df.Network_ID),'1');
            Network_parcels_index = df.Parietal_ID == NETWORK_LOOP;
            Network_parcels = df.Parcels(Network_parcels_index,:);
            
            for i = 1:length(Network_parcels)
                x = Network_parcels{i,1};
                x_vectors = [x_vectors; x];
                Network_marged_parcels{NETWORK_LOOP,1} = x_vectors;
            end
        end
        % ~~~~~~ END of A SUBJECT 

    EachHemNetworks{SUB_LOOP,1} = Network_marged_parcels;
end
% ~~~~~ END of ALL SUBJECTS


% ~~~~~ New data table similar to ALLSUB_parcels_table

% Change the name
ALLSUB_ALFF_ParietalLobule = EachHemNetworks;

% get the NETWORK_NAMES
uniqueNetworkArrays = uniqueOrderPreserved(NetworkNames)';
characterToAdd_LH = 'LH_'; % Specify the character to add
characterToAdd_RH = 'RH_'; % Specify the character to add
uniqueNetworkArrays_LH = cellfun(@(x) [characterToAdd_LH, x], uniqueNetworkArrays, 'UniformOutput', false);
uniqueNetworkArrays_RH = cellfun(@(x) [characterToAdd_RH, x], uniqueNetworkArrays, 'UniformOutput', false);
uniqueNetworkArrays = [uniqueNetworkArrays_LH;uniqueNetworkArrays_RH];
% modify the network name for parietal lobule
Int_Netname = ["LH_Control A", "LH_Control B", "LH_Salience/VenAttn A", "LH_Salience/VenAttn B"... 
    "LH_Dorsal Attention A", "LH_Dorsal Attention B", "RH_Control A", "RH_Control B", "RH_Salience/VenAttn A", "RH_Salience/VenAttn B"... 
    "RH_Dorsal Attention A", "RH_Dorsal Attention B"];
idx = ismember(uniqueNetworkArrays, Int_Netname);
uniqueNetworkArrays = uniqueNetworkArrays(idx,:);

%get the NETWORK_ID
uniqueIDArrays = 1:12;

for SUB_LOOP = 1:length(ALLSUB_ALFF_ParietalLobule)
    % make the table for each SUB
    x = array2table(EachHemNetworks{SUB_LOOP,1},"VariableNames",{'Parcels'});
    x.NETWORK_NAMES = uniqueNetworkArrays;
    x.NETWORK_ID = uniqueIDArrays';
    % Modify the data
    ALLSUB_ALFF_ParietalLobule{SUB_LOOP,1} = x;
end




%% ~~~~~~~~~~ dfALFF
% clearvars -except ALLSUB_dALFF_EachHemNetworks ALLSUB_ALFF_EachHemNetworks;
% Download the data
load('ALLSUB_dALFF_200parcels.mat')
% Merge the parcels into Functional Networks per hemisphere

% Functional networks csv file
FNcsv_path = fullfile(atlas_path,'Schaefer2018_200Parcels_Kong2022_parietallobule.csv');
FNcsv_data = readmatrix(FNcsv_path);
FNcsv_data_table = readtable(FNcsv_path);

%% ~~~~~~ Assign the index for functional networks for each hemisphere, from 1:34
% ~~~~~ Create new data structure from 55(SUB)*1 cell inside 34(Each network)*3(Parcels, NETWORK_Names, Network_ID) 
% Only stepsize2
% Transfrom from cell to table
NetworkNames = FNcsv_data_table.NETWORKNAME;
ALLSUB_parcels_table = {};
for SUB_LOOP = 1:length(ALLSUB_parcel)
    Sub_data = ALLSUB_parcel{SUB_LOOP,1};
    Sub_table = array2table(Sub_data, 'VariableNames', {'Parcels', 'Mean','NETWORK_Names'});
    Sub_table.Parietal_ID = FNcsv_data_table.ParietalID;
    ALLSUB_parcels_table{SUB_LOOP,1} = Sub_table;
end

%% ~~~~~ Merge percels for each hemisphere Network
% ~~~~~ ALL SUBJECTS
EachHemNetworks = [];
for SUB_LOOP = 1:length(ALLSUB_parcels_table)

    df = ALLSUB_parcels_table{SUB_LOOP, 1};
    x_vectors = [];
        % ~~~~~~ A SUBJECT 
        for NETWORK_LOOP = 1:12
            % OneNetwork_parcels_index = strcmp(string(df.Network_ID),'1');
            Network_parcels_index = df.Parietal_ID == NETWORK_LOOP;
            Network_parcels = df.Parcels(Network_parcels_index,:);
            
            for i = 1:length(Network_parcels)
                x = Network_parcels{i,1};
                x_vectors = [x_vectors; x];
                Network_marged_parcels{NETWORK_LOOP,1} = x_vectors;
            end
        end
        % ~~~~~~ END of A SUBJECT 

    EachHemNetworks{SUB_LOOP,1} = Network_marged_parcels;
end
% ~~~~~ END of ALL SUBJECTS


% ~~~~~ New data table similar to ALLSUB_parcels_table

% Change the name
ALLSUB_dALFF_ParietalLobule = EachHemNetworks;

% get the NETWORK_NAMES
uniqueNetworkArrays = uniqueOrderPreserved(NetworkNames)';
characterToAdd_LH = 'LH_'; % Specify the character to add
characterToAdd_RH = 'RH_'; % Specify the character to add
uniqueNetworkArrays_LH = cellfun(@(x) [characterToAdd_LH, x], uniqueNetworkArrays, 'UniformOutput', false);
uniqueNetworkArrays_RH = cellfun(@(x) [characterToAdd_RH, x], uniqueNetworkArrays, 'UniformOutput', false);
uniqueNetworkArrays = [uniqueNetworkArrays_LH;uniqueNetworkArrays_RH];
% modify the network name for parietal lobule
Int_Netname = ["LH_Control A", "LH_Control B", "LH_Salience/VenAttn A", "LH_Salience/VenAttn B"... 
    "LH_Dorsal Attention A", "LH_Dorsal Attention B", "RH_Control A", "RH_Control B", "RH_Salience/VenAttn A", "RH_Salience/VenAttn B"... 
    "RH_Dorsal Attention A", "RH_Dorsal Attention B"];
idx = ismember(uniqueNetworkArrays, Int_Netname);
uniqueNetworkArrays = uniqueNetworkArrays(idx,:);

%get the NETWORK_ID
uniqueIDArrays = 1:12;

for SUB_LOOP = 1:length(ALLSUB_dALFF_ParietalLobule)
    % make the table for each SUB
    x = array2table(EachHemNetworks{SUB_LOOP,1},"VariableNames",{'Parcels'});
    x.NETWORK_NAMES = uniqueNetworkArrays;
    x.NETWORK_ID = uniqueIDArrays';
    % Modify the data
    ALLSUB_dALFF_ParietalLobule{SUB_LOOP,1} = x;
end

%% Merge with the excel file
excel_info_path = fullfile(info_path,"lnif_prep_status.xlsx");
excel_info_data = readtable(excel_info_path);
excel_info_data(58:end,:) = [];
str_exclusion = {'ToF_02_54','DeTe_93'};
exclusion_index = ~ismember(excel_info_data.ID, str_exclusion);
ALLSUB_lists = excel_info_data(exclusion_index,:);

% Add fALFF and dfALFF score
ALLSUB_lists.ALFF = ALLSUB_ALFF_ParietalLobule;
ALLSUB_lists.dALFF = ALLSUB_dALFF_ParietalLobule;

df_name = fullfile(code_path,'ALLSUB_dALFF&ALFF_ParietalLobule.mat');
save(df_name, 'ALLSUB_lists');




%% ~~~~~~~~~~~~ Download the data
clear all
% ~~~~~~~~~~ specify the path
dataset_root = fullfile(pwd); % /Volumes/mumi-ext-001/mumi-data/Project-Glioma-FrequencyFluctuations/Codes
addpath(genpath(fullfile(pwd)))
addpath(genpath(fullfile(dataset_root,'Tools/DPABI_V7.0_230110/DPABI_V7.0_230110/Subfunctions')));

% path for the folders
dataset_root = dataset_root(1:end-6);
code_path   = fullfile(dataset_root,'Codes');
data_path   = fullfile(dataset_root,'Data');
info_path   = fullfile(dataset_root,'Infos');
atlas_path  = fullfile(dataset_root,'Atlas');
result_path  = fullfile(dataset_root,'Results_ALLSUB');

% extract all the files in the data folder
data_list   = dir(data_path);
data_list   = data_list(4:end); % first three rows are unnecessary

% Load the data
load('ALLSUB_dALFF&ALFF_ParietalLobule.mat') 


%% Extract the interested patients from 55 total

str_Loc = {'frontal','fronto-insulare'};
str_Lateral = {'LEFT','RIGHT'};
Loc_idx = ismember(ALLSUB_lists.Tumor_loc, str_Loc);
LEFT_idx = strcmp(ALLSUB_lists.Tumor_lateralization,str_Lateral(1));
RIGHT_idx = strcmp(ALLSUB_lists.Tumor_lateralization,str_Lateral(2));
combined_idx_LEFT = Loc_idx & LEFT_idx;
combined_idx_RIGHT = Loc_idx & RIGHT_idx;

df_LEFT = ALLSUB_lists(combined_idx_LEFT,:);
df_RIGHT = ALLSUB_lists(combined_idx_RIGHT,:);

%% Make a dALFF table for R
% LEFT frontal patients
Col_1 = [];
Col_2 = [];
Col_3 = [];
Col_4 = [];
Col_5 = [];
for SUB_LOOP = 1:size(df_LEFT,1)
    Mean_dALFF_LEFT = cellfun(@mean,df_LEFT.dALFF{SUB_LOOP,1}{:,"Parcels"});
    Col_1 = [Col_1; Mean_dALFF_LEFT];
    Col_2 = [Col_2; df_LEFT.dALFF{SUB_LOOP,1}{:,"NETWORK_NAMES"}];
    Col_3 = [Col_3; repmat({df_LEFT.ID{SUB_LOOP,1}},12,1)];
    Col_4 = [Col_4; repmat({df_LEFT.Tumor_lateralization{SUB_LOOP,1}},12,1)];
    for NETWORK_LOOP = 1:length(df_LEFT.dALFF{SUB_LOOP,1}{:,"NETWORK_NAMES"}) %1:34
        if strncmp(df_LEFT.dALFF{SUB_LOOP,1}{NETWORK_LOOP,"NETWORK_NAMES"},'L',1)
            Col_5 = [Col_5; {'Ipsilesional'}];
        else
            Col_5 = [Col_5; {'Contralesional'}];
        end
    end
end

% RIGHT frontal patients
for SUB_LOOP = 1:size(df_RIGHT,1)
    Mean_dALFF_RIGHT = cellfun(@mean,df_RIGHT.dALFF{SUB_LOOP,1}{:,"Parcels"});
    Col_1 = [Col_1; Mean_dALFF_RIGHT];
    Col_2 = [Col_2; df_RIGHT.dALFF{SUB_LOOP,1}{:,"NETWORK_NAMES"}];
    Col_3 = [Col_3; repmat({df_RIGHT.ID{SUB_LOOP,1}},12,1)];
    Col_4 = [Col_4; repmat({df_RIGHT.Tumor_lateralization{SUB_LOOP,1}},12,1)];
    for NETWORK_LOOP = 1:length(df_RIGHT.dALFF{SUB_LOOP,1}{:,"NETWORK_NAMES"}) %1:34
        if strncmp(df_RIGHT.dALFF{SUB_LOOP,1}{NETWORK_LOOP,"NETWORK_NAMES"},'R',1)
            Col_5 = [Col_5; {'Ipsilesional'}];
        else
            Col_5 = [Col_5; {'Contralesional'}];
        end
    end
end

% Table
R_table = table(Col_1, Col_2, Col_3, Col_4,Col_5, 'VariableNames', {'Mean_dALFF', 'NETWORK_Names','SUB_Names','Lateralization','Location'});
df_name = fullfile(dataset_root,'CSV_file/dALFF_frontal_parietallobule.csv');
writetable(R_table, df_name);  



%% Make ALFF table for R
% LEFT frontal patients
Col_1 = [];
Col_2 = [];
Col_3 = [];
Col_4 = [];
Col_5 = [];
for SUB_LOOP = 1:size(df_LEFT,1)
    Mean_ALFF_LEFT = cellfun(@mean,df_LEFT.ALFF{SUB_LOOP,1}{:,"Parcels"});
    Col_1 = [Col_1; Mean_ALFF_LEFT];
    Col_2 = [Col_2; df_LEFT.ALFF{SUB_LOOP,1}{:,"NETWORK_NAMES"}];
    Col_3 = [Col_3; repmat({df_LEFT.ID{SUB_LOOP,1}},12,1)];
    Col_4 = [Col_4; repmat({df_LEFT.Tumor_lateralization{SUB_LOOP,1}},12,1)];
    for NETWORK_LOOP = 1:length(df_LEFT.ALFF{SUB_LOOP,1}{:,"NETWORK_NAMES"}) %1:34
        if strncmp(df_LEFT.ALFF{SUB_LOOP,1}{NETWORK_LOOP,"NETWORK_NAMES"},'L',1)
            Col_5 = [Col_5; {'Ipsilesional'}];
        else
            Col_5 = [Col_5; {'Contralesional'}];
        end
    end
end

% RIGHT frontal patients
for SUB_LOOP = 1:size(df_RIGHT,1)
    Mean_ALFF_RIGHT = cellfun(@mean,df_RIGHT.ALFF{SUB_LOOP,1}{:,"Parcels"});
    Col_1 = [Col_1; Mean_ALFF_RIGHT];
    Col_2 = [Col_2; df_RIGHT.ALFF{SUB_LOOP,1}{:,"NETWORK_NAMES"}];
    Col_3 = [Col_3; repmat({df_RIGHT.ID{SUB_LOOP,1}},12,1)];
    Col_4 = [Col_4; repmat({df_RIGHT.Tumor_lateralization{SUB_LOOP,1}},12,1)];
    for NETWORK_LOOP = 1:length(df_RIGHT.ALFF{SUB_LOOP,1}{:,"NETWORK_NAMES"}) %1:34
        if strncmp(df_RIGHT.ALFF{SUB_LOOP,1}{NETWORK_LOOP,"NETWORK_NAMES"},'R',1)
            Col_5 = [Col_5; {'Ipsilesional'}];
        else
            Col_5 = [Col_5; {'Contralesional'}];
        end
    end
end

% Table
R_table = table(Col_1, Col_2, Col_3, Col_4,Col_5, 'VariableNames', {'Mean_ALFF', 'NETWORK_Names','SUB_Names','Lateralization','Location'});
df_name = fullfile(dataset_root,'CSV_file/ALFF_frontal_parietallobule.csv');
writetable(R_table, df_name); 

%% T-statistics

Networkname_lists = unique(R_table.NETWORK_Names);
var1 = zeros(0, 1);
var2 = cell(0, 1); 
tstat_garage = table(var1,var2);
ranksum_garage = table(var1,var2);
for NETWORK_LOOP = 1:17

    Ipsi_idx = ismember(R_table.Location,'Ipsilesional');
    Contra_idx = ismember(R_table.Location,'Contralesional');
    str_Networks = [Networkname_lists(NETWORK_LOOP),Networkname_lists(NETWORK_LOOP+17)];
    Network_idx = ismember(R_table.NETWORK_Names,str_Networks);
    
    ispi_combined_idx = Ipsi_idx & Network_idx;
    contra_combined_idx = Contra_idx & Network_idx;
    
    group_ipsi = R_table.Mean_ALFF(ispi_combined_idx);
    group_contra = R_table.Mean_ALFF(contra_combined_idx);
    
    % Perform paired t-test
    [h, p, ci, stats] = ttest(group_contra,group_ipsi); % paired-ttest, default is two-tailed and alpha = .05
    tstat_garage(NETWORK_LOOP,1) = {h};
    tstat_garage{NETWORK_LOOP,2} = Networkname_lists(NETWORK_LOOP);
    tstat_garage(NETWORK_LOOP,3) = {p};
    tstat_garage(NETWORK_LOOP,4) = {ci};
    tstat_garage{NETWORK_LOOP,5} = stats;

    % Perform ranksum test (Wilcoxon test)
    [p,h,stats] = ranksum(group_contra,group_ipsi);
    ranksum_garage(NETWORK_LOOP,1) = {h};
    ranksum_garage{NETWORK_LOOP,2} = Networkname_lists(NETWORK_LOOP);
    ranksum_garage(NETWORK_LOOP,3) = {p};
    ranksum_garage{NETWORK_LOOP,4} = stats;
end







%% ~~~~~~~~ Each time course of dALFF

% ~~~~ get the HGG index, read info.excel file
excel_info_path = fullfile(info_path,"lnif_prep_status.xlsx");
excel_info_data = readtable(excel_info_path);
excel_info_data(58:end,:) = [];
str_exclusion = {'ToF_02_54','DeTe_93'};
exclusion_index = ~ismember(excel_info_data.ID, str_exclusion);
data_list = data_list(exclusion_index,:);
excel_info_data = excel_info_data(exclusion_index,:);

% Functional networks template
parcellation_path = fullfile(atlas_path, 'Schaefer2018_200Parcels_Kong2022_17Networks_order_FSLMNI152_2mm.nii.gz');
parcellation_data = niftiread(parcellation_path);
% % Functional networks csv file
FNcsv_path = fullfile(atlas_path,'Schaefer2018_200Parcels_Kong2022_parietallobule.csv');
FNcsv_data_table = readtable(FNcsv_path);

cogsco_path = fullfile(csv_path,'CognitiveScore_PL.csv');
df_cogsco = readtable(cogsco_path);
% list of interested patient
listpatient = unique(df_cogsco.Patient);
% remove FaLa_45
idx = ismember(listpatient,'FaLa_45');
listpatient = listpatient(~idx);
% extract 19 patients from 55 
idx = ismember(excel_info_data.ID,listpatient);
df_timecourse = excel_info_data(idx,:);


%% read normalized_4D_dALFF.nii

for SUB_LOOP = 1:size(df_timecourse,1)

    SUBJ = char(df_timecourse.ID(SUB_LOOP));
    % SUBJ = data_list(SUB_LOOP).name;
    subj_folder = fullfile(data_path,SUBJ);
    anat_folder = dir(fullfile(subj_folder,'*.anat'));
    anat_path = fullfile(subj_folder,anat_folder.name);
    func_path = fullfile(subj_folder,'func','sdc_unc');
    result_folder = fullfile(result_path,SUBJ); 
    
    % ~~~~~~~ import dALFF map, FN template,network csv file    
    % import dALFF (4D map)
    step_path = fullfile(result_folder,'dALFF','Stepsize_2','normalized_4D_dALFF.nii.gz');
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
        idx = FNcsv_data_table.ParietalID>0;
        garage_parcel = cell2table(garage_parcel(idx,:));
        garage_parcel.Var4 = FNcsv_data_table.ParietalID(idx);

        % Sum the parcels
        % each network
        N_all = [];
        for i = 1:12
        idx = garage_parcel.Var4==i;
        x = garage_parcel.Var1(idx);
        x_all=[];
            % ~~~ merge parcels
            for ii = 1:sum(idx)      
                s = cell2mat(x(ii));
                x_all = [x_all; s];
            end
            %~~~~ merge parcels

        N_all = [N_all; mean(x_all)];
        end % end of network loop

        % clean the data from 26 to 12 network
        % Find unique rows
        [b,m] = unique(garage_parcel.Var4, 'rows', 'stable');
        Ex_garage_parcel = garage_parcel(m,3:4);
        Ex_garage_parcel.Var5 = N_all;

        ALLWINDOW_parcel{WINDOW_LOOP,1} = Ex_garage_parcel;
    end

    df_timecourse.dALFF(SUB_LOOP) = {ALLWINDOW_parcel};
    disp(sprintf('%s is done', SUBJ))
    
end %SUB_LOOP end
%
filename = fullfile(code_path,'SUB_dALFF_200parcels_ParietalLobule');
save(filename, 'df_timecourse');

%% ~~~~~ Visualize timecourse

%import data
load('SUB_dALFF_200parcels_ParietalLobule.mat')
% % Functional networks csv file
FNcsv_path = fullfile(atlas_path,'Schaefer2018_200Parcels_Kong2022_parietallobule.csv');
FNcsv_data_table = readtable(FNcsv_path);


%% Clean the data
df = {};
x_all = [];

% SUB_LOOP = 19
for SUB_LOOP = 1:size(df_timecourse,1)

    
    % WINDOW LOOP = 124
    for WINDOW_LOOP = 1:size(df_timecourse.dALFF{1,1},1)

        % extract a data into x variable(Table)
        x = df_timecourse.dALFF{SUB_LOOP,1}{WINDOW_LOOP,1};
        allVars = 1:width(x);
        x = renamevars(x, allVars, ["Network","ID","Mean"]);
        
        % remove unnesarry strings
        net_vec = ["Control A","Control B","Salience_VenAttnA","Salience_VenAttnB","DorsalAttnA","DorsalAttnB"...
            , "Control A","Control B","Salience_VenAttnA","Salience_VenAttnB","DorsalAttnA","DorsalAttnB"];
        x.Network = net_vec';

        % create a location vector
        cont_loc = repmat("contralesional",6,1);
        ipsi_loc = repmat("ipsilesional",6,1);
        glio_loc = df_timecourse.Tumor_lateralization;
        if glio_loc == "RIGHT"
            x.Location = [cont_loc; ipsi_loc];
        else
            x.Location = [ipsi_loc; cont_loc];
        end
        
        % create window vector
        rep_vec = repmat(WINDOW_LOOP,height(x),1);
        x.Window = rep_vec;
        % create SUB vector
        rep_vec = repmat(string(df_timecourse.ID(SUB_LOOP)),height(x),1);
        x.Patient = rep_vec;

        % create a dataframe
        x_all = [x_all; x];

    end
end

% download csv.file
filename = fullfile(csv_path,"6_Network_parietal.csv");
writetable(x_all, filename);


%% Clean the data for not separating ipsi and contra
df = {};
x_all = [];
xx = table();

% SUB_LOOP = 19
for SUB_LOOP = 1:size(df_timecourse,1)

    
    % WINDOW LOOP = 124
    for WINDOW_LOOP = 1:size(df_timecourse.dALFF{1,1},1)

        % extract a data into x variable(Table)
        x = df_timecourse.dALFF{SUB_LOOP,1}{WINDOW_LOOP,1};
        allVars = 1:width(x);
        x = renamevars(x, allVars, ["Network","ID","Mean"]);
        mean_vec = [];
        for i = 1:6
            m = mean(x.Mean(x.ID==i & i+6));
            mean_vec = [mean_vec;m];
        end

        xx.ID = [1:6]';
        xx.Mean = mean_vec;
        % remove unnesarry strings
        net_vec = ["Control A","Control B","Salience_VenAttnA","Salience_VenAttnB","DorsalAttnA","DorsalAttnB"];
        xx.Network = net_vec';
        % create window vector
        rep_vec = repmat(WINDOW_LOOP,height(xx),1);
        xx.Window = rep_vec;
        % create SUB vector
        rep_vec = repmat(string(df_timecourse.ID(SUB_LOOP)),height(xx),1);
        xx.Patient = rep_vec;

        % create a dataframe
        x_all = [x_all; xx];

    end
end

% download csv.file
filename = fullfile(csv_path,"6_Network_parietal_WhoBra.csv");
writetable(x_all, filename);



%% Plot
% SUB_LOOP
SUB_LOOP = 1;

% WINDOW_LOOP
WINDOW_LOOP = 1;

% NETWORK_LOOP
NETWORK_LOOP = 1;

df{SUB_LOOP,1}{WINDOW_LOOP,1}.Mean(NETWORK_LOOP)


NETWORK_idx = strcmp(res_table.Net_names,NETWORK_df(ii));
res_timeseries = timeseries(res_table.res(NETWORK_idx));
res_timeseries.Name = 'fdALFF';
res_timeseries.TimeInfo.Units = 'windows';

% time series plot
figure(1)    
plot(res_timeseries)
% plot(contra,'-b','color','r')
% grid on
title(sprintf('%s',SUB(1)))
xlabel('Time window')
ylabel('fdALFF (Ipsi - Contra)')
xlim([-1,124])
ylim([-0.3,0.3])
fontsize(gca, 18,'points')
line(xlim(), [0,0], 'LineWidth', 1, 'Color', 'k');

hold on

txt_ipsi = 'Ipsilesional dominant';
txt_contra = 'Contralesional dominant';
text(0,0.28,txt_ipsi,"FontSize",13)
text(0,-0.28,txt_contra,"FontSize",13)
txt_SalienceA = sprintf('―%s',NETWORK_df(ii-1));
txt_SalienceB = sprintf('―%s',NETWORK_df(ii));
text(85,0.28,txt_SalienceA,"FontSize",13,"Color",'b')
text(85,0.25,txt_SalienceB,"FontSize",13,"Color",'r')
hold off
%     filename = sprintf('%s_timecourse.png',SUB(1));
%     fig_file = fullfile(dest_dir, filename);
%     saveas(gcf,fig_file)

% Swarmchart plot
x = categorical(res_table.Net_names);
y = res_table.res;
swarmchart(x,y,20,'filled')
set(gca, 'FontSize', 15);
ylabel('ipsi - contra for each time window','FontSize',18);
ylim([-0.3,0.3])
line(xlim(), [0,0], 'LineWidth', 1, 'Color', 'k');
    
%     dest_dir = sprintf("/Users/muku/Desktop/fdALFF");
%     mkdir(dest_dir);
%     filename = sprintf('%s.png',SUB(1));
%     fig_file = fullfile(dest_dir, filename);
%     saveas(gcf,fig_file)

