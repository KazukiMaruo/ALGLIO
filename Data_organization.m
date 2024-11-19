% Statistical Analysis
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

% extract all the files in the data folder
data_list   = dir(data_path);
data_list   = data_list(4:end); % first three rows are unnecessary



%% ~~~~~ dALFF
% Download the data
load('ALLSUB_fdALFF_200parcels.mat')

%% Merge the parcels into Functional Networks per hemisphere

% ~~~~~~ Assign the index for functional networks for each hemisphere, from 1:34
networkID = double(FNcsv_data_table.NETWORKID);
networkID(101:end) = networkID(101:end)+17;
% Update ALLSUB_parcel
for ss = 1:length(ALLSUB_parcel) % 1:the number of patients, eg., 1:55
    for STEP_LOOP = 1:2 % step2 and step3
        for ii = 1:200
            ALLSUB_parcel{ss,STEP_LOOP}{ii,4} = networkID(ii);
        end
    end
end


% ~~~~~ Create new data structure from 55(SUB)*1 cell inside 34(Each network)*3(Parcels, NETWORK_Names, Network_ID) 
% Only stepsize2
% Transfrom from cell to table
NetworkNames = FNcsv_data_table.NETWORKNAME;
ALLSUB_parcels_table = {};
for SUB_LOOP = 1:length(ALLSUB_parcel)
    Sub_data = ALLSUB_parcel{SUB_LOOP,1};
    Sub_table = array2table(Sub_data, 'VariableNames', {'Parcels', 'Mean','NETWORK_Names','Network_ID' });
    Sub_table.NETWORK_Names = NetworkNames;
    ALLSUB_parcels_table{SUB_LOOP,1} = Sub_table;
end

% ~~~~~ Merge percels for each hemisphere Network
% ~~~~~ ALL SUBJECTS
EachHemNetworks = [];
for SUB_LOOP = 1:length(ALLSUB_parcels_table)

    df = ALLSUB_parcels_table{SUB_LOOP, 1};

    x_vectors = [];
        % ~~~~~~ A SUBJECT 
        for NETWORK_LOOP = 1:34
            % OneNetwork_parcels_index = strcmp(string(df.Network_ID),'1');
            Network_parcels_index = cell2mat(df.Network_ID) == NETWORK_LOOP;
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



% ~~~~~~ New data table similar to ALLSUB_parcels_table

% Change the name
ALLSUB_dALFF_EachHemNetworks = EachHemNetworks;

% get the NETWORK_NAMES
uniqueNetworkArrays = uniqueOrderPreserved(NetworkNames)';
characterToAdd_LH = 'LH_'; % Specify the character to add
characterToAdd_RH = 'RH_'; % Specify the character to add
uniqueNetworkArrays_LH = cellfun(@(x) [characterToAdd_LH, x], uniqueNetworkArrays, 'UniformOutput', false);
uniqueNetworkArrays_RH = cellfun(@(x) [characterToAdd_RH, x], uniqueNetworkArrays, 'UniformOutput', false);
uniqueNetworkArrays = [uniqueNetworkArrays_LH;uniqueNetworkArrays_RH];

%get the NETWORK_ID
uniqueIDArrays = uniqueOrderPreserved(ALLSUB_parcels_table{1,1}.Network_ID)';

for SUB_LOOP = 1:length(ALLSUB_dALFF_EachHemNetworks)
    % make the table for each SUB
    x = array2table(EachHemNetworks{SUB_LOOP,1},"VariableNames",{'Parcels'});
    x.NETWORK_NAMES = uniqueNetworkArrays;
    x.NETWORK_ID = uniqueIDArrays;
    % Modify the data
    ALLSUB_dALFF_EachHemNetworks{SUB_LOOP,1} = x;
end



%% ~~~~~~~~~~ ALFF

% clearvars -except ALLSUB_dALFF_EachHemNetworks ALLSUB_ALFF_EachHemNetworks;
% Download the data
load('ALLSUB_fALFF_200parcels.mat')
% Merge the parcels into Functional Networks per hemisphere

% ~~~~~~ Assign the index for functional networks for each hemisphere, from 1:34
networkID = double(FNcsv_data_table.NETWORKID);
networkID(101:end) = networkID(101:end)+17;
% Update ALLSUB_parcel
for ss = 1:length(ALLSUB_parcel) % 1:the number of patients, eg., 1:55
    for ii = 1:200
        ALLSUB_parcel{ss,1}{ii,4} = networkID(ii);
    end
end


% ~~~~~ Create new data structure from 55(SUB)*1 cell inside 34(Each network)*3(Parcels, NETWORK_Names, Network_ID) 
% Only stepsize2
% Transfrom from cell to table
NetworkNames = FNcsv_data_table.NETWORKNAME;
ALLSUB_parcels_table = {};
for SUB_LOOP = 1:length(ALLSUB_parcel)
    Sub_data = ALLSUB_parcel{SUB_LOOP,1};
    Sub_table = array2table(Sub_data, 'VariableNames', {'Parcels', 'Mean','NETWORK_Names','Network_ID' });
    Sub_table.NETWORK_Names = NetworkNames;
    ALLSUB_parcels_table{SUB_LOOP,1} = Sub_table;
end

%% ~~~~~ Merge percels for each hemisphere Network
% ~~~~~ ALL SUBJECTS
EachHemNetworks = [];
for SUB_LOOP = 1:length(ALLSUB_parcels_table)

    df = ALLSUB_parcels_table{SUB_LOOP, 1};
    x_vectors = [];
        % ~~~~~~ A SUBJECT 
        for NETWORK_LOOP = 1:34
            % OneNetwork_parcels_index = strcmp(string(df.Network_ID),'1');
            Network_parcels_index = cell2mat(df.Network_ID) == NETWORK_LOOP;
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
ALLSUB_ALFF_EachHemNetworks = EachHemNetworks;

% get the NETWORK_NAMES
uniqueNetworkArrays = uniqueOrderPreserved(NetworkNames)';
characterToAdd_LH = 'LH_'; % Specify the character to add
characterToAdd_RH = 'RH_'; % Specify the character to add
uniqueNetworkArrays_LH = cellfun(@(x) [characterToAdd_LH, x], uniqueNetworkArrays, 'UniformOutput', false);
uniqueNetworkArrays_RH = cellfun(@(x) [characterToAdd_RH, x], uniqueNetworkArrays, 'UniformOutput', false);
uniqueNetworkArrays = [uniqueNetworkArrays_LH;uniqueNetworkArrays_RH];

%get the NETWORK_ID
uniqueIDArrays = uniqueOrderPreserved(ALLSUB_parcels_table{1,1}.Network_ID)';

for SUB_LOOP = 1:length(ALLSUB_ALFF_EachHemNetworks)
    % make the table for each SUB
    x = array2table(EachHemNetworks{SUB_LOOP,1},"VariableNames",{'Parcels'});
    x.NETWORK_NAMES = uniqueNetworkArrays;
    x.NETWORK_ID = uniqueIDArrays;
    % Modify the data
    ALLSUB_ALFF_EachHemNetworks{SUB_LOOP,1} = x;
end




%% Merge with the excel file
excel_info_path = fullfile(info_path,"lnif_prep_status.xlsx");
excel_info_data = readtable(excel_info_path);
excel_info_data(58:end,:) = [];
str_exclusion = {'ToF_02_54','DeTe_93'};
exclusion_index = ~ismember(excel_info_data.ID, str_exclusion);
ALLSUB_lists = excel_info_data(exclusion_index,:);

% Add dALFF score
ALLSUB_lists.dALFF = ALLSUB_dALFF_EachHemNetworks;
ALLSUB_lists.ALFF = ALLSUB_ALFF_EachHemNetworks;

save('ALLSUB_fdALFF&ALFF.mat', 'ALLSUB_lists')
