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

%% Download the data

load('ALLSUB_fdALFF&ALFF.mat') 

%% Extract the interested patients from 55 total

% clearvars -except ALLSUB_lists;
% Extract left frontal and fronto-insulare
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
    Col_3 = [Col_3; repmat({df_LEFT.ID{SUB_LOOP,1}},34,1)];
    Col_4 = [Col_4; repmat({df_LEFT.Tumor_lateralization{SUB_LOOP,1}},34,1)];
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
    Col_3 = [Col_3; repmat({df_RIGHT.ID{SUB_LOOP,1}},34,1)];
    Col_4 = [Col_4; repmat({df_RIGHT.Tumor_lateralization{SUB_LOOP,1}},34,1)];
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
% writetable(R_table, 'fdALFF_frontal.csv');  



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
    Col_3 = [Col_3; repmat({df_LEFT.ID{SUB_LOOP,1}},34,1)];
    Col_4 = [Col_4; repmat({df_LEFT.Tumor_lateralization{SUB_LOOP,1}},34,1)];
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
    Col_3 = [Col_3; repmat({df_RIGHT.ID{SUB_LOOP,1}},34,1)];
    Col_4 = [Col_4; repmat({df_RIGHT.Tumor_lateralization{SUB_LOOP,1}},34,1)];
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
% writetable(R_table, 'fALFF_frontal.csv');  

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
