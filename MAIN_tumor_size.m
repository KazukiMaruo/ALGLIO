clc
close all
clear

% Dataset folders
dataset_root = fullfile(pwd,'..');
data_path = fullfile(dataset_root,'DATA_PROCESSED');
data_list = dir(fullfile(dataset_root,'DATA_PROCESSED','*set*'));
results_folder = fullfile(dataset_root,'RESULTS');
figures_folder = fullfile(dataset_root,'FIGURES');

% Load packages
addpath(genpath(fullfile(dataset_root,'..','USERS','Manuela_Moretto','TOOLS','MATLAB_UTILS','NIfTI_20140122')))
addpath(genpath(fullfile(pwd,'..','..','USERS','Manuela_Moretto','TOOLS','MATLAB_UTILS','github_repo')))

%% COMPUTE LESION SIZE IN T1 SPACE
% Analyse single subject tumor masks
tm_mask_name    = 'lesionmask.nii.gz';
brain_mask_name = 'T1_biascorr_optiBET_brain_mask.nii.gz';
    
pos = 0;
for sets = 1:length(data_list)
    batch = data_list(sets).name;
    subj_list = dir(fullfile(data_path,batch,'sub-*'));
    
    for subjs = 1:length(subj_list)
        pos  = pos + 1;
        subj = subj_list(subjs).name;
 
        subject_func_path = fullfile(data_path,batch,subj,'func','sdc_unc');
        tmp_path = dir(fullfile(subject_func_path,'*init.mat'));
        SUBJ = tmp_path.name;
        limi = strfind(string(SUBJ),'M00');
        SUBJ = SUBJ(1:limi-2);
        disp(['Working on ' SUBJ])
        results{pos,1} = SUBJ;
        
        % Load tumor mask
        subject_anat_path = dir(fullfile(data_path,batch,subj,'sub-*ses-01_T1w.anat'));
        tm_mask_path = fullfile(data_path,batch,subj,subject_anat_path(1).name,tm_mask_name);
        tm_mask_file = load_untouch_nii(tm_mask_path);
        tm_mask = double(tm_mask_file.img);
        
        % Load brain mask
        brain_mask_path = fullfile(data_path,batch,subj,subject_anat_path(1).name,brain_mask_name);
        
        % Fill holes in brain mask
        system(['fslmaths ' brain_mask_path ' -fillh ' fullfile(data_path,batch,subj,subject_anat_path(1).name,'T1_biascorr_optiBET_brain_mask_filled.nii.gz')]);
        fill_brain_mask_path = fullfile(data_path,batch,subj,subject_anat_path(1).name,'T1_biascorr_optiBET_brain_mask_filled.nii.gz');
        brain_mask_file = load_untouch_nii(fill_brain_mask_path);
        brain_mask = double(brain_mask_file.img);

        % N voxels tumor mask
        results{pos,2} = nnz(tm_mask);
        
        % Volume tumor mask
        vox_vol = tm_mask_file.hdr.dime.pixdim(2)*tm_mask_file.hdr.dime.pixdim(3)*tm_mask_file.hdr.dime.pixdim(4);
        results{pos,3} = nnz(tm_mask)*vox_vol;
        
        % Perc Normalized volume tumor mask
        results{pos,4} = 100*nnz(tm_mask)/nnz(brain_mask);
        
    end
end
% save(fullfile(results_folder,'tumor_extension.mat'),'results');
save(fullfile(results_folder,'tumor_extension_20230926.mat'),'results');

%% DIVIDE BETWEEN HGG AND LGG
% load(fullfile(results_folder,'tumor_extension.mat'),'results');
load(fullfile(results_folder,'tumor_extension_20230926.mat'),'results');
load(fullfile(results_folder,'Patients_info.mat'),'all_info','pt_surgery','tumor_hemi','tumor_loc','tumor_grade','tumor_stad','tumor_idh','pt_names','idx_tm_mask')

pos_hgg = 1;
pos_lgg = 1;
for ss=1:size(results,1)
    SUBJ = results{ss,1};
    if SUBJ=="BuE_07_62"
        SUBJ="BuE_07_52";
    end
    idx_pt = find(strcmp(SUBJ,[pt_names]));
    
    if tumor_grade{idx_pt}=="HGG"
        idx_pt_hgg(pos_hgg,1) = ss;
        pos_hgg = pos_hgg+1;
    else
        idx_pt_lgg(pos_lgg,1) = ss;
        pos_lgg = pos_lgg+1;
    end
end

% Divide results in HGG and LGG subgroups
results_hgg = results(idx_pt_hgg,:);
results_lgg = results(idx_pt_lgg,:);
avg_FD_hgg  = mean(cell2mat(results_hgg(:,4)));
avg_FD_lgg  = mean(cell2mat(results_lgg(:,4)));
disp(['Average % normalized tumor volume - LGG=' num2str(avg_FD_lgg)])
disp(['Average % normalized tumor volume - HGG=' num2str(avg_FD_hgg)])

% Statistics
[P,H] = ranksum(cell2mat(results_lgg(:,4)),cell2mat(results_hgg(:,4)));

%save(fullfile(results_folder,'tumor_extension.mat'),'results_hgg','results_lgg','-append');
save(fullfile(results_folder,'tumor_extension_20230926.mat'),'results_hgg','results_lgg','-append');

%% Figures - stem
load(fullfile(results_folder,'tumor_extension_20230926.mat'));

figure
stem(cell2mat(results(:,4)),'filled','Color',[0 0 1])
xlim([0,length(results(:,4))])
xticks(1:length(results(:,4)))
xticklabels(strrep(results(:,1),'_',' '))
set(gca,'FontSize',8,'Fontweight','bold')
xlabel('Patient ID','FontSize',14,'Fontweight','bold')
title('% tumor volume with respect of brain volume','FontSize',14,'Fontweight','bold')
view(90,90)
set(gcf,'Position',get(0,'Screensize'));
set(gcf,'Color',[1 1 1]);
%export_fig(gcf,fullfile(figures_folder,['Tumor_extension.png']))
export_fig(gcf,fullfile(figures_folder,['Tumor_extension_20230926.png']))
close

% Divide in LGG and HGG
figure
subplot(121)
stem(cell2mat(results_lgg(:,4)),'filled','Color',[0 0 1])
xlim([0,length(results_lgg(:,4))])
xticks(1:length(results_lgg(:,4)))
xticklabels(strrep(results_lgg(:,1),'_',' '))
set(gca,'FontSize',8,'Fontweight','bold')
xlabel('Patient ID','FontSize',14,'Fontweight','bold')
title(['% tumor volume with respect of brain volume-' num2str(length(idx_pt_lgg)) 'LGG'],'FontSize',14,'Fontweight','bold')
view(90,90)
ylim([0 10])
subplot(122)
stem(cell2mat(results_hgg(:,4)),'filled','Color',[0 0 1])
xlim([0,length(results_hgg(:,4))])
xticks(1:length(results_hgg(:,4)))
xticklabels(strrep(results_hgg(:,1),'_',' '))
set(gca,'FontSize',8,'Fontweight','bold')
xlabel('Patient ID','FontSize',14,'Fontweight','bold')
title(['% tumor volume with respect of brain volume-' num2str(length(idx_pt_hgg)) 'HGG'],'FontSize',14,'Fontweight','bold')
view(90,90)
ylim([0 10])
set(gcf,'Position',get(0,'Screensize'));
set(gcf,'Color',[1 1 1]);
%export_fig(gcf,fullfile(figures_folder,['Tumor_extension_HGG_LGG.png']))
export_fig(gcf,fullfile(figures_folder,['Tumor_extension_HGG_LGG_20230926.png']))
close