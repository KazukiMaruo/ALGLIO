#!/bin/bash
PATH=/mnt/storage/tier2/MUMI-EXT-001/mumi-data/Project-Glioma-FrequencyFluctuations/Codes/
FSLDIR=/opt/fsl
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH
. ${FSLDIR}/etc/fslconf/fsl.sh

names="AgFa_58 AmL_03_93 AmM_04_54 ArAn_59 ArMa_64 AvVi_46 BoFl_52 BoR_07_66 BuE_07_62 BuG_03_45 CaEl_66 CaL_01_47 CeA_10_92 CoN_02_97 DaM_11_81 DeMa_56 DeTe_58 DrGa_90 FaLa_45 FeM_05_79 FiGa_53 FoR_05_46 GaJ_06_78 GuMa_78 IoPa_61 LoGi_45 MaAl_87 MaL_07_68 MaR_02_64 MaR_09_83 MaVi_40 MeA_12_67 MeEd_38 MeIt_63 MiF_11_75 MoI_12_64 MuRa_67 NiGi_62 PaRo_62 PiA_06_76 PiAn_52 PiS_07_70 PrVi_50 RiCa_52 RuGi_49 SaGi_67 ScMa_01 SeMa_66 SpPi_66 SuI_04_68 TeAl_83 ToCh_01 YeH_08_85 ZeAl_96 ZoE_03_65"
for SUB_names in $names; do
    #input_folder_step2="/Volumes/mumi-data/Project-Glioma-FrequencyFluctuations/Results_ALLSUB/$SUB_names/dALFF/Stepsize_2"
    input_folder_step2="/mnt/storage/tier2/MUMI-EXT-001/mumi-data/Project-Glioma-FrequencyFluctuations/Results_ALLSUB/$SUB_names/dALFF/Stepsize_2"
    cd "$input_folder_step2"
# compute mean, std, meadian for stepsize_2
    fslmaths $input_folder_step2/normalized_4D_fdALFF.nii.gz -Tmean $input_folder_step2/nmz_avg_fdALFF.nii.gz
    fslmaths $input_folder_step2/normalized_4D_fdALFF.nii.gz -Tstd $input_folder_step2/nmz_std_fdALFF.nii.gz
    fslmaths $input_folder_step2/normalized_4D_fdALFF.nii.gz -Tmedian $input_folder_step2/nmz_median_fdALFF.nii.gz

# compute mean, std, meadian for stepsize_3
    #input_folder_step2="/Volumes/mumi-data/Project-Glioma-FrequencyFluctuations/Results_ALLSUB/$SUB_names/dALFF/Stepsize_2"
    input_folder_step3="/mnt/storage/tier2/MUMI-EXT-001/mumi-data/Project-Glioma-FrequencyFluctuations/Results_ALLSUB/$SUB_names/dALFF/Stepsize_3"
    cd "$input_folder_step3"
# compute mean, std, meadian for stepsize_2
    fslmaths $input_folder_step3/normalized_4D_fdALFF.nii.gz -Tmean $input_folder_step3/nmz_avg_fdALFF.nii.gz
    fslmaths $input_folder_step3/normalized_4D_fdALFF.nii.gz -Tstd $input_folder_step3/nmz_std_fdALFF.nii.gz
    fslmaths $input_folder_step3/normalized_4D_fdALFF.nii.gz -Tmedian $input_folder_step3/nmz_median_fdALFF.nii.gz
    
done