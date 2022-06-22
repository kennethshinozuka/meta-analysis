#!/bin/bash

cd /Users/kshinozuka/Documents/Oxford/Research/Rebirth/Literature/Meta-Analysis/Code\ Review/Data/FC/Parcellation\ Consistency/

declare -A parcels_lookup
num_rows=2      
num_columns=7  

parcels_lookup[1,1]="brainnetome"
parcels_lookup[1,3]="craddock"
parcels_lookup[1,4]="harvard_oxford"
parcels_lookup[1,5]="power"
parcels_lookup[1,6]="raichle"
parcels_lookup[1,7]="schaefer"
parcels_lookup[1,8]="shen"
parcels_lookup[2,1]=246
parcels_lookup[2,2]=200
parcels_lookup[2,3]=132
parcels_lookup[2,4]=27
parcels_lookup[2,5]=36
parcels_lookup[2,6]=200
parcels_lookup[2,7]=268

for (( i = 0; i <= ${num_columns}; i++ )); do
    for (( j = 0; j <= ${parcels_lookup[2,${i}]}; j++ )); do
        l=`echo "${j} - 0.5" | bc`
        u=`echo "${j} + 0.5" | bc`
        parcellation=${parcels_lookup[1,${i}]}
        parcellation_filename="${parcellation}_atlas.nii"
        fslmaths ${parcellation_filename} -thr ${l} -uthr ${u} -bin "${parcellation}_mask_$( printf '%d' $j )";
    done
done