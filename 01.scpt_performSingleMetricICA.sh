
# script to perform ICA on single metric maps

cd ../singleMetricICA/

for m in alff area  dc  ec  falff  lgi  meancurv  reho  reho2  sulc  thickness  volume
do
	melodic -i forICA_${m}.nii.gz --outdir=${m}.ica --mask=allOne.nii.gz --nobet --no_mm --Oorig -v
done

