#INDIR=/g100_work/OGS_prodC/MIT/V1M-prod/archive/current/products
INDIR=/g100_work/OGS_prodC/MIT/V1M-prod/archive/
FILE_CHLA=20220317_h-OGS--PFTC-MITgcmBFM-pilot8-b20220317_fc-v01.nc
MASKFILE=/g100_work/OGS_prodC/MIT/V1/prod/wrkdir/POSTPROC/meshmask.nc
INPUT_CHLA_DIR=$INDIR
# SAT validation:
mkdir -p Chla_SAT/Tserie
SAT_DAILY_DIR=/g100_work/OGS_prodC/MIT/V1/prod/inpdir/SAT/


mkdir -p ./Chla_SAT/Tserie/ ./Chla_SAT/BIAS_RMSD/
mkdir -p ./Chla_SAT/Tserie/everywhere ./Chla_SAT/BIAS_RMSD/everywhere
mkdir -p ./Chla_SAT/Tserie/open_sea ./Chla_SAT/BIAS_RMSD/open_sea
mkdir -p ./Chla_SAT/Tserie/coast ./Chla_SAT/BIAS_RMSD/coast


python ScMYvalidation_plan_STD_CORR_valid_NAdr.py -s $SAT_DAILY_DIR -i $INPUT_CHLA_DIR -m $MASKFILE -c everywhere -o export_data_ScMYValidation_plan_everywhere.pkl -t_start "20230101" -t_end "20240110"

python ScMYvalidation_plan_STD_CORR_valid_NAdr.py -s $SAT_DAILY_DIR -i $INPUT_CHLA_DIR -m $MASKFILE -c open_sea -o export_data_ScMYValidation_plan_open_sea.pkl -t_start "20230101" -t_end "20240110"

python ScMYvalidation_plan_STD_CORR_valid_NAdr.py -s $SAT_DAILY_DIR -i $INPUT_CHLA_DIR -m $MASKFILE -c coast -o export_data_ScMYValidation_plan_coast.pkl -t_start "20230101" -t_end "20240110"


python plot_timeseries_STD_NAdr.py -i export_data_ScMYValidation_plan_everywhere.pkl  -o ./Chla_SAT/Tserie/everywhere -v "chl"

python plot_timeseries_STD_NAdr.py -i export_data_ScMYValidation_plan_open_sea.pkl  -o ./Chla_SAT/Tserie/open_sea -v "chl"

python plot_timeseries_STD_NAdr.py -i export_data_ScMYValidation_plan_coast.pkl  -o ./Chla_SAT/Tserie/coast -v "chl"


python plot_timeseries_RMS_CORR_NAdr.py -i export_data_ScMYValidation_plan_everywhere.pkl -o ./Chla_SAT/BIAS_RMSD/everywhere #-v "chl"

python plot_timeseries_RMS_CORR_NAdr.py -i export_data_ScMYValidation_plan_coast.pkl -o ./Chla_SAT/BIAS_RMSD/coast #-v "chl"

python plot_timeseries_RMS_CORR_NAdr.py -i export_data_ScMYValidation_plan_open_sea.pkl -o ./Chla_SAT/BIAS_RMSD/open_sea #-v "chl"

#echo date +%F-%T > timestamp.txt

date --utc --iso-8601=seconds > ./Chla_SAT/timestamp.txt
