#!/bin/bash -x
function Initial_Conditions_ERA5() {
#-----------------------------------------------------------------------------#
#                                   DIMNT/INPE                                #
#-----------------------------------------------------------------------------#
#BOP
#
# !SCRIPT: Initial_Conditions_ERA5
#
# !DESCRIPTION:
#      
#     Program to donwload the initial boundary 
#     condition from  ERA5  to MPAS model.
#
#        ./Initial_Conditions_ERA5.bash   
#           o IC_HOME   : Path to save the IC  
#           o EXP_IBC   : Forcing : ERA5
#           o LABELI    : Initial : date 2015030600
#           o LABELF    : Final date date 2015030600
#           o LABELF    : Hours Step to download Boundary Conditios
#           o LAT_INI   : Initial latitude
#           o LAT_INI   : Final latitude
#           o LON_INI   : Initial longitude
#           o LON_FIN   : Final longitude
#
# !REVISION HISTORY:
#
# !REMARKS:
#
#EOP
#-----------------------------------------------------------------------------!
#EOC

function usage(){
   sed -n '/^# !CALLING SEQUENCE:/,/^# !/{p}' ./scripts/initial_conditions_era5.sh | head -n -1
}

if [ $# -ne 9 ]; then
   usage
fi

############################################
# Args
#############################################

EXP_IBC=${1}
ERA5_DIR=${2}
LABELI=${3}
LABELF=${4}
ERA5_BCHOURSSTEP=${5}
LAT_INI=${6}
LAT_FIN=${7}
LON_INI=${8}
LON_FIN=${9}

#####################################################################

SCRDIR=${DIR_HOME}/run/scripts

#####################################################################
start_date=${LABELI:0:4}-${LABELI:4:2}-${LABELI:6:2}_${LABELI:8:2}:00:00
end_date=${LABELF:0:4}-${LABELF:4:2}-${LABELF:6:2}_${LABELF:8:2}:00:00

#####################################################################
#ERA5
#####################################################################
if [ ${EXP_IBC} = "ERA5" ]; then

	download="no"

	BNDDIR=${ERA5_DIR}

	if [ ! -d ${BNDDIR} ]; then
		mkdir ${BNDDIR} 
	fi

	if [ ! -e ${BNDDIR}/e5.oper.*.*.${LABELI}.grib ]; then
		download="yes"
	fi

	if [ ! -e ${BNDDIR}/e5.oper.*.*.${LABELF}.grib ]; then
		download="yes"
	fi

	#Checking the folder and the if exists the #data 
	if [ ${download} = "yes" ]; then

	   mkdir ${BNDDIR}/python

	   #copy function of python for the main_download program
	   cp ${SCRDIR}/python/dfunctions.py ${BNDDIR}/python

	   #####################################################################
	   
	   #Program to download directly ERA5	
	   sed -e "s,#datein#,${LABELI},g;\
		s,#dateout#,${LABELF},g; \
		s,#era5_data#,${ERA5_DIR},g;\
		s,#lat_init#,${LAT_INI},g;\
		s,#lat_fin#,${LAT_FIN},g;\
		s,#lon_init#,${LON_INI},g;\
		s,#lon_fin#,${LON_FIN},g;\
		s,#nh#,${ERA5_BCHOURSSTEP},g" \
		${SCRDIR}/python/main_download.py > \
	       	${BNDDIR}/python/main_download_era5_${LABELI}_${LABELF}.py

	   echo "Condicao de contorno inexistente !"
	   echo "Sera baixada do ERA5"
	   echo "$0 ${LABELI}"

	   cd  ${BNDDIR}/python


	   #********************************************************************
	   echo  "RUN: python main_download_era5_${LABELI}_${LABELF}.py"
	   #
	   python main_download_era5_${LABELI}_${LABELF}.py 
           # 
	   #********************************************************************
	
	else

	   echo "Condicao de contorno existente"
	   echo "Na pasta${BNDDIR}"

	fi

fi

#GFS
####################################################################

if [ ${EXP_IBC} = "GFS" ]; then

	OPERDIR=/oper/dados/ioper/tempo/${EXP}
	BNDDIR=$OPERDIR/0p25/brutos/${LABELI:0:4}/${LABELI:4:2}/${LABELI:6:2}/${LABELI:8:2}
fi

start=${LABELI}
end=${LABELF}
current="$start"



while [[ "${current}" -le "${end}" ]]
do

    #*************************************************
    #BC
    if [ ! e5.oper.*.*.$current.grib ]; then

    	echo -e "${RED}==>${NC}A condicao inicial  $current nao  foi baixanda" 
	echo -e "${RED}==>${NC}corretamente na pasta do experimento."
	echo -e "${RED}==>${NC}Verifique antes de rodar o modelo "

    	#*************************************************
    	formatted_date=$(date -d "${current:0:8} ${current:8}:00:00" "+%Y-%m-%d %H")
    fi
    current=$(date -d "$formatted_date + ${ERA5_BCHOURSSTEP} hour" "+%Y%m%d%H")

done

echo -e "${GREEN}==>${NC}Finished initial_boudary_condition.bash script" 

}
