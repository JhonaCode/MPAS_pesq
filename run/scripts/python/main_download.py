'''
**************************
Program Download the data of ERA5 to run de MPAS 

#################################################
To run:

Define :
    INITIAL DATA  
        datain ='1995010100'
    FINAL DATA  
        datain ='1995010123'
    Path to save the data 
        path= '/pwd'

#################################################
Before run the first time :

    see: https://cds.climate.copernicus.eu/api-how-to

    1) Go to  CDS registration page and make the sing in 
    2) Get the uid and api-key in your login information, CLick in your name upper rigth conner 
    3) Install pytho lib of era5: pip install cdsapi 

#################################################
Data:25-04-24
Created by: Jhonatan A. A. Manco
**************************
python 3.9

'''

import  os,sys

import  dfunctions as df 

import subprocess, sys


"""
path='/pesq/dados/bamc/jhonatan.aguirre/DATA/ERA5/MPAS_data'
datein ='2014090100'
dateout='2014090200'
datetd=df.dates2download(datein,dateout,nhpull=1)
lat=[15,-60]
lon=[-90,-30]
"""

#path to save the datain 
path='#era5_data#'

#Inicial day of BC
datein ='#datein#'

#Final day of BC
dateout='#dateout#'

nhour=#nh#

lat=[#lat_init#,#lat_fin#]
lon=[#lon_init#,#lon_fin#]

#Determine which dates to download 
#using the datein and the dateout and the number 
#of hour betweent these dates.
datetd,hours=df.dates2download(datein,dateout,nhpull=nhour)


#Directory
dire='./'
dire2='%s'%(path)

if  os.path.exists(dire2):

    print('-----------------------------------------')
    print('The BC for %s_%s already exist'%(datein,dateout))
    print('The BC will be read for the data folder')
    print('or download in %s if not'%(path))
    print('-----------------------------------------')

#Directory 2
# Check if the directory exists
if not os.path.exists(dire2):
    # If it doesn't exist, create it
    os.makedirs(dire2)


k=0
for date in datetd:

    #Level variables
    df.download_lev('temperature'        ,'t','129',date,dire,lat,lon)
    df.download_lev('geopotential'       ,'z','130',date,dire,lat,lon)
    df.download_lev('u_component_of_wind','u','131',date,dire,lat,lon)
    df.download_lev('v_component_of_wind','v','132',date,dire,lat,lon)
    df.download_lev('relative_humidity'  ,'r','157',date,dire,lat,lon)
    
    #Surface variables
    df.download_sfc('sea_ice_cover'                 ,'ci'   ,'031',date,dire,lat,lon)
    df.download_sfc('snow_density'                  ,'rsn'  ,'033',date,dire,lat,lon)
    df.download_sfc('sea_surface_temperature'       ,'sstk' ,'034',date,dire,lat,lon)
    df.download_sfc('volumetric_soil_water_layer_1' ,'swvl1','039',date,dire,lat,lon)
    df.download_sfc('volumetric_soil_water_layer_2' ,'swvl2','040',date,dire,lat,lon)
    df.download_sfc('volumetric_soil_water_layer_3' ,'swvl3','041',date,dire,lat,lon)
    df.download_sfc('volumetric_soil_water_layer_4' ,'swvl4','042',date,dire,lat,lon)
    df.download_sfc('surface_pressure'              ,'sp'   ,'134',date,dire,lat,lon)
    df.download_sfc('soil_temperature_level_1'      ,'stl1' ,'139',date,dire,lat,lon)
    df.download_sfc('soil_temperature_level_2'      ,'stl2' ,'170',date,dire,lat,lon)
    df.download_sfc('soil_temperature_level_3'      ,'stl3' ,'183',date,dire,lat,lon)
    df.download_sfc('soil_temperature_level_4'      ,'stl4' ,'236',date,dire,lat,lon)
    df.download_sfc('snow_depth'                    ,'sd'   ,'141',date,dire,lat,lon)
    df.download_sfc('mean_sea_level_pressure'       ,'msl'  ,'151',date,dire,lat,lon)
    df.download_sfc('10m_u_component_of_wind'       ,'10u'  ,'165',date,dire,lat,lon)
    df.download_sfc('10m_v_component_of_wind'       ,'10v'  ,'166',date,dire,lat,lon)
    df.download_sfc('2m_dewpoint_temperature'       ,'2d'   ,'167',date,dire,lat,lon)
    df.download_sfc('2m_temperature'                ,'2t'   ,'168',date,dire,lat,lon)
    df.download_sfc('skin_temperature'              ,'skt'  ,'235',date,dire,lat,lon)

    #invariant
    df.download_sfc('geopotential'                  ,'z'    ,'129',date,dire,lat,lon)
    df.download_sfc('land_sea_mask'                 ,'lsm'  ,'172',date,dire,lat,lon)


    #The files are download in the runing folder, after that they are joint 
    #with the comand cat and move to the dataforder
    subprocess.run('cat  *%s.grib > e5.oper.f0%s.ll025sc.%s.grib'%(date,hours[k],date), shell = True, executable="/bin/bash")
    subprocess.run('echo Moving to  %s'%(dire2), shell = True, executable="/bin/bash")
    subprocess.run('mv   e5.oper.f0%s.ll025sc.%s.grib %s'%(hours[k],date,dire2), shell = True, executable="/bin/bash")
    subprocess.run('rm   *.grib', shell = True, executable="/bin/bash")

    k+=1
