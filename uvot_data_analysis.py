#!/usr/bin/env python
# -*- coding utf-8 -*-
#Created by ypc

import glob
import os
import os.path as path
from tqdm import tqdm
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.table import Table

def check_ext(filepaths):
	'''searching all image file absolute path and classifying these as single extension and multiple extension file paths, respectively
	'''

	sin_img_ext_paths = []
	multi_img_ext_paths = []
	for filepath in tqdm(filepaths):
		filename = path.split(filepath)[1]
		print('\n****************************************************************************************************************************')
		print('checking the number of image extenson for '+filepath + '.')
		mainfits = fits.open(filepath)
		num_image_ext = []
		for extension in mainfits[0:]:
			try:
				keyword = extension.header['XTENSION']
				num_image_ext.append(keyword)
			except KeyError:
				print("Keyword 'XTENSION' not found in " + str(extension) + '.')
		if num_image_ext.count('IMAGE') > 1:
			multi_img_ext_paths.append(filepath)
			print(filename + ' contains multiple image extensions.')
		else:
			sin_img_ext_paths.append(filepath)
			print(filename + ' has only a single image extensions.')
	return [sin_img_ext_paths,multi_img_ext_paths]
	
def sum_multi_ext(mul_ext_path,uvotimsumfile):
	'''Wrapper for running ``UVOTIMSUM`` on all Multiple Image Extension file.
	   After running the wrapper, the multiple extension files will become single extension file.
	   
	Returns:
            uvotimsumfile(str): The absolute path of the output file obtained by running the'uvot' tool
	'''
	
	uvotimsum = 'uvotimsum infile=%s outfile=%s clobber=yes' %(mul_ext_path,uvotimsumfile)
	os.system(uvotimsum)
	return uvotimsumfile
	
def measure_flux(sin_ext_path,uvotsourcefile):
	'''Wrapper for running ``UVOTSOURCE`` on all Single Image Extension
	'''
	
	uvotsource = 'uvotsource image=%s srcreg=src.reg bkgreg=bkg.reg sigma=3 outfile=%s clobber=YES' %(sin_ext_path,uvotsourcefile)
	os.system(uvotsource)

def get_observation_time(sin_ext_path):
	'''getting the start and end time from the header of image file.
	
	Returns:
            (float): the middle of the observation time in ISOt format
	'''
	
	mainfits = fits.open(sin_ext_path)
	obs_start = Time(mainfits[0].header['DATE-OBS'],format='isot')
	obs_end = Time(mainfits[0].header['DATE-END'],format='isot')
	return obs_start + (obs_end - obs_start)/2

def get_observation_data(sin_ext_path,uvotsourcefile,obsID):
	''' Parse the output of ``UVOTSOURCE`` to extract essential photometry information.
	
	returns:
	list: observing band, observing time (MJD), magnitude and 1$\sigma$ error in UVOT system, 
	flux and 1$\sigma$ error in erg/s/cm$^{2}$/$\AA$, flux and 1$\sigma$ error in milliJanskies
	'''
	obsid = obsID
	obstime = get_observation_time(sin_ext_path)
	data = fits.getdata(uvotsourcefile)
	band = data['FILTER'][0]
#	ct_rate = data['CORR_RATE'][0]
#	ct_rate_err = data['CORR_RATE_ERR'][0]
#	mag = data['MAG'][0]
#	magerr = data['MAG_ERR'][0]
#	flux_AA = data['FLUX_AA'][0]
#	flux_AA_err = data['FLUX_AA_ERR'][0]
	flux_mjy = data['FLUX_HZ'][0]
	flux_mjy_err = data['FLUX_HZ_ERR'][0]
#	return [obsid,band,obstime.mjd,ct_rate,ct_rate_err,mag,magerr,flux_AA,flux_AA_err,flux_mjy,flux_mjy_err]
#	return [obsid,band,obstime.mjd,flux_AA,flux_AA_err]
	return [obsid,band,obstime.mjd,flux_mjy,flux_mjy_err]
#--------------------------------------------------------------------------------------------------------------------

filepaths = glob.glob('/home/pengcheng/Work/data/J1753/00*/uvot/image/*sk.img')
#filepaths = glob.glob('/home/pengcheng/pythonwork/test/00*/uvot/image/*sk.img')  #It's for testing

#Listing all absolute paths of Single Image Extension file and Multiple Image Extension file.
print('###############################################')
print('         creating Absolute Path List')
print('###############################################\n')
ls = check_ext(filepaths)
sin_ext_paths= ls[0]    #Absolute path list of Single Image Extension file
mul_ext_paths = ls[1]   #Absolute path list of Multiple Image Extension file

#Running 'UVOTIMSUM' tool on the Multiple Image Extension file and 
#adding the absolute path of the output file to the absolute path list of the Single Image Extension file.
print('\n####################################################################')
print('  Running ''UVOTIMSUM'' tool on the Multiple Image Extension file')
print('####################################################################\n')
for mul_ext_path in tqdm(mul_ext_paths):
	uvotimsumfile = '%s.comb' %mul_ext_path
	summed = sum_multi_ext(mul_ext_path,uvotimsumfile)
	sin_ext_paths.append(summed)
'''
#creating Astropy table for storing the photometry information
ptab = Table(names=('obsID','Filter','MJD','CountRate','CountRateErr','Mag','MagErr','FluxDensity','FluxDensityErr','FluxDensitymJy','FluxDensitymJyErr'),
                     dtype=('S11','S4','f8','f8','f8','f8','f8','f8','f8','f8','f8'))
#defining units for each column
ptab['MJD'].unit = u.d
ptab['CountRate'].unit = u.ct/u.second
ptab['CountRateErr'].unit = u.ct/u.second
ptab['Mag'].unit = u.mag
ptab['MagErr'].unit = u.mag
ptab['FluxDensity'].unit = u.erg/u.cm/u.cm/u.second/u.AA
ptab['FluxDensityErr'].unit = u.erg/u.cm/u.cm/u.second/u.AA
ptab['FluxDensitymJy'].unit = u.mJy
ptab['FluxDensitymJyErr'].unit = u.mJy

#1. Runing 'UVOTSOURCE' tool on all Single Image Extension files for calculating flux
#2. extracting photometry information from all 'uvotsourcefile' and adding photometry information to Astropy table row by row.
#3. sorting the table information by 'filter'
print('\n##############################################################################################################')
print('  Running ''UVOTSOURCE'' Tool on the Image File AND Creating a table for storing the photometry information')
print('##############################################################################################################\n')
for sin_ext_path in tqdm(sin_ext_paths):
	dirpath,filename = path.split(sin_ext_path)
	base,extn = filename.split('_')
	obsID = base[2:-3]
	band = base[-2:]
	uvotsourcefile = path.join(dirpath,'uvotsource_%s_%s.fits' %(obsID,band))
	measure_flux(sin_ext_path,uvotsourcefile)
	objphot = get_observation_data(sin_ext_path,uvotsourcefile,obsID)
	ptab.add_row(objphot)

sort_ptab = ptab.group_by(['Filter','MJD'])
sort_ptab.write('photometry.fits',overwrite=True)
'''

#1. Runing 'UVOTSOURCE' tool on all Single Image Extension files for calculating flux density.

print('\n##############################################################################################################')
print('  Running ''UVOTSOURCE'' Tool on the Image File AND Collecting the photometry information')
print('##############################################################################################################\n')

v_flux_density = open('v_flux_density.dat', 'w')
b_flux_density = open('b_flux_density.dat', 'w')
u_flux_density = open('u_flux_density.dat', 'w')
uvw1_flux_density = open('uvw1_flux_density.dat', 'w')
uvm2_flux_density = open('uvm2_flux_density.dat', 'w')
uvw2_flux_density = open('uvw2_flux_density.dat', 'w')

for sin_ext_path in tqdm(sin_ext_paths):
	dirpath,filename = path.split(sin_ext_path)
	base,extn = filename.split('_')
	obsID = base[2:-3]
	band = base[-2:]
	uvotsourcefile = path.join(dirpath,'uvotsource_%s_%s.fits' %(obsID,band))
	measure_flux(sin_ext_path,uvotsourcefile)
	objphot = get_observation_data(sin_ext_path,uvotsourcefile,obsID)
	if band == 'vv':
		print(objphot[0],objphot[1],objphot[2],objphot[3],objphot[4],file=v_flux_density)
	if band == 'bb':
		print(objphot[0],objphot[1],objphot[2],objphot[3],objphot[4],file=b_flux_density)
	if band == 'uu':
		print(objphot[0],objphot[1],objphot[2],objphot[3],objphot[4],file=u_flux_density)
	if band == 'w1':
		print(objphot[0],objphot[1],objphot[2],objphot[3],objphot[4],file=uvw1_flux_density)
	if band == 'm2':
		print(objphot[0],objphot[1],objphot[2],objphot[3],objphot[4],file=uvm2_flux_density)
	if band == 'w2':
		print(objphot[0],objphot[1],objphot[2],objphot[3],objphot[4],file=uvw2_flux_density)

v_flux_density.close()
b_flux_density.close()
u_flux_density.close()
uvw1_flux_density.close()
uvm2_flux_density.close()
uvw2_flux_density.close()













































	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
				
