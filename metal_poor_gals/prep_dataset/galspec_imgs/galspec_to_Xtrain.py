# Imports
#import bz2
import datetime
import numpy as np
import pdb
import os
import urllib.request
import warnings
warnings.filterwarnings('ignore')

from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u, wcs

# Load the galSpec query results
#galspec_final = pd.read_csv('../galspec_sf_final')
start_time = datetime.datetime.now()

galspec_final = np.load('galspec_final_Xy.npy', allow_pickle=True)
ra = galspec_final[:,1]
dec = galspec_final[:,2]
uimg = galspec_final[:,7]

Ngal_galspec = len(galspec_final) # number of systems in galSpec

pix_sz = 16 # half of the pixel size of the image
galspec_flux = np.zeros((Ngal_galspec, 2*pix_sz, 2*pix_sz, 5)) # empty array to store 32x32 fluxes

for i in range(Ngal_galspec):
	start = datetime.datetime.now()

	# Save name of the u-band .fits file
	#url = galspec_final['uimg'][i]
	#u_frame = galspec_final['uimg'][i].split('/')[-1] #.strip('.bz2')
	url = uimg[i]
	u_frame = uimg[i].split('/')[-1]
	print ('Working on galaxy', i, ':', u_frame)

	# Download SDSS ugriz .fits frames if they do not exist and unzip
	if not os.path.exists(u_frame):
		print ('\t Downloading SDSS ugriz .fits frames')
		urllib.request.urlretrieve(url, os.getcwd() + '/' + u_frame)
		urllib.request.urlretrieve(url.replace('-u-', '-g-'), os.getcwd() + '/' + u_frame.replace('-u-', '-g-'))
		urllib.request.urlretrieve(url.replace('-u-', '-r-'), os.getcwd() + '/' + u_frame.replace('-u-', '-r-'))
		urllib.request.urlretrieve(url.replace('-u-', '-i-'), os.getcwd() + '/' + u_frame.replace('-u-', '-i-'))
		urllib.request.urlretrieve(url.replace('-u-', '-z-'), os.getcwd() + '/' + u_frame.replace('-u-', '-z-'))

		#print ('\t Unzipping .fits files')
		#for item in os.listdir(os.getcwd()):
		#	if item.endswith('.bz2'):
		#		zipfile = bz2.BZ2File(item)
		#		data = zipfile.read()
		#		open(item[:-4], 'wb').write(data)
		#		os.remove(item)

	else:
		print ('\t SDSS ugriz .fits frames already exist for this system')

	# Open ugriz .fits files, replacing letter 'u' in u_frame for the other bands
	print ('\t Loading .fits files and coords')
	uband = fits.open(u_frame)
	gband = fits.open(u_frame.replace('-u-', '-g-'))
	rband = fits.open(u_frame.replace('-u-', '-r-'))
	iband = fits.open(u_frame.replace('-u-', '-i-'))
	zband = fits.open(u_frame.replace('-u-', '-z-'))

	#coords = SkyCoord(galspec_final['ra'][i], galspec_final['dec'][i], unit=(u.deg, u.deg))
	coords = SkyCoord(ra[i], dec[i], unit=(u.deg, u.deg))

	# Calculate the center of the galaxy in (x, y) pixel values given the RA, DEC
	print ('\t Calculating central pixels and cropping to 32x32 box')
	ux, uy = wcs.utils.skycoord_to_pixel(coords, wcs.WCS(uband[0].header))
	gx, gy = wcs.utils.skycoord_to_pixel(coords, wcs.WCS(gband[0].header))
	rx, ry = wcs.utils.skycoord_to_pixel(coords, wcs.WCS(rband[0].header))
	ix, iy = wcs.utils.skycoord_to_pixel(coords, wcs.WCS(iband[0].header))
	zx, zy = wcs.utils.skycoord_to_pixel(coords, wcs.WCS(zband[0].header))

	#print (ux, uy)
	#print (gx, gy)
	#print (rx, ry)
	#print (ix, iy)
	#print (zx, zy)

	# x and y range defining the 32x32 pixel image
	uxmin, uxmax = int(np.round(ux))-pix_sz, int(np.round(ux, 0))+pix_sz
	uymin, uymax = int(np.round(uy, 0))-pix_sz, int(np.round(uy, 0))+pix_sz

	gxmin, gxmax = int(np.round(gx, 0))-pix_sz, int(np.round(gx, 0))+pix_sz
	gymin, gymax = int(np.round(gy, 0))-pix_sz, int(np.round(gy, 0))+pix_sz

	rxmin, rxmax = int(np.round(rx, 0))-pix_sz, int(np.round(rx, 0))+pix_sz
	rymin, rymax = int(np.round(ry, 0))-pix_sz, int(np.round(ry, 0))+pix_sz

	ixmin, ixmax = int(np.round(ix, 0))-pix_sz, int(np.round(ix, 0))+pix_sz
	iymin, iymax = int(np.round(iy, 0))-pix_sz, int(np.round(iy, 0))+pix_sz

	zxmin, zxmax = int(np.round(zx, 0))-pix_sz, int(np.round(zx, 0))+pix_sz
	zymin, zymax = int(np.round(zy, 0))-pix_sz, int(np.round(zy, 0))+pix_sz

	# Save flux values into giant array
	print ('\t Saving flux values for this system')
	galspec_flux[i, :, :, 0] = uband[0].data[uymin:uymax, uxmin:uxmax]
	galspec_flux[i, :, :, 1] = gband[0].data[gymin:gymax, gxmin:gxmax]
	galspec_flux[i, :, :, 2] = rband[0].data[rymin:rymax, rxmin:rxmax]
	galspec_flux[i, :, :, 3] = iband[0].data[iymin:iymax, ixmin:ixmax]
	galspec_flux[i, :, :, 4] = zband[0].data[zymin:zymax, zxmin:zxmax]

	# Delete SDSS ugriz .fits files if the next system doesn't use the same files
	if uimg[i] == uimg[i+1]:
		print ('\t Next system in same SDSS field -- not deleting files')

	else:
		print ('\t Deleting SDSS .fits files')
		os.remove(u_frame)
		os.remove(u_frame.replace('-u-', '-g-'))
		os.remove(u_frame.replace('-u-', '-r-'))
		os.remove(u_frame.replace('-u-', '-i-'))
		os.remove(u_frame.replace('-u-', '-z-'))

	if i % 1000 == 0 and i != 0:
		print ('\t Saving binary flux file for 1000 systems')
		np.save('galspec_flux_thru_{0:d}'.format(i), galspec_flux)
	
	end = datetime.datetime.now()
	print ('\t Took', (end-start).seconds, 'seconds \n')

print ('Done with all galaxies..saving final binary galspec_flux file..!')
np.save('galspec_flux_all', galspec_flux)

print ('Took', datetime.datetime.now()-start_time)