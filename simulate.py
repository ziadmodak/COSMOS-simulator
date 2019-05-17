#import string
import numpy as np
import glob
#import random
import os

from astropy.io import fits



def model_simulate(ref_file, outfile, flux, source_size_maj, source_size_min, point_src):
    """ Generate a model CASA image file based on the properties of an input image from 
    the COSMOS ALMA Archive.
    
    Args: 
    ref_file: The input Archive file to represent
    outfile: The name given to the output files like (model image and pointing file)
    flux: Flux (in mJy) of the model
    source_size_maj: Major axis size of the Gaussian source
    source_size_min: Minor axis size of the Gaussian source
    point_src: When True makes a point source model instead of a Gaussian source when False 

    Returns:
    hdu: The header of the input Archive file
    im_direction: The central sky coordinates of the generated image

    """

    # Get the header of the input file from the archive and determine the direction
    # of the image to be simulated
    hdu = fits.getheader(ref_file)
    im_ra = au.deg2radec(ra=hdu['CRVAL1'], prec=3, delimiter=' ')[:13]
    im_dec = au.deg2radec(dec=hdu['CRVAL2'], prec=3, delimiter=' ')[13:].replace(':','.')
    im_direction = im_ra + im_dec
    
    # Generate an empty component list
    cl.done()
    
    # Add a component (point source or Gaussian) at the central pixel and at the frequency 
    # of the original .fits image
    if point_src == True:
        cl.addcomponent(dir=im_direction, flux=flux, fluxunit='mJy', freq=str(hdu['CRVAL3'])+'Hz', 
    	    shape='point')
    
    elif point_src == False:
    	cl.addcomponent(dir=im_direction, flux=flux, fluxunit='mJy', freq=str(hdu['CRVAL3'])+'Hz', 
    	    shape='Gaussian', majoraxis=str(source_size_maj)+'arcsec', 
    	    minoraxis=str(source_size_min)+'arcsec', positionangle='0.0deg')
    
    # Generate and export a .fits file using CASA Toolkits
    ia.fromshape(outfile+'.im', [hdu['NAXIS1'],hdu['NAXIS2'],hdu['NAXIS3'],hdu['NAXIS4']], 
    	overwrite=True)
    csys = ia.coordsys()
    csys.setunits(['rad','rad','','Hz'])
    cell_rad = qa.convert(qa.quantity(hdu['CDELT2'], hdu['CUNIT2']),'rad')['value']
    csys.setincrement([-cell_rad,cell_rad],'direction')
    csys.setreferencevalue([qa.convert(im_ra,'rad')['value'],qa.convert(im_dec,'rad')['value']], 
    	type='direction')
    csys.setreferencevalue(hdu['CRVAL3'],'spectral')
    csys.setincrement('10GHz','spectral')
    ia.setcoordsys(csys.torecord())
    ia.setbrightnessunit('Jy/pixel')
    ia.modify(cl.torecord(),subtract=False)
    exportfits(imagename=outfile+'.im', fitsimage=outfile+'.fits', overwrite=True)

    # Generate a pointing file for simobserve() to use
    f = open(outfile+'_ptg.txt', 'w')
    f.write(''.join(('J2000 ', im_direction)))
    f.close()

    # Write the name of the input file into the generated .fits file's header
    fits.setval(filename=outfile+'.fits', keyword='OBSERVER', value=ref_file)
    return hdu, im_direction



def model_simobserve(ref_file, t_int, direction, inp_seed):
    """ Simulate ALMA observations MS from the model file with a beamsize as close as possible
    to the beamsize of the input file from the Archive.

    Args: 
    ref_file: The ALMA archive .fits file to get the beamsize
    t_int: The integration time (in sec) of the simulated ALMA observation
    direction: The central sky direction of the generated MS

    Returns:
    beam_size: The synthesized beam size of the input Archive image
    cell_size: The size of a pixel in the MS

    """

    # Get the beamsize, calculate the cellsize and calculate scan times to cover the map
    beam_size = '%0.3f' %(fits.getval(filename=ref_file, keyword='BMAJ') * 3600)+'arcsec'
    cell_size = str(float('%0.3s' %beam_size) / 8)+'arcsec'
    fwhm_rad = 1.13*3e8/header['CRVAL3']/12
    fwhm = fwhm_rad*180/np.pi*3600

    # Use CASA task simobserve() to simulate ALMA observations
    default(simobserve)
    simobserve(project='sim_'+str(ct),
               skymodel='simulations/model_'+str(ct)+'.fits',
               integration='10s',
               maptype='ALMA',
               indirection = direction,
               incell='',
               inwidth = '7.5GHz',
               incenter='',
               setpointings=False,
               ptgfile='simulations/model_'+str(ct)+'_ptg.txt',
               direction = direction,
               obsmode = 'int',
               antennalist = 'alma;'+beam_size,
               totaltime = t_int+'s', 
               mapsize = str(fwhm)+'arcsec',
               thermalnoise = 'tsys-atm',
               user_pwv = 0.5,
               seed=inp_seed)
    return beam_size, cell_size



def model_clean(t_int, dir, thresh):
	""" Generate a CLEANed image from the simulated MS.

	Args:
	t_int: The integration time (in sec)  ## Remove this, NOT NEEDED
	dir: The central sky direction of the image
	thresh: The threshold (in mJy) of the image to CLEAN down to

	Returns:
	None  

	"""
	default(tclean)
	tclean(vis='sim_'+str(ct)+'/sim_'+str(ct)+'.alma_'+str(beam_size)+'.noisy.ms',
		   imagename='sim_'+str(ct)+'/dirty'+str(ct),
		   imsize = header['NAXIS1'],
		   cell = str(header['CDELT2']*3600)+'arcsec',
		   phasecenter = dir,
		   niter = 0,
		   threshold = thresh+'mJy',
		   gridder='standard',
		   pbcor=True)

	return None

	

# Genearte a single array of Archive files, required integration time to reach the noise level
# in them and the noise level
path='/vol/arc2/archive2/ziad/data/'
setpath=path+'set'
files = glob.glob(path+'*.image.fits')
tint = np.loadtxt(path+'t_x_pntg.out')
sig = np.loadtxt(path+'sigma_x.out')
file_data = np.vstack((files, tint.round(2), (5*sig*1e-3).round(4))).T


# Lists of source sizes and fluxes to be simulated
src_size = [0, 0.2, 0.5, 1.0]
fluxes = [0.5, 0.2, 0.1, 0.07, 0.05, 0.03, 0.02, 0.01, 0.005]
setnum = ['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19']

# Loop over the size and flux lists and generate 100 simulations in the matrix of sizes and fluxes

#### CHECK THE END OF THE FILE FOR PARALLEL CODE IN A HACKY WAY
#np.random.seed(1909765)
#sim_seeds = np.random.randint(0,999999,size=72000)
#set_seed = np.random.randint(0,9999999, size=20)

sim_seed_ct = 0
for x in setnum:
	ct = 0
	os.chdir(setpath+x)
	os.system('rm -rf simulations')
	os.system('rm -rf sim_*')
	os.system('mkdir simulations')
	for i in src_size:
		for j in fluxes:
			np.random.seed(set_seed[int(x)])
			for k in range(0,100):
				filename, t_int, sigma = file_data[np.random.randint(file_data.shape[0]), :]
				#filename, t_int, sigma = file_data[1569, :] ## Test line
				if i == 0:
					header, im_direction = model_simulate(ref_file=filename,
						outfile='simulations/model_'+str(ct), flux=j,
						source_size_maj=i, source_size_min=i, point_src=True)
				else:
					header, im_direction = model_simulate(ref_file=filename,
						outfile='simulations/model_'+str(ct), flux=j,
						source_size_maj=i, source_size_min=i, point_src=False)

				beam_size, cell_size = model_simobserve(ref_file=filename,
					t_int=t_int, direction=im_direction, inp_seed=sim_seeds[sim_seed_ct])
				model_clean(t_int=t_int, dir=im_direction, thresh=sigma)
				f = open('sim_'+str(ct)+'/log.out', 'w')
				f.write(''.join(('File: ', filename, '\n', 'T_int (s): ',
					t_int, '\n', 'Threshold (mJy): ', sigma, '\n', 'Model Source Size (arcsec): ', str(i),
					'\n', 'Model Source Flux (mJy): ', str(j), '\n')))
				f.close()
				ct += 1
				sim_seed_ct+=1
	
os.chdir(path)


