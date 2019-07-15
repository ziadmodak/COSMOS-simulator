import numpy as np
import glob
import os

from astropy.io import fits



def model_simulate(ref_file, outfile, flux, point_src, source_size_maj, source_size_min):
    """ Generate a model CASA image based on the properties of an input image from 
    the COSMOS ALMA Archive.
    
    Args: 
    ref_file: The input Archive file to represent
    outfile: The name given to the output files like model image and pointing file
    flux: Flux (in mJy) of the model
    point_src: When True, makes a point source model instead of a Gaussian source when False
    source_size_maj: Major axis size of the Gaussian source. Not used when point_src = True
    source_size_min: Minor axis size of the Gaussian source. Not used when point_src = True 

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
    
    # Initialise an empty component list
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
    inp_seed: Seed for random number generator used to add noise to the visibilities

    Returns:
    beam_size: The synthesized beam size of the input Archive image
    cell_size: The size of a pixel in the MS

    """

    # Get the beamsize, cellsize, and map size
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



def model_clean(sky_dir, thresh):
	""" Generate a dirty image from the simulated MS.

	Args:
	sky_dir: The central sky direction of the image
	thresh: The threshold (in mJy) of the image to CLEAN down to

	Returns:
	None  

	"""
	default(tclean)
	tclean(vis='sim_'+str(ct)+'/sim_'+str(ct)+'.alma_'+str(beam_size)+'.noisy.ms',
		   imagename='sim_'+str(ct)+'/dirty'+str(ct),
		   imsize = header['NAXIS1'],
		   cell = str(header['CDELT2']*3600)+'arcsec',
		   phasecenter = sky_dir,
		   niter = 0,
		   threshold = thresh+'mJy',
		   gridder='standard',
		   pbcor=True)

	return None
