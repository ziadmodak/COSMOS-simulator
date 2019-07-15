# Generate a single array of Archive files, required integration time to reach the noise level
# in them and the noise level
path='data/'
files = glob.glob(path+'*.image.fits')
tint = np.array([])
sig = np.array([])
file_data = np.vstack((files, tint.round(2), (5*sig*1e-3).round(4))).T


# Lists of source sizes and fluxes to be simulated
src_size = [0, 0.2, 0.5, 1.0] # FWHM arcsec
fluxes = [0.5, 0.2, 0.1, 0.07, 0.05, 0.03, 0.02, 0.01, 0.005] #mJy

# Loop over the size and flux lists and generate 100 simulations in the matrix of sizes and fluxes
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
			t_int=t_int, direction = im_direction, inp_seed=sim_seeds[sim_seed_ct])
			model_clean(t_int=t_int, dir=im_direction, thresh=sigma)
			f = open('sim_'+str(ct)+'/log.out', 'w')
			f.write(''.join(('File: ', filename, '\n', 'T_int (s): ',
				t_int, '\n', 'Threshold (mJy): ', sigma, '\n', 'Model Source Size (arcsec): ', str(i),
				'\n', 'Model Source Flux (mJy): ', str(j), '\n')))
			f.close()
			ct += 1
			sim_seed_ct+=1
	
os.chdir(path)

