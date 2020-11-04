from __future__ import print_function
import matplotlib
matplotlib.use('agg')
import numpy as np
import pandas as pd
from astropy.io import fits
#import fitsio
from scipy import interpolate
import glob
import math
import os

import h5py
import bisect
from matplotlib import rc
import matplotlib.pyplot as plt


def fit_to_kasen(kasen_dirname, wavelength_array, df_obs_new, t, z,zerr):

    names = glob.glob(kasen_dirname+'/knova_*.h5')
    #print(kasen_dirname)
    #print(names)
    #print ('this are the names of kn templates')
    bestSum2 = 9.999e+99
    bestName = ''
    bestz=0
    c = 2.99e10   # speed of light in cm/sec...
    
    # Loop over all models at redshift z...

    z_min=z-zerr
    if z_min <0.001:
        z_min=0.001
   
    z_max=z+zerr
    zs=np.arange(z_min,z_max+0.001,0.001)
 
    for z_test in zs:
        for name in names:
        
            # open model file
            fin    = h5py.File(name,'r')

            # frequency in Hz
            nu    = np.array(fin['nu'],dtype='d')
            # array of time in seconds
            times = np.array(fin['time'])
            # covert time to days
            times = times/3600.0/24.0

            # specific luminosity (ergs/s/Hz)
            # this is a 2D array, Lnu[times][nu]
            Lnu_all   = np.array(fin['Lnu'],dtype='d')

            # index corresponding to t
            it = bisect.bisect(times,t)
            # spectrum at this epoch
            Lnu = Lnu_all[it,:]

            # if you want thing in Flambda (ergs/s/Angstrom)
            lam0  = c/nu*1e8
            lam   = lam0*(1+z_test)
            Llam = Lnu*nu**2.0/c/1e8

            df_model = pd.DataFrame({'LAMBDA0':lam0, 'LAMBDA':lam, 'Llam':Llam})

            spec_flux_model = interpolate.interp1d(df_model.LAMBDA, df_model.Llam,bounds_error=False, fill_value=0.,kind='linear')
            spec_flux_model_array = spec_flux_model(wavelength_array)

            df_model_new = pd.DataFrame({'LAMBDA':wavelength_array, 'Llam':spec_flux_model_array})

            norm = df_model_new['Llam'].median()

            df_model_new['normLlam'] = df_model_new['Llam']/norm

            merged = pd.merge(df_obs_new, df_model_new, on='LAMBDA')
        
            merged['delta_flux'] = merged['normFLUX'] - merged['normLlam']
            mask=[]
            mask = ((merged.normFLUX.isnull()) | (merged.normLlam.isnull()))
            mask = (mask == False)
        
            if np.any(mask):

            
                sum2 = np.sum(np.square(merged[mask].delta_flux)) / np.sum(mask)

                if sum2 < bestSum2:
                    bestSum2 = sum2
                    bestName = name
                    bestz=z_test
                    #print(type(np.sum(np.square(merged[mask].delta_flux))))
                    #print (type(np.sum(mask)))  
                    
        #else:
        #    print ('Empty files!')
            
       
    return bestName, bestSum2,bestz



def plot_obs_and_bestfit_model(obsname, kasen_dirname, name, wavelength_array, df_obs_new, t, z, sum2,plot_name='kn_fit.png', z_best=0, tex = True, pos=[0.7,0.2]):

    c = 2.99e10   # speed of light in cm/sec...
    if z_best==0:
       z_best=z

    basename = name
    name = kasen_dirname+'/'+basename

    fin    = h5py.File(name,'r')

    # frequency in Hz
    nu    = np.array(fin['nu'],dtype='d')
    # array of time in seconds
    times = np.array(fin['time'])
    # covert time to days
    times = times/3600.0/24.0

    # specific luminosity (ergs/s/Hz)
    # this is a 2D array, Lnu[times][nu]
    Lnu_all   = np.array(fin['Lnu'],dtype='d')

    # index corresponding to t
    it = bisect.bisect(times,t)
    # spectrum at this epoch
    Lnu = Lnu_all[it,:]

    # if you want thing in Flambda (ergs/s/Angstrom)
    lam0  = c/nu*1e8
    lam   = lam0*(1+z_best)
    Llam = Lnu*nu**2.0/c/1e8

    df_model = pd.DataFrame({'LAMBDA0':lam0, 'LAMBDA':lam, 'Llam':Llam})

    spec_flux_model = interpolate.interp1d(df_model.LAMBDA, df_model.Llam,bounds_error=False, fill_value=0.,kind='linear')
    spec_flux_model_array = spec_flux_model(wavelength_array)

    df_model_new = pd.DataFrame({'LAMBDA':wavelength_array, 'Llam':spec_flux_model_array})

    norm = df_model_new['Llam'].median()

    df_model_new['normLlam'] = df_model_new['Llam']/norm
    if (tex == True):
        rc('text', usetex=True)
    else:
        rc('text', usetex=False)
    ax = df_model_new.plot('LAMBDA','normLlam',grid=True)
    df_obs_new.plot('LAMBDA', 'normFLUX', ax=ax)
    title = """ $MSE$=%.2f, $z$=%.3f, $z_b$=%.3f""" % (sum2,z,z_best)
    #plt.figtext(pos[0], pos[1], title)
    plt.title(title)
    plt.savefig(plot_name)

    return 0








######================================ here code starts









kasen_dirname='/home/cleciobom/lib/Kasen_Kilonova_Models_2017/systematic_kilonova_model_grid'#'/home/cleciobom/lib/Kasen_Kilonova_Models_2017/kilonova_models'#


# Redshift and 1sigma error in redshift based on LIGO luminosity distance, 1sigma error in luminosity distance...
z = 0.059
sigma_z = 0.011  # not currently used..


#inputFile_obs_list = ['2019noq.flm', 
#                      'AT2019npw-3.flm', 
#                      'AT2019ntp.flm', 
#                      'AT2019num.flm', 
#                      '2019ntn.flm', 
#                      'AT2019ntr.flm', 
#                      'AT2019omx.flm' ]

# Time [in days] past the merger event when the spectrum was observed...
#t_dict = {'2019noq.flm': 6., 
#          'AT2019npw-3.flm': 12., 
#          'AT2019ntp.flm': 17., 
#          'AT2019num.flm': 12., 
#          '2019ntn.flm': 6., 
#          'AT2019ntr.flm': 14., 
#          'AT2019omx.flm': 14.} 


inputFile_obs_list = ['2019noq.flm', 
                      'AT2019npw-3.flm', 
                      'AT2019ntp.flm', 
                      'AT2019num.flm', 
                      'desgw-190814c.flm', 
                      '2019ntn.flm', 
                      'AT2019ntr.flm', 
                      'AT2019omx.flm', 
                      'desgw-190814d.flm']
#                      'AT2019nte.flm', 

inputFile_dir='/home/cleciobom/lib/S190814bv_UCSC/'

#inputFile_obs_list=[inputFile_dir+inputFile_obs_list[s] for s in inputFile_obs_list] 

# Time [in days] past the merger event when the spectrum was observed...
t_dict = {'2019noq.flm': 6.49, 
          'AT2019npw-3.flm': 12.45, 
          'AT2019ntp.flm': 17.41, 
          'AT2019num.flm': 12.45, 
          'desgw-190814c.flm': 2.44, 
          '2019ntn.flm': 6.49, 
          'AT2019ntr.flm': 14.45, 
          'AT2019omx.flm': 14.45, 
          'desgw-190814d.flm': 2.44}
#          'AT2019nte.flm': ???,


for inputFile_obs in inputFile_obs_list:
    
    t = t_dict[inputFile_obs]

    obsname = os.path.basename(inputFile_obs)
    print(inputFile_dir+inputFile_obs) 
    
    df_obs = pd.read_csv(inputFile_dir+inputFile_obs, sep='\s+', header = None, names = ['LAMBDA', 'FLUX'])
    
    lambda_lo = df_obs.LAMBDA.min()
    lambda_hi = df_obs.LAMBDA.max()

    wavelength_array = np.arange(lambda_lo, lambda_hi, 10.)

    spec_flux_obs = interpolate.interp1d(df_obs.LAMBDA, df_obs.FLUX,bounds_error=False, fill_value=0.,kind='linear')
    spec_flux_obs_array = spec_flux_obs(wavelength_array)

    df_obs_new = pd.DataFrame({'LAMBDA':wavelength_array, 'FLUX':spec_flux_obs_array})

    norm = df_obs_new['FLUX'].median()
    df_obs['normFLUX'] = df_obs['FLUX']/norm
    df_obs_new['normFLUX'] = df_obs_new['FLUX']/norm

    

    name, sum2,z_best = fit_to_kasen(kasen_dirname, wavelength_array, df_obs_new, t, z, sigma_z)
    
    name = os.path.basename(name)
    
    print (name, sum2)
    
    status = plot_obs_and_bestfit_model(obsname, kasen_dirname, name, wavelength_array, df_obs_new, t, z, sum2,plot_name=obsname+'_'+name+'.png',z_best=z_best)


