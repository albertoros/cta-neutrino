import numpy as np
from ebltable.tau_from_model import OptDepth
from random import uniform
import xml_generator as xml
from astropy.io import fits

tau =  OptDepth.readmodel(model = 'dominguez')

input_model='NeutrinoAlerts_10000_1e4_transient_100s_MD2014SFR_SC_2.13.out'

imax = 10000

gam = 2.13

ep = 100.

ttrans = 100.

hdr = fits.Header()
hdr['EXTNAME'] = 'Time profile'
hdr['MJDREFI'] = '59000'
hdr['MJDREFF'] = '5.0000000000E-01'
hdr['TIMEUNIT'] = 's'
hdr['TIMESYS'] = 'TT'
hdr['TIMEREF'] = 'LOCAL'

declination,redshift,A = np.loadtxt(input_model, skiprows=11, unpack=True)

fsrcpoint=open('nu_src_point.dat', 'w')

for i in xrange(0, imax):
    z = redshift[i]
    if z < 4.:
        delem = uniform(0.,ttrans)
        delal = uniform(20.,80.)
        delrp = uniform(20.,50.)
        tobs = ttrans * (1. + z)
        delobs = delem * (1. + z)
        tsigstart = delobs + delal + delrp
        if tsigstart < tobs:
            lib,doc = xml.CreateLib()
            LCfile = 'LC_nu_'+str(i+1)+'.fits'
            ta=np.empty(500)
            na=np.empty(500)
            for j in xrange(0, 500):
                ta[j]=j
                if ta[j] < tsigstart:
                    na[j] = 0.
                elif tsigstart < ta[j] < tobs:
                    na[j] = 1.
                else:
                    na[j] = 0.
            time = fits.Column(name='TIME', array=ta, format='1D', unit='s')
            norm = fits.Column(name='NORM', array=na, format='1D')
            t = fits.BinTableHDU.from_columns([time, norm],header=hdr)
            t.writeto(LCfile,overwrite=True) 
            ra = uniform(0.,360.)
            dec = declination[i]
            pointing = str(i+1)+' '+str(ra)+' '+str(dec)+'\r\n'
            fsrcpoint.write(pointing)
            ETeV = np.logspace(-2,2.5,45)
            EMeV = ETeV * 1e6
            if z < 0.01:
                atten = 1.
            else:
                atten = np.exp(-1. * tau.opt_depth(z,ETeV))
            prefac = A[i] * 1e-13
            spec = prefac * (ETeV / ep) ** (-gam)
            specebl = spec * atten
            sourcename = 'nu'+str(i+1)
            Filefunction = 'spec_nu_ebl_'+str(i+1)+'.dat'
            np.savetxt(Filefunction, np.column_stack([EMeV,specebl + 1.e-300]))
            speci = xml.addFileFunction(lib, sourcename, type = "PointSource", filefun=Filefunction, flux_free=1, flux_value=1., flux_scale=1., flux_max=100000000.0, flux_min=0.0)
            spatial = xml.AddPointLike(doc,ra,dec)
            temporal = xml.AddLCTrans(doc, LCfile, 1.)
            speci.appendChild(spatial)
            speci.appendChild(temporal)
            lib.appendChild(speci)
    
            bkg = xml.addCTAIrfBackground(lib)
            lib.appendChild(bkg)

            open('nu_sources'+str(i+1)+'.xml', 'w').write(doc.toprettyxml('  '))
            
fsrcpoint.close()
