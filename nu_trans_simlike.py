import gammalib
import ctools
import cscripts
import numpy as np
from ebltable.tau_from_model import OptDepth
from random import uniform
import xml_generator as xml
from astropy.io import fits
from xml.dom import minidom

tau =  OptDepth.readmodel(model = 'dominguez')

input_model='NeutrinoAlerts_10000_1e4_transient_100s_MD2014SFR_SC_2.13.out'

imin = 0
imax = 10000

gam = 2.13

ep = 100.

ttrans = 100.

tobscta = 600.

debug = True
edisp = True

caldb='prod3b-v1'
irf='North_z20_average_30m'

hdr = fits.Header()
hdr['EXTNAME'] = 'Time profile'
hdr['MJDREFI'] = '59000'
hdr['MJDREFF'] = '5.0000000000E-01'
hdr['TIMEUNIT'] = 's'
hdr['TIMESYS'] = 'TT'
hdr['TIMEREF'] = 'LOCAL'

declination,redshift,A = np.loadtxt(input_model, skiprows=11, unpack=True)

realsrc=open('nu_src_ts_'+str(int(ttrans))+'s_'+irf+'_'+str(int(tobscta))+'s_'+str(imin+1)+'-'+str(imax)+'.dat', 'w')
lowrealsrc=open('nu_src_low_ts_'+str(int(ttrans))+'s_'+irf+'_'+str(int(tobscta))+'s_'+str(imin+1)+'-'+str(imax)+'.dat', 'w')
fakesrc=open('nu_src_fake_'+str(int(ttrans))+'s_'+irf+'_'+str(int(tobscta))+'s_'+str(imin+1)+'-'+str(imax)+'.dat', 'w')

for i in xrange(imin, imax):
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
            tsig = tobs - tsigstart
            ra = uniform(0.,360.)
            dec = declination[i]
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

            open('nu_sources_'+str(i+1)+'.xml', 'w').write(doc.toprettyxml('  '))
                        
            foutmodel='nu_'+str(i+1)+'_ts.xml'
            
            sim = ctools.ctobssim()
            sim['inmodel']   = 'nu_sources_'+str(i+1)+'.xml'
            sim['caldb']     = caldb
            sim['irf']       = irf
            sim['outevents'] = 'events_nu_'+str(ttrans)+'s_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
            sim['ra']        = ra
            sim['dec']       = dec
            sim['rad']       = 5.0
            sim['tmin']      = '2020-05-31T12:00:00'
            sim['tmax']      = '2020-05-31T12:10:00'
            sim['emin']      = 0.02
            sim['emax']      = 199.0
            sim['maxrate']   = 1.0e9
            sim['debug']     = debug
            sim['edisp']     = edisp
            sim.execute()

            like = ctools.ctlike()
            like['inobs']     = 'events_nu_'+str(ttrans)+'s_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
            like['caldb']     = caldb
            like['irf']       = irf
            like['inmodel']   = 'nu_sources_'+str(i+1)+'.xml'
            like['outmodel']  = foutmodel
            like['debug']     = debug
            like['edisp']     = edisp
            like.execute()
            
            outfile = minidom.parse(foutmodel)
            srcs = outfile.getElementsByTagName('source')
            ts = float(srcs[0].attributes['ts'].value)
            srcsp = outfile.getElementsByTagName('parameter')
            normsp = float(srcsp[0].attributes['value'].value)
            normsp_error = float(srcsp[0].attributes['error'].value)                
                            
            if ts >= 25.:
                if normsp > 2. or normsp < 0.5:
                    fake = str(i+1)+' '+str(ts)+' '+str(normsp)+' '+str(normsp_error)+' '+str(ra)+' '+str(dec)+' '+str(tsig)+'\n'
                    fakesrc.write(fake)
                else:
                    real_nu = str(i+1)+' '+str(ts)+' '+str(normsp)+' '+str(normsp_error)+' '+str(ra)+' '+str(dec)+' '+str(tsig)+'\n'
                    realsrc.write(real_nu)
            else:
                lowreal_nu = str(i+1)+' '+str(ts)+' '+str(normsp)+' '+str(normsp_error)+' '+str(ra)+' '+str(dec)+' '+str(tsig)+'\n'
                lowrealsrc.write(lowreal_nu)

realsrc.close()
lowrealsrc.close()
fakesrc.close()
