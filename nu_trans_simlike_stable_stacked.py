import gammalib
import ctools
import cscripts
import numpy as np
from ebltable.tau_from_model import OptDepth
from random import randint, uniform
import xml_generator as xml
from astropy.io import fits
from xml.dom import minidom

tau =  OptDepth.readmodel(model = 'dominguez')

input_model='3e-9_all.out.alert'

imin = 0

gam = 2.19

ep = 100.

tobscta = 600.

debug = True
edisp = True

caldb='prod3b-v1'
irf='North_z20_average_30m'

declination,redshift,A = np.loadtxt(input_model, unpack=True)

imax = len(redshift)

nusrcts=open('nu_src_ts_'+irf+'_'+str(int(tobscta))+'s_'+str(imin+1)+'-'+str(imax)+'.dat', 'w')

for i in xrange(imin, imax):
    z = redshift[i]
    if z < 4.:
        lib,doc = xml.CreateLib()
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
        speci.appendChild(spatial)

        lib.appendChild(speci)
    
        bkg = xml.addCTAIrfBackground(lib)
        lib.appendChild(bkg)

        open('nu_sources_'+str(i+1)+'.xml', 'w').write(doc.toprettyxml('  '))
                        
        foutmodel='nu_'+str(i+1)+'_ts.xml'
        
        nuseed = randint(1, 1000000000)
            
        sim = ctools.ctobssim()
        sim['inmodel']   = 'nu_sources_'+str(i+1)+'.xml'
        sim['caldb']     = caldb
        sim['irf']       = irf
        sim['outevents'] = 'events_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        sim['ra']        = ra
        sim['dec']       = dec
        sim['rad']       = 5.0
        sim['tmin']      = '2020-05-31T12:00:00'
        sim['tmax']      = '2020-05-31T12:10:00'
        sim['emin']      = 0.02
        sim['emax']      = 199.0
        sim['maxrate']   = 1.0e9
        sim['seed']      = nuseed
        sim['debug']     = debug
        sim['edisp']     = edisp
        sim.execute()
        
        binning = ctools.ctbin()
        binning['inobs']    = 'events_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        binning['xref']     = ra
        binning['yref']     = dec
        binning['coordsys'] = 'CEL'
        binning['proj']     = 'CAR'
        binning['binsz']    = 0.02
        binning['nxpix']    = 500
        binning['nypix']    = 500
        binning['ebinalg']  = 'LOG'
        binning['emin']     = 0.02
        binning['emax']     = 199.0
        binning['enumbins'] = 40
        binning['outcube']  = 'cntcube_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        binning['debug']    = debug
        binning.execute()
       
        expcube = ctools.ctexpcube()
        expcube['inobs']   = 'events_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        expcube['caldb']   = caldb
        expcube['irf']     = irf
        expcube['incube']  = 'cntcube_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        expcube['outcube'] = 'expcube_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        expcube['debug']   = debug
        expcube.execute()
        
        psfcube = ctools.ctpsfcube()
        psfcube['inobs']    = 'events_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        psfcube['caldb']    = caldb
        psfcube['irf']      = irf
        psfcube['incube']   = 'NONE'
        psfcube['ebinalg']  = 'LOG'
        psfcube['emin']     = 0.02
        psfcube['emax']     = 199.0
        psfcube['enumbins'] = 40
        psfcube['nxpix']    = 12
        psfcube['nypix']    = 12
        psfcube['binsz']    = 1.0
        psfcube['coordsys'] = 'CEL'
        psfcube['proj']     = 'CAR'
        psfcube['xref']     = ra
        psfcube['yref']     = dec
        psfcube['outcube']  = 'psfcube_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        psfcube['debug']    = debug
        psfcube.execute()

        bkgcube = ctools.ctbkgcube()
        bkgcube['inobs']    = 'events_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        bkgcube['caldb']    = caldb
        bkgcube['irf']      = irf
        bkgcube['inmodel']  = 'nu_sources_'+str(i+1)+'.xml'
        bkgcube['incube']   = 'cntcube_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        bkgcube['outcube']  = 'bkgcube_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        bkgcube['outmodel'] = 'stacked_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.xml'
        bkgcube['debug']    = debug
        bkgcube.execute()
        
        edispcube = ctools.ctedispcube()
        edispcube['inobs']    = 'events_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        edispcube['caldb']    = caldb
        edispcube['irf']      = irf
        edispcube['incube']   = 'NONE'
        edispcube['ebinalg']  = 'LOG'
        edispcube['emin']     = 0.02
        edispcube['emax']     = 199.0
        edispcube['enumbins'] = 40
        edispcube['nxpix']    = 12
        edispcube['nypix']    = 12
        edispcube['binsz']    = 1.0
        edispcube['coordsys'] = 'CEL'
        edispcube['proj']     = 'CAR'
        edispcube['xref']     = ra
        edispcube['yref']     = dec
        edispcube['outcube']  = 'edispcube_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        edispcube['debug']    = debug
        edispcube.execute()
       
        like = ctools.ctlike()
        like['inobs']     = 'cntcube_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        like['edisp']     = edisp
        like['edispcube'] = 'edispcube_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        like['expcube']   = 'expcube_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        like['psfcube']   = 'psfcube_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        like['bkgcube']   = 'bkgcube_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        like['inmodel']   = 'stacked_nu_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.xml'
        like['outmodel']  = foutmodel
        like['debug']     = debug
        
        like.execute()
            
        outfile = minidom.parse(foutmodel)
        srcs = outfile.getElementsByTagName('source')
        ts = float(srcs[0].attributes['ts'].value)
        srcsp = outfile.getElementsByTagName('parameter')
        normsp = float(srcsp[0].attributes['value'].value)
        normsp_error = float(srcsp[0].attributes['error'].value)                
                            
        real_nu = str(i+1)+' '+str(ts)+' '+str(normsp)+' '+str(normsp_error)+' '+str(ra)+' '+str(dec)+' '+str(nuseed)+'\n'
        nusrcts.write(real_nu)
        
nusrcts.close()
