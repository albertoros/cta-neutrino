import gammalib
import ctools
import cscripts
import numpy as np
from ebltable.tau_from_model import OptDepth
from random import randint, uniform
import xml_generator as xml
from astropy.io import fits
from xml.dom import minidom
import argparse

tau =  OptDepth.readmodel(model = 'dominguez')

parser = argparse.ArgumentParser()
parser.add_argument('-alert', action='store', dest='alertfile',
                        default='3e-9_all.out.alert', help='File with alerts')
parser.add_argument('--irf', action='store', dest='irf',
                        default='North_z20_average_30m', help='IRF')
parser.add_argument('--obs', action='store', dest='tobs',
                        type=float, default=600.,
                        help='Observation duration time in [s]')
parser.add_argument('--inter', action='store', dest='interaction',
                        default='no',
                        help='Interaction type: pp (proton-proton), pph (proton-photon), txs (TXS-like sources), no (no scaling)')
options = parser.parse_args()

input_model = options.alertfile

imin = 0

gam = 2.19

ep = 100.

tobscta = options.tobs

debug = True
edisp = True

caldb='prod3b-v1'
irf=options.irf

declination,redshift,A = np.loadtxt(input_model, unpack=True)

# flux scaling according to intearction type pp, p-gamma or no scaling
#Ref https://arxiv.org/pdf/1805.11112.pdf
if options.interaction == 'no':
    A_prefix = 1.0
if options.interaction == 'pp':
    A_prefix = np.pow(2.,-gam-1)
if options.interaction == 'pph':
    A_prefix = np.pow(2.,-gam)

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
        if options.interaction == 'txs': # reference: https://arxiv.org/abs/1811.07439
            prefac = A[i] * 1e-13
            spec = prefac * (ETeV / ep) ** (-2) * exp(-0.1*(z+1)/ETeV - ETeV/(20.*(z+1)))
        else:
            prefac = A[i] * A_prefix * 1e-13
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
        emin=0.02
        emax=199.0
        binz=0.1
        x_pixs=100
        y_pixs=100
        n_bins=40

        sim = ctools.ctobssim()
        sim['inmodel']   = 'nu_sources_'+str(i+1)+'.xml'
        sim['caldb']     = caldb
        sim['irf']       = irf
        sim['outevents'] = 'events_nu_'+'s_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        sim['ra']        = ra
        sim['dec']       = dec
        sim['rad']       = 5.0
        sim['tmin']      = '2020-05-31T12:00:00'
        sim['tmax']      = '2020-05-31T12:10:00'
        sim['emin']      = emin
        sim['emax']      = emax
        sim['maxrate']   = 1.0e9
        sim['debug']     = debug
        sim['edisp']     = edisp
        sim.execute()
        
        c_bin = ctools.ctbin()
        c_bin['inobs']       = 'events_nu_'+'s_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        c_bin['xref']        = ra
        c_bin['yref']        = dec
        c_bin['proj']        = 'CAR'
        c_bin['coordsys']    = 'CEL'
        c_bin['binsz']       = binz
        c_bin['nxpix']       = x_pixs
        c_bin['nypix']       = y_pixs
        c_bin['ebinalg']     = 'LOG'
        c_bin['emin']        = emin
        c_bin['emax']        = emax
        c_bin['enumbins']    = n_bins
        c_bin['outcube']     = 'cube_'+str(i+1)+'.fits'
        c_bin['debug']       = debug
        c_bin.execute()
        
        exp_cube = ctools.ctexpcube()
        exp_cube['inobs']       ='events_nu_'+'s_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        exp_cube['incube']      = 'cube_'+str(i+1)+'.fits'
        exp_cube['caldb']       = caldb
        exp_cube['irf']         = irf
        exp_cube['xref']        = ra
        exp_cube['yref']        = dec
        exp_cube['proj']        = 'CAR'
        exp_cube['coordsys']    = 'CEL'
        exp_cube['binsz']       = binz
        exp_cube['nxpix']       = x_pixs
        exp_cube['nypix']       = y_pixs
        exp_cube['ebinalg']     = 'LOG'
        exp_cube['emin']        = emin
        exp_cube['emax']        = emax
        exp_cube['enumbins']    = n_bins
        exp_cube['outcube']     = 'exp_cube_'+str(i+1)+'.fits'
        exp_cube['debug']       = debug
        exp_cube.execute()
        
        psf_cube = ctools.ctpsfcube()
        psf_cube['inobs']       ='events_nu_'+'s_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        psf_cube['incube']      = 'NONE'
        psf_cube['caldb']       = caldb
        psf_cube['irf']         = irf
        psf_cube['xref']        = ra
        psf_cube['yref']        = dec
        psf_cube['proj']        = 'CAR'
        psf_cube['coordsys']    = 'CEL'
        psf_cube['binsz']       = binz
        psf_cube['nxpix']       = x_pixs
        psf_cube['nypix']       = y_pixs
        psf_cube['ebinalg']     = 'LOG'
        psf_cube['emin']        = emin
        psf_cube['emax']        = emax
        psf_cube['enumbins']    = n_bins
        psf_cube['outcube']     = 'psf_cube_'+str(i+1)+'.fits'
        psf_cube['debug']       = debug
        psf_cube.execute()
        
        bkg_cube = ctools.ctbkgcube()
        bkg_cube['inobs']       ='events_nu_'+'s_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        bkg_cube['incube']      = 'cube_'+str(i+1)+'.fits'
        bkg_cube['inmodel']   = 'nu_sources_'+str(i+1)+'.xml'
        bkg_cube['caldb']       = caldb
        bkg_cube['irf']         = irf
        bkg_cube['xref']        = ra
        bkg_cube['yref']        = dec
        bkg_cube['proj']        = 'CAR'
        bkg_cube['coordsys']    = 'CEL'
        bkg_cube['binsz']       = binz
        bkg_cube['nxpix']       = x_pixs
        bkg_cube['nypix']       = y_pixs
        bkg_cube['ebinalg']     = 'LOG'
        bkg_cube['emin']        = emin
        bkg_cube['emax']        = emax
        bkg_cube['enumbins']    = n_bins
        bkg_cube['outcube']     = 'bkg_cube_'+str(i+1)+'.fits'
        bkg_cube['outmodel']     = 'bkg_cube_'+str(i+1)+'.xml'
        bkg_cube['debug']       = debug
        bkg_cube.execute()
        
        edisp_cube = ctools.ctedispcube()
        edisp_cube['inobs']     = 'events_nu_'+'s_'+irf+'_'+str(int(tobscta))+'s_'+str(i+1)+'.fits'
        edisp_cube['caldb']     = caldb
        edisp_cube['irf']       = irf
        edisp_cube['incube']    = 'NONE'
        edisp_cube['xref']      = ra
        edisp_cube['yref']      = dec
        edisp_cube['proj']      = 'CAR'
        edisp_cube['coordsys']  = 'CEL'
        edisp_cube['binsz']       = binz
        edisp_cube['nxpix']       = x_pixs
        edisp_cube['nypix']       = y_pixs
        edisp_cube['ebinalg']     = 'LOG'
        edisp_cube['emin']        = emin
        edisp_cube['emax']        = emax
        edisp_cube['enumbins']    = n_bins
        edisp_cube['outcube']     = 'edisp_cube_'+str(i+1)+'.fits'
        edisp_cube['debug']       = debug
        edisp_cube.execute()

        like = ctools.ctlike()
        like['inobs']     = 'cube_'+str(i+1)+'.fits'
        like['expcube']   = 'exp_cube_'+str(i+1)+'.fits'
        like['psfcube']   = 'psf_cube_'+str(i+1)+'.fits'
        like['edispcube']   = 'edisp_cube_'+str(i+1)+'.fits'
        like['bkgcube']   = 'edisp_cube_'+str(i+1)+'.fits'
        like['inmodel']   = 'bkg_cube_'+str(i+1)+'.xml'
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
        
        real_nu = str(i+1)+' '+str(ts)+' '+str(normsp)+' '+str(normsp_error)+' '+str(ra)+' '+str(dec)+' '+str(nuseed)+'\n'
        nusrcts.write(real_nu)
        
nusrcts.close()
