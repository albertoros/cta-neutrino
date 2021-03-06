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

input_model= options.alertfile

imin = 0
imax = 10000

gam = 2.19

ep = 100.

tobscta = 600.

debug = True
edisp = True

caldb='prod3b-v1'
irf='North_z20_average_30m'

declination,redshift,A = np.loadtxt(input_model, unpack=True)

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
        
        nuseed = randint(1, 1000000000)
            
        sim = ctools.ctobssim()
        sim['inmodel']   = 'nu_sources_'+str(i+1)+'.xml'
        sim['caldb']     = caldb
        sim['irf']       = irf
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
        sim.run()

        like = ctools.ctlike(sim.obs())
        like['debug']     = debug
        like['edisp']     = edisp
        like.run()
            
        nuts = like.obs().models()[sourcename].ts()
        nunormsp = like.obs().models()[sourcename].spectral()['Normalization'].value()
        nunormsp_error = like.obs().models()[sourcename].spectral()['Normalization'].error()
              
        
        real_nu = str(i+1)+' '+str(nuts)+' '+str(nunormsp)+' '+str(nunormsp_error)+' '+str(ra)+' '+str(dec)+' '+str(nuseed)+'\n'
        nusrcts.write(real_nu)
        
nusrcts.close()
