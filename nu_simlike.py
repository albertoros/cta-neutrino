import gammalib
import ctools
import cscripts
import numpy as np
from random import uniform

debug = True

nuid0,ra0,dec0 = np.loadtxt('nu_src_point.dat', unpack=True)
nuid = nuid0.astype(int)

caldb0='prod3b-v1'
irf0='North_z20_average_30m'

imax = np.size(nuid)

for i in xrange(0, imax):
    sim = ctools.ctobssim()
    sim['inmodel']   = 'nu_sources'+str(nuid[i])+'.xml'
    sim['caldb']     = caldb0
    sim['irf']       = irf0
    sim['outevents'] = 'events_nu_'+irf0+'_'+str(nuid[i])+'.fits'
    sim['ra']        = ra0[i]
    sim['dec']       = dec0[i]
    sim['rad']       = 2.0
    sim['tmin']      = '2020-05-31T12:00:00'
    sim['tmax']      = '2020-05-31T12:10:00'
    sim['emin']      = 0.02
    sim['emax']      = 199.0
    sim['maxrate']   = 1.0e9
    sim['debug']     = debug
    sim['edisp']     = True
    sim.execute()

    like = ctools.ctlike()
    like['inobs']     = 'events_nu_'+irf0+'_'+str(nuid[i])+'.fits'
    like['caldb']     = caldb0
    like['irf']       = irf0
    like['edisp']     = True
    like['inmodel']   = 'nu_sources'+str(nuid[i])+'.xml'
    like['outmodel']  = 'nu_'+str(nuid[i])+'_ts.xml'
    like["debug"]     = debug
    like.execute()
