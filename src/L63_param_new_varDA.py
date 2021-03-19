#!/usr/bin/env python

###############################################################
# < next few lines under version control, D O  N O T  E D I T >
# $Date$
# $Revision$
# $Author$
# $Id$
###############################################################

###############################################################
# L63_param_varDA.py - parameters for variational DA on L63
###############################################################

###############################################################
__author__    = "Rahul Mahajan"
__email__     = "rahul.mahajan@nasa.gov"
__copyright__ = "Copyright 2012, NASA / GSFC / GMAO"
__license__   = "GPL"
__status__    = "Prototype"
###############################################################

###############################################################
import numpy as np
from module_Lorenz import Lorenz
from module_DA import DataAssim, VarDataAssim
from module_IO import Container
###############################################################

# insure the same sequence of random numbers EVERY TIME
np.random.seed(0)

model = Lorenz()
DA    = DataAssim()
varDA = VarDataAssim()

# Initialize Lorenz model class
Name = 'L63'                      # model name
Ndof = 3                          # model degrees of freedom
Par  = [10.0, 28.0, 8.0/3.0]      # model parameters F, dF
dt   = 1.0e-2                     # model time-step
model.init(Name=Name,Ndof=Ndof,Par=Par,dt=dt)

nassim = 200                 # no. of assimilation cycles
ntimes = 0.25                # do assimilation every ntimes non-dimensional time units
maxouter = 1
DA.init(nassim=nassim,ntimes=ntimes,maxouter=maxouter)

Q = np.ones(model.Ndof)         # model error covariance ( covariance model is white for now )
Q = np.diag(Q) * 0.0

H = np.ones(model.Ndof)         # obs operator ( eye(Ndof) gives identity obs )
H = np.diag(H)

R = np.ones(model.Ndof)         # observation error covariance
R = 2 * R
R = np.diag(R)

update                  = 1              # DA method (1= 3Dvar; 2= 4Dvar)
precondition            = 1              # precondition before minimization (0= None; 1= sqrtB; 2= FullB)
maxouter                = 1              # no. of outer loops
maxiter                 = 1000           # maximum iterations        
tol                     = 1e-4           # tolerance to end the variational minimization iteration
inflate                 = True           # inflate [ > 1.0 ] / deflate [ < 1.0 ] static covariance
infl_fac                = 1.85           # inflate static covariance
infl_adp                = True           # inflate adaptively (cuts inflation as a function  of OL)
localize                = 1              # localization (0= None, 1= Gaspari-Cohn, 2= Boxcar, 3= Ramped)
cov_cutoff              = 0.0625         # normalized covariance cutoff = cutoff / ( 2*normalized_dist )
cov_trunc               = model.Ndof     # truncate localization matrix (cov_trunc <= model.Ndof)

if   ( update == 1 ):
    window              = 0.00           # length of the assimilation window
    offset              = 1.00           # time offset: forecast from analysis to background time
    nobstimes           = 1              # no. of evenly spaced obs. times in the window
elif ( update == 2 ):
    window              = DA.ntimes      # length of the 4Dvar assimilation window
    offset              = 0.5            # time offset: forecast from analysis to background time
    nobstimes           = 5              # no. of evenly spaced obs. times in the window
varDA.init(model,DA,\
           update=update,precondition=precondition,\
           maxouter=maxouter,tol=tol,\
           inflate=inflate,infl_fac=infl_fac,infl_adp=infl_adp,\
           localize=localize,cov_cutoff=cov_cutoff,cov_trunc=cov_trunc,\
           window=window,offset=offset,nobstimes=nobstimes)


    
filename   = model.Name + '_varDA_diag.nc4'
attributes = {'model'       : model.Name,
              'sigma'       : model.Par[0],
              'rho'         : model.Par[1],
              'beta'        : model.Par[2],
              'ntimes'      : DA.ntimes,
              'dt'          : model.dt,
              'Vupdate'     : varDA.update,
              'maxouter'    : varDA.maxouter,
              'precondition': varDA.precondition,
              'Vlocalize'   : varDA.localization.localize,
              'Vcov_cutoff' : varDA.localization.cov_cutoff,
              'Vcov_trunc'  : varDA.localization.cov_trunc,
              'maxiter'     : varDA.minimization.maxiter,
              'tol'         : varDA.minimization.tol}
diag_file = Container(filename=filename,attributes=attributes)
        
if ( varDA.update == 2 ):
        diag_file.attributes.update({'offset'    : varDA.fdvar.offset,
                                     'window'    : varDA.fdvar.window,
                                     'nobstimes' : varDA.fdvar.nobstimes})

# restart conditions
time     = None              # None == default | -N...-1 0 1...N
filename = ""
restart = Container(time=time,filename=filename)
