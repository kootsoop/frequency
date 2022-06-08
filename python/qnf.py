#
# QNF - The Quinn - Fernandes frequency estimator.
#
#          Inputs:  signal - T x N matrix where
#                            T = data length
#                            N = number of signals
#                            (i.e. N signals in columns).
#
#          Outputs: est - N Quinn-Fernandes frequency estimates.
#
# [1]  B.G. Quinn & J.M. Fernandes, "A fast technique for the estimation of frequency,"
#      Biometrika, Vol. 78(3), pp. 489--497, 1991.

# $Id: qnf.m 1.1 2000/06/07 18:57:16 PeterK Exp PeterK $

# File: qnf.m
#
# Copyright (C) 1993 CRC for Robust & Adaptive Systems
# 
# This software provided under the GNU General Public Licence, which
# is available from the website from which this file originated. For
# a copy, write to the Free Software Foundation, Inc., 675 Mass Ave, 
# Cambridge, MA 02139, USA.

import numpy as np
import scipy.signal as signal

def qnf(sig):
    if type(sig) is not np.ndarray:
        raise Exception("signal must be an np.ndarray")
        
    #
    # Initializations
    #
    shape=sig.shape;
    t = shape[0]
    if t < 4:
        raise Exception("signal must be longer than 1 point")
        
    if (len(shape)==1):
        ns = 1
        sig = np.array([sig]).T
    else:
        ns = shape[1]
        
    xb = np.mean(sig, axis=0)
    
    if (len(xb.shape) == 1):
        xbm = np.multiply(np.ones([t,1]),xb)
    else:
        xbm = np.matmul(np.ones([t,1]),xb)

    sig=np.subtract(sig,xbm)
    t3 = t+1
    
    y = np.fft.fft(sig, n=2*t, axis=0)
    
    z=np.multiply(y,np.conj(y))
    z=z[2:t3,]

    [m,j]= z[2:t-1,].max(axis=0), z[2:t-1,].argmax(axis=0) 

    j=j+1 # TODO: Not needed because of zero indexing?

    a=2*np.cos(np.pi*j/t)
    y=y[1:2*t:2]
    
    #
    # Quinn-Fernandes method
    #
    b=[1];
    nm=t-1;
    for jjj in [1,2,3,4,5,6,7]:
        for q in np.arange(ns):
            c=[1,-a[q],1]
            y[:,q] = signal.lfilter(b,c,sig[:,q])

        v = np.sum(np.divide(np.multiply(sig[2:t,],y[1:nm,]),np.sum(np.multiply(y[1:nm,],y[1:nm,]))))
        a = np.add(a,2*v);

    return np.real(np.arccos(a/2))

    # Author: SJS 1992; Adapted from code within ttinpie.m (author PJK)
    #
    # Based on: P.J. Kootsookos, S.J. Searle and B.G. Quinn, 
    # "Frequency Estimation Algorithms," CRC for Robust and 
    # Adaptive Systems Internal Report, June 1993.

