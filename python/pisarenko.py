#
# PISARENKO Find the Pisarenko frequency estimates of the signals 
#          in each column of signal.
#
#          omegahat = pisarenko(signal)
#
# [1] V. F.  Pisarenko, "On the Estimation of Spectra by 
#     Means of Non-linear Functions of the Covariance Matrix,"
#     Geophys. J. Roy. astr. Soc., Vol. 28, pp. 511-531, 1972. 
#
# [2] V. F. Pisarenko, ``The Retrieval of Harmonics from a 
#     Covariance Function," Geophys. J. Roy. astr. Soc., 
#     Vol.   33, pp. 347-366, 1973. 
#
# $Id: pisarenko.m,v 1.2 1999/01/09 11:10:35 PeterK Exp PeterK $

# File: pisarenko.m
#
# Copyright (C) 1999 Peter J. Kootsookos
# 
# This software provided under the GNU General Public Licence, which
# is available from the website from which this file originated. For
# a copy, write to the Free Software Foundation, Inc., 675 Mass Ave, 
# Cambridge, MA 02139, USA.

import numpy as np
import scipy.signal as signal

def pisarenko(sig):

    if type(sig) is not np.ndarray:
        raise Exception("signal must be an np.ndarray")
        
    #
    # Initializations
    #
    eps = 1
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
      
    Rss1 = np.sum(sig[1:t-1,:]*sig[0:t-2,:], axis=0)
    Rss2 = np.sum(sig[2:t-1,:]*sig[0:t-3,:], axis=0)
    
    alpha = ( Rss2 + np.sqrt(np.power(Rss2,2) + 8*np.power(Rss1,2)) ) / ( 4*Rss1 + eps ) ;
    
    if np.any(np.abs(alpha) > 1):
        alpha[np.any(np.abs(alpha) > 1)] = 1

    return np.real(np.arccos(alpha))

# Author: Peter J. Kootsookos (p.kootsookos@ieee.org)
#
# Based on: P.J. Kootsookos, S.J. Searle and B.G. Quinn, 
# "Frequency Estimation Algorithms," CRC for Robust and 
# Adaptive Systems Internal Report, June 1993.
#