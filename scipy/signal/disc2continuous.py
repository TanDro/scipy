"""
Discrete to continuous transformations for state-space and transfer function.
"""
from __future__ import division, print_function, absolute_import

# Author: Jens Wurm <wurm.jens@gmail.com>
# October 01, 2015

import numpy as np
from scipy import linalg 

from .ltisys import tf2ss, ss2tf, zpk2ss, ss2zpk

__all__ = ['disc2continuous']

def disc2continuous(sys, dt, method="zoh", beta=None):
    """
    Transform a discrete to a continuous state-space system.
    Parameters
    ----------
    sys : a tuple describing the system.
        The following gives the number of elements in the tuple and
        the interpretation:
           * 2: (num, den)
           * 3: (zeros, poles, gain)
           * 4: (A, B, C, D)
    dt : float
        The discretization time step.
    method : {"impulse", "matched", "tustin", "foh", "zoh"}, optional
        Which method to use:
           * tustin: bilinear tustin method
           * zoh: zero-order hold (default)
    Returns
    -------
    sysd : tuple containing the continuous system
        Based on the input type, the output will be of the form
        * (num, den)   for transfer function input
        * (zeros, poles, gain)   for zeros-poles-gain input
        * (A, B, C, D) for state-space system input
    Notes
    -----
    By default, the routine uses a Zero-Order Hold (zoh) method to perform
    the transformation. 
	
	"""
    if len(sys) == 2:
        sysd = disc2continuous(tf2ss(sys[0], sys[1]), dt, method=method)
        return ss2tf(sysd[0], sysd[1], sysd[2], sysd[3]) + (dt,)
    elif len(sys) == 3:
        sysd = disc2continuous(zpk2ss(sys[0], sys[1], sys[2]), dt, method=method)
        return ss2zpk(sysd[0], sysd[1], sysd[2], sysd[3]) + (dt,)
    elif len(sys) == 4:
        a, b, c, d = sys
    else:
		raise ValueError("First argument must either be a tuple of 2 (tf), "
                         "3 (zpk), or 4 (ss) arrays.")

    if method == 'bilinear' or method == 'tustin':
        sigma = 2 / dt
        ac = linalg.solve(np.eye(Ad.shape[0]) + Ad, sigma * (Ad - np.eye(Ad.shape[0])))
        bc = linalg.solve(
    elif method == 'zoh':
		# Build an logarithmic matrix
		lm_upper = np.hstack((a, b))

		# Need to stack zeros under the a and b matrices
		lm_lower = np.hstack((np.zeros((b.shape[1], a.shape[0])),
					  np.ones((b.shape[1], b.shape[1])) ))

		lm = np.vstack((lm_upper, lm_lower))
		ms = (1./dt) * linalg.logm(lm)

		# Dispose of the lower rows
			ms = ms[:a.shape[0], :]

		ac = ms[:, 0:a.shape[1]]
		bc = ms[:, a.shape[1]:]

		cc = c
		dc = d
		
    else:
		raise ValueError("Unknown transformation method '%s'" % method)

    return ac, bc, cc, dc
