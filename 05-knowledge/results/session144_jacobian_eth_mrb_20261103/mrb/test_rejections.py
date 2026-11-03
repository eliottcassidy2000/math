#!/usr/bin/env python3
from fractions import Fraction
from mpmath import mp
from check import check_mean_formula, check_shift_sign, check_abel_counterterm, require
mp.dps=70
Q=mp.euler**2/2+mp.pi**2/12
cases=[
 ('even unit mask incorrectly converges',lambda:check_mean_formula([1,0],Fraction(0))),
 ('odd mask loses parity sidecar',lambda:check_mean_formula([1,0,0],Fraction(-1,3))),
 ('perturbation sign reversed',lambda:check_shift_sign(1,-1)),
 ('Abel logarithmic counterterm omitted',lambda:check_abel_counterterm(mp.mpf(0),Q)),
 ('Abel finite shift reversed',lambda:check_abel_counterterm(-mp.euler,-Q)),
]
for name,callback in cases:
    try:callback()
    except ValueError:print('REJECTED:',name)
    else:raise RuntimeError('Hostile mutation escaped: '+name)
print('PASS: all five hostile mutations rejected by explicit exceptions.')
