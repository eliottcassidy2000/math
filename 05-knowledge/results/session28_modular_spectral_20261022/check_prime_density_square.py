"""Exact square certificate for decimal prime-density moment obstruction.
Counts are externally cited, not all recomputed here. No RH assumption.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
pi=[0,4,25,168,1229,9592,78498,664579,5761455,50847534,455052511,4118054813]
d=[Q(pi[k]-pi[k-1],9*10**(k-1)) for k in range(1,len(pi))]
f=[0,-1,3,-1,-3,2]
sq=[sum(f[i]*f[n-i] for i in range(len(f)) if 0<=n-i<len(f)) for n in range(11)]
v=sum(a*b for a,b in zip(sq,d))
assert sq==[0,0,1,-6,11,0,-21,18,5,-12,4]
assert v==Q(-2447633,22500000000)
# Positive control, uniform measure on [0,1].
control=sum(Q(a,n+1) for n,a in enumerate(sq))
assert control>0
out={'prime_pi_at_powers10':pi,'densities':list(map(str,d)), 'moment_index':'m_n=delta_(n+1), n>=0','polynomial':'x*(x-1)^2*(x+1)*(2*x-1)','coefficients':f,'square_coefficients':sq,'negative_square_value':str(v),'uniform_measure_positive_control':str(control),'uses_only_density_indices':[n+1 for n,a in enumerate(sq) if a], 'scope':'No positive scalar Hamburger moment representation of the full native shifted density sequence; not a claim about RH or transformed observables.'}
print(json.dumps(out,indent=2));Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
