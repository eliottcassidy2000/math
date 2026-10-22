#!/usr/bin/env python3
"""Independent exact-rational and interval checks for the global theorem."""
from fractions import Fraction as Q
from pathlib import Path
import json
from mpmath import iv
iv.dps=100

# e>sum_{k=0}^3 1/k!=8/3 (omitted tail is strictly positive).
assert Q(1)+Q(1)+Q(1,2)+Q(1,6)==Q(8,3)
# For 0<t<1, exp(t)<1+t+t^2/(2*(1-t)), since n!>=2 for n>=2,
# with strict inequality already at n=3.
t=Q(3,8)
upper=1+t+t*t/(2*(1-t))
assert upper==Q(119,80)<Q(3,2)
assert Q(5,3)**2<3
# sqrt(3)>5/3, pi<22/7 imply sinc(pi/3)>35/44.
sinc_lower=3*Q(5,3)/(2*Q(22,7))
assert sinc_lower==Q(35,44)
c_lower=2*Q(2,3)*sinc_lower
assert c_lower==Q(35,33)>1
eta_upper=Q(4,9)/c_lower
assert eta_upper==Q(44,105)<1
# For every integer m>=3, 35m^2-66m-66>0: value at 3 is 51,
# and its forward difference is 70m-31>=179>0.
assert 35*3**2-66*3-66==51>0
assert 70*3-31==179>0
# Therefore 2(m+1)<(35/33)m^2<c*m^2 for all such m.

c=3*iv.sqrt(3)*iv.exp(-1/iv.exp(iv.mpf(1)))/iv.pi
eta=4/(9*c)
assert bool(c.a>iv.mpf(35)/33)
assert bool(eta.b<iv.mpf(44)/105)
assert bool(iv.pi.b<iv.mpf(22)/7)
assert bool(iv.exp(iv.mpf(3)/8).b<iv.mpf(119)/80)

result={'status':'PASS','c_interval':str(c),'eta_interval':str(eta),
 'exact_c_lower_bound':str(c_lower),'exact_eta_upper_bound':str(eta_upper),
 'exp_3_over_8_upper_bound':str(upper),
 'all_integer_growth_check':{'polynomial':'35*m^2-66*m-66','value_at_3':51,
                            'forward_difference':'70*m-31','minimum_difference_for_m_ge_3':179},
 'edge_cases':{'ties':'Either minimizer gives the same absolute residual; the recursive sign records the choice.',
               'exact_termination':'Stop at zero. The finite signed identity holds, and no next index is asserted.'}}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
