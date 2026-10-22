from fractions import Fraction as F
from math import lcm,log,isclose,gcd
import json
r=lambda x:x-1/x
q=lambda x:min(r(x),F(0))
assert q(F(2))==0 and -q(F(1,2))==F(3,2)
for x in [F(1,4),F(1,2),F(1),F(2),F(3)]: assert r(1/x)==-r(x)
for e in [.1,.01,.0001]: assert isclose(log(2)-log(e)+(1-e*e)/2+log(e),log(2)+(1-e*e)/2)
S=list(range(1,13))+[182]
d=lambda x:min(x%1,1-x%1)
assert min(d(F(v*3,40)) for v in S)==F(3,40)>F(1,14)
# Explicit compressed primitive covering family, seven far runners.
# For ANY integer Q>=14, L=lcm(1,...,Q): S={1,...,6,L,2L,...,7L}.
# All q<=Q divide L, so runner L is at an integer at every a/q.
# Thus these times cannot be lonely. This proves no universal fixed-Q census.
Q=50;L=lcm(*range(1,Q+1));S=list(range(1,7))+[k*L for k in range(1,8)]
assert len(S)==len(set(S))==13 and S[-1]<=13*S[-2]
assert all(L%q==0 for q in range(1,Q+1))
assert all(min(d(F(v*a,q)) for v in S)==0 for q in range(1,Q+1) for a in range(q))
assert all(any(v%q==0 for v in S) for q in range(2,15))
assert all(v>7*6 for v in S[6:])
assert all(all(n%d for d in range(2,int(n**.5)+1)) for n in [7,13,19,31])
# Exact maximum proof: eight circular points 0,Lt,...,7Lt have a gap <=1/8,
# whose pair difference is jLt, 1<=j<=7. Thus M(S)<=1/8 for all real t.
# The following exact witnesses prove the reverse bound in tested instances;
# the symbolic argument in README proves it for all Q>=14.
witness_checks=[]
for Qw in [14,20,50,100]:
    Lw=lcm(*range(1,Qw+1))
    assert Lw%8==0 and Lw>=6
    tw=F(Lw+1,8*Lw)
    assert gcd(Lw+1,8*Lw)==1 and tw.denominator==8*Lw
    fast=[(j*Lw*tw)%1 for j in range(1,8)]
    slow=[i*tw for i in range(1,7)]
    assert fast==[F(j,8) for j in range(1,8)]
    assert all(F(1,8)<=x<=F(7,8) for x in slow)
    assert min(d(x) for x in fast+slow)==F(1,8)>F(1,14)
    witness_checks.append({'Q':Qw,'L':Lw,'witness':str(tw),'reduced_denominator':tw.denominator,'exact_min_distance':'1/8'})
print(json.dumps({'parity_counterexample':{'x':'1/2','q(2)':'0','-q(1/2)':'3/2'},'untruncated_parity':'PASS','cutoff_sum':'PASS','deep_well_at_3/40':'3/40 > 1/14','compressed_lcm_test':{'Q':Q,'L':L,'runner_count':len(S),'primitive':True,'covering_2_to_14':True,'seven_far_relative_to_core_max_6':True,'all_rational_times_denominator_at_most_Q_blocked':True},'exceptional_table_primality_arithmetic_only':'PASS','exact_maximum_witness_checks':witness_checks,'universal_upper_bound':'1/8 by eight-point circular-gap proof'},indent=2))
