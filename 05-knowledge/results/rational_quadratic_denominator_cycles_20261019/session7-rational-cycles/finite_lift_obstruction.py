#!/usr/bin/env python3
"""Verify every period-six residue mod43^2 fails period-six lift mod43^3."""
import json
p=43;q=p*p;r=p**3
inv4q=pow(4,-1,q);inv4r=pow(4,-1,r)
def iterate(x,n,mod):
 inv4=pow(4,-1,mod)
 for _ in range(n):x=(x*x-29)*inv4%mod
 return x
period6=[x for x in range(q) if iterate(x,6,q)==x and all(iterate(x,k,q)!=x for k in (1,2,3))]
assert len(period6)==126
lifts=[x+k*q for x in period6 for k in range(p)]
assert len(lifts)==5418
bad=[x for x in lifts if iterate(x,6,r)==x]
assert not bad
assert all(iterate(x,258,r)==x for x in lifts)
assert all(iterate(x,129,r)!=x and iterate(x,86,r)!=x and iterate(x,6,r)!=x for x in lifts)
print(json.dumps({'map':'T(a)=(a^2-29)/4','p':p,'period6_mod_p2_points':len(period6),'period6_mod_p2_cycles':len(period6)//6,'lifts_mod_p3':len(lifts),'period6_lifts_mod_p3':len(bad),'all_lifts_exact_period':258,'cycles_mod_p3_over_old_period6':len(lifts)//258,'status':'PASS'},indent=2))
