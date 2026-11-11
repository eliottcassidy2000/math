#!/usr/bin/env python3
"""Structural controls only; no A000568 exact-cache queries."""
import json, math
from fractions import Fraction
from pathlib import Path
from predict import parts,core,forecast,tail_log,logweight
import mpmath as mp
ROOT=Path(__file__).resolve().parent

def orbit_count(p,fixed):
    perm=[]; start=0
    for j in p+(1,)*fixed:
        perm.extend(list(range(start+1,start+j))+[start]);start+=j
    remaining={(i,j) for i in range(start) for j in range(i+1,start)}
    count=0
    while remaining:
        pair=next(iter(remaining)); count+=1
        while pair in remaining:
            remaining.remove(pair)
            pair=tuple(sorted((perm[pair[0]],perm[pair[1]])))
    return count

checks=0
for m in range(3,20):
    for p in parts(m):
        mm,r,e,z=core(p)
        for fixed in [0,1,7]:
            n=m+fixed; actual=orbit_count(p,fixed)
            assert actual==e+fixed*r+fixed*(fixed-1)//2
            D=n*(n-1)//2-actual
            assert 3*D>=m*(2*n-m-1)
            checks+=1
# Coefficients in falling-factorial basis for defect levels 1..6.
coeff=[]
for k in range(1,7):
    by_m={}; inventory=[]
    for m in range(2*k+1,3*k+1):
        for p in parts(m):
            mm,r,e,z=core(p)
            if m-r!=2*k:continue
            power=m*(m-r)-m*(m-1)//2+e
            c=Fraction(2**power,z)
            by_m[m]=by_m.get(m,Fraction())+c
            inventory.append({'cycles':p,'falling_degree':m,'coefficient':str(c)})
    assert by_m[3*k]==Fraction(2**(3*k*k+k),3**k*math.factorial(k))
    coeff.append({'defect':2*k,'k':k,'falling_factorial_coefficients':{str(m):str(c) for m,c in by_m.items()},'cycle_inventory':inventory})
(ROOT/'defect_coefficients.json').write_text(json.dumps(coeff,indent=2)+'\n')
# Independent Stirling residual and adaptive envelope consistency.
rows=[]
for n in [30,100,300,1000,10**6,10**100]:
    r=forecast(n)
    assert mp.mpf(r['stirling_vs_loggamma_abs_discrepancy']) <= mp.mpf(r['stirling_log10_absolute_error_bound'])*mp.mpf('1.001') or mp.mpf(r['stirling_vs_loggamma_abs_discrepancy'])<mp.mpf('1e-70')
    rows.append({'n':str(n),'adaptive_cutoff':r['adaptive_support_cutoff'],'achieved':r['adaptive_bound_meets_request']})
report={'status':'PASS','pair_orbit_and_deficit_checks':checks,'defect_levels':6,'numerical_controls':rows,'scope':'structural identities and precision controls; no exact sequence checking; not interval certification'}
(ROOT/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
