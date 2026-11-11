#!/usr/bin/env python3
"""Numerical diagnostic for Proposition U1; no directed interval claim."""
import json
from pathlib import Path
import mpmath as mp

def ratio(L, dps):
    with mp.workdps(dps):
        L=mp.mpf(L)
        lo=mp.exp(-L/2)
        f=lambda u: u**6*mp.besselk(1,u)*(L/(L+mp.log(u)))**2/384
        # Integration is in log coordinate at the lower endpoint.
        g=lambda v: f(mp.exp(v))*mp.exp(v)
        value=mp.quad(g,[-L/2,-4,0,2,mp.log(100)])
        low_bound=L**2*mp.exp(-3*L)/2304
        # For u>=100 and L>0, ratio factor <=1. This positive integral
        # bounds the omitted upper tail (evaluated numerically, not certified).
        upper=mp.quad(lambda u:u**6*mp.besselk(1,u)/384,[100,mp.inf])
        return value,low_bound,upper

def main():
    rows=[]
    for L in [10,20,40,80]:
        a,low,upper=ratio(L,35)
        b,_,_=ratio(L,50)
        assert 0<a<1 and abs(a-b)<mp.mpf('1e-25')
        rows.append({'L':L,'ratio':mp.nstr(b,24),'precision_difference':mp.nstr(abs(a-b),6),
            'omitted_lower_bound':mp.nstr(low,8),'omitted_upper_numerical_bound':mp.nstr(upper,8)})
    assert all(mp.mpf(rows[i]['ratio'])<mp.mpf(rows[i+1]['ratio']) for i in range(len(rows)-1))
    result={'status':'PASS_NUMERICAL_DIAGNOSTIC_NOT_INTERVAL_CERTIFICATE',
        'analytic_bessel_moment':384,'dps_comparison':[35,50],'rows':rows,
        'infrared_controls':{'gap_only':'sigma=0','gapless_continuum':'sigma(ds)=1_(0,1) ds','isolated_atom':'sigma=delta_(1)'}}
    Path(__file__).with_name('OUTPUT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
