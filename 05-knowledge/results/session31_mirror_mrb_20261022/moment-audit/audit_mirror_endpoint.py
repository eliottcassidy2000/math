#!/usr/bin/env python3
"""Independent mirror enclosure proof using only rational logarithm bounds.
This does not reuse the author's Padé construction or numerical Lambert W.
Endpoint-model integrations remain separately identified diagnostics.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
from hashlib import sha256
import json, sys
import mpmath as mp

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'tree-moments'
if not SOURCE.is_dir():
    SOURCE=HERE.parent/'session31-tree-moments'
mp.mp.dps=100
certpath=Path(sys.argv[1]) if len(sys.argv)>1 else SOURCE/'mirror_endpoint_verification.json'
j=json.loads(certpath.read_text())

def log_bounds(q,K=100):
    assert q>1
    t=(q-1)/(q+1)
    lo=2*sum((t**(2*r+1)/(2*r+1) for r in range(K)),F(0))
    tail=2*t**(2*K+1)/((2*K+1)*(1-t*t))
    return lo,lo+tail

ll,lu=log_bounds(F(3))
assert ll==F(j['log3_lower_exact'])
assert lu==F(j['log3_upper_exact'])
assert lu-ll==F(j['log3_tail_exact'])
zl,zu=ll/3,lu/3
b=F(2719,1000)
elow=sum((F(1,factorial(k)) for k in range(11)),F(0))
eupper=elow+F(1,10*factorial(10))
assert b*zu<1 and eupper<b
report={'status':'PASS','method':'Independent exact rational atanh bounds plus monotonicity of log(r)/r on (1,e); no Lambert-W numerics and no Padé solver needed.',
        'source_sha256':{p.name:sha256(p.read_bytes()).hexdigest() for p in SOURCE.iterdir() if p.is_file()},
        'mirror_endpoint_checks':[]}
for row in j['mirror_certificates']:
    lower=F(row['lower_decimal_certified']); upper=F(row['upper_decimal_certified'])
    assert 1<lower<elow
    log_lower_lo,log_lower_hi=log_bounds(lower)
    assert log_lower_hi/lower<zl
    if upper>eupper:
        upper_method='upper endpoint exceeds certified upper bound for e'
    else:
        assert upper<elow
        log_upper_lo,log_upper_hi=log_bounds(upper)
        assert log_upper_lo/upper>zu
        upper_method='log(upper)/upper is greater than the certified upper bound for log(3)/3'
    report['mirror_endpoint_checks'].append({'N':row['N'],'lower_decimal_certified':row['lower_decimal_certified'],'upper_decimal_certified':row['upper_decimal_certified'],'upper_method':upper_method,'independent_strict_rational_proof':'PASS'})

report['endpoint_model_diagnostics']=[]
for p in [mp.mpf('.01'),mp.mpf('.1'),mp.mpf('.5'),mp.mpf(1),mp.mpf(2)]:
    integrated=mp.quad(lambda t:mp.exp(-p*t)*t/(1+t),[0,1,10,100,mp.inf])
    closed=1/p-mp.exp(p)*mp.e1(p)
    assert closed>0 and abs(closed-integrated)<mp.mpf('1e-85')
    report['endpoint_model_diagnostics'].append({'p':str(p),'integral_absolute_error':mp.nstr(abs(closed-integrated),8),'positive':True})
report['endpoint_asymptotic_diagnostics']=[]
for k in [10**3,10**6,10**12,10**20]:
    k=mp.mpf(k); p=1/k
    h=1/p-mp.exp(p)*mp.e1(p)
    truth=mp.exp(mp.log(k)/k)
    model=1/(p*p*mp.gamma(p)*h)
    scaled=(truth-model)*k*k/mp.log(k)**2
    assert scaled<0
    assert abs(scaled+mp.mpf('.5'))<1/mp.log(k)
    report['endpoint_asymptotic_diagnostics'].append({'k':str(k),'signed_error_times_k2_over_log2k':mp.nstr(scaled,40)})
(HERE/'independent_mirror_endpoint.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS: every published mirror interval independently certified using exact rational logarithms.')
print('PASS: continuous endpoint-model integral and signed leading asymptotic diagnostics.')
print('Maximum N:',report['mirror_endpoint_checks'][-1]['N'])
print('Certified final mirror interval:')
print(report['mirror_endpoint_checks'][-1]['lower_decimal_certified'])
print(report['mirror_endpoint_checks'][-1]['upper_decimal_certified'])
