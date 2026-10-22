#!/usr/bin/env python3
"""Exact rational certificate at the transcendental mirror point log(3)/3.
Endpoint-model checks are numerical diagnostics, not error certificates.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import json, time
from verify_tree_moments import a, dec, mpq, pade_gauss, val, mp

HERE=Path(__file__).resolve().parent
mp.mp.dps=100
K=100
half=F(1,2)
log_lower=sum((2*half**(2*j+1)/(2*j+1) for j in range(K)),F(0))
log_tail=2*half**(2*K+1)/((2*K+1)*(1-half*half))
log_upper=log_lower+log_tail
zl,zu=log_lower/3,log_upper/3
b=F(2719,1000)
e_upper=sum((F(1,factorial(j)) for j in range(11)),F(0))+F(1,10*factorial(10))
assert e_upper<b and b*zu<1
actual_z=mp.log(3)/3
actual_R=-mp.lambertw(-actual_z)/actual_z
assert mpq(zl)<actual_z<mpq(zu)
cert={'status':'PASS','input':'z=log(3)/3',
      'certification_scope':'All enclosure endpoints are computed with exact fractions, then rounded outward. log(3) is enclosed by its positive atanh series and a proved geometric tail. Lambert W is used only as a numerical diagnostic.',
      'log3_terms':K,'log3_lower_exact':str(log_lower),'log3_upper_exact':str(log_upper),
      'log3_tail_exact':str(log_tail),'b':str(b),'e_upper':str(e_upper),
      'b_times_z_upper_approx':mp.nstr(mpq(b*zu),50),
      'input_width_approx':mp.nstr(mpq(zu-zl),12)}
rows=[]
for n in [4,8,16,32,64]:
    started=time.monotonic()
    moments=[a(j) for j in range(2*n+1)]
    transformed=[b*moments[j]-moments[j+1] for j in range(2*n)]
    lower=val(pade_gauss(moments,n),zl)
    upper=(1-zu*val(pade_gauss(transformed,n),zu))/(1-b*zu)
    assert lower<upper
    # Diagnostic only; certified directions are proved independently.
    assert mpq(lower)<actual_R<mpq(upper)
    raw_lower=sum((a(j)*zl**j for j in range(2*n)),F(0))
    raw_upper=sum((a(j)*zu**j for j in range(2*n)),F(0))+a(2*n)*zu**(2*n)/(1-e_upper*zu)
    assert raw_lower<lower<upper<raw_upper
    lower_out=dec(lower,60)
    upper_out=dec(upper,60,upper=True)
    # Independently check the emitted finite decimals are actually outward.
    assert F(lower_out)<=lower and F(upper_out)>=upper
    row={'N':n,'lower_decimal_certified':lower_out,'upper_decimal_certified':upper_out,
         'bracket_width_approx':mp.nstr(mpq(upper-lower),16),
         'raw_series_degree':2*n-1,
         'raw_series_width_approx':mp.nstr(mpq(raw_upper-raw_lower),16),
         'width_improvement_factor_approx':mp.nstr(mpq((raw_upper-raw_lower)/(upper-lower)),16),
         'moment_indices_used':f'0..{2*n}'}
    rows.append(row)
    print(f"N={n:2d}: exact-rational width={row['bracket_width_approx']}; raw degree {2*n-1} width={row['raw_series_width_approx']}; factor={row['width_improvement_factor_approx']}; seconds={time.monotonic()-started:.2f}",flush=True)
cert['mirror_certificates']=rows

# Positive but nonnormalized endpoint reference measure on (0,1).
def M0(p):
    return 1/p-mp.exp(p)*mp.e1(p)

def M(p):
    return mp.exp((p-1)*mp.log(p))/mp.gamma(1+p)

checks=[]
for k in [10,100,1000,10000,1000000]:
    p=mp.mpf(1)/k
    model=1/(p*p*mp.gamma(p)*M0(p))
    truth=mp.exp(mp.log(k)/k)
    error=model-truth
    assert error>0
    checks.append({'k':k,'M0':mp.nstr(M0(p),40),
                   'kernel_error':mp.nstr(error,20),
                   'error_scaled_by_k2_over_log2k':mp.nstr(error*k*k/mp.log(k)**2,20),
                   'one_node_bad_kernel_over_k':mp.nstr(k/mp.gamma(p),20)})
cert['endpoint_model_diagnostics']=checks
cert['reference_value_diagnostic']=mp.nstr(actual_R,80)
(HERE/'mirror_endpoint_verification.json').write_text(json.dumps(cert,indent=2)+'\n')
print('PASS: directed rational mirror-input certificate and endpoint diagnostics.',flush=True)
print('R(log(3)/3), numerical diagnostic:',mp.nstr(actual_R,80),flush=True)
