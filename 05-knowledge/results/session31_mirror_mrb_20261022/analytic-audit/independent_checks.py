#!/usr/bin/env python3
"""Independent hostile audit: direct-binomial Euler weights, radius 2x/3.

Does not import, modify, or execute the submitted source scripts.
Only mpmath and the Python standard library are needed.
"""
from pathlib import Path
from math import comb
import json, re
from mpmath import iv, mp

mp.dps = 170
iv.dps = 180
BASE = Path(__file__).resolve().parents[1]
M, N = 384, 128
alpha = iv.mpf(2)/3
L = iv.log((1+alpha)*M)+iv.pi/4
b = L/((1-alpha)*M)
base_bound = iv.factorial(N)/(1-alpha)/(2*alpha)**N * (
    L/iv.mpf(M)**(N+1)+(L/N+1/iv.mpf(N)**2)/iv.mpf(M)**N)
# Flatten the entire forward-difference triangle into exact integer weights.
weights = [sum(comb(j,k)*2**(N-j-1) for j in range(k,N)) for k in range(N)]
assert all(0 < w < 2**N for w in weights)

def h(n):
    n = iv.mpf(n)
    return iv.log(n)/n

def amplitude(n):
    n = iv.mpf(n)
    return 2/iv.pi * iv.exp(-h(n))*iv.sin(iv.pi/n)

def enclosure(t=iv.mpf(1), derivative=False, abs_t_bound=iv.mpf(1)):
    def value(n):
        hn = h(n)
        return hn*iv.exp(t*hn) if derivative else iv.exp(t*hn)-1
    partial = sum(((-1)**n*value(n) for n in range(1,M)),iv.mpf(0))
    smooth_tail = sum(((-1)**k*weights[k]*value(M+k) for k in range(N)),iv.mpf(0))/2**N
    remainder = base_bound*iv.exp(abs_t_bound*b)*(1 if derivative else abs_t_bound)
    return partial+(-1)**M*smooth_tail+iv.mpf([-1,1])*remainder, remainder

def striv(x): return str(x)

def midpoint(x): return mp.mpf(x._mpi_[0])+ (mp.mpf(x._mpi_[1])-mp.mpf(x._mpi_[0]))/2

C,B = enclosure()
assert bool(B.b < iv.mpf('1.6e-132'))
recorded_mrb = BASE/'mrb/certify_mrb.out'
if not recorded_mrb.is_file():
    recorded_mrb = BASE/'session31-mrb/certify_mrb.out'
text = recorded_mrb.read_text()
old = re.search(r'^C = \[([^,]+), ([^\]]+)\]',text,re.M)
old_lo,old_hi=iv.mpf(old[1]),iv.mpf(old[2])
assert bool(C.a > old_lo.b) and bool(C.b < old_hi.a)
D = amplitude(8)-C
q4 = iv.exp(-4*iv.pi*iv.sqrt(2))
assert bool(D.a > 0) and bool((q4-D).a > 0)
# theta_3+theta_4-2 = 4*sum_{j>=1}q^(4j^2) proves theta tail>D.

indices = [8,105719004,77401040908104410,24845630402060641702669818431121505]
rows=[]
r=C
sgn=1
reconstruction=iv.mpf(0)
for i,k in enumerate(indices):
    a=amplitude(k)
    lower_margin=r-(a+amplitude(k+1))/2
    upper_margin=(amplitude(k-1)+a)/2-r
    assert bool(lower_margin.a>0) and bool(upper_margin.a>0)
    signed=r-a
    sign=-1 if i==0 else 1
    assert bool((-signed if sign<0 else signed).a>0)
    # Bracket continuous inverse independently, without root finding.
    lo=k if sign<0 else k-1
    assert bool((amplitude(lo)-r).a>0)
    assert bool((r-amplitude(lo+1)).a>0)
    reconstruction += sgn*a
    sgn *= sign
    r = -signed if sign<0 else signed
    assert bool((1/iv.mpf(lo)**2-r).a>0)
    rows.append({'step':i+1,'index':str(k),'continuous_inverse_bracket':[str(lo),str(lo+1)],
                 'lower_nearest_midpoint_margin':striv(lower_margin),
                 'upper_nearest_midpoint_margin':striv(upper_margin),
                 'signed_residual':striv(signed),'next_absolute_residual':striv(r)})
# Interval arithmetic contains zero in the exact signed reconstruction identity.
recon_error=C-reconstruction-sgn*r
assert bool(recon_error.a<=0) and bool(recon_error.b>=0)

# Certify deformation sign change independently of both nsum and quadrature.
tlo=iv.mpf('0.999999'); thi=iv.mpf('1.000001')
def A8t(t): return 2/iv.pi*iv.exp(-t*h(8))*iv.sin(iv.pi*t/8)
def A8dt(t):
    return 2/iv.pi*iv.exp(-t*h(8))*(iv.pi/8*iv.cos(iv.pi*t/8)-h(8)*iv.sin(iv.pi*t/8))
Cl,_=enclosure(tlo,abs_t_bound=thi)
Ch,_=enclosure(thi,abs_t_bound=thi)
Dl=A8t(tlo)-Cl; Dh=A8t(thi)-Ch
assert bool(Dl.a>0) and bool(Dh.b<0)
tinterval=iv.mpf(['0.999999','1.000001'])
Cd,_=enclosure(tinterval,derivative=True,abs_t_bound=thi)
Dd=A8dt(tinterval)-Cd
assert bool(Dd.b<0)
Cd1,_=enclosure(derivative=True)
Dd1=A8dt(iv.mpf(1))-Cd1
assert bool(Dd1.a<iv.mpf('-0.0896013973074147415174576834930579928139735691707304'))
assert bool(Dd1.b>iv.mpf('-0.0896013973074147415174576834930579928139735691707305'))
# Tight enclosure of the numerical zero: two exact-decimal endpoint certificates.
zlo=iv.mpf('1.000000211135851595066050802731745766168538790567401655840')
zhi=iv.mpf('1.000000211135851595066050802731745766168538790567401655842')
Zcl,_=enclosure(zlo,abs_t_bound=thi); Zch,_=enclosure(zhi,abs_t_bound=thi)
Zdl=A8t(zlo)-Zcl; Zdh=A8t(zhi)-Zch
assert bool(Zdl.a>0) and bool(Zdh.b<0)

# Numerical adversarial controls use formulas derived in the audit, not source code.
A=lambda x: 2/mp.pi*mp.exp(-mp.log(x)/x)*mp.sin(mp.pi/x)
Cmid=midpoint(C)
residual_scale=[]
for x in [mp.mpf('100'),mp.mpf('10000'),mp.mpf('100000000')]:
    rr=A(x); X=2/rr; lx=mp.log(X)
    approx=X-lx+(lx-lx*lx/2-mp.pi**2/6)/X
    residual_scale.append({'x':str(x),'scaled_inverse_remainder':mp.nstr((x-approx)*X**2/lx**3,40)})
# Nearest amplitude really is not nearest x: reciprocal-like convexity shifts cutoff.
x=mp.mpf('100.499'); rr=A(x)
assert abs(rr-A(101))<abs(rr-A(100)) and mp.nint(x)==100
# Direct finite sums with explicit geometric-polynomial truncation independent of AP.
perturbations=[]
for K in range(1,9):
    exact=(-1)**(K+1)*mp.factorial(K)/(mp.e+1)**(K+1)
    direct=mp.fsum((-1)**n*mp.exp(-n)*mp.fprod(n-j for j in range(1,K+1)) for n in range(K+1,1000))
    assert abs(exact-direct)<mp.mpf('1e-160')
    perturbations.append({'K':K,'exact_shift':mp.nstr(exact,55),'absolute_difference':mp.nstr(abs(exact-direct),8)})

out={'status':'PASS: independent analytic derivations and interval certificates',
     'method':{'M':M,'N':N,'cauchy_radius':'2*x/3','iv_decimal_precision':iv.dps,
               'finite_sum':'direct exact-integer binomial weights; no forward-difference recursion',
               'remainder':striv(B)},
     'C':striv(C),'defect':striv(D),'q4_minus_defect':striv(q4-D),
     'nested_strictly_inside_submitted_C_interval':True,
     'certified_nearest_indices':rows,
     'deformation':{'D8_0.999999':striv(Dl),'D8_1.000001':striv(Dh),
                    'D8_prime_at_1':striv(Dd1),
                    'D8_prime_over_full_bracket':striv(Dd),
                    'unique_real_zero_bracket':[str(zlo),str(zhi)],
                    'zero_lower_endpoint_sign':striv(Zdl),'zero_upper_endpoint_sign':striv(Zdh)},
     'inverse_asymptotic_numerical_controls':residual_scale,
     'nearest_not_rounding_example':{'x':'100.499','rounded_x':100,'nearest_amplitude_index':101},
     'perturbation_numerical_controls':perturbations}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS: C interval is strictly nested inside submitted interval.')
print('Independent remainder bound:',B)
print('PASS: all four greedy indices are now rigorously interval-certified.')
print('PASS: D8 has exactly one real zero in [0.999999,1.000001].')
print('Certified zero bracket: [1.000000211135851595066050802731745766168538790567401655840,')
print('                         1.000000211135851595066050802731745766168538790567401655842].')
print('PASS: exact deformation, branch-shift, theta sign, inverse and rounding controls.')
