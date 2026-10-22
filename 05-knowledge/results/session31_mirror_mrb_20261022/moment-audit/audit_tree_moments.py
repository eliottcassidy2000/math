#!/usr/bin/env python3
"""Independent Session31 analytic-packet replay.
Exact rational identities certify the submitted rational arithmetic conditional
on the audited analytic proof. mpmath integrations are diagnostics only.
Does not import the author's program or mutate Git/repository files.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
from hashlib import sha256
import json
import mpmath as mp

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'tree-moments'
if not SOURCE.is_dir():
    SOURCE=HERE.parent/'session31-tree-moments'
mp.mp.dps=90

def a(n):
    return F((n+1)**n, factorial(n+1))
def ev(p,z):
    return sum((c*z**j for j,c in enumerate(p)),F(0))
def mq(q):
    return mp.mpf(q.numerator)/q.denominator

author=json.loads((SOURCE/'verification.json').read_text())
report={'status':'PASS','scope':'Independent exact arithmetic replay plus separately nonrigorous numerical diagnostics; see INDEPENDENT_MOMENT_AUDIT.md for analytic proof audit.',
        'source_sha256':{p.name:sha256(p.read_bytes()).hexdigest() for p in SOURCE.iterdir() if p.is_file()},
        'exact_pade_replays':[]}
b=F(2719,1000)
ecap=sum((F(1,factorial(j)) for j in range(11)),F(0))+F(1,10*factorial(10))
assert ecap==F(4697191,1728000)<b<3
z=F(1,3)
for row in author['gauss_rational_enclosures_R_one_third']:
    n=row['N']
    for suffix in ['', '_b']:
        p=list(map(F,row['P'+suffix])); q=list(map(F,row['Q'+suffix]))
        moments=[a(k) if not suffix else b*a(k)-a(k+1) for k in range(2*n)]
        assert len(p)==n and len(q)==n+1 and q[0]==1
        # q is the reversal of the monic x-variable orthogonal polynomial.
        orth=q[::-1]
        for k in range(n):
            assert sum((orth[j]*moments[k+j] for j in range(n+1)),F(0))==0
        for k in range(2*n):
            conv=sum((q[j]*moments[k-j] for j in range(min(k,n)+1)),F(0))
            assert conv==(p[k] if k<n else 0)
        assert ev(q,z)>0
    p,q=map(lambda key:list(map(F,row[key])),['P','Q'])
    pb,qb=map(lambda key:list(map(F,row[key])),['P_b','Q_b'])
    lo=ev(p,z)/ev(q,z)
    hi=(1-z*ev(pb,z)/ev(qb,z))/(1-b*z)
    assert lo==F(row['lower_exact'])<F(row['upper_exact'])==hi
    assert hi-lo==F(row['width_exact'])
    assert F(row['lower_decimal_outward'])<=lo
    assert F(row['upper_decimal_outward'])>=hi
    assert lo-F(row['lower_decimal_outward'])<F(1,10**80)
    assert F(row['upper_decimal_outward'])-hi<F(1,10**80)
    series=sum((a(k)*z**k for k in range(2*n)),F(0))
    tail=a(2*n)*z**(2*n)/(1-ecap*z)
    assert series<=lo and hi<series+tail
    report['exact_pade_replays'].append({'nodes':n,'orthogonality_and_series_match':'PASS','decimal_outward_rounding':'PASS','width':str(hi-lo),'width_diagnostic':mp.nstr(mq(hi-lo),25)})

# Exact denominator hostile: z is below 1/e, since z*ecap<1,
# but above 1/b. The formula advertised as an upper bound becomes negative.
zh=F(1839,5000)
assert zh*ecap<1<zh*b
c0=b-1; c1=b-F(3,2)
gb=c0*c0/(c0-zh*c1)
uh=(1-zh*gb)/(1-b*zh)
assert uh<0
report['exact_denominator_hostile']={'z':str(zh),'z_less_than_1_over_e':'certified by z*ecap<1','z_greater_than_1_over_b':True,'false_upper_bound_N1':str(uh),'reason':'U<0<R, so the condition bz<1 cannot be dropped.'}

# Original density integration, not the author's integration-by-parts kernel.
def lx(v):
    if not v: return mp.mpf(1)
    if v==mp.pi: return mp.ninf
    return v*mp.cot(v)+mp.log(mp.sin(v)/v)
def D(v):
    if not v: return mp.mpf(0)
    return v*v+(1-v*mp.cot(v))**2
mesh=[0,mp.pi/4,mp.pi/2,3*mp.pi/4,mp.pi*mp.mpf('.95'),mp.pi*mp.mpf('.995'),mp.pi]
def directM(s):
    p=s+1
    def f(v):
        if v==mp.pi or v==0: return mp.mpf(0)
        return mp.exp(p*lx(v))*D(v)/mp.pi
    return mp.quad(f,mesh)
def formulaM(s):
    p=s+1
    return mp.exp((p-2)*mp.log(p))/mp.gamma(p)
report['original_density_mellin_diagnostics']=[]
for s in [mp.mpf('-.99'),mp.mpf('-.9'),mp.mpf('-.5'),mp.mpf(0),mp.mpf(1),mp.mpf(5),mp.mpf('-.7')+mp.mpf('.2')*1j,mp.mpf('.2')+mp.mpf('.7')*1j]:
    actual=directM(s); target=formulaM(s)
    err=abs(actual-target)/max(1,abs(target))
    assert err<mp.mpf('1e-65'),(s,err)
    report['original_density_mellin_diagnostics'].append({'s':str(s),'relative_or_absolute_error':mp.nstr(err,8)})

report['complex_reconstruction_diagnostics']=[]
for w in [mp.mpf('.5'),mp.mpf(1),mp.mpf(7),mp.mpf(50),1+1j,2+3j,mp.mpf('.5')+mp.mpf('.75')*1j]:
    w=mp.mpc(w)
    p=1/w
    from_measure=w*w/(mp.gamma(p)*directM(p-1))
    target=mp.exp(mp.log(w)/w)
    err=abs(from_measure-target)
    assert err<mp.mpf('1e-65'),(w,err)
    report['complex_reconstruction_diagnostics'].append({'w':str(w),'absolute_error':mp.nstr(err,8),'principal_log_reciprocal_error':mp.nstr(abs(mp.log(p)+mp.log(w)),6)})

report['boundary_and_endpoint_diagnostics']=[]
for v in [mp.mpf('.01'),mp.mpf('.7'),mp.mpf('1.9'),mp.mpf('2.9')]:
    x=mp.exp(lx(v)); eta=x*mp.mpf('1e-35')
    up=-mp.lambertw(-1/(x+1j*eta))
    down=-mp.lambertw(-1/(x-1j*eta))
    target=v*mp.cot(v)-1j*v
    assert up.imag<0<down.imag
    err=abs(up-target)
    assert err<mp.mpf('1e-29')
    report['boundary_and_endpoint_diagnostics'].append({'v':str(v),'upper_bank_error':mp.nstr(err,8),'lower_bank_conjugacy_error':mp.nstr(abs(down-up.conjugate()),8)})
for d in [mp.mpf('1e-3'),mp.mpf('1e-5')]:
    upper=(mp.e-mp.exp(lx(d)))/(mp.e*d*d/2)
    lower=mp.exp(lx(mp.pi-d))/(mp.e*d/mp.pi*mp.exp(-mp.pi/d))
    assert abs(upper-1)<d*d
    assert abs(lower-1)<2*d
    report['boundary_and_endpoint_diagnostics'].append({'endpoint_parameter':str(d),'upper_square_root_ratio':mp.nstr(upper,20),'lower_essential_decay_ratio':mp.nstr(lower,20)})

# Strict inequalities on a range not restricted to the submitted z=1/3.
report['strict_remainder_diagnostics']=[]
for N in [0,1,3,8]:
    for zz in [mp.mpf('.1'),mp.mpf('1')/3,mp.mpf('.367'),mp.mpf('-2')]:
        value=-mp.lambertw(-zz)/zz
        partial=sum(mq(a(k))*zz**k for k in range(N+1))
        residual=value-partial
        r=mq(a(N+2)/a(N+1))
        scale=mq(a(N+1))*abs(zz)**(N+1)
        if zz>0:
            low=scale/(1-r*zz)
            high=scale*(1-(mp.e-r)*zz)/(1-mp.e*zz)
            assert low<residual<high
        else:
            assert mp.sign(residual)==(-1)**(N+1)
            assert scale/(1+r*abs(zz))<abs(residual)<scale
        report['strict_remainder_diagnostics'].append({'N':N,'z':str(zz),'status':'PASS'})

(HERE/'independent_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS: independent exact replay of all submitted Gauss/Padé coefficients, bounds, and 80-place outward rounding.')
print('PASS: exact counterexample when bz<1 is omitted.')
print('PASS: original-density real/complex Mellin integration and full complex kernel reconstruction diagnostics.')
print('PASS: branch-bank signs, both endpoint asymptotics, Jensen/chord and negative-axis remainder diagnostics.')
print('Largest audited N:',report['exact_pade_replays'][-1]['nodes'])
print('Largest-N exact enclosure width, diagnostic decimal:',report['exact_pade_replays'][-1]['width_diagnostic'])
