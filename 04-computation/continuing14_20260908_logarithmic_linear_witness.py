"""Exact controls for the analytic all-degree linear-witness obstruction."""
from pathlib import Path
import hashlib
import json
import sys
import sympy as s

sys.stdout.reconfigure(newline='\n')

x, t, v, r, b = s.symbols('x t v r b')
gates = 0


def check(ok, name):
    global gates
    gates += 1
    if not bool(ok):
        raise RuntimeError(name)


def eq(a, bb, name):
    check(s.cancel(a-bb) == 0, name)


def jac(a, bb):
    return s.diff(a, x)*s.diff(bb, t)-s.diff(a, t)*s.diff(bb, x)


def infinity(a):
    return s.expand(a.subs({x:1/r, t:-r**2-r**4*b}, simultaneous=True))


def min_r(a):
    terms = s.expand(a).as_ordered_terms()
    return min(int(term.as_powers_dict().get(r, 0)) for term in terms)


def source_smooth(a, name):
    # This is a separate ideal-of-partials control, not the logarithmic test.
    ideal = s.groebner([s.diff(a, x), s.diff(a, t)], t, x)
    check(list(ideal) == [1], name)


records = {'positive': [], 'negative': [], 'mixed': []}

# Independent finite positive-root universe: 3 roots x 3 B x 4 P degrees.
for alpha in [-1, 0, 2]:
    for BB in [s.Integer(3), x+4, x*x+4]:
        for n in range(4):
            PP = s.sympify(s.prod(v+7+2*j for j in range(n)))
            AA = x-alpha
            HH = AA*t+BB
            FF = s.expand(AA*PP.subs(v, HH))
            eq(jac(FF, HH), FF, 'positive bracket')
            check(PP.subs(v, BB.subs(x, alpha)) != 0, 'positive root condition')
            eq(s.diff(FF, x).subs(x, alpha), PP.subs(v, BB.subs(x, alpha)), 'positive jet')
            pole = -min_r(infinity(FF))
            check(pole > 0, 'positive global obstruction')
            records['positive'].append([alpha, str(BB), n, pole])

# Negative-root universe: 3 roots x 4 C x 3 R degrees.
for beta in [-1, 0, 2]:
    for CC in [s.Integer(0), s.Integer(1), x+2, x*x+1]:
        for m in range(3):
            b0 = s.Integer(3)
            RR = s.sympify(s.prod(v+7+2*j for j in range(m)))
            LL = -t+CC
            HH = b0+(x-beta)*LL
            FF = s.expand(LL*RR.subs(v, HH))
            eq(jac(FF, HH), FF, 'negative bracket')
            eq(HH.subs(x, beta), b0, 'negative H root')
            eq(s.diff(LL, t), -1, 'negative simple component')
            FF_inf = infinity(FF)
            actual_global = min_r(FF_inf) >= 0
            expected_global = CC == 0 or (not CC.has(x) and m == 0)
            check(actual_global == expected_global, 'negative complete globality')
            if actual_global:
                eq(s.diff(FF_inf, r).subs(r, 0), 0, 'whole D normal derivative')
                eq(s.diff(FF_inf, b).subs(r, 0), 0, 'whole D tangent derivative')
                if CC == 0:
                    eq(s.expand(FF_inf).coeff(r, 2), RR.subs(v, b0), 'negative exact boundary order2')
            records['negative'].append([beta, str(CC), m, bool(actual_global), min_r(FF_inf)])

# Mixed-root universe: 3 distinct root pairs x 2 C x 3 R degrees.
for alpha, beta in [(1, 0), (-1, 2), (2, -1)]:
    for CC in [s.Integer(1), x*x+1]:
        for m in range(3):
            b0 = s.Integer(3)
            RR = s.sympify(s.prod(v+31+2*j for j in range(m)))
            AA = (x-alpha)*(x-beta)/s.Integer(alpha-beta)
            BB = b0+(x-beta)*CC
            HH = AA*t+BB
            LL = (x-alpha)*t/s.Integer(alpha-beta)+CC
            FF = s.expand((x-alpha)*LL*RR.subs(v, HH))
            PP = (v-b0)*RR
            eq(FF, (x-alpha)/(x-beta)*PP.subs(v, HH), 'mixed cancellation')
            eq(jac(FF, HH), FF, 'mixed bracket')
            check(PP.subs(v, BB.subs(x, alpha)) != 0, 'mixed positive root')
            eq(PP.subs(v, BB.subs(x, beta)), 0, 'mixed negative root')
            eq(s.diff(FF, x).subs(x, alpha), PP.subs(v, BB.subs(x, alpha))/(alpha-beta), 'mixed vertical jet')
            check(min_r(infinity(FF)) < 0, 'mixed global obstruction')
            records['mixed'].append([alpha, beta, str(CC), m])

# Fully symbolic separation and residue normalizations.
aa, bb, c0, c1, c2 = s.symbols('aa bb c0 c1 c2')
for AA, YY in [(x-aa, x-aa), (-(x-bb), 1/(x-bb)), ((x-aa)*(x-bb)/(aa-bb), (x-aa)/(x-bb))]:
    eq(s.diff(YY, x)/YY, 1/AA, 'all-parameter logarithmic derivative')
    eq(jac(YY*(c0+c1*v+c2*v*v).subs(v, AA*t+x*x+1), AA*t+x*x+1), YY*(c0+c1*v+c2*v*v).subs(v, AA*t+x*x+1), 'all-parameter rational eigenfunction')

# Independent smoothness controls, including nonlinear B and C.
smooth_samples = [
    x, x*((x*t+1)**2+1),
    (-t+x*x+1)*(3+(x-2)*(-t+x*x+1)+7),
    (x-1)*((x-1)*t+1),
    t*(1+x*t), -t+x*x,
]
for i, FF in enumerate(smooth_samples):
    source_smooth(s.expand(FF), 'source ideal control '+str(i))

# Genuine nonzero order-one unit: paid rational field and two labelled jets.
FF = t*(1+x*t)
HH = -x*t
GG = -x/(1+x*t)
eq(jac(FF, HH), FF, 'order1 logarithmic identity')
eq(jac(FF, GG), 1, 'order1 rational mate')
eq(t, FF/(1-HH), 'order1 field inverse t')
eq(x, -HH*(1-HH)/FF, 'order1 field inverse x')
eq(HH.subs(t, 0), 0, 'regular labelled component')
eq(HH.subs(t, -1/x), 1, 'pole labelled component')
eq(jac(-t+x*x, x), 1, 'coordinate unit-zero boundary')
FF_inf = infinity(FF)
check(min_r(FF_inf) == 2, 'order1 hostile boundary multiplicity')
eq(s.diff(FF_inf, r).subs(r, 0), 0, 'order1 hostile normal critical')
eq(s.diff(FF_inf, b).subs(r, 0), 0, 'order1 hostile tangent critical')

# A-degree3 survives logarithmic exactness but not source submersion.
AA = x*(x*x-1)/2
YY = (x*x-1)/(x*x)
HH = AA*t
FF = s.expand((x*x-1)**3*t*t/4)
eq(s.diff(YY, x)/YY, 1/AA, 'cubic logarithmic hostile')
eq(YY*HH**2, FF, 'cubic polynomial cancellation')
eq(jac(FF, HH), FF, 'cubic bracket hostile')
residues = [s.cancel(1/s.diff(AA, x).subs(x, z)) for z in [-1, 0, 1]]
check(residues == [1, -2, 1], 'cubic forbidden residue')
eq(s.diff(FF, x).subs(t, 0), 0, 'cubic critical x derivative')
eq(s.diff(FF, t).subs(t, 0), 0, 'cubic critical t derivative')

# Newton reconstruction is checked independently for monic degrees 2..6.
newton_degrees = []
for n in range(2, 7):
    roots = list(range(1, n+1))
    powers = {j:sum(z**j for z in roots) for j in range(1, n+1)}
    elementary = [s.Integer(1)]
    for k in range(1, n+1):
        elementary.append(s.cancel(sum((-1)**(j-1)*elementary[k-j]*powers[j] for j in range(1, k+1))/k))
    reconstructed = sum((-1)**k*elementary[k]*x**(n-k) for k in range(n+1))
    eq(reconstructed, s.prod(x-z for z in roots), 'Newton reconstruction')
    check(2*n-2 >= n, 'available all-degree moment count')
    newton_degrees.append(n)

source = Path(__file__).resolve()
out_dir = source.parent
if source.parent.name == '04-computation':
    out_dir = source.parent.parent/'05-knowledge'/'results'
certificate = {
    'status':'FINITE-EXACT controls; analytic proof supplies unbounded quantifiers',
    'scope':'Original source-linear logarithmic witness; no bound on first-function degree',
    'gates':gates,
    'declared_universe':{'positive_cases':36, 'negative_cases':36, 'mixed_cases':18, 'newton_degrees':newton_degrees},
    'records':records,
    'source_smooth_groebner_controls':len(smooth_samples),
    'hostile_order1_F':'t*(1+x*t)',
    'hostile_order1_boundary':'entire D critical',
    'cubic_A_residues':[str(z) for z in residues],
}
target = out_dir/(source.stem+'_certificate.json')
target.write_bytes((json.dumps(certificate, indent=2, sort_keys=True)+'\n').encode())
print('PASS: logarithmic witness obstruction exact controls')
print('normal forms: 36 positive, 36 negative, 18 mixed')
print('separate source ideal controls: '+str(len(smooth_samples)))
print('gates: '+str(gates))
print('certificate SHA256: '+hashlib.sha256(target.read_bytes()).hexdigest())
