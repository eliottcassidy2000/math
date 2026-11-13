"""Exact finite controls; numerical illustrations are separately labelled."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
from math import comb, expm1, fsum, gamma, lgamma, log, log1p, pi
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
checks = 0

def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


def tt(n):
    return [sum(1 << j for j in range(i+1,n)) for i in range(n)]


def gate(t):
    n = len(t)
    # Old core -> b -> c -> old core.
    return [row | (1 << n) for row in t] + [1 << (n+1), (1 << n)-1]


def h_dp(t):
    n = len(t)
    states = [dict() for _ in range(1 << n)]
    for i in range(n):
        states[1 << i][i] = 1
    for mask, lastcounts in enumerate(states):
        for last, count in lastcounts.items():
            choices = t[last] & ~mask
            while choices:
                bit = choices & -choices
                nxt = bit.bit_length()-1
                newmask = mask | bit
                states[newmask][nxt] = states[newmask].get(nxt,0)+count
                choices -= bit
    return sum(states[-1].values())


def path_words(t):
    n = len(t)
    for word in permutations(range(n)):
        if all(t[u] & (1 << v) for u,v in zip(word,word[1:])):
            yield word


def gate_poly(p):
    # q[P(q+1)+qP(q)-P(q-1)], exact integer coefficients.
    out = [0] * (len(p)+2)
    for k,c in enumerate(p):
        out[k+2] += c
        for j in range(k):
            if (k-j) % 2:
                out[j+1] += 2*c*comb(k,j)
    return out

# Independent direct path enumeration gives the entire multiaffine polynomial.
carrier=[]
for N in range(2,9):
    core = tt(N-1)
    t = gate(core)
    subset_counts = Counter()
    words = list(path_words(t))
    for word in words:
        cut = word.index(N-1)
        mask = sum(1 << v for v in word[:cut] if v < N-1)
        subset_counts[mask] += 1
    expected = Counter({mask:1 for mask in range(1 << (N-1))})
    expected[(1 << (N-1))-1] += 1
    require(subset_counts == expected, 'weighted subset polynomial')
    require(len(words) == h_dp(t) == 2**(N-1)+1, 'carrier H')
    carrier.append({'N':N,'vertices':N+1,'H':len(words),'monomials':len(subset_counts)})

# Direct finite towers versus a separate polynomial implementation, then all
# finite depths through 30 of the exact identity and its proved estimates.
C3=[2,4,1]
polys = [[0,1], [0,0,0,1], [0,2,0,1]]
tournaments = [[0],tt(3),C3]
towers=[]
for depth in range(31):
    values = [sum(p) for p in polys]
    one, trans, cyclic = values
    require(cyclic-trans == 2*one, 'absolute gate defect')
    require(one >= 3**depth, 'gate growth')
    require((depth+1)*one <= trans, 'relative gate forgetting')
    if depth <= 3:
        for i,t in enumerate(tournaments):
            require(h_dp(t) == values[i], 'gate polynomial versus subset DP')
            require(sum(1 for _ in path_words(t)) == values[i], 'gate versus permutations')
        tournaments=[gate(t) for t in tournaments]
    towers.append({'r':depth,'H_K1':one,'H_TT3':trans,'H_C3':cyclic,'defect':2*one})
    polys=[gate_poly(p) for p in polys]

dyadic_controls=[]
for j in range(5):
    r=2**j
    row=towers[r]
    eps=Fraction(row['defect'],row['H_TT3'])
    require(0<eps<=Fraction(2,r+1),'summable relative error control')
    require(row['defect']>=2*3**r,'dyadic absolute defect control')
    dyadic_controls.append({'j':j,'r':r,'relative_error':str(eps),'absolute_defect':row['defect']})
require(sum((Fraction(2,2**j+1) for j in range(20)),Fraction())<4,'relative error envelope prefix')

# Exact Taylor coefficients of log(1+s*a), s=+/-1, split into alternating
# artanh and the even logarithmic part. This checks signs and all first 30
# coefficients; the infinite remainder proof is in PROOF.md.
for s in [-1,1]:
    for k in range(1,31):
        direct=Fraction((-1)**(k+1)*s**k,k)
        split=Fraction(s,k) if k%2 else Fraction(-1,k)
        require(direct == split, 'odd/even logarithm split')

# Rational Machin certificate: arctan alternating truncations, no float pi.
def atan_interval(q, terms=25):
    value=sum((Fraction((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(terms)),Fraction())
    next_term=Fraction((-1)**terms,(2*terms+1)*q**(2*terms+1))
    return min(value,value+next_term),max(value,value+next_term)
a,b=atan_interval(5)
c,d=atan_interval(239)
pi_lo, pi_hi=16*a-4*d,16*b-4*c
require(Fraction(3)<pi_lo<pi_hi<Fraction(22,7),'pi bounds')
# Volume coefficient: V_n=coefficient*pi**power.
volume_coeff=[Fraction(1),Fraction(2)]
for n in range(2,31):
    volume_coeff.append(Fraction(2,n)*volume_coeff[n-2])
intervals=[]
for n,c in enumerate(volume_coeff):
    k=n//2
    intervals.append((c*pi_lo**k,c*pi_hi**k))
require(intervals[0] == (1,1),'V0 convention')
require(all(intervals[n][0]>1 for n in range(1,13)),'strict superlevel prefix')
require(all(intervals[n][1]<1 for n in range(13,31)),'post-cutoff finite controls')
require(all(intervals[5][0]>intervals[n][1] for n in range(31) if n != 5),'dimension-five maximum controls')
require(volume_coeff[10]/volume_coeff[3] == Fraction(1,160),'3,10 volume ratio')
require(volume_coeff[11]/volume_coeff[4] == Fraction(128,10395),'4,11 volume ratio')
require(pi_hi**4 < 160,'first gap-seven ratio below one')
require(128*pi_hi**3 < 10395,'second gap-seven ratio below one')
require(Fraction(160)*128 < Fraction(10395)*pi_lo,'gap-seven radius order')

# Proof-preserving boundary probes for the claimed renormalization depth.
thresholds=[]
for alpha in [Fraction(1,8),Fraction(1,6),Fraction(1,4),Fraction(1,3),Fraction(1,2),Fraction(3,4),Fraction(1)]:
    J=(Fraction(1,2)/alpha).__floor__()
    require(2*(J+1)*alpha>1,'tail exponent exceeds one')
    require(J==0 or 2*J*alpha<=1,'last included divergent counterterm')
    thresholds.append({'alpha':str(alpha),'J':J,'tail_exponent':str(2*(J+1)*alpha)})

# Numerical illustrations only; no use as infinite-series certificates.
illustrations=[]
for alpha in [.25,.5,.75,1.0]:
    logterms=[]
    fterms=[]
    for n in range(2,100001):
        f=expm1(log(n)/n)
        a=f**alpha
        logterms.append(log1p(a if n%2==0 else -a))
        fterms.append(f)
        if n in [100,1000,10000,100000]:
            logP=fsum(logterms)
            row={'alpha':alpha,'N':n,'log_P':logP}
            if alpha == .5:
                row['critical_log_scaled']=logP+log(n)**2/4
                row['exact_counterterm_residual']=logP+fsum(fterms)/2
            if alpha < .5:
                row['leading_log_ratio']=logP/(-n**(1-2*alpha)*log(n)**(2*alpha)/(2*(1-2*alpha)))
            illustrations.append(row)

sources=[
'01-canon/theorems/THM-4900-gate-and-doubling-in-the-marked-resolvent-golden-towers-and-tower-ocf.md',
'05-knowledge/results/session127_gsma_norm_collapse_20261031/volume/THEORY.md',
'05-knowledge/results/session162_mrb_square_20261103/mrb_weighted/RESULT.md',
'05-knowledge/results/session144_jacobian_eth_mrb_20261103/mrb/RESEARCH.md',
'05-knowledge/results/microcosm_gate_doubling_resolvent_20261106.md',
]
pins={p:sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode('utf-8')).hexdigest() for p in sources}
receipt={
'status':'FINITE-EXACT controls plus separately labelled numerical illustrations',
'explicit_checks':checks,
'universe':{'carrier_N':list(range(2,9)),'direct_gate_depths':list(range(4)),'polynomial_gate_depths':list(range(31)),'volume_dimensions':list(range(31)),'log_coefficients':list(range(1,31))},
'filters':'No isomorphism quotient; labelled explicit tournaments only; no inference to all tournaments.',
'controls':{'positive':'whole marked subset polynomial; direct H versus polynomial gate operator','hostile':'ratio tends to one but exact defect grows; signed partition is distinct from H; parameter boundary retained'},
'carrier':carrier,
'gate_towers':towers,
'dyadic_summable_relative_hostile':dyadic_controls,
'pi_rational_interval':{'lower':str(pi_lo),'upper':str(pi_hi),'method':'Machin, 25 alternating terms per arctan'},
'unit_volume_threshold':{'strict_integer_superlevel':list(range(1,13)),'nonstrict_integer_superlevel':list(range(13)),'finite_tail_check_through':30,'all_tail_proof':'two-step recurrence / PROOF.md C5'},
'renormalization_boundary_controls':thresholds,
'NUMERICAL_NOT_CERTIFIED':{'MRB_products':illustrations,'gap_seven_equality_radii':[(160/pi**4)**(1/7),(10395/(128*pi**3))**(1/7)],'ball_scaled_radius_critical':1/(2*pi*2.718281828459045)**.5},
'source_pins_text_lf_sha256':pins,
}
(HERE/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'checks':checks,'carrier_cases':len(carrier),'tower_depths':len(towers),'status':'PASS','receipt':str(HERE/'verification.json')}))
