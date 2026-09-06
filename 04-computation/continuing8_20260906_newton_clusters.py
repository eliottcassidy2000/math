"""Effective separated-cluster Newton-circuit law.

All numerical controls use rational coefficient convolution. Universal steps
are symbolic identities and a finite constructive threshold for any profile.
No repository imports; raw stdout and JSON are UTF-8 LF.
"""
from fractions import Fraction as F
from math import comb, ceil, prod, factorial
from pathlib import Path
import hashlib
import json
import sys
import sympy as S

sys.stdout.reconfigure(newline="\n")
GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise ArithmeticError(label)


def sign(a):
    return (a > 0) - (a < 0)


def changes(word):
    w = [a for a in word if a]
    return sum(a != b for a, b in zip(w, w[1:]))


def e_row(roots):
    e = [F(1)]
    for root in roots:
        e.append(F(0))
        for k in range(len(e)-1, 0, -1):
            e[k] += root*e[k-1]
    return e


def ratios(e):
    d = len(e)-1
    out = [None]
    for k in range(1, d):
        out.append(e[k]**2/e[k-1]/e[k+1]
                   * F(comb(d, k-1)*comb(d, k+1), comb(d, k)**2))
    return out


def circuits(R):
    return {k: R[k]/R[k-1] for k in range(2, len(R))}


def packet(profile):
    need(len(profile) >= 2 and all(m >= 2 for m in profile), "cluster multiplicity scope")
    d = sum(profile)
    boundaries, top_coeff = [], [1]
    a = 0
    for m in profile:
        top_coeff.extend(comb(m, l) for l in range(1, m+1))
        a += m
        boundaries.append(a)
    boundaries.pop()
    kappa = ratios(list(map(F, top_coeff)))
    circuit = circuits(kappa)
    good = []
    good += [circuit[k] for k in range(2, boundaries[0])]
    good += [1/circuit[k] for k in range(boundaries[-1]+2, d)]
    for a, b in zip(boundaries, boundaries[1:]):
        good += [circuit[k+1]/circuit[k] for k in range(a+2, b-1)]
    need(all(q > 1 for q in good), "strict factorial-band margin")
    eta = min(good, default=F(2))
    M = 1
    while F(M+1, M)**8 >= eta:
        M *= 2
    tau = F(M+1, M)
    spike_budget = max([1/circuit[b] for b in boundaries]
                       + [circuit[b+1] for b in boundaries])
    T = max(3*M*2**d, ceil(tau**4*spike_budget)+1)
    eps = F(1, 6*M*d)
    need(tau**8 < eta, "strict interior error budget")
    need(T > tau**4*spike_budget, "both signs at every cluster boundary")
    need((1+eps)**d*(1+F(2**d, T)) < tau, "uniform coefficient error bound")
    return d, boundaries, top_coeff, kappa, circuit, M, tau, T, eps


# Exact factorial interior identity, and the positivity behind its convexity.
l, A, B, m = S.symbols("l A B m", positive=True)
f = (l+1)*(A+l)/(l*(A+l+1))
phi1 = -1/(l*(l+1))+1/((l+A)*(l+A+1))
phi2 = 1/l**2-1/(l+1)**2-1/(l+A)**2+1/(l+A+1)**2
need(S.simplify(S.diff(S.log(f), l)-phi1) == 0, "log factorial quotient derivative")
need(S.simplify(S.diff(phi1, l)-phi2) == 0, "log factorial quotient curvature")
den = l**2*(l+1)**2*(l+A)**2*(l+A+1)**2
num = S.Poly(S.cancel(phi2*den/A), A, l)
need(all(c > 0 for c in num.coeffs()), "strict positive symbolic curvature numerator")
interior = (l+1)*(m-l+1)*(A+l)*(B+m-l)/(l*(m-l)*(A+l+1)*(B+m-l+1))
need(S.cancel(interior-f*(m-l+1)*(B+m-l)/((m-l)*(B+m-l+1))) == 0,
     "limiting band ratio is the product of two factorial quotients")

profiles = [(2,2), (2,3), (3,2), (3,5), (2,2,2), (3,3,3),
            (2,5,3), (4,5,3), (6,2,4), (2,7,3), (5,4,2,3),
            (3,3,3,3), (2,3,4,2,3), (2,2,2,2,2), (2,2,2,2,2,2)]
records, ties = [], 0
for profile in profiles:
    d, bounds, tops, kap, ideal, M, tau, T, eps = packet(profile)
    record = {"profile": profile, "M": M, "tau": str(tau), "separation": T,
              "relative_cluster_width": str(eps), "boundaries": bounds,
              "kappa": [str(q) for q in kap[1:]], "controls": []}
    for mode in ("coalesced_geometric", "unequal_gaps_spread"):
        bases = [F(1)]*len(profile)
        for j in range(len(profile)-2, -1, -1):
            bases[j] = bases[j+1]*T*(1 if mode == "coalesced_geometric" else j+2)
        roots = []
        for base, count in zip(bases, profile):
            roots.extend(base*(1 if mode == "coalesced_geometric" else 1+eps*F(i, count-1))
                         for i in range(count))
        e = e_row(roots)
        R = ratios(e)
        c = circuits(R)
        word = [sign(q-1) for q in c.values()]
        need(changes(word) == 2*len(profile)-3, "complete zero-filtered cluster law")
        for b in bounds:
            need(c[b] > 1 and c[b+1] < 1, "actual boundary sign spike pair")
        need(all(c[k] > 1 for k in range(2, bounds[0]+1)), "whole first circuit block positive")
        need(all(c[k] < 1 for k in range(bounds[-1]+1, d)), "whole final circuit block negative")
        for a, b in zip(bounds, bounds[1:]):
            need(changes([sign(c[k]-1) for k in range(a+1, b+1)]) == 1,
                 "exactly one inter-band sign change")
            need(all(c[k+1] > c[k] for k in range(a+2, b-1)), "strict bounded-band circuit ordering")
        # Literal coefficient comparison with the dominant binomial row.
        k, top_scale = 0, F(1)
        need(e[0] == 1, "constant leading coefficient")
        for base, count in zip(bases, profile):
            for ll in range(1, count+1):
                k += 1
                leading = top_scale*base**ll*comb(count, ll)
                need(1 <= e[k]/leading < tau, "uniform coefficient envelope on actual roots")
            top_scale *= base**count
        rev = circuits(ratios(e_row([1/root for root in roots])))
        need(all(rev[k] == 1/c[d+1-k] for k in c), "exact reversal negates and reverses circuits")
        if all(mm == 2 for mm in profile):
            need(changes(word) == d-3 and all(word), "open maximal-alternation cone")
        ties += sum(q == 1 for q in c.values())
        record["controls"].append({"mode": mode, "bases": list(map(str, bases)),
                                    "roots": list(map(str, roots)), "sign_word": word})
    records.append(record)
    print("PROFILE", profile, "T", T, "M", M,
          "words", ["".join("+" if a > 0 else "-" if a < 0 else "0" for a in q["sign_word"])
                    for q in record["controls"]])
need(ties > 0, "zero-filtering exercised by exact reciprocal ties")

# Canonical hostiles delimit the theorem rather than broadening the census.
hostiles = []
for label, roots, expected in [
    ("singleton_bands_1_2_1", [F(100), F(10), F(10), F(1)], 1),
    ("one_narrow_band_need_not_have_zero_changes", [F(1), F(1), F(1001,1000), F(1001,1000)], 1),
    ("THM3004_two_end_classifier", [F(1),F(1),F(3),F(3),F(8)], 2),
]:
    c = circuits(ratios(e_row(roots)))
    word = [sign(q-1) for q in c.values()]
    need(changes(word) == expected, "named multiplicity/radius/classifier hostile")
    hostiles.append({"label": label, "roots": list(map(str, roots)), "word": word})
print("HOSTILES", [(h["label"], h["word"]) for h in hostiles])

# Literal factorial-row consumer: translate each PF factor to make a narrow
# block, dilate the blocks independently, and only then multiply. This changes
# the original coefficient rows and does not infer their original circuit law.
n = S.symbols("n")
profile = (2,3,4)
d, bounds, tops, kap, ideal, M, tau, T, eps = packet(profile)
bases = [F(T**2), F(T), F(1)]
target = S.Integer(1)
factor_controls = []
for degree, base in zip(profile, bases):
    factor = S.Poly(sum(S.Rational(factorial(degree)*comb(degree,j),factorial(j))*n**j
                        for j in range(degree+1)),n)
    bound = 1+max(abs(v) for v in factor.all_coeffs()[1:])
    need(factor.LC() == 1 and factor.count_roots(-bound,0) == degree,
         "independent literal Laguerre factor root scope")
    shift = F(bound)/eps
    scale = base/shift
    moved = S.Poly(S.expand(S.Rational(scale)**degree
                           * factor.as_expr().subs(n,n/S.Rational(scale)+S.Rational(shift))),n)
    need(moved.LC() == 1 and all(v>0 for v in moved.all_coeffs()),
         "positive monic rational transported factorial row")
    target *= moved.as_expr()
    factor_controls.append({"degree":degree,"original_descending_coefficients":list(map(str,factor.all_coeffs())),
                            "Cauchy_bound":str(bound),"shift":str(shift),"scale":str(scale),"base":str(base)})
target = S.Poly(S.expand(target),n)
e = list(map(F,target.all_coeffs()))
c = circuits(ratios(e))
word = [sign(q-1) for q in c.values()]
need(changes(word) == 3, "transported factorial rows give the complete three-block law")
need(all(c[b]>1 and c[b+1]<1 for b in bounds), "transported factorial boundary spike pairs")
k, top_scale = 0,F(1)
for base,count in zip(bases,profile):
    for ll in range(1,count+1):
        k += 1
        need(1 <= e[k]/(top_scale*base**ll*comb(count,ll)) < tau,
             "coefficient-only transported factorial envelope")
    top_scale *= base**count
transport = {"profile":profile,"factors":factor_controls,"target_descending_coefficients":list(map(str,e)),
             "sign_word":word,"T":T,"epsilon":str(eps)}
print("FACTORIAL_TRANSPORT",profile,"T",T,"word",word)

result = {"scope": "effective 2K-3 law for fixed multiplicities>=2, K>=2, narrow well-separated positive root blocks",
          "curvature_numerator": str(num.as_expr()), "profiles": records, "hostiles": hostiles,
          "factorial_transport":transport}
here = Path(__file__).resolve().parent
dest = here.parent/"05-knowledge/results" if here.name == "04-computation" else here
path = dest/(Path(__file__).stem+"_certificate.json")
path.write_bytes((json.dumps(result,indent=2,sort_keys=True)+"\n").encode("utf-8"))
print("CERTIFICATE_SHA256",hashlib.sha256(path.read_bytes()).hexdigest())
print("PASS",GATES,"always-active exact gates; actual LF")
