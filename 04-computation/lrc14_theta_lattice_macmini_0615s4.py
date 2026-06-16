#!/usr/bin/env python3
"""
lrc14_theta_lattice_macmini_0615s4.py  (mac-mini-2026-06-15-S4, THM-515)

THE THETA-FUNCTION / LONELY-MEASURE FORM OF THE LRC(14) SINGULAR SERIES.

Claim (THM-515): with h(0)=6/7, h(t)=-sin(πt/7)/(πt) for t≠0, and Λ={t∈Z^13:Σt_i v_i=0}
the relation lattice of S={v_1..v_13},
        L(S) = Σ_{t∈Λ} ∏_{i=1}^{13} h(t_i)        (the singular series, THM-501 form)
             = ∫_0^1 ∏_{i=1}^{13} 1_safe(v_i τ) dτ  (the LONELY-SET MEASURE),
because ĥ(θ)=Σ_t h(t)e(tθ)=1-1_danger(θ)=1_safe(θ) ≥ 0 (h positive-definite, Bochner),
and Poisson summation over Λ collapses the lattice sum to the orbit average = lonely
measure. Recovers THM-501's L=lim D(q,S)/q and makes L≥0 STRUCTURAL.

This script:
 (1) verifies ĥ = 1_safe numerically (the positive-definite identity);
 (2) for several S, compares the THETA partial sum Σ_{t∈Λ,|t|≤B}∏h  vs the LONELY
     MEASURE D(Q,S)/Q at large Q -- confirming the theta form & its short-relation
     dominance;
 (3) correlates L with RELATION-LATTICE geometry: λ_1 = shortest relation (Σ|t_i|),
     and the count of short relations, across generic / near-tight / AP-core S;
 (4) probes the infimum & extremizer via the lonely measure.
"""
import sys, itertools, math
from math import gcd, sin, pi

sys.stdout.reconfigure(line_buffering=True)

def s(t):
    return sin(pi*t/7)/(pi*t) if t != 0 else 1/7

def h(t):
    return 6/7 if t == 0 else -s(t)

def lonely_measure(S, Q):
    rad = Q // 14
    cnt = 0
    for a in range(Q):
        ok = True
        for v in S:
            r = (v*a) % Q
            if r <= rad or r >= Q - rad:
                ok = False; break
        if ok: cnt += 1
    return cnt / Q

def short_relations(S, B):
    """all integer relations Σ t_i v_i = 0, max|t_i|<=B, not all zero, by support<=4 brute force
    (captures the dominant short vectors of the relation lattice)."""
    n = len(S); rels = []
    for sz in range(2, 5):
        for T in itertools.combinations(range(n), sz):
            # coefficients in [-B,B]\{0} on T
            for coeffs in itertools.product([c for c in range(-B, B+1) if c != 0], repeat=sz):
                if sum(c*S[i] for c, i in zip(coeffs, T)) == 0:
                    # primitive (gcd of coeffs = 1) to avoid double counting multiples within box
                    if gcd_list(coeffs) == 1:
                        t = [0]*n
                        for c, i in zip(coeffs, T): t[i] = c
                        rels.append(tuple(t))
    return rels

def gcd_list(xs):
    g = 0
    for x in xs: g = gcd(g, x)
    return g

def theta_partial(S, B):
    """Σ_{t∈Λ, support<=4, |t|<=B} ∏ h(t_i), grouped: main term + short-relation corrections.
    NOTE: includes ALL integer multiples k*rel within the |t|<=B box automatically via the
    coefficient enumeration when not requiring primitivity -- here we enumerate the full box."""
    n = len(S); total = (6/7)**n  # t=0 main term
    seen = set()
    for sz in range(2, 5):
        for T in itertools.combinations(range(n), sz):
            for coeffs in itertools.product([c for c in range(-B, B+1) if c != 0], repeat=sz):
                if sum(c*S[i] for c, i in zip(coeffs, T)) == 0:
                    key = tuple(sorted(zip(T, coeffs)))
                    if key in seen: continue
                    seen.add(key)
                    term = (6/7)**(n-sz)
                    for c in coeffs: term *= h(c)
                    total += term
    return total

def lam1(S, B):
    rels = short_relations(S, B)
    if not rels: return None, 0
    norms = [sum(abs(x) for x in r) for r in rels]
    return min(norms), len(rels)

print("="*74)
print("(1) POSITIVE-DEFINITE IDENTITY:  ĥ(θ) = Σ_t h(t) e(tθ) = 1_safe(θ) ≥ 0")
print("="*74)
# evaluate the partial Fourier sum at sample thetas; danger band |θ|<=1/14
for theta in [0.0, 0.02, 0.05, 0.0714, 0.1, 0.25, 0.5]:
    val = h(0) + sum(2*h(t)*math.cos(2*pi*t*theta) for t in range(1, 400))  # even -> 2cos
    band = "DANGER(|θ|≤1/14=.0714)" if min(theta, 1-theta) <= 1/14 + 1e-9 else "safe"
    print(f"   θ={theta:.4f}: ĥ≈{val:+.4f}  (expect {0.0 if min(theta,1-theta)<=1/14 else 1.0}) [{band}]")

print("\n" + "="*74)
print("(2) THETA partial sum vs LONELY MEASURE D(Q)/Q  (theta form check)")
print("="*74)
configs = {
  "generic 14*{1,2,4,8,..} (Sidon-ish)": sorted(set([14*x for x in (1,2,4,8,16,32,64,128,256,512,1024,2048,4096)])),
  "evader 7*{1..12}∪{13}": sorted(set([7*i for i in range(1,13)] + [13])),
  "near-tight {1..12}∪{14}": sorted(set(list(range(1,13)) + [14])),
  "interior-drop {1..13}\\{6}∪{56}": sorted(set([x for x in range(1,14) if x!=6] + [56])),
}
for name, S in configs.items():
    if len(S) != 13:
        print(f"   {name}: |S|={len(S)} (skip)"); continue
    Lm = lonely_measure(S, 14000)
    th3 = theta_partial(S, 3); th5 = theta_partial(S, 5)
    l1, nshort = lam1(S, 5)
    print(f"   {name}:")
    print(f"      lonely measure D(14000)/14000 = {Lm:.5f}   (main (6/7)^13={ (6/7)**13:.5f})")
    print(f"      theta partial |t|≤3: {th3:.5f}; |t|≤5: {th5:.5f}; λ_1(Σ|t|)={l1}, #short rels(B=5)={nshort}")

print("\n" + "="*74)
print("(3) LATTICE-GEOMETRY PREDICTOR: does λ_1(Λ) track L?  (lonely measure, more S)")
print("="*74)
import random
random.seed(6154)
rows = []
# dilated AP cores d*{1..12}∪{14m} + random primitive mult-of-14 sets
for d in [1,3,5]:
    for m in [1,13,17]:
        core = [d*i for i in range(1,13)]; rr = 14*m
        if rr in core: continue
        S = sorted(set(core+[rr]))
        if len(S)!=13 or gcd_list(S)!=1 or not any(v%14==0 for v in S): continue
        rows.append((f"d{d}m{m}", S))
for _ in range(8):
    base = sorted(random.sample(range(1,45),12)); S = sorted(set(base+[random.choice([14,28,42])]))
    if len(S)==13 and gcd_list(S)==1: rows.append(("rand", S))
print(f"   {'tag':10s} {'lonely L':>9s} {'λ_1(Σ|t|)':>9s} {'#short(B=4)':>11s}")
data = []
for tag, S in rows:
    Lm = lonely_measure(S, 11000); l1, ns = lam1(S, 4)
    data.append((Lm, l1, ns, tag));
    print(f"   {tag:10s} {Lm:9.5f} {str(l1):>9s} {ns:11d}")
# crude correlation: do small-lambda1 / many-short sets have small L?
data.sort()
print(f"   --> lowest-L configs: {[(round(L,4),tag,'λ1='+str(l1),'short='+str(ns)) for L,l1,ns,tag in data[:4]]}")
print(f"   --> highest-L configs: {[(round(L,4),tag,'λ1='+str(l1)) for L,l1,ns,tag in data[-3:]]}")
print("\nDONE.")
