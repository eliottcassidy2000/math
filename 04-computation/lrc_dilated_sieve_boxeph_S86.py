#!/usr/bin/env python3
"""
The DILATED-AP SIEVE: a provable M>=1/13 witness for the compact minimizers (boxeph-S86)
========================================================================================
LEMMA (explicit witness, elementary): if V contains d*{1,...,12} (a dilated AP of length
12) and every OTHER element w satisfies dist(w, 13d*Z) >= d, then M(V) >= 1/13, via
t* = 1/(13d):  the AP part d*i lands on i/13 (all >= 1/13), and each extra w has
||w/(13d)|| = dist(w,13d*Z)/(13d) >= d/(13d) = 1/13.

This proves M>=1/13 for the conjectured MINIMIZERS 2*{1..12}u{13}, d*{1..12}u{killer},
i.e. the exact families that defeat every single-runner perturbation tool.  It is the
dilation-aware analogue of the empty-circle sieve (sieve_frac).  Here we verify the
lemma and measure how much of the M<1/13-adjacent / compact-covering stratum it covers.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, itertools

def norm(x):
    r=x%1; return min(r,1-r)
def exact_M(V):
    best=F(0); qs=set([13,14])
    for i in range(len(V)):
        for j in range(i,len(V)):
            s=V[i]+V[j]
            for dd in range(1,s+1):
                if s%dd==0: qs.add(dd)
    for q in qs:
        for a in range(1,q):
            if gcd(a,q)==1:
                m=min(min((v*a)%q,q-(v*a)%q) for v in V); c=F(m,q)
                if c>best: best=c
    return best
def cov(V,n=14): return all(any(v%b==0 for v in V) for b in range(2,n+1))
def prim(V): return reduce(gcd,V)==1

def dist_to(w, mod):
    r = w % mod
    return min(r, mod-r)

def dilated_sieve_certifies(V):
    """Return (d, t*) if V contains d*{1..12} and all extras w have dist(w,13d*Z)>=d."""
    Vs = set(V)
    for d in range(1, max(V)//12 + 2):
        ap = set(d*i for i in range(1,13))
        if ap <= Vs:
            extras = [w for w in V if w not in ap]
            if all(dist_to(w, 13*d) >= d for w in extras):
                return d, F(1,13*d)
    return None

# ---- verify the lemma on the conjectured minimizers ----
print("="*84)
print("DILATED-AP SIEVE LEMMA  (V >= d*{1..12} + safe extras => M>=1/13 at t=1/(13d))")
print("="*84)
mins = [
    ("2*{1..12}u{13}",   sorted([2*i for i in range(1,13)]+[13])),
    ("3*{1..12}u{13}",   sorted([3*i for i in range(1,13)]+[13])),
    ("2*{1..12}u{39}",   sorted([2*i for i in range(1,13)]+[39])),
    ("5*{1..12}u{13}",   sorted([5*i for i in range(1,13)]+[13])),
    ("2*{1..12}u{65}",   sorted([2*i for i in range(1,13)]+[65])),
]
for name,V in mins:
    if len(set(V))!=13:
        print(f"  {name}: not 13 distinct"); continue
    cert = dilated_sieve_certifies(V)
    M = exact_M(V)
    if cert:
        d,t = cert
        # double-check witness
        mg = min(norm(F(v)*t) for v in V)
        print(f"  {name:16s} cover={cov(V)} prim={prim(V)}  CERT d={d} t*={t} min_gap={mg} (>=1/13:{mg>=F(1,13)})  M={M}")
    else:
        print(f"  {name:16s} NOT certified by dilated sieve;  M={M}")

# ---- coverage: how many compact covering families with small M does it certify? ----
print("\n" + "="*84)
print("COVERAGE on structured compact covering families (dilated-AP + killer)")
print("="*84)
random.seed(41)
tested=0; cert_cnt=0; small_cert=0; small_total=0; missed_small=[]
def consider(V):
    global tested,cert_cnt,small_cert,small_total
    V=sorted(set(V))
    if len(V)!=13 or not prim(V) or not cov(V): return
    if F(max(V),sorted(V)[-2])>=13: return  # compact only
    tested+=1
    M=exact_M(V)
    c=dilated_sieve_certifies(V)
    if c: cert_cnt+=1
    if M<=F(1,13):
        small_total+=1
        if c: small_cert+=1
        elif len(missed_small)<8: missed_small.append((M,V))
# dilated AP + killer
for d in range(1,9):
    ap=[d*i for i in range(1,13)]
    for w in range(1,400):
        consider(ap+[w])
# random compact covering
for _ in range(4000):
    consider(random.sample(range(1,60),13))

print(f"compact covering tested: {tested}")
print(f"  dilated-sieve certifies (M>=1/13): {cert_cnt}/{tested}")
print(f"  among the TIGHTEST (M<=1/13): certified {small_cert}/{small_total}")
if missed_small:
    print("  M<=1/13 families NOT certified by dilated sieve (need generalization):")
    for M,V in missed_small: print(f"     M={M} V={V}")
else:
    print("  => EVERY compact covering family with M<=1/13 is certified by the dilated sieve!")
print("\nREADING: the dilated sieve PROVES M>=1/13 exactly for the minimizer stratum")
print("(the families that defeat single-runner perturbation). If it certifies ALL")
print("M<=1/13 compact families, the compact floor = dilated-sieve + (M>1/13 for the rest).")
