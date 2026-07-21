#!/usr/bin/env python3
"""large_n_circulant_patterns_boxeph_S197.py -- boxeph-2026-07-21-S197

n>=7 is where full enumeration dies (A000568: 456, 6880, 191536 at n=7,8,9). But the CIRCULANT
(rotational) tournaments stay enumerable at every odd n (only 2^{(n-1)/2} connection sets), so they
are the window into large-n structure. This computes circulant families at n=7,9,11,13,19 and reports
the patterns: count of iso classes, c3, eigenvalue spread var(|lambda|^2) (kps's GIT scalar), the
Re=-1/2 universality, and which family is Paley (flat, doubly regular) vs interval vs max-spread.
"""
import numpy as np
from itertools import permutations
from math import pi, comb

def circ(n, C):
    A=[[0]*n for _ in range(n)]
    for i in range(n):
        for k in C: A[i][(i+k)%n]=1
    return A
def num_c3(A,n):
    # c3 = C(n,3) - sum C(score,2)
    sc=[sum(A[i]) for i in range(n)]
    return comb(n,3)-sum(comb(s,2) for s in sc)
def canon_small(A,n):
    # canonical code only feasible n<=8; for larger use rotation+multiplier orbit reduction
    best=None
    for p in permutations(range(n)):
        c=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or c<best: best=c
    return best
def circ_orbit_code(n,C):
    """canonical connection-set under rotation-invariance + multiplier group (fast, exact for circulants)."""
    best=None
    for a in range(1,n):
        if np.gcd(a,n)!=1: continue
        Cp=frozenset((a*k)%n for k in C)
        # tournament requires exactly one of {k,n-k}; multiplier preserves this
        key=tuple(sorted(Cp))
        if best is None or key<best: best=key
    return best

print("n : #circ-connsets  #circ-iso(multiplier-orbits)  c3-range        var(|lam|^2)-range   Paley?")
for n in range(7,20,2):
    pairs=[(k,n-k) for k in range(1,(n+1)//2)]
    seen=set(); data=[]
    for choice in range(1<<len(pairs)):
        C=set(a if (choice>>t&1) else b for t,(a,b) in enumerate(pairs))
        code=circ_orbit_code(n,C)
        if code in seen: continue
        seen.add(code)
        A=circ(n,list(C))
        w=np.exp(2j*pi/n)
        lam=[sum(w**((j*k)%n) for k in C) for j in range(n)]
        mods=[abs(lam[j])**2 for j in range(1,n)]
        var=float(np.var(mods))
        remax=max(abs(lam[j].real+0.5) for j in range(1,n))   # deviation from Re=-1/2
        data.append((num_c3(A,n), var, remax, tuple(sorted(C))))
    c3s=sorted(set(d[0] for d in data)); vars_=sorted(set(round(d[1],2) for d in data))
    paley=[d for d in data if d[1]<1e-6]  # zero spread = Gauss sum flat = Paley (n=3 mod 4)
    reok=max(d[2] for d in data)
    print("%2d :  %4d            %4d                        [%d..%d]        [%.2f..%.2f]        %s  (Re=-1/2 dev max %.1e)"
          %(n, 1<<len(pairs), len(seen), min(c3s),max(c3s), min(vars_),max(vars_),
            "YES" if paley else "no ", reok))

print()
print("PALEY-p specifics (p = 3 mod 4): the doubly-regular flat-spectrum circulant")
for p in (7,11,19,23):
    QR=sorted({(k*k)%p for k in range(1,p)})
    A=circ(p,QR)
    c3=num_c3(A,p)
    # Paley is doubly regular: every pair of vertices has (p-3)/4 common out-neighbors
    w=np.exp(2j*pi/p); lam=[sum(w**((j*k)%p) for k in QR) for j in range(p)]
    modsq=sorted(set(round(abs(lam[j])**2,3) for j in range(1,p)))
    # c3 closed form for Paley: p(p-1)(p+1)/24
    c3_formula=p*(p-1)*(p+1)//24
    print("  Paley-%2d: c3=%d (formula p(p^2-1)/24=%d, match=%s); |lam_nonPerron|^2=%s (=p? %s); doubly-reg lambda=(p-3)/4=%d"
          %(p,c3,c3_formula,c3==c3_formula, modsq, modsq==[float(p)], (p-3)//4))

print()
print("KNOWN SEQUENCES reaching n>=7 (for extrapolation):")
A000568=[1,1,1,2,4,12,56,456,6880,191536,9733056]  # n=0..10 iso tournaments
strong_iso=[None,1,0,1,1,6,35,353,6008,178133]      # n=1..9 strongly-connected iso (OEIS A051337-ish)
print("  A000568 iso tournaments n=0..10:", A000568)
print("  strong-connected iso n=1..9   :", strong_iso[1:])
for n in range(3,9):
    tot=A000568[n]; st=strong_iso[n] if n<len(strong_iso) and strong_iso[n] else None
    if st: print("    n=%d: total=%d strong=%d reducible=%d  strong-fraction=%.3f" % (n,tot,st,tot-st,st/tot))
