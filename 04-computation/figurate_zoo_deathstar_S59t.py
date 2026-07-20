#!/usr/bin/env python3
"""
death-star-2026-07-20-S59t (HYP-8175) -- the figurate zoo: continuations of the
triangular numbers, staggered at slopes 0/1/2, combined by SUM and PRODUCT, all
identified.  Plus the n*2^x+1 staggerings and the relational reading.
"""
from math import comb, log
from fractions import Fraction as Fr

# ---------- a small sequence identifier ----------
KNOWN = {
 "2^n":            [2**n for n in range(12)],
 "Fibonacci":      [1,1,2,3,5,8,13,21,34,55,89,144],
 "Lucas":          [2,1,3,4,7,11,18,29,47,76,123,199],
 "Moser A000127":  [1,2,4,8,16,31,57,99,163,256,386,562],
 "Padovan":        [1,1,1,2,2,3,4,5,7,9,12,16],
 "Narayan cows":   [1,1,1,2,3,4,6,9,13,19,28,41],
 "tribonacci":     [0,1,1,2,4,7,13,24,44,81,149,274],
 "factorial":      [1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800],
 "Catalan":        [1,1,2,5,14,42,132,429,1430,4862,16796,58786],
 "squares":        [n*n for n in range(12)],
 "triangular":     [n*(n+1)//2 for n in range(12)],
}
def fit_recurrence(seq, maxord=4):
    """find shortest linear recurrence with small integer coeffs (|c|<=3)."""
    import itertools
    seq=list(seq)
    for order in range(1,maxord+1):
        if len(seq)<2*order+1: break
        for coeffs in itertools.product(range(-3,4),repeat=order):
            if all(c==0 for c in coeffs): continue
            ok=all(seq[i]==sum(coeffs[j]*seq[i-1-j] for j in range(order)) for i in range(order,len(seq)))
            if ok:
                return coeffs
    return None
def growth(seq):
    seq=[abs(x) for x in seq if x!=0]
    if len(seq)<4: return None
    rs=[seq[i+1]/seq[i] for i in range(len(seq)-3,len(seq)-1) if seq[i]]
    return sum(rs)/len(rs) if rs else None
def identify(seq):
    seq=list(seq)
    for name,k in KNOWN.items():
        m=min(len(seq),len(k))
        if m>=5 and seq[:m]==k[:m]: return name+" (exact)"
        # shifted match
        for sh in range(1,4):
            if m-sh>=5 and seq[:m-sh]==k[sh:m]: return f"{name} (shift {sh})"
    rec=fit_recurrence(seq)
    g=growth(seq)
    tag=f"rec{rec}" if rec else "no-small-rec"
    if g: tag+=f" growth~{g:.4f}"
    return tag

# ---------- figurate families as triangles T(row m>=1, col c>=0) ----------
def simplicial(m,c):   # polyhedral: c-dim simplex number = C(m+c-1, c+1)?? use C(m-1+c, c+1)
    return comb(m-1+c, c+1) if m>=1 else 0
def pascal(m,c):       # C(m,c)
    return comb(m,c)
def polygonal(m,c):    # opus-S317: i + (j-1)*C(i,2), row i>=1, col j>=1 -> shift col
    i=m; j=c+1
    return i + (j-1)*comb(i,2)
def powersum(m,c):     # Sum_{k=1}^{m} k^c
    return sum(k**c for k in range(1,m+1))
def centered(m,c):     # centered (c+2)-gonal at index m: 1 + (c+2)*T_{m-1}
    k=c+2
    return 1 + k*(m*(m-1)//2)
def kary(m,c):         # C(m, c+1) = number of (c+1)-subsets = (c+1)-ary relations
    return comb(m, c+1)

FAMILIES={"Pascal":pascal,"simplicial":simplicial,"polygonal":polygonal,
          "power-sum":powersum,"centered":centered,"k-ary(C(m,c))":kary}

def diag_sum(fam, m, s):   # slope-s: sum_c fam(m - s*c, c)
    tot=0
    for c in range(0, m+2):
        r=m - s*c
        if r < (1 if fam in (polygonal,powersum,centered) else 0): break
        v=fam(r,c)
        if v==0 and c>0 and s>0: 
            tot+=v; continue
        tot+=v
    return tot
def diag_prod(fam, m, s):
    prod=1; cnt=0
    for c in range(0, m+2):
        r=m - s*c
        if r < 1: break
        v=fam(r,c)
        if v>0:
            prod*=v; cnt+=1
    return prod if cnt else 0

print("=== FIGURATE ZOO: diagonal SUMS at slopes 0,1,2 (each identified) ===")
NR=11
for name,fam in FAMILIES.items():
    print(f"\n{name}:")
    for s in (0,1,2):
        seq=[diag_sum(fam,m,s) for m in range(1,NR)]
        print(f"  slope {s} SUM : {seq[:9]}   [{identify(seq)}]")

print("\n=== PRODUCTS (slope 1) — growth (log) ===")
for name,fam in FAMILIES.items():
    seq=[diag_prod(fam,m,1) for m in range(1,9)]
    logs=[round(log(x),2) if x>0 else 0 for x in seq]
    print(f"  {name:16s} slope1 PROD: {seq[:6]}  logs {logs}")
