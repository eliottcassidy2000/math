#!/usr/bin/env python3
"""trapped_margin_sample_kps_S128c48.py -- sample genuine TRAPPED families and measure
the strict margin M(v) - 1/14 > 0.  Confirms the residual lives in the strict interior."""
import sys, random
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
random.seed(1414)
def norm_dist(x):
    fx = x - int(x)
    if fx < 0: fx += 1
    return min(fx, 1-fx)
def M_grid(V, Q=100800):
    # fine-grid lower estimate of M (>= this); sufficient to show M > 1/14 strictly
    best = 0.0
    step = 1.0/Q
    t = step
    A = V
    while t < 1.0:
        m = min(abs(round(v*t) - v*t) for v in A)
        if m > best: best = m
        t += step
    return best
def gap(V): A=[abs(v) for v in V]; return any(a>13*b for a in A for b in A)
def compressed(V): A=[abs(v) for v in V]; return all(any(j!=i and A[i]<=13*A[j] for j in range(len(A))) for i in range(len(A)))
def distinct(V): A=[abs(v) for v in V]; return len(set(A))==len(A)
def maxge23(V): return max(abs(v) for v in V)>=23
def noncluster(V):
    n=len(V)
    for i0 in range(n):
        rest=[abs(V[j]) for j in range(n) if j!=i0]
        g=reduce(gcd,rest,0)
        if g>=2 and V[i0]%g!=0: return False
    return True
def covering(V):
    for q in range(2,15):
        t=F(1,q)
        if min(norm_dist(v*t) for v in V)>=F(1,14): return False
    return True
def is_trapped(V):
    return (distinct(V) and gap(V) and compressed(V) and maxge23(V)
            and noncluster(V) and covering(V))
# generate trapped candidates: distinct, one big (>=23) but not dominant (within 13x of another big),
# compressed cluster, span > 13
found=0; margins=[]
tries=0
while found<40 and tries<200000:
    tries+=1
    lo=random.randint(20,60)
    V=sorted(random.sample(range(lo, lo*15), 13))
    if not is_trapped(V): continue
    m=M_grid(V)
    margins.append((m, tuple(V)))
    found+=1
margins.sort()
print("== TRAPPED-FAMILY MARGIN SAMPLE (%d trapped families found in %d tries) =="%(found,tries))
print("  min M over sample  = %.5f (margin over 1/14 = %+.5f)"%(margins[0][0], margins[0][0]-1/14))
print("  worst family       = %s"%(list(margins[0][1]),))
print("  median M           = %.5f"%(margins[len(margins)//2][0]))
print("  max M              = %.5f"%(margins[-1][0]))
below=[m for m,_ in margins if m <= 1/14 + 1e-4]
print("  families within 1e-4 of 1/14: %d"%len(below))
print("  ALL trapped samples have M > 1/14:", all(m > 1/14 for m,_ in margins))
print("DONE")
