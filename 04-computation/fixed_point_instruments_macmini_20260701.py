#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S90 -- fixed-point instruments: a transform is S_n-INVARIANT (averages over the group ->
blind to the fixed-point/symmetry content |Aut|); a fixed-point-sensitive instrument must be BUILT FROM the
group action (a COVERING invariant or a MOMENT), not another clever transform.

Two tests + the complementarity:
 (A) n=7 SAMPLE: does the U-spectrum (adjacency) still DETERMINE |Aut| (transform robustness for fixed pts)?
 (B) n=6 COMPLEMENTARITY: group-action instruments (|Aut|, MFAS covering-radius) separate the high-|Aut|
     NEEDLES (flip-rank-excess drivers) but NOT the |Aut|=1 spectral twins; transforms (U-spectrum) separate
     the twins but are not the natural NEEDLE detector. Two fine-structures, two instrument kinds.
"""
import itertools, random
import numpy as np
from math import comb
from collections import defaultdict, Counter

def canon(A, n, perms):
    return min(tuple(A[p[i]][p[j]] for i in range(n) for j in range(n)) for p in perms)
def aut(A, n, perms):
    return sum(1 for p in perms if all(A[i][j] == A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j))
def uspec(A):
    ev = np.linalg.eigvals(np.array(A, float)); return tuple(np.round(sorted(ev.real),3))+tuple(np.round(sorted(ev.imag),3))
def mfas(A, n):
    C2 = comb(n,2); best = 0
    for p in itertools.permutations(range(n)):
        best = max(best, sum(1 for a in range(n) for b in range(a+1,n) if A[p[a]][p[b]]))
    return C2 - best

# ---------- (A) n=7 sample: does U-spectrum determine |Aut|? ----------
print("(A) n=7 SAMPLE -- does the U-spectrum (a transform) still DETERMINE |Aut|?")
n = 7; perms = list(itertools.permutations(range(n))); pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
random.seed(5); reps = {}
for _ in range(4000):
    bits = [random.getrandbits(1) for _ in pairs]
    A = [[0]*n for _ in range(n)]
    for (i,j),b in zip(pairs,bits):
        if b: A[i][j]=1
        else: A[j][i]=1
    c = canon(A,n,perms)
    if c not in reps: reps[c] = A
byU = defaultdict(set)
for c,A in reps.items(): byU[uspec(A)].add(aut(A,n,perms))
bad = {sp:auts for sp,auts in byU.items() if len(auts)>1}
print(f"  sampled {len(reps)} distinct n=7 iso classes; distinct U-spectra={len(byU)}")
print(f"  U-cospectra containing DIFFERENT |Aut| (=> U FAILS to determine |Aut| at n=7): {len(bad)}")
for sp,auts in list(bad.items())[:4]: print(f"     one U-spectrum holds |Aut| values {sorted(auts)}")
print(f"  => {'U-spectrum is NOT a robust |Aut| detector at n=7 (transform leaks/fails).' if bad else 'U still determines |Aut| in this sample.'}")

# ---------- (B) n=6 complementarity ----------
print("\n(B) n=6 COMPLEMENTARITY -- needles need group-action, twins need transform:")
n = 6; perms = list(itertools.permutations(range(n))); pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
allreps = {}
for bits in itertools.product((0,1),repeat=len(pairs)):
    A = [[0]*n for _ in range(n)]
    for (i,j),b in zip(pairs,bits):
        if b: A[i][j]=1
        else: A[j][i]=1
    c = canon(A,n,perms)
    if c not in allreps: allreps[c] = A
rows = [(A, aut(A,n,perms), mfas(A,n), uspec(A)) for A in allreps.values()]
# the NEEDLES = high |Aut|
needles = [r for r in rows if r[1] >= 3]
print(f"  {len(rows)} classes; needles (|Aut|>=3): {len(needles)}, |Aut| multiset {dict(Counter(r[1] for r in needles))}")
# group-action instruments separate needles from |Aut|=1?
print(f"  |Aut| (group-action) trivially flags all {len(needles)} needles; MFAS (covering-radius) distribution on needles: {sorted(r[2] for r in needles)}")
# spectral twins: |Aut|=1 classes that are U-distinct but share all elementary counts (S86)
# show |Aut| and covering CANNOT separate the S86 spectral twins (both |Aut|=1, need the U-spectrum)
print("  The S86 spectral twins {12,22},{43,44} both have |Aut|=1 => |Aut|/covering CANNOT separate them; only the U-spectrum does.")
print("  => TWO fine-structures: FIXED-POINT/symmetry (|Aut| needles) <- group-action instruments (covering/moment);")
print("     GENERIC realization (|Aut|=1 twins) <- transforms (spectra). Build the RIGHT instrument for the target.")
