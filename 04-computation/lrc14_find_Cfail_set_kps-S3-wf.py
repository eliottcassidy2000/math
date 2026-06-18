"""
Pinpoint the covering S3 k=2 set(s) where C(S) FAILS for all removed runners, and compute
exact M(S) to decide: criterion-gap (M>=1/14) vs TRUE LRC14 counterexample (M<1/14).

Same enumeration as the boundary sweep (11-subset P of {1..13} + pair {a,b}, 14<=a<b<=80),
but records EVERY C-failure with full exact data. Fast because Mval is only called on failures.
kind-pasteur-S3-wf
"""
from fractions import Fraction as F
from math import gcd
import itertools

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b - a for a, b in sc]
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        ws.append((sc[0][1]) + (1 - sc[-1][0]))
    return max(ws)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def Mval(S):
    b = F(0); arg = None
    for t in cand(S):
        v = min(nrm(x*t) for x in S)
        if v > b: b = v; arg = t
    return b, arg

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def gcd_all(S):
    g = 0
    for x in S: g = gcd(g, x)
    return g

def case_of(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13*Vmin: return 'S2'
    return 'S3'

def best_C(S):
    best = F(0); arg = None
    for v in S:
        A = [x for x in S if x != v]
        r = Wwidth(A) * 7 * v
        if r > best: best = r; arg = v
    return best, arg

def main(BCAP=80):
    fails = []
    n_total = 0
    for P in itertools.combinations(range(1, 14), 11):
        P = list(P)
        for a in range(14, BCAP+1):
            for b in range(a+1, BCAP+1):
                S = sorted(set(P) | {a, b})
                if len(S) != 13: continue
                if gcd_all(S) != 1: continue
                if not is_covering(S): continue
                if case_of(S) != 'S3': continue
                n_total += 1
                bc, arg = best_C(S)
                if bc <= 1:
                    m, tau = Mval(S)
                    fails.append((S, bc, m, tau))
                    import sys
                    print(f"  >>> C-FAIL FOUND: S={S} margin={float(bc):.4f} M={m}={float(m):.6f} "
                          f"(>=1/14? {m>=F(1,14)}) {'**TRUE COUNTEREXAMPLE**' if m<F(1,14) else '(criterion gap)'}",
                          flush=True)
    print(f"total covering S3 k=2 sets (BCAP={BCAP}): {n_total}")
    print(f"C-failures (C(S) false for ALL removed runners): {len(fails)}")
    truecex = 0
    for S, bc, m, tau in fails:
        is_cex = m < F(1, 14)
        if is_cex: truecex += 1
        print(f"  S={S}")
        print(f"     best-C-margin={bc}={float(bc):.4f} (<1)  M(S)={m}={float(m):.6f} at tau={tau}")
        print(f"     >=1/14? {m>=F(1,14)}   {'**TRUE LRC14 COUNTEREXAMPLE**' if is_cex else '(criterion gap only; M still safe)'}")
        # also show per-runner ratios
        rr = []
        for v in S:
            A=[x for x in S if x!=v]
            rr.append((v, float(Wwidth(A)*7*v)))
        print("     per-v ratio:", ", ".join(f"{v}:{r:.3f}" for v,r in rr))
    print(f"\nTRUE LRC14 counterexamples (M<1/14): {truecex}")

if __name__ == '__main__':
    main()
