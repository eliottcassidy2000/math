"""
Targeted adversarial sweep of the DANGEROUS S3 sub-regime: k=2 (exactly two large
runners) and near the S2/S3 boundary (Vmax just >= 13*Vmin). This is where the criterion
C(S) is most constrained and where the known exception [1,2,3,5,7,8,9,10,11,12,13,27,28]
lives (Vmax=28=13*Vmin+... actually Vmin=1 so 13*Vmin=13, but degenerate L'={27}).

We EXHAUSTIVELY enumerate, for each 11-subset P of {1..13} and each pair (a,b) of large
speeds with 14<=a<b<=BCAP:
   S = P u {a,b}, require covering, primitive, case S3.
Then:
   - exact M(S); flag M<1/14 (true counterexample).
   - C(S) via ANY removed runner; flag total failures.
   - record min over-board C margin and which v achieves it.

k=2 means W(S\{v}) for v large leaves only ONE large runner in A -> the large part gives a
COARSE comb (period 1/a or 1/b), the hardest case for finding a wide safe arc. So if C(S)
ever fails, it fails here.

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
    b = F(0)
    for t in cand(S):
        v = min(nrm(x*t) for x in S)
        if v > b: b = v
    return b

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

def main(BCAP=120):
    print(f"BOUNDARY k=2 exhaustive sweep, large speeds in [14,{BCAP}]")
    n_total = 0; n_Mbelow = 0; n_Cfail = 0
    Mbelows = []; Cfails = []
    minM = F(1); worstM = None
    minmargin = F(100); worstmargin = None
    # enumerate 11-subsets of 1..13 (there are C(13,11)=78)
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
                m = Mval(S)
                if m < minM: minM = m; worstM = S
                if m < F(1, 14):
                    n_Mbelow += 1; Mbelows.append((S, m))
                bc, arg = best_C(S)
                if bc <= 1:
                    n_Cfail += 1; Cfails.append((S, m, bc))
                if bc < minmargin: minmargin = bc; worstmargin = (S, arg)
    print(f"  total covering S3 k=2 sets: {n_total}")
    print(f"  #M<1/14 (TRUE COUNTEREXAMPLES): {n_Mbelow}")
    for S, m in Mbelows[:30]:
        print("     M<1/14:", S, "M=", m, float(m))
    print(f"  #C(S) total-failures: {n_Cfail}")
    for S, m, bc in Cfails[:30]:
        print("     C FAILS:", S, "M=", m, float(m), "bestmargin=", float(bc))
    print(f"  min M = {minM} = {float(minM):.5f} at {worstM}")
    print(f"  min best-C-margin = {minmargin} = {float(minmargin):.4f} at {worstmargin}")

if __name__ == '__main__':
    main()
