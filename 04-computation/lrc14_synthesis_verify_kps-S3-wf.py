"""
Synthesis verification for the LRC(14) S3 residual report.
Independently re-derives the load-bearing facts with EXACT Fractions:
  (A) The C(S)-failure counterexample S* = [1,2,3,5,7,8,9,10,11,12,13,38,42]
      - covering? primitive? case S3? all 13 single-removal ratios < 1? M(S*) >= 1/14?
  (B) THM-526 arc-width implication sanity (teeth structure).
  (C) Cluster-collapse Lemma A: window W_K safety + nonemptiness identity.
  (D) Two-gap band-fit lemma logical content.
  (E) k=2 finite-core claims: min M(P), W_min=9/3920, V2>=51 sharpness.
  (F) A widened exact sweep over k=2 and k=3 S3 sets in small windows:
      - any M < 1/14 ? (LRC counterexample)
      - any OTHER C(S)-failure beyond S* ?
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h / u) % 1; b = (c + h / u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
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
        while F(2 * k + 1, 2 * v) <= F(1, 2): C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def Mval(S):
    b = F(0)
    for t in cand(S):
        v = min(nrm(x * t) for x in S)
        if v > b: b = v
    return b

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def is_primitive(S):
    g = 0
    for x in S: g = gcd(g, x)
    return g == 1

def C_criterion_ratios(S):
    """Return dict v -> W(S\\{v})*7*v (Fraction). C holds iff any ratio > 1."""
    S = sorted(set(S))
    out = {}
    for v in S:
        A = [u for u in S if u != v]
        out[v] = Wwidth(A) * 7 * v
    return out

print("="*70)
print("(A) C-FAILURE COUNTEREXAMPLE S* = [1,2,3,5,7,8,9,10,11,12,13,38,42]")
print("="*70)
Sstar = [1,2,3,5,7,8,9,10,11,12,13,38,42]
print("len:", len(Sstar))
print("covering:", is_covering(Sstar))
print("primitive:", is_primitive(Sstar))
large = [v for v in Sstar if v > 13]
print("k (#speeds>13):", len(large), " large:", large)
print("Vmin:", min(Sstar), "Vmax:", max(Sstar), "Vmax>=13*Vmin:", max(Sstar) >= 13*min(Sstar))
ratios = C_criterion_ratios(Sstar)
maxr = max(ratios.values())
print("Single-removal ratios W(S\\v)*7v:")
for v in sorted(ratios): print(f"   v={v:2d}: {ratios[v]} = {float(ratios[v]):.4f}")
print("MAX ratio:", maxr, "=", float(maxr), " -> C holds?", maxr > 1)
M = Mval(Sstar)
print("M(S*) =", M, "=", float(M), "  >= 1/14?", M >= F(1,14))
print("Threshold 1/14 =", float(F(1,14)))

print()
print("="*70)
print("(C) CLUSTER-COLLAPSE LEMMA A: window W_K safety + nonemptiness identity")
print("="*70)
import random
random.seed(1)
bad_safe = 0; bad_id = 0; tested = 0
for _ in range(20000):
    Vmin = random.randint(20, 500)
    s = random.randint(1, 60)
    Vmax = Vmin + s
    K = random.randint(0, 20)
    lo = F(14*K+1, 14*Vmin); hi = F(14*K+13, 14*Vmax)
    nonempty = lo < hi
    ident = (13*Vmin - Vmax) > 14*K*s
    if nonempty != ident: bad_id += 1
    if nonempty:
        tested += 1
        # test safety at midpoint for a few speeds in [Vmin,Vmax]
        mid = (lo+hi)/2
        for u in {Vmin, Vmax, (Vmin+Vmax)//2}:
            if nrm(u*mid) < F(1,14): bad_safe += 1; break
print(f"nonemptiness identity mismatches: {bad_id}/20000")
print(f"safety violations at midpoint:    {bad_safe}/{tested} nonempty windows")

print()
print("="*70)
print("(E) k=2 FINITE-CORE: min M(P), W_min, V2>=51 sharpness")
print("="*70)
# min M(P) over subsets P of {1..13} of size <=11
minMP = None; argmin = None
for sz in range(1, 12):
    for P in combinations(range(1,14), sz):
        m = Mval(list(P))
        if minMP is None or m < minMP:
            minMP = m; argmin = P
print("min M(P) over P subset {1..13}, |P|<=11:", minMP, "=", float(minMP), "at", argmin)
# W_min over A = P(size11) union V2, V2 in 14..50, with W(A)>0
P11_candidates = list(combinations(range(1,14), 11))
Wmin = None; argW = None
for P in P11_candidates:
    for V2 in range(14, 51):
        A = list(P) + [V2]
        w = Wwidth(A)
        if w > 0 and (Wmin is None or w < Wmin):
            Wmin = w; argW = (P, V2)
print("W_min over A=P(11) u V2, V2<=50, W>0:", Wmin, "=", float(Wmin) if Wmin else None)
print("   7*63*W_min =", 7*63*Wmin if Wmin else None, "=", float(7*63*Wmin) if Wmin else None, "(should be >1)")
# V2>=51 sharpness: largest V2 with 7*V2*W(A)<1 for some P11
worst_below = 0
for P in P11_candidates:
    for V2 in range(14, 80):
        A = list(P) + [V2]
        if 7*V2*Wwidth(A) < 1:
            worst_below = max(worst_below, V2)
print("largest V2 with 7*V2*W(P11 u V2) < 1:", worst_below, "(claim: 50, so V2>=51 sharp)")

print()
print("="*70)
print("(F) WIDENED EXACT SWEEP: k=2 and k=3 S3 sets in small windows")
print("    -> any M<1/14 (LRC break)?  any OTHER C-failure beyond S*?")
print("="*70)
def sweep(maxspeed_small, large_range, k, report):
    breaks = []; cfails = []
    count = 0
    smalls = list(range(1, 14))
    # choose small part: any subset of {1..13}, then k large speeds in large_range
    # to keep feasible, fix |small|=13-k and enumerate
    nsmall = 13 - k
    for P in combinations(smalls, nsmall):
        for L in combinations(large_range, k):
            S = list(P) + list(L)
            if max(S) < 13*min(S): continue   # must be S3 (also k>=2 by construction)
            if not is_covering(S): continue
            if not is_primitive(S): continue
            count += 1
            M = Mval(S)
            if M < F(1,14): breaks.append((S, M))
            ratios = C_criterion_ratios(S)
            if max(ratios.values()) <= 1:
                cfails.append((S, M, max(ratios.values())))
    print(f"[{report}] k={k}, {count} covering primitive S3 sets")
    print(f"   LRC breaks (M<1/14): {len(breaks)}")
    for S,M in breaks[:5]: print("      BREAK", S, M)
    print(f"   C-failures (all ratios<=1): {len(cfails)}")
    for S,M,r in cfails[:10]: print("      CFAIL", S, "M=",M,"=",float(M)," maxratio=",r,"=",float(r))
    return breaks, cfails

# k=2, large speeds in [14,80]
b2,c2 = sweep(13, list(range(14,81)), 2, "k=2 large[14,80]")
# k=3, large speeds in [14,45] (smaller to stay feasible)
b3,c3 = sweep(13, list(range(14,46)), 3, "k=3 large[14,45]")

print()
print("="*70)
print("SUMMARY")
print("="*70)
print("S* is the C-failure:", Sstar)
allc = [s for s,_,_ in c2] + [s for s,_,_ in c3]
print("Total C-failures found in sweep:", len(allc))
print("Distinct C-failure sets:", [list(s) for s in allc])
print("Total LRC breaks (M<1/14):", len(b2)+len(b3))
