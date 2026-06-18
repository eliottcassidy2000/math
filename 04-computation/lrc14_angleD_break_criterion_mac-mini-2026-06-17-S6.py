#!/usr/bin/env python3
"""
lrc14_angleD_break_criterion — mac-mini-2026-06-17-S6  (ANGLE D: honesty / break C)

GOAL: try HARD to find a covering 13-set where criterion C FAILS:
   C(S) = [ EXISTS v in S:  W(S\\{v}) > 1/(7v) ]
W(A) = widest arc of the level-1/14 SAFE set of A (complement of the danger teeth
       ||u tau|| < 1/14 for u in A; each u contributes u teeth of full width 1/(7u)).

C(S) => M(S) >= 1/14  (generalized arc-width lemma, THM-526 extended).

We attack from every direction:
  - maximally clustered covering sets (all runners large, in a tight window)
  - covering sets with two/three EQUAL-largest runners
  - scaled tight cores (kind-pasteur THM-522 scale-invariance: M(cS)=M(S))
  - engineered "every-removal-tiny" sets
  - greedy local descent to MINIMIZE the best criterion margin

If C FAILS on any S: we then check whether M(S) >= 1/14 still holds.
  * C-fail with M>=1/14  => criterion sufficient-not-necessary; needs strengthening.
  * C-fail with M<1/14   => would REFUTE LRC(14). Triple-check with exact M.

CAREFUL ROBUSTNESS NOTE on Wsafe: if the danger teeth COVER the whole circle
(total danger measure >= 1 with no gap), there is NO safe arc and W should be 0.
The S5 Wsafe could return a spurious value in that degenerate case; we harden it.
"""
from fractions import Fraction as F
import random

C = F(1, 14)

# ---------- exact arc-width machinery (hardened) ----------
def darcs(v, c=C):
    hw = F(c, v)
    return [(F(k, v) - hw, F(k, v) + hw) for k in range(v)]

def wrapU(iv):
    """Union of intervals on the circle R/Z, returned as sorted disjoint [lo,hi) in [0,1)."""
    o = []
    for lo, hi in iv:
        s = lo - (lo % 1)          # integer floor offset
        a = lo - s; b = hi - s
        if b <= 1:
            o.append((a, b))
        else:
            o.append((a, F(1))); o.append((F(0), b - 1))
    o = sorted(o); r = []; cl, ch = o[0]
    for lo, hi in o[1:]:
        if lo <= ch:
            ch = ch if ch > hi else hi
        else:
            r.append((cl, ch)); cl, ch = lo, hi
    r.append((cl, ch))
    return r

def Wsafe(A, c=C):
    """Widest SAFE arc = widest gap between danger components on the circle.
    Returns exact Fraction. 0 if danger covers the whole circle."""
    dz = []
    for v in set(A):
        dz += darcs(v, c)
    if not dz:
        return F(1)
    dz = wrapU(dz)
    # total danger measure; if it wraps to cover the circle, gaps are the safe arcs.
    best = F(0)
    n = len(dz)
    for i in range(n):
        hi = dz[i][1]
        # next component start, wrapping +1 for the last
        lo = dz[(i + 1) % n][0] + (1 if i == n - 1 else 0)
        gap = lo - hi
        if gap > best:
            best = gap
    # gap can't exceed 0 if danger fully covers; max with 0 for safety
    return best if best > 0 else F(0)

def covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def crit(S):
    """Return (holds, list of (v,W,thr,W>thr), best_margin_tuple)."""
    res = []
    best = None
    for v in sorted(set(S)):
        A = [u for u in S if u != v]
        W = Wsafe(A); thr = F(1, 7 * v); m = W - thr
        ok = W > thr
        res.append((v, W, thr, ok, m))
        if best is None or m > best[4]:
            best = (v, W, thr, ok, m)
    holds = any(r[3] for r in res)
    return holds, res, best

# ---------- exact M ----------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); Cc = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            Cc.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        Cc.add(F(k, d)); k += 1
    Cc.add(F(1, 2)); return Cc
def M(S): return max(g(S, t) for t in cand(S))

# ===================================================================
# ATTACK 1: maximally-clustered covering sets (all runners large)
# ===================================================================
def clustered_covering(N, win, rng):
    used = set(); S = []
    qs = list(range(2, 15)); rng.shuffle(qs)
    for q in qs:
        cands = [x for x in range(N, N + win + 1) if x % q == 0 and x not in used]
        if not cands:
            return None
        x = rng.choice(cands); used.add(x); S.append(x)
    S = sorted(set(S))
    return S if len(S) == 13 and covering(S) else None

print("=" * 78)
print("ANGLE D — BREAKING criterion C: search for a covering 13-set with C(S) FALSE")
print("=" * 78)

rng = random.Random(12345)
ntot = 0; cfail = 0; mbreak = 0
worst = (F(99), None, None)
cfail_examples = []

print("\n[ATTACK 1] maximally-clustered covering sets (all runners large, tight window)")
for _ in range(40000):
    N = rng.choice([100, 200, 400, 800, 1500, 3000, 7000, 20000, 60000])
    win = rng.choice([14, 20, 28, 40, 60, 90, 140, 220])
    S = clustered_covering(N, win, rng)
    if S is None:
        continue
    ntot += 1
    holds, res, best = crit(S)
    if not holds:
        cfail += 1
        Mv = M(S)
        cfail_examples.append((S, best, Mv))
        if Mv < C:
            mbreak += 1
    else:
        if best[4] < worst[0]:
            worst = (best[4], S, best[0])
print(f"  tested {ntot} clustered covering 13-sets")
print(f"  C FAILED: {cfail}   (of those M<1/14: {mbreak})")
print(f"  tightest successful margin W-1/(7v) = {float(worst[0]):.8f} at v={worst[2]}")
print(f"     S = {worst[1]}")

# ===================================================================
# ATTACK 2: equal-largest runners (engineer two/three tied tops)
# ===================================================================
print("\n[ATTACK 2] covering sets with TWO/THREE equal-largest runners")
# pick a base small core that already covers, then duplicate a large multiple twice
# Need 13 DISTINCT runners; equal-largest means two near-equal large values.
def equal_top_covering(rng):
    # cover 2..14 with small runners using <=11 distinct, then add 2 large equalish
    drop = rng.sample([1,2,3,4,5,6,7,8,9,10,11,12,13], rng.choice([1,2,3]))
    base = [v for v in range(1, 14) if v not in drop]
    # large multiples chosen so they (with base) still cover, near-equal magnitude
    k = rng.choice([84,168,252,420,840,1260])
    a = k * rng.randint(1,3)
    b = a + 14 * rng.choice([0,0,0,0,0])  # try exact-equal-ish via lcm shift; keep distinct
    # build a near-equal pair that are distinct multiples covering what base misses
    # choose two distinct large multiples of lcm(missing) close together
    missing = [q for q in range(2,15) if not any(v%q==0 for v in base)]
    from math import lcm
    L = lcm(*missing) if missing else 1
    base_mult = rng.randint(1, 20)
    x1 = L * base_mult
    x2 = L * (base_mult + 1)
    S = sorted(set(base + [x1, x2]))
    # pad/trim to 13 distinct
    if len(S) < 13:
        extra_pool = [v for v in range(1,14) if v not in S]
        for e in extra_pool:
            if len(S) >= 13: break
            S = sorted(set(S+[e]))
    if len(S) != 13 or not covering(S):
        return None
    return S

for _ in range(40000):
    S = equal_top_covering(rng)
    if S is None: continue
    ntot += 1
    holds, res, best = crit(S)
    if not holds:
        cfail += 1
        Mv = M(S)
        cfail_examples.append((S, best, Mv))
        if Mv < C: mbreak += 1
    else:
        if best[4] < worst[0]:
            worst = (best[4], S, best[0])
print(f"  cumulative tested {ntot}; C failures so far {cfail}")
print(f"  tightest successful margin {float(worst[0]):.8f} at v={worst[2]}  S={worst[1]}")

# ===================================================================
# ATTACK 3: scaled tight cores (scale-invariance of M; arc-width thr shrinks ~1/v)
# ===================================================================
print("\n[ATTACK 3] scaled covering cores cS (M(cS)=M(S); thresholds 1/(7cv) shrink)")
base_cores = [
    [1,2,3,4,5,6,7,8,9,10,11,13,84],
    [1,2,3,4,5,7,8,9,10,11,12,13,84],
    [2,3,4,5,6,7,9,10,11,12,13,168,858],
    [1,2,3,4,5,6,7,8,9,10,11,12,13][:12]+[91],  # placeholder, fixed below
]
# fix the 4th to be covering 13-set
base_cores[3] = [1,2,3,4,5,6,7,8,9,10,11,13,91]
for core in base_cores:
    if len(set(core))!=13 or not covering(core):
        continue
    for c in [1,2,3,5,7,11,13,1001]:
        S = sorted(set(x*c for x in core))
        if len(S)!=13 or not covering(S):
            continue
        ntot += 1
        holds,res,best = crit(S)
        if not holds:
            cfail += 1; Mv=M(S); cfail_examples.append((S,best,Mv))
            if Mv<C: mbreak+=1
        else:
            if best[4] < worst[0]:
                worst=(best[4],S,best[0])

# ===================================================================
# ATTACK 4: greedy local descent — MINIMIZE the best criterion margin
# ===================================================================
print("\n[ATTACK 4] greedy local descent to MINIMIZE best margin over covering 13-sets")
def best_margin(S):
    holds,res,best = crit(S)
    return best[4], holds
def neighbors(S, rng):
    """perturb one large runner to a nearby multiple keeping covering."""
    S=sorted(set(S)); out=[]
    for i,v in enumerate(S):
        if v < 50:  # keep small core
            continue
        for dv in (-14,-7,7,14,-28,28,-1,1,-2,2):
            w = v+dv
            if w<=0: continue
            T = sorted(set(S[:i]+S[i+1:]+[w]))
            if len(T)==13 and covering(T):
                out.append(T)
    rng.shuffle(out)
    return out[:30]
# seed from the tightest found so far
descent_worst = worst
for restart in range(40):
    if worst[1] is None: break
    S = list(worst[1])
    for step in range(200):
        cm,_ = best_margin(S)
        improved=False
        for T in neighbors(S,rng):
            tm,holds = best_margin(T)
            ntot+=1
            if not holds:
                cfail+=1; Mv=M(T); cfail_examples.append((T,crit(T)[2],Mv))
                if Mv<C: mbreak+=1
            if tm < cm:
                S=T; cm=tm; improved=True
                if tm < descent_worst[0]:
                    descent_worst=(tm,T,crit(T)[2][0])
                break
        if not improved:
            break
print(f"  descent tightest margin = {float(descent_worst[0]):.8f}")
print(f"     at S = {descent_worst[1]}")

# ===================================================================
# VERDICT
# ===================================================================
print("\n" + "=" * 78)
print("VERDICT")
print("=" * 78)
print(f"  total covering 13-sets tested: {ntot}")
print(f"  criterion C FAILURES (no v with W>1/(7v)): {cfail}")
print(f"  of those, ACTUAL LRC breaks M(S)<1/14: {mbreak}")
gm = min(worst[0], descent_worst[0])
print(f"  global tightest successful margin: {float(gm):.8f}")
if cfail == 0:
    print("  ==> C NEVER FAILED. Strong evidence C is universal on covering 13-sets.")
    print("      (Sufficient: proving C universal proves LRC(14).)")
else:
    print(f"  ==> C failed {cfail} times. Examples (S, best_margin, M):")
    for S,best,Mv in cfail_examples[:5]:
        print(f"      S={S}  best_margin={float(best[4]):.6f}  M={Mv}={float(Mv):.5f}  "
              f"{'*** LRC BREAK ***' if Mv<C else '(C-only fail, M ok)'}")
