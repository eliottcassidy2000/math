#!/usr/bin/env python3
"""
lrc14_angleD_break_criterion — mac-mini-2026-06-17-S6  (ANGLE D: honesty / break C)

GOAL: try HARD to find a covering 13-set where criterion C FAILS:
   C(S) = [ EXISTS v in S:  W(S\\{v}) > 1/(7v) ]
W(A) = widest arc of the level-1/14 SAFE set of A (complement of the danger teeth
       ||u tau|| < 1/14 for u in A; each u contributes u teeth of full width 1/(7u)).

C(S) => M(S) >= 1/14   (generalized arc-width lemma, THM-526 extended).

SCALE-INVARIANCE (proved below, ATTACK 0): the criterion margin for cS via cv equals
the margin for S via v divided by c. So scaling NEVER flips C; we may restrict runner
MAGNITUDES to a bounded window and still probe every geometric configuration. This is
the key efficiency win (Wsafe cost ~ sum of runner sizes, so we keep runners moderate).

Attacks: (1) clustered-large covering sets (bounded magnitude), (2) two/three equal-largest
runners, (3) private-obligation sole-cover sporadics (codex HYP-2579), (4) greedy descent
that MINIMIZES the best margin.

If C FAILS: recompute M exactly. C-fail with M>=1/14 => criterion needs strengthening.
C-fail with M<1/14 => would REFUTE LRC(14) (triple-checked).
"""
from fractions import Fraction as F
from math import lcm
import random

C = F(1, 14)

# ---------- efficient exact arc-width on the circle ----------
def Wsafe(A, c=C):
    """Widest SAFE arc (exact). Build sorted danger-interval endpoints, sweep gaps.
    Each u in A gives u teeth centered k/u, half-width c/u. Returns 0 if fully covered."""
    iv = []
    for u in set(A):
        hw = F(c, u)
        for k in range(u):
            ctr = F(k, u)
            iv.append((ctr - hw, ctr + hw))
    if not iv:
        return F(1)
    # normalize each into [0,1), splitting wrap-arounds
    norm = []
    for lo, hi in iv:
        shift = lo - (lo % 1)
        a = lo - shift; b = hi - shift
        if b <= 1:
            norm.append((a, b))
        else:
            norm.append((a, F(1))); norm.append((F(0), b - 1))
    norm.sort()
    # union
    merged = []
    cl, ch = norm[0]
    for lo, hi in norm[1:]:
        if lo <= ch:
            if hi > ch: ch = hi
        else:
            merged.append((cl, ch)); cl, ch = lo, hi
    merged.append((cl, ch))
    # widest gap (safe arc) including the wrap gap
    best = F(0); n = len(merged)
    for i in range(n):
        hi = merged[i][1]
        lo = merged[(i + 1) % n][0] + (1 if i == n - 1 else 0)
        gap = lo - hi
        if gap > best: best = gap
    return best if best > 0 else F(0)

def covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def crit(S):
    res = []; best = None
    for v in sorted(set(S)):
        A = [u for u in S if u != v]
        W = Wsafe(A); thr = F(1, 7 * v); m = W - thr
        res.append((v, W, thr, W > thr, m))
        if best is None or m > best[4]: best = (v, W, thr, W > thr, m)
    return any(r[3] for r in res), res, best

# ---------- exact M (guarded for large runners) ----------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); Cc = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1,2): Cc.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1,2): Cc.add(F(k, d)); k += 1
    Cc.add(F(1,2)); return Cc
def M(S): return max(g(S, t) for t in cand(S))

# ===================================================================
print("="*78)
print("ANGLE D — BREAKING criterion C: search for a covering 13-set with C(S) FALSE")
print("="*78)

# ---------- ATTACK 0: scale-invariance check ----------
print("\n[ATTACK 0] scale-invariance: margin(cS via cv) = margin(S via v)/c (cannot flip C)")
core = [1,2,3,4,5,6,7,8,9,10,11,13,84]
_,_,b1 = crit(core)
ok_scale = True
for c in [2,3,5,7,30]:
    S = sorted(x*c for x in core)
    _,_,b2 = crit(S)
    if b1[4] != b2[4]*c: ok_scale = False
    print(f"  c={c:3d}: margin(S)={b1[4]}  margin(cS)={b2[4]}  c*margin(cS)={c*b2[4]}  match={b1[4]==c*b2[4]}")
print(f"  scale-invariance EXACT: {ok_scale}  => scaled cores can NEVER break C.")

rng = random.Random(12345)
ntot = 0; cfail = 0; mbreak = 0
worst = (F(99), None, None)
cfail_examples = []

# ---------- ATTACK 1: clustered-large (bounded magnitude) ----------
def clustered_covering(N, win, rng):
    used = set(); S = []
    qs = list(range(2, 15)); rng.shuffle(qs)
    for q in qs:
        cands = [x for x in range(N, N + win + 1) if x % q == 0 and x not in used]
        if not cands: return None
        x = rng.choice(cands); used.add(x); S.append(x)
    S = sorted(set(S))
    return S if len(S) == 13 and covering(S) else None

print("\n[ATTACK 1] clustered covering sets, all runners large (bounded window <=1200)")
for _ in range(30000):
    N = rng.choice([60, 100, 200, 400, 700, 1000])
    win = rng.choice([14, 20, 28, 40, 60, 90, 140, 220])
    S = clustered_covering(N, win, rng)
    if S is None: continue
    ntot += 1
    holds, res, best = crit(S)
    if not holds:
        cfail += 1; Mv = M(S); cfail_examples.append((S, best, Mv))
        if Mv < C: mbreak += 1
    elif best[4] < worst[0]:
        worst = (best[4], S, best[0])
print(f"  tested {ntot}; C failures {cfail}; tightest margin {float(worst[0]):.8f} at v={worst[2]}")
print(f"     S={worst[1]}")

# ---------- ATTACK 2: equal-largest runners ----------
print("\n[ATTACK 2] two near-equal large runners (sole cover of missing moduli)")
def equal_top_covering(rng):
    drop = rng.sample(list(range(1,14)), rng.choice([1,2,3]))
    base = [v for v in range(1, 14) if v not in drop]
    missing = [q for q in range(2,15) if not any(v%q==0 for v in base)]
    L = lcm(*missing) if missing else 1
    k = rng.randint(1, 8)
    S = sorted(set(base + [L*k, L*(k+1)]))
    fill = [v for v in range(1,14) if v not in S]; rng.shuffle(fill)
    for v in fill:
        if len(S) >= 13: break
        S = sorted(set(S+[v]))
    if len(S)!=13 or not covering(S): return None
    return S
for _ in range(30000):
    S = equal_top_covering(rng)
    if S is None: continue
    ntot += 1
    holds, res, best = crit(S)
    if not holds:
        cfail += 1; Mv = M(S); cfail_examples.append((S, best, Mv))
        if Mv < C: mbreak += 1
    elif best[4] < worst[0]:
        worst = (best[4], S, best[0])
print(f"  cumulative tested {ntot}; C failures {cfail}; tightest margin {float(worst[0]):.8f} at v={worst[2]}")
print(f"     S={worst[1]}")

# ---------- ATTACK 3: private-obligation sole-cover sporadics ----------
print("\n[ATTACK 3] private-obligation hard cores (sole-cover large incommensurate)")
def build_private_obligation(rng):
    privset = sorted(set(rng.choice([
        [8,9,11,13],[8,11,13],[9,11,13],[11,13],[7,8,9,11,13],
        [8,9,10,11,13],[8,9,11,12,13],[5,7,8,9,11,13],
    ])))
    other = [q for q in range(2,15) if q not in privset]
    pool = list(range(1,14)); rng.shuffle(pool)
    core=[]; need=set(other)
    for x in pool:
        covs = {q for q in other if x%q==0}
        if covs & need: core.append(x); need -= covs
        if not need: break
    if need: return None
    larges=[]; mode = rng.choice(['one_per','combined','mixed'])
    if mode=='one_per':
        for q in privset:
            larges.append(q*rng.choice([1,2,3,5,7,11])*rng.choice([1,1,13,17,19]))
    elif mode=='combined':
        L = lcm(*privset); larges=[L*rng.choice([1,2,3,5,7]), L*rng.choice([11,13,17])]
    else:
        rng.shuffle(privset); h=len(privset)//2 or 1
        g1,g2 = privset[:h], privset[h:]
        if g1: larges.append(lcm(*g1)*rng.choice([1,3,5,7,11]))
        if g2: larges.append(lcm(*g2)*rng.choice([1,3,5,7,11]))
    S = sorted(set(core+larges))
    fill = [v for v in range(1,14) if v not in S]; rng.shuffle(fill)
    for v in fill:
        if len(S)>=13: break
        S = sorted(set(S+[v]))
    if len(S)!=13 or not covering(S): return None
    # keep runners bounded so M/Wsafe stay fast (scale-invariance justifies)
    if max(S) > 3000: return None
    return S
for _ in range(60000):
    S = build_private_obligation(rng)
    if S is None: continue
    ntot += 1
    holds, res, best = crit(S)
    if not holds:
        cfail += 1; Mv = M(S); cfail_examples.append((S, best, Mv))
        if Mv < C: mbreak += 1
    elif best[4] < worst[0]:
        worst = (best[4], S, best[0])
print(f"  cumulative tested {ntot}; C failures {cfail}; tightest margin {float(worst[0]):.8f} at v={worst[2]}")
print(f"     S={worst[1]}")

# ---------- ATTACK 4: greedy descent minimizing best margin ----------
print("\n[ATTACK 4] greedy descent: minimize best criterion margin (bounded runners)")
def best_margin(S):
    holds,_,best = crit(S); return best[4], holds
def neighbors(S, rng):
    S=sorted(set(S)); out=[]
    for i,v in enumerate(S):
        if v < 50: continue
        for dv in (-14,-7,7,14,-28,28,-1,1,-2,2,-42,42):
            w = v+dv
            if w<=0 or w>3000: continue
            T = sorted(set(S[:i]+S[i+1:]+[w]))
            if len(T)==13 and covering(T): out.append(T)
    rng.shuffle(out); return out[:24]
descent_worst = worst
for restart in range(15):
    if worst[1] is None: break
    S = list(worst[1]); cm,_ = best_margin(S)
    for step in range(40):
        improved=False
        for T in neighbors(S,rng):
            tm,holds = best_margin(T); ntot+=1
            if not holds:
                cfail+=1; Mv=M(T); cfail_examples.append((T,crit(T)[2],Mv))
                if Mv<C: mbreak+=1
            if tm < cm:
                S=T; cm=tm; improved=True
                if tm < descent_worst[0]: descent_worst=(tm,T,crit(T)[2][0])
                break
        if not improved: break
print(f"  descent tightest margin = {float(descent_worst[0]):.8f}  at S={descent_worst[1]}")

# ---------- VERDICT ----------
print("\n" + "="*78); print("VERDICT"); print("="*78)
print(f"  total covering 13-sets tested: {ntot}")
print(f"  criterion C FAILURES (no v with W>1/(7v)): {cfail}")
print(f"  of those, ACTUAL LRC breaks M(S)<1/14: {mbreak}")
gm = min(worst[0], descent_worst[0])
print(f"  global tightest successful margin: {float(gm):.8f}")
if cfail == 0:
    print("  ==> C NEVER FAILED. Strong evidence C is universal on covering 13-sets.")
else:
    print(f"  ==> C failed {cfail} times. Examples (S, best_margin, M):")
    for S,best,Mv in cfail_examples[:6]:
        print(f"      S={S}  best_margin={float(best[4]):.6f}  M={Mv}={float(Mv):.5f}  "
              f"{'*** LRC BREAK ***' if Mv<C else '(C-only fail, M ok)'}")
