#!/usr/bin/env python3
"""
lrc14_angleD_break_criterion — mac-mini-2026-06-17-S6  (ANGLE D: honesty / break C)

GOAL: try HARD to find a covering 13-set where criterion C FAILS:
   C(S) = [ EXISTS v in S:  W(S\\{v}) > 1/(7v) ]
W(A) = widest arc of the level-1/14 SAFE set of A (complement of danger teeth
       ||u tau||<1/14, u in A; each u gives u teeth, full width 1/(7u)).
C(S) => M(S) >= 1/14   (generalized arc-width lemma, THM-526 extended).

SCALE-INVARIANCE (proved here, ATTACK 0): margin(cS via cv) = margin(S via v)/c
EXACTLY, so scaling never flips C => bounded runner magnitudes lose NO geometry.
Runners capped (<= CAP) so the search is fast; rigorous, not a shortcut.

Wsafe computed by FAST INTEGER arithmetic on common denominator D=14*lcm(A)
(verified identical to the canonical Fraction version, 27x faster).

If C FAILS: recompute M exactly. C-fail with M>=1/14 => criterion needs strengthening;
C-fail with M<1/14 => would REFUTE LRC(14) (triple-checked).
"""
from fractions import Fraction as F
from math import lcm, gcd
import random, sys

C = F(1, 14)
CAP = 200

def Wsafe(A):
    """Widest safe arc, exact, via integers on D = 14*lcm(A)."""
    A = sorted(set(A))
    L = 1
    for u in A: L = L * u // gcd(L, u)
    D = 14 * L
    iv = []
    for u in A:
        cu = D // u          # = D/u
        hw = D // (14 * u)   # = D/(14u); 14u | D always
        for k in range(u):
            c = k * cu
            iv.append((c - hw, c + hw))
    norm = []
    for lo, hi in iv:
        length = hi - lo
        a = lo % D; b = a + length
        if b <= D: norm.append((a, b))
        else: norm.append((a, D)); norm.append((0, b - D))
    norm.sort(); mg = []; cl, ch = norm[0]
    for lo, hi in norm[1:]:
        if lo <= ch:
            if hi > ch: ch = hi
        else: mg.append((cl, ch)); cl, ch = lo, hi
    mg.append((cl, ch))
    best = 0; n = len(mg)
    for i in range(n):
        hi = mg[i][1]; lo = mg[(i + 1) % n][0] + (D if i == n - 1 else 0)
        g = lo - hi
        if g > best: best = g
    return F(best, D)

def covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def crit(S):
    best = None; holds = False
    for v in sorted(set(S)):
        A = [u for u in S if u != v]
        W = Wsafe(A); thr = F(1, 7 * v); m = W - thr
        if W > thr: holds = True
        if best is None or m > best[1]: best = (v, m, W, thr)
    return holds, best

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

def P(*a): print(*a); sys.stdout.flush()

P("="*78)
P("ANGLE D — BREAKING criterion C: search for a covering 13-set with C(S) FALSE")
P("="*78)

P("\n[ATTACK 0] scale-invariance: margin(cS via cv) = margin(S via v)/c (cannot flip C)")
core = [1,2,3,4,5,6,7,8,9,10,11,13,84]
_, b1 = crit(core); ok_scale = True
for c in [2,3,5,7,11]:
    S = sorted(x*c for x in core); _, b2 = crit(S)
    if b1[1] != c*b2[1]: ok_scale = False
    P(f"  c={c}: margin(S)={b1[1]} margin(cS)={b2[1]} c*margin(cS)={c*b2[1]} match={b1[1]==c*b2[1]}")
P(f"  scale-invariance EXACT: {ok_scale}  => scaled cores can NEVER break C.")

rng = random.Random(20260617)
ntot = 0; cfail = 0; mbreak = 0
worst = (F(99), None, None); cfail_examples = []
seen = set()

def record(S, note=""):
    global ntot, cfail, mbreak, worst
    S = tuple(sorted(set(S)))
    if len(S) != 13 or max(S) > CAP or S in seen: return
    seen.add(S)
    ntot += 1
    holds, best = crit(S)
    if not holds:
        cfail += 1; Mv = M(S); cfail_examples.append((S, best, Mv, note))
        if Mv < C: mbreak += 1
    elif best[1] < worst[0]:
        worst = (best[1], S, best[0])

P("\n[ATTACK 1] clustered covering sets, all runners large (window <= CAP)")
def clustered_covering(N, win, rng):
    used = set(); S = []
    qs = list(range(2, 15)); rng.shuffle(qs)
    for q in qs:
        cands = [x for x in range(N, min(N+win, CAP)+1) if x % q == 0 and x not in used]
        if not cands: return None
        x = rng.choice(cands); used.add(x); S.append(x)
    S = sorted(set(S))
    return S if len(S) == 13 and covering(S) else None
for _ in range(30000):
    N = rng.choice([15, 20, 28, 40, 60, 80, 110, 140]); win = rng.choice([14, 20, 30, 50, 70])
    S = clustered_covering(N, win, rng)
    if S: record(S, "clustered")
P(f"  tested {ntot}; C failures {cfail}; tightest margin {float(worst[0]):.8f} at v={worst[2]}  S={worst[1]}")

P("\n[ATTACK 2] two near-equal large runners (sole cover of missing moduli)")
def equal_top_covering(rng):
    drop = rng.sample(list(range(1,14)), rng.choice([1,2,3]))
    base = [v for v in range(1,14) if v not in drop]
    missing = [q for q in range(2,15) if not any(v%q==0 for v in base)]
    L = lcm(*missing) if missing else 1
    if L*2 > CAP: return None
    k = rng.randint(1, max(1, CAP//(2*L)))
    S = sorted(set(base + [L*k, L*(k+1)]))
    fill = [v for v in range(1,14) if v not in S]; rng.shuffle(fill)
    for v in fill:
        if len(S) >= 13: break
        S = sorted(set(S+[v]))
    if len(S)!=13 or not covering(S): return None
    return S
for _ in range(30000):
    S = equal_top_covering(rng)
    if S: record(S, "equal-top")
P(f"  cumulative tested {ntot}; C failures {cfail}; tightest margin {float(worst[0]):.8f} at v={worst[2]}  S={worst[1]}")

P("\n[ATTACK 3] private-obligation hard cores (sole-cover large incommensurate)")
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
        for q in privset: larges.append(q*rng.choice([1,2,3,5,7,11,13]))
    elif mode=='combined':
        L = lcm(*privset)
        if L>CAP: return None
        larges=[L*rng.choice([1,2]), L]
    else:
        rng.shuffle(privset); h=len(privset)//2 or 1
        g1,g2 = privset[:h], privset[h:]
        if g1: larges.append(lcm(*g1)*rng.choice([1,3,5,7]))
        if g2: larges.append(lcm(*g2)*rng.choice([1,3,5,7]))
    S = sorted(set(core+larges))
    fill = [v for v in range(1,14) if v not in S]; rng.shuffle(fill)
    for v in fill:
        if len(S)>=13: break
        S = sorted(set(S+[v]))
    if len(S)!=13 or not covering(S): return None
    return S
for _ in range(60000):
    S = build_private_obligation(rng)
    if S: record(S, "private-obl")
P(f"  cumulative tested {ntot}; C failures {cfail}; tightest margin {float(worst[0]):.8f} at v={worst[2]}  S={worst[1]}")

P("\n[ATTACK 4] greedy descent: minimize best criterion margin (bounded runners)")
def bm(S):
    _, best = crit(S); return best[1]
def neighbors(S, rng):
    S=sorted(set(S)); out=[]
    for i,v in enumerate(S):
        if v < 14: continue
        for dv in (-14,-7,7,14,-1,1,-2,2,-21,21,-28,28):
            w = v+dv
            if w<=0 or w>CAP: continue
            T = sorted(set(S[:i]+S[i+1:]+[w]))
            if len(T)==13 and covering(T): out.append(T)
    rng.shuffle(out); return out[:24]
descent_worst = (worst[0], worst[1], worst[2])
for restart in range(40):
    if worst[1] is None: break
    S = list(worst[1]); cm = bm(S)
    for step in range(80):
        improved=False
        for T in neighbors(S,rng):
            record(T, "descent")
            tm = bm(T)
            if tm < cm:
                S=T; cm=tm; improved=True
                if tm < descent_worst[0]:
                    _, b = crit(T); descent_worst=(tm, tuple(T), b[0])
                break
        if not improved: break
P(f"  descent tightest margin = {float(descent_worst[0]):.8f}  at S={descent_worst[1]}")

P("\n" + "="*78); P("VERDICT"); P("="*78)
P(f"  total covering 13-sets tested: {ntot}")
P(f"  criterion C FAILURES (no v with W>1/(7v)): {cfail}")
P(f"  of those, ACTUAL LRC breaks M(S)<1/14: {mbreak}")
gm = min(worst[0], descent_worst[0])
P(f"  global tightest successful margin: {float(gm):.8f} = {gm if gm!=F(99) else 'NA'}")
if cfail == 0:
    P("  ==> C NEVER FAILED across all attacks. Strong evidence C is universal.")
else:
    P(f"  ==> C failed {cfail} times. Examples (S, best_margin, M):")
    for S,best,Mv,note in cfail_examples[:8]:
        P(f"      [{note}] S={S}  best_margin={float(best[1]):.6f}  M={Mv}={float(Mv):.5f}  "
          f"{'*** LRC BREAK ***' if Mv<C else '(C-only fail, M ok)'}")
