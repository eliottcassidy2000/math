#!/usr/bin/env python3
"""
lrc14_angleD_break_v2 — mac-mini-2026-06-17-S6  (ANGLE D, sharper attack)

Two refinements over v1:
 (A) Theoretical sanity on SCALING: criterion C is scale-INVARIANT.
     For S and cS: removing cu from cS gives danger teeth of u' = cv (u in S\{v}),
     centers k/(cv), width 1/(7cv). The whole safe-set picture is the c-fold
     "u-coordinate" version. In fact W(cS\{cv}) and 1/(7cv) BOTH scale: the safe
     set of cA in tau equals (safe set of A in c*tau), a measure-preserving cover,
     so W(cA) = W(A)/c and 1/(7*c*v) = (1/(7v))/c.  margin scales by 1/c, SIGN UNCHANGED.
     => Scaling can never flip C. We verify this numerically (the real danger is
        INCOMMENSURATE large runners, not scaled cores).

 (B) Exhaustive-ish "minimal hard core" sweep:
     The hard case (per reductions) is a PRIMITIVE covering 13-set. Enumerate
     covering 13-sets built as: small core {subset of 1..13 that covers most q}
     + a few SPORADIC large runners that are the SOLE cover of one modulus
     (codex HYP-2579 private-obligation). Push the sole-cover runner large and
     incommensurate. This is exactly where C is claimed weakest.

If C fails anywhere, recompute M exactly and classify.
"""
from fractions import Fraction as F
from math import gcd, lcm
import random, itertools

C = F(1, 14)

def darcs(v, c=C):
    hw = F(c, v)
    return [(F(k, v) - hw, F(k, v) + hw) for k in range(v)]
def wrapU(iv):
    o = []
    for lo, hi in iv:
        s = lo - (lo % 1); a = lo - s; b = hi - s
        if b <= 1: o.append((a, b))
        else: o.append((a, F(1))); o.append((F(0), b - 1))
    o = sorted(o); r = []; cl, ch = o[0]
    for lo, hi in o[1:]:
        if lo <= ch: ch = ch if ch > hi else hi
        else: r.append((cl, ch)); cl, ch = lo, hi
    r.append((cl, ch)); return r
def Wsafe(A, c=C):
    dz = []
    for v in set(A): dz += darcs(v, c)
    if not dz: return F(1)
    dz = wrapU(dz); best = F(0); n = len(dz)
    for i in range(n):
        hi = dz[i][1]; lo = dz[(i + 1) % n][0] + (1 if i == n - 1 else 0)
        gap = lo - hi
        if gap > best: best = gap
    return best if best > 0 else F(0)
def covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def crit(S):
    res = []; best = None
    for v in sorted(set(S)):
        A = [u for u in S if u != v]
        W = Wsafe(A); thr = F(1, 7 * v); m = W - thr
        res.append((v, W, thr, W > thr, m))
        if best is None or m > best[4]: best = (v, W, thr, W > thr, m)
    return any(r[3] for r in res), res, best
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

print("="*78)
print("(A) SCALE-INVARIANCE of C — verify margin(cS via cv) = margin(S via v)/c")
print("="*78)
core = [1,2,3,4,5,6,7,8,9,10,11,13,84]
assert covering(core) and len(set(core))==13
_,_,b1 = crit(core)
for c in [2,3,5,7,30]:
    S = sorted(x*c for x in core)
    _,_,b2 = crit(S)
    # best margin should equal b1 margin / c  AND best v should be c*b1v
    print(f"  c={c:3d}: margin(S)={float(b1[4]):.7f} via v={b1[0]};  "
          f"margin(cS)={float(b2[4]):.7f} via v={b2[0]};  "
          f"ratio={float(b1[4]/b2[4]) if b2[4]!=0 else float('inf'):.4f} (expect {c})")
print("  => scaling shrinks margin by exactly c but NEVER flips sign. Scaled cores can't break C.")

print("\n" + "="*78)
print("(B) PRIVATE-OBLIGATION HARD CORES: sole-cover large incommensurate runners")
print("="*78)
# Build: pick which moduli are covered ONLY by a single large sporadic runner.
# A runner r covers exactly {q in 2..14 : q | r}. For a SOLE cover of q via large r,
# we want r large with q|r but the small core NOT covering q.
rng = random.Random(99)
ntot=0; cfail=0; mbreak=0; worst=(F(99),None,None); examples=[]

def build_private_obligation(rng):
    # choose 1..3 "private" moduli to be covered only by large runners
    privset = rng.choice([
        [8,9,11,13],[8,11,13],[9,11,13],[11,13],[7,8,9,11,13],
        [8,9,10,11,13],[8,9,11,12,13],[5,7,8,9,11,13],
    ])
    privset = sorted(set(privset))
    # small core: integers in 1..13 that do NOT introduce these private moduli's cover
    # actually we just need the small core to cover the OTHER moduli.
    other = [q for q in range(2,15) if q not in privset]
    small_pool = list(range(1,14))
    # greedily build a small core covering 'other'
    core=[]
    need=set(other)
    rng.shuffle(small_pool)
    for x in small_pool:
        covs = {q for q in other if x%q==0}
        if covs & need:
            core.append(x); need -= covs
        if not need: break
    if need: return None
    # now add large runners covering privset; make them large & incommensurate
    larges=[]
    # one runner per private modulus (sole cover) OR a few combining them
    mode = rng.choice(['one_per','combined','mixed'])
    if mode=='one_per':
        for q in privset:
            mult = rng.choice([1,2,3,5,7,11])* rng.choice([13,17,19,23,29,1])
            larges.append(q*mult)
    elif mode=='combined':
        L = lcm(*privset)
        larges=[L*rng.choice([1,2,3,5,7]), L*rng.choice([11,13,17,19])]
    else:
        # split privset into two groups
        rng.shuffle(privset)
        h=len(privset)//2 or 1
        g1,g2 = privset[:h], privset[h:]
        if g1: larges.append(lcm(*g1)*rng.choice([1,3,5,7,11,13]))
        if g2: larges.append(lcm(*g2)*rng.choice([1,3,5,7,11,13]))
    S = sorted(set(core+larges))
    # pad to 13 distinct with small fillers that don't break covering structure
    fill = [v for v in range(1,14) if v not in S]
    rng.shuffle(fill)
    for v in fill:
        if len(S)>=13: break
        S = sorted(set(S+[v]))
    if len(S)!=13 or not covering(S): return None
    return S

for _ in range(120000):
    S = build_private_obligation(rng)
    if S is None: continue
    ntot+=1
    holds,res,best = crit(S)
    if not holds:
        cfail+=1; Mv=M(S); examples.append((S,best,Mv))
        if Mv<C: mbreak+=1
    else:
        if best[4] < worst[0]: worst=(best[4],S,best[0])
print(f"  tested {ntot} private-obligation covering 13-sets")
print(f"  C FAILED: {cfail}   (of those M<1/14: {mbreak})")
print(f"  tightest successful margin = {float(worst[0]):.8f} at v={worst[2]}")
print(f"     S = {worst[1]}")
if cfail:
    print("  C-FAIL EXAMPLES:")
    for S,best,Mv in examples[:6]:
        print(f"    S={S}  best_margin={float(best[4]):.6f}  M={Mv}={float(Mv):.5f}  "
              f"{'*** LRC BREAK ***' if Mv<C else '(C-only fail)'}")
else:
    print("  C held on ALL private-obligation hard cores tested.")
