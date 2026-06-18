#!/usr/bin/env python3
"""
INDEPENDENT adversarial probe #2 of the cluster-collapse-threedistance advance (LRC14 S3).

Focus: the CRUX. Does C(S) (best_margin = max_v W(S\{v})*7*v  > 1) hold for EVERY
covering+primitive S3 13-set?  If even one fails, S3 is NOT closed by this criterion.

Strategy differences vs the first script:
  1. STRUCTURE-AWARE generation: a covering set must contain a multiple of every q in 2..14.
     The 'hard' large moduli (11,13,14) force specific large speeds. We build the small part
     to cover 2..13 cheaply, then build a tight cluster whose members supply the remaining
     large divisibility (esp. 14, and 11/13 if not in small part). This yields MANY more valid
     covering S3 sets than rejection sampling.
  2. ADVERSARIAL near-tight clusters: push Vmax/Vmin toward 13 (the S2/S3 boundary) and toward
     small spreads; sweep offset patterns including ones designed to minimize the widest gap
     after deleting the cluster teeth (three-distance worst case).
  3. LARGE V0: test V0 up to ~360360 to probe the carry-phase limit regime (where margin_limit->1).
  4. Independent recompute of M(S) on every candidate where best_margin is small, to check the
     ACTUAL LRC inequality M(S)>=1/14 directly (not just the sufficient criterion C).
  5. Re-derive the v=max counterexample existence and the cluster-collapse formula exactly.
"""
from fractions import Fraction as F
from math import gcd, lcm
import random, itertools

H = F(1,14)

def nrm(x):
    r = x - int(x)
    if r < 0: r += 1
    return r if r <= F(1,2) else 1 - r

def safe_components(A, h=H):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u)
            a = (c - h/u) % 1
            b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else:
                iv.append((a, F(1))); iv.append((F(0), b))
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
        while F(2*k+1, 2*v) <= F(1,2):
            C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1,2):
                        C.add(F(k, d)); k += 1
    C.add(F(1,2)); return C

def Mval(S):
    b = F(0)
    for t in cand(S):
        v = min(nrm(x*t) for x in S)
        if v > b: b = v
    return b

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2,15))

def is_primitive(S):
    g = 0
    for v in S: g = gcd(g, v)
    return g == 1

def best_margin(S):
    S = sorted(set(S)); best = F(-1); arg = None
    for v in S:
        m = Wwidth([u for u in S if u != v]) * 7 * v
        if m > best: best = m; arg = v
    return best, arg

def case_of(S):
    S=sorted(set(S)); Vmin=S[0]; Vmax=S[-1]
    k=sum(1 for v in S if v>13)
    if k<=1: return "S1"
    if Vmax < 13*Vmin: return "S2"
    return "S3"

# ============================================================
# PART A: exact re-derivation of cluster-collapse formula + threshold
print("="*70); print("PART A: cluster-collapse formula, exact"); print("="*70)
def cc_formula(V0,Vmax):
    lo=F(1,14*V0); hi=F(13,14*Vmax)
    w=hi-lo; m=w*7*Vmax; f=(F(13)-F(Vmax,V0))/2
    return m, f, (m==f)
allok=True
for V0 in [14,15,29,100,360360]:
    for Vmax in [V0+1, V0+14, 2*V0, 10*V0, 11*V0-1, 11*V0, 12*V0]:
        if 13*V0<=Vmax: continue
        m,f,ok=cc_formula(V0,Vmax)
        allok = allok and ok
print("formula width*7*Vmax == (13-Vmax/V0)/2 for all tested:", allok)
print("threshold >1 iff Vmax<11*V0:",
      all(((cc_formula(V0,Vmax)[1]>1)==(Vmax<11*V0)) for V0 in [14,100,1000] for Vmax in range(V0+1,13*V0)))

# ============================================================
# PART B: structure-aware covering S3 generator
print(); print("="*70)
print("PART B: structure-aware adversarial hunt for C(S) failures (best_margin<=1)")
print("="*70)

def gen_covering_S3(rng, V0range, spread_max):
    """
    Build small part covering 2..13 (using a 7-9 element subset of {1..13} that hits 2..13),
    then a tight cluster near V0 that supplies a multiple of 14 (and any missing of 11,13).
    Returns a covering primitive S3 set or None.
    """
    # small part: must include enough to cover 2..13. Easy choice: pick from {1..13} to hit each.
    # Greedy: include 12 (covers 2,3,4,6,12), 8 (covers 8), 9(9), 10(2,5,10), 11(11), 13(13), 7(7).
    # That's {7,8,9,10,11,12,13} covers 2,3,4,5,6,7,8,9,10,11,12,13. size 7.
    base_small = [7,8,9,10,11,12,13]  # covers all of 2..13
    # randomly augment small with a few more from {1..6} to vary
    extra_pool = [1,2,3,4,5,6]
    nextra = rng.randint(0,4)
    small = sorted(set(base_small + rng.sample(extra_pool, nextra)))
    # cluster size = 13 - len(small) - (we'll just fill to 13 total with cluster speeds)
    csize = 13 - len(small)
    if csize < 2: return None
    V0 = rng.randint(*V0range)
    spread = rng.randint(2, spread_max)
    # cluster MUST contain a multiple of 14 (since small part has none).
    # Find a multiple of 14 in [V0, V0+spread]; if none, slide V0.
    m14 = ((V0 + 13)//14)*14
    if m14 > V0 + spread:
        V0 = m14  # anchor cluster start at the multiple of 14
    cl = set([m14 if m14<=V0+spread else ((V0+13)//14)*14])
    cl = set()
    # guarantee a mult of 14 inside [V0,V0+spread]
    cand14 = [x for x in range(V0, V0+spread+1) if x%14==0]
    if not cand14:
        return None
    cl.add(rng.choice(cand14))
    while len(cl) < csize:
        cl.add(V0 + rng.randint(0, spread))
    L = sorted(cl)
    S = sorted(set(small) | set(L))
    if len(S) != 13: return None
    if not is_primitive(S): return None
    if not is_covering(S): return None
    if case_of(S) != "S3": return None
    return S

rng = random.Random(424242)
worst = F(10**9); worst_set=None
nC_fail=0; nM_fail=0; tested=0
cfail_ex=[]; mfail_ex=[]
configs = [
    ((14,200), 45),
    ((14,500), 45),
    ((14,2000), 30),
    ((100,5000), 45),
    ((1000,50000), 45),
    ((10000,360360), 60),
    ((14,60), 14),   # near-tight S3 boundary
    ((14,40), 20),
]
for cfg in configs:
    for _ in range(40000):
        S = gen_covering_S3(rng, *cfg)
        if S is None: continue
        tested += 1
        bm, arg = best_margin(S)
        if bm < worst: worst = bm; worst_set = S[:]
        if bm <= 1:
            nC_fail += 1
            M = Mval(S)
            if len(cfail_ex) < 12: cfail_ex.append((S[:], bm, arg, M))
            if M < H:
                nM_fail += 1
                if len(mfail_ex) < 12: mfail_ex.append((S[:], M))

print(f"valid covering primitive S3 sets tested: {tested}")
print(f"C(S) failures (best_margin<=1): {nC_fail}")
print(f"TRUE M<1/14 counterexamples: {nM_fail}")
print(f"worst (min) best_margin: {worst} = {float(worst):.6f} at {worst_set}")
for s,b,a,M in cfail_ex[:12]:
    print("   C-FAIL set:", s, "bm=",b,float(b),"arg=",a,"M=",M,float(M))
for s,M in mfail_ex[:12]:
    print("   !!! M<1/14:", s, M, float(M))

# ============================================================
# PART C: directly minimize best_margin via targeted near-tight construction.
# The arc-width lemma fires through the LARGE cluster. The danger is a set where the
# widest arc of S\{v} is just under 1/(7v) for the best v. We hand-craft adversarial
# offset patterns to shrink the widest gap (pack cluster teeth to cover the safe band).
print(); print("="*70)
print("PART C: hand-crafted adversarial near-tight S3 (minimize best_margin)")
print("="*70)
# Use small part that forces tight covering, cluster with offsets spread to break the wide arc.
def build(small, V0, offsets):
    S = sorted(set(list(small) + [V0+o for o in offsets]))
    return S

# try to drive best_margin below 1.33 (their realized floor) and ideally to <=1
small_opts = [
    [7,8,9,10,11,12,13],
    [5,7,8,9,10,11,12,13],
    [6,7,8,9,10,11,12,13],
    [3,5,7,8,9,10,11,12,13],
    [1,7,8,9,10,11,12,13],
]
best_seen = F(10**9); best_seen_set=None
rng2 = random.Random(13)
for small in small_opts:
    csize = 13 - len(small)
    if csize < 2: continue
    for _ in range(60000):
        V0 = rng2.choice([14,28,42,56,140,280,1400,14014,140140,360360]) + rng2.randint(0,13)
        spread = rng2.randint(csize, 50)
        # offsets must include a multiple-of-14 anchor
        cand14 = [x for x in range(0, spread+1) if (V0+x)%14==0]
        if not cand14: continue
        offs = set([rng2.choice(cand14)])
        while len(offs) < csize:
            offs.add(rng2.randint(0, spread))
        S = build(small, V0, sorted(offs))
        if len(S) != 13: continue
        if not is_primitive(S) or not is_covering(S): continue
        if case_of(S) != "S3": continue
        bm,arg = best_margin(S)
        if bm < best_seen:
            best_seen = bm; best_seen_set = (S[:], arg, Mval(S))
print(f"min best_margin found in Part C: {best_seen} = {float(best_seen):.6f}")
if best_seen_set:
    s,a,M = best_seen_set
    print(f"   set={s} arg={a} M={M}={float(M):.6f}  M>=1/14? {M>=H}")

# ============================================================
# PART D: the v=max determinism refutation + the specific claimed numbers
print(); print("="*70)
print("PART D: v=max determinism (claimed REFUTED) + claimed constants")
print("="*70)
print("claimed v=max counterexample margin 242015/245024 =", float(F(242015,245024)), "(<1)")
print("claimed M for it 123/989 =", float(F(123,989)), " >=1/14?", F(123,989)>=H)
# Search for S3 set where v=max margin<=1 but C(S) holds via other v (v=max determinism FALSE)
rng3 = random.Random(7)
vmax_fail = 0; ex=[]
for _ in range(200000):
    S = gen_covering_S3(rng3, (14,3000), 45)
    if S is None: continue
    vmax = max(S)
    mvmax = Wwidth([u for u in S if u!=vmax])*7*vmax
    if mvmax <= 1:
        bm,arg = best_margin(S)
        vmax_fail += 1
        if len(ex)<6: ex.append((S[:], mvmax, bm, arg, Mval(S)))
print(f"S3 sets where v=max margin<=1 (determinism FALSE): {vmax_fail}")
for s,mv,bm,arg,M in ex:
    print(f"   S={s}")
    print(f"     vmax_margin={float(mv):.5f}  best_margin(other v ok)={float(bm):.5f} arg={arg}  M={M}={float(M):.6f}")
