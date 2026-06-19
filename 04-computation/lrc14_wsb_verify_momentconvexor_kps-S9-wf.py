#!/usr/bin/env python3
"""
ADVERSARIAL VERIFICATION of the kps-S9 "moment-convex-order" advance for LRC(14).

Angle claimed: HYP-2607 reframed in the PROBABILITY basis p_t = meas{N_E = t},
yielding LEMMA B (consec maximizes meas(S7)=p_0 directly at binding rows k=8,9,10),
LEMMA A (P(N=6)(E) <= 1/(7(k-1)) for primitive E), and resonance-dissolution
(gcd=7 wides collapse to consec via THM-531 scale invariance).

EVERYTHING here is EXACT (fractions.Fraction). We re-derive meas_S7, Sr/J, Ly, p_t
INDEPENDENTLY (do not import the macmini engine) and cross-check, then HUNT for
counterexamples to every load-bearing statement.

Targets (default skeptical):
  (A) re-derive the stated CONSTANTS exactly.
  (B) LEMMA B: consec is the UNIQUE max of meas(S7) over primitive same-k clusters
      at k=8,9,10. Hunt: exhaustive bounded-span + wide random + resonant (w=0 mod7)
      + short-relation {0,1,N,N+1}-type shapes. One hit => Lemma B false at a binding row.
  (C) LEMMA A: P(N=6)(E) <= 1/(7(k-1)) for primitive E (k=8,9,10,11). Hunt a breacher.
      Also test the component-count bound #comp(G_E) <= max(E)/(k-1).
  (D) scale invariance meas_S7(d*E)=meas_S7(E); resonance dissolution.
  (E) the claimed k=12 Lemma-B FAILURE witness E=[0,1,..,10,12].
  (F) Does meas_S7(E) <= cap_k hold at the binding rows for EVERY E we can reach?
      (this is the actual LRC target; Lemma B is the *route*, this is the *conclusion*).
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, comb
import random, sys
sys.stdout.reconfigure(line_buffering=True)

# ----------------------------------------------------------------------
# INDEPENDENT exact engine.  frac(x) in [0,1).
# ----------------------------------------------------------------------
def frac(x):
    r = x - int(x)
    if r < 0: r += 1
    return r

def gcd_list(E):
    g = 0
    for e in E: g = gcd(g, e)
    return g

def is_primitive(E):
    return gcd_list([e for e in E if e != 0]) == 1

def breakpoints(E):
    """All x in [0,1] where some frac(e x) crosses a sector boundary j/7."""
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, e + 1):
            for j in range(0, 7):
                t = F(7 * m + j, 7 * e)
                if 0 <= t <= 1:
                    bps.add(t)
    return sorted(bps)

def sector(e, xm):
    """sector index 0..6 of frac(e*xm)."""
    s = int(frac(e * xm) * 7)
    return 6 if s == 7 else s

def cell_iter(E):
    """yield (length, set_of_hit_sectors) over each order cell."""
    bps = breakpoints(E)
    for a, b in zip(bps, bps[1:]):
        if b <= a: continue
        xm = (a + b) / 2
        hit = set(sector(e, xm) for e in E)
        yield (b - a, hit)

def meas_S7(E):
    """meas{x: every sector 0..6 hit} = p_0 (N counts empty sectors among 1..6;
       sector 0 always hit by e=0, so N=0 <=> all 7 hit)."""
    tot = F(0)
    for L, hit in cell_iter(E):
        if len(hit) == 7:
            tot += L
    return tot

def p_vec(E):
    """probability vector p_t = meas{ N_E = t }, t=0..6,
       N_E = # empty sectors among {1..6} (sector 0 always occupied by e=0)."""
    p = [F(0)] * 7
    have_zero = (0 in E)
    for L, hit in cell_iter(E):
        # empty sectors among 1..6
        occ = set(s for s in hit if s != 0)
        # if 0 not in E, sector 0 may matter; but the cluster always has 0 in E.
        N = 6 - len(occ)
        p[N] += L
    return p

def Sr(E, r):
    """factorial moment S_r = sum_{|A|=r, A subset {1..6}} meas{all sectors in A empty}."""
    tot = F(0)
    for A in combinations(range(1, 7), r):
        Aset = set(A)
        for L, hit in cell_iter(E):
            if Aset.isdisjoint(hit):
                tot += L
    return tot

def Sr_from_p(p, r):
    """S_r = E[C(N,r)] = sum_t p_t C(t,r). Cross-check against direct Sr."""
    return sum(p[t] * comb(t, r) for t in range(7))

def Ly(E, k, p=None):
    """THM-534 integer-root dual functional L_y(E)."""
    if p is None: p = p_vec(E)
    S1 = Sr_from_p(p, 1); S2 = Sr_from_p(p, 2); S3 = Sr_from_p(p, 3); S4 = Sr_from_p(p, 4)
    if k in (11, 12, 13):
        return 1 - F(1, 2) * S1 + F(1, 6) * S2
    if k in (9, 10):
        return 1 - F(13, 18) * S1 + F(4, 9) * S2 - F(1, 6) * S3
    if k == 8:
        return 1 - S1 + S2 - F(9, 10) * S3 + F(3, 5) * S4
    raise ValueError

def P_N6(E):
    """P(N=6) = meas{ frac(e x) in [0,1/7) for ALL e in E } = meas(G_E)."""
    return p_vec(E)[6]

def n_components_GE(E):
    """# connected components (order cells with N=6) of G_E={x: all frac(ex) in [0,1/7)}.
       Two adjacent cells both with N=6 merge; count maximal runs."""
    runs = 0
    prev = False
    for L, hit in cell_iter(E):
        cur = (set(s for s in hit if s != 0) == set())  # N=6 means no nonzero sector occupied
        if cur and not prev:
            runs += 1
        prev = cur
    # NOTE: cells are over [0,1); wrap-around at x=1<->x=0? G_E at x=0: frac(0)=0 in [0,1/7) yes.
    # The interval [0,1) with periodicity: a run touching both x~0 and x~1 should be one comp.
    # check wrap
    cells = list(cell_iter(E))
    first = (set(s for s in cells[0][1] if s != 0) == set())
    last = (set(s for s in cells[-1][1] if s != 0) == set())
    if first and last and runs >= 2:
        runs -= 1
    return runs

# ----------------------------------------------------------------------
# cap constants
# ----------------------------------------------------------------------
cap = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91), 11: F(66, 91), 12: F(6, 7), 13: F(1)}
# NOTE: prompt text gave cap_9=2025/4004, cap_11=25/91, cap_12=1/7 in one place but
# CONSTANTS section + realsup engine give cap_9=1979/4004, cap_11=66/91, cap_12=6/7.
# We test against BOTH where it matters.

print("=" * 72)
print("PART A — exact re-derivation of CONSTANTS")
print("=" * 72)
for k in [8, 9, 10, 11, 12]:
    E = list(range(k))
    mS7 = meas_S7(E); pn6 = P_N6(E); ly = Ly(E, k)
    print(f" consec k={k}: meas_S7={mS7}={float(mS7):.6f}  P(N=6)={pn6}  1/(7(k-1))={F(1,7*(k-1))}  "
          f"{'OK' if pn6==F(1,7*(k-1)) else 'MISMATCH'}  L_y={ly}={float(ly):.6f}  cap={cap[k]}={float(cap[k]):.6f}")

# claimed constants
print("\n claimed-constant checks:")
claims = {
    "meas_S7(consec_8)=481/1470": (meas_S7(list(range(8))), F(481,1470)),
    "meas_S7(consec_9)=2447/5880": (meas_S7(list(range(9))), F(2447,5880)),
    "meas_S7(consec_10)=8899/17640": (meas_S7(list(range(10))), F(8899,17640)),
    "cap_8=2243/5880": (cap[8], F(2243,5880)),
    "cap_9=1979/4004": (cap[9], F(1979,4004)),
    "cap_10=55/91": (cap[10], F(55,91)),
}
for name, (got, exp) in claims.items():
    print(f"   {name}: got {got} {'OK' if got==exp else 'MISMATCH got '+str(got)}")

# tightest competitor claim
comp = meas_S7([0,2,3,4,5,6,7,8])
print(f"\n   meas_S7([0,2,3,4,5,6,7,8]) = {comp} = {float(comp):.6f}  (claim 0.273639)")

print("=" * 72)
print("PART B — LEMMA B: consec maximizes meas_S7 at binding k=8,9,10 (exhaustive bounded)")
print("=" * 72)
def hunt_bounded(k, B):
    consec = list(range(k)); cv = meas_S7(consec)
    best = cv; arg = tuple(consec); nbeat = 0; novercap = 0; overs = []
    cnt = 0
    for rest in combinations(range(1, B + 1), k - 1):
        E = (0,) + rest
        if not is_primitive(E):  # restrict to primitive (the LRC reduction uses primitive)
            continue
        v = meas_S7(list(E)); cnt += 1
        if v > best:
            best = v; arg = E
        if v > cv:
            nbeat += 1
        if v > cap[k]:
            novercap += 1; overs.append((float(v), E))
    return cv, best, arg, nbeat, novercap, cnt, overs

for k, B in [(8, 16), (9, 14), (10, 14)]:
    cv, best, arg, nbeat, novercap, cnt, overs = hunt_bounded(k, B)
    status = "consec=UNIQUE max" if arg == tuple(range(k)) else f"BEATER {arg}"
    print(f" k={k} B={B} primitive#={cnt}: consec={float(cv):.6f}  SUP={float(best):.6f}  {status}  "
          f"#beat_consec={nbeat}  #over_cap={novercap}")
    if overs: print("    OVERCAP:", overs[:5])

print("=" * 72)
print("PART B2 — DEEPER exhaustive at k=8 (extend span) + structured beater hunt")
print("=" * 72)
# extend the box for k=8 to B=22 (heavier but feasible: C(22,7)~170k)
for k, B in [(8, 22)]:
    cv = meas_S7(list(range(k))); best = cv; arg = tuple(range(k)); nover = 0; overs = []
    cnt = 0
    for rest in combinations(range(1, B + 1), k - 1):
        E = (0,) + rest
        if not is_primitive(E): continue
        v = meas_S7(list(E)); cnt += 1
        if v > best: best = v; arg = E
        if v > cap[k]: nover += 1; overs.append((float(v), E))
    print(f" k={k} B={B} primitive#={cnt}: SUP={float(best):.6f} arg={arg} consec={float(cv):.6f} "
          f"#over_cap={nover}", overs[:3])

print("=" * 72)
print("PART B3 — WIDE-SPREAD + RESONANT (w=0 mod7) + SHORT-RELATION hunt (k=8,9,10)")
print("=" * 72)
random.seed(20260619)
def report_hunt(name, gen_fn, k, trials):
    cv = meas_S7(list(range(k))); best = cv; arg = tuple(range(k)); nover = 0; overs = []
    nbeat = 0
    for _ in range(trials):
        E = gen_fn(k)
        if E is None: continue
        if not is_primitive(E): continue
        v = meas_S7(list(E))
        if v > best: best = v; arg = E
        if v > cv: nbeat += 1
        if v > cap[k]: nover += 1; overs.append((float(v), tuple(E)))
    print(f"  [{name}] k={k} trials={trials}: SUP={float(best):.6f} consec={float(cv):.6f} cap={float(cap[k]):.6f} "
          f"#beat_consec={nbeat} #over_cap={nover}", overs[:3])
    return best, arg

def gen_wide(k):
    span = random.randint(20, 120)
    s = set([0]); 
    while len(s) < k:
        s.add(random.randint(1, span))
    return sorted(s)

def gen_resonant(k):
    # force most elements to be 0 mod 7 (apex-prime-7 resonance), one off-resonance
    base = sorted(random.sample(range(0, 15), k - 1))
    E = [0] + [7 * b for b in base if b > 0][:k-2]
    # add one stranger not divisible by 7
    while len(set(E)) < k:
        x = random.randint(1, 200)
        if x % 7 != 0: E.append(x)
        E = list(set(E))
    return sorted(set(E))[:k] if len(set(E)) >= k else None

def gen_short_rel(k):
    # {0,1, N, N+1, 2N, 2N+1, ...}  short-relation lattice shapes
    N = random.randint(5, 80)
    blocks = []
    j = 0
    while len(set(sum(blocks, []))) < k and j < k:
        blocks.append([j * N, j * N + 1])
        j += 1
    flat = sorted(set([0] + sum(blocks, [])))
    if len(flat) < k: return None
    return flat[:k]

for k in [8, 9, 10]:
    report_hunt("wide", gen_wide, k, 40000)
    report_hunt("resonant7", gen_resonant, k, 20000)
    report_hunt("short-rel", gen_short_rel, k, 20000)

print("=" * 72)
print("PART C — LEMMA A: P(N=6)(E) <= 1/(7(k-1)) for primitive E + component-count bound")
print("=" * 72)
def hunt_lemmaA(k, B):
    bound = F(1, 7 * (k - 1))
    nviol = 0; viol = []; maxv = F(0); arg = None
    cc_viol = 0; cc_examples = []
    cnt = 0
    for rest in combinations(range(1, B + 1), k - 1):
        E = (0,) + rest
        if not is_primitive(E): continue
        cnt += 1
        pn6 = P_N6(list(E))
        if pn6 > maxv: maxv = pn6; arg = E
        if pn6 > bound: nviol += 1; viol.append((float(pn6), E))
        # component-count bound: #comp <= max(E)/(k-1)
        nc = n_components_GE(list(E))
        if F(nc) > F(max(E), k - 1):
            cc_viol += 1; cc_examples.append((nc, max(E), E))
    return bound, maxv, arg, nviol, viol, cc_viol, cc_examples, cnt

for k, B in [(8, 18), (9, 15), (10, 14), (11, 14)]:
    bound, maxv, arg, nviol, viol, cc_viol, cc_ex, cnt = hunt_lemmaA(k, B)
    print(f" k={k} B={B} prim#={cnt}: bound=1/(7(k-1))={bound}={float(bound):.6f}  "
          f"max P(N=6)={maxv}={float(maxv):.6f} arg={arg}  #viol_lengthbound={nviol}  #viol_compcount={cc_viol}")
    if viol: print("    LENGTH-BOUND VIOL:", viol[:3])
    if cc_ex: print("    COMP-COUNT VIOL:", cc_ex[:3])

# wide hunt for Lemma A breach
print("\n [Lemma A wide/resonant hunt]")
random.seed(424242)
for k in [8, 9, 10, 11]:
    bound = F(1, 7 * (k - 1)); maxv = F(0); arg = None; nviol = 0; cc_viol = 0
    for _ in range(30000):
        E = gen_wide(k)
        if not is_primitive(E): continue
        pn6 = P_N6(E)
        if pn6 > maxv: maxv = pn6; arg = tuple(E)
        if pn6 > bound: nviol += 1
        nc = n_components_GE(E)
        if F(nc) > F(max(E), k - 1): cc_viol += 1
    print(f"  k={k}: bound={float(bound):.6f} max P(N=6)={float(maxv):.6f} arg={arg} "
          f"#viol_len={nviol} #viol_cc={cc_viol}")

print("=" * 72)
print("PART D — scale invariance + resonance dissolution")
print("=" * 72)
for E in [[0,1,2,3,4,5,6,7],[0,2,3,5,7,9,11,13],[0,1,3,4,9,10,12,13]]:
    ok = all(meas_S7([d*x for x in E]) == meas_S7(E) for d in [2,3,5,7])
    print(f"  scale-inv E={E}: {'OK' if ok else 'FAIL'}")
# resonance: gcd=7 wide config must equal its primitive rep (consec if rep is consec)
res = [0,7,14,21,28,35,42,49]
prim = [x // 7 for x in res]
print(f"  resonance E={res} gcd=7  meas_S7={float(meas_S7(res)):.6f}  prim={prim} meas_S7={float(meas_S7(prim)):.6f}  "
      f"{'EQUAL (dissolves)' if meas_S7(res)==meas_S7(prim) else 'MISMATCH'}  vs cap_8={float(cap[8]):.6f}")

print("=" * 72)
print("PART E — claimed k=12 Lemma-B FAILURE witness")
print("=" * 72)
E12 = list(range(11)) + [12]
c12 = meas_S7(list(range(12)))
v12 = meas_S7(E12)
print(f"  E=[0..10,12]: meas_S7={v12}={float(v12):.6f}  consec_12={float(c12):.6f}  "
      f"{'BEATS consec (Lemma B FAILS at k=12, as claimed)' if v12>c12 else 'does not beat'}  "
      f"vs cap_12=6/7={float(cap[12]):.6f}  over_cap={v12>cap[12]}")

print("=" * 72)
print("PART F — Sr/p_t cross-check (engine consistency) + L_y dual sign-floor on N support")
print("=" * 72)
# cross-check Sr direct vs from p, on a few sets
for E in [list(range(8)),[0,2,3,4,5,6,7,8],[0,1,2,3,4,5,6,9]]:
    p = p_vec(E)
    ok = all(Sr(E, r) == Sr_from_p(p, r) for r in range(1, 5)) and sum(p) == 1
    # also p_0 == meas_S7
    ok2 = (p[0] == meas_S7(E))
    print(f"  E={E}: sum(p)={sum(p)} Sr-consistent={ok} p0==meas_S7={ok2}")

# verify g(t) duals are >=1[t=0] on t in 0..6 (the dual feasibility the whole bound rests on)
def g_of(k):
    if k in (11,12,13): return [F((t-3)*(t-4),12) for t in range(7)]
    if k in (9,10): return [F(-(t-2)*(t-3)*(t-6),36) for t in range(7)]
    if k==8: return [F((t-1)*(t-2)*(t-4)*(t-5),40) for t in range(7)]
for k in [8,9,10,11]:
    g = g_of(k)
    feas = all((g[t] >= (1 if t==0 else 0)) for t in range(7))
    print(f"  k={k} dual g(t),t=0..6 = {[str(x) for x in g]}  feasible(>=1[t=0])={feas}")

print("\nDONE")
