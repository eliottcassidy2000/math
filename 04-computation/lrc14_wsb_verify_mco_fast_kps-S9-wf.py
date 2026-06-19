#!/usr/bin/env python3
"""
FAST adversarial verification (companion to lrc14_wsb_verify_momentconvexor_kps-S9-wf.py).

Strategy: a FLOAT screen of meas_S7 (cheap) to find any candidate that beats consec
or exceeds cap, then EXACT (Fraction) confirmation only on the survivors. This lets us
push the search span and trial counts far higher than the all-exact engine.

We re-verify the float engine against the exact one on consec to calibrate, then hunt:
  B2  : exhaustive bounded primitive sweep, k=8 to span 30, k=9 to 24, k=10 to 20 (FLOAT screen)
  B3  : wide / resonant(w=0 mod7) / short-relation random (FLOAT screen, EXACT confirm)
  C   : Lemma A  P(N=6) <= 1/(7(k-1))  + component-count bound  (EXACT, cheaper than full p_vec)
  E   : k=12 Lemma-B failure witness (EXACT)
All candidate beaters / over-caps are RE-CONFIRMED exactly before being reported.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, comb
import random, sys
sys.stdout.reconfigure(line_buffering=True)

# ---------------- FLOAT engine (fast screen) ----------------
def meas_S7_float(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    if not Enz: return 0.0
    bps = {0.0, 1.0}
    for e in Enz:
        for m in range(e):
            for j in range(7):
                bps.add((7 * m + j) / (7.0 * e))
    bps = sorted(b for b in bps if 0.0 <= b <= 1.0)
    tot = 0.0
    for a, b in zip(bps, bps[1:]):
        if b <= a: continue
        xm = 0.5 * (a + b)
        hit = set()
        for e in E:
            fx = (e * xm) % 1.0
            s = int(fx * 7)
            if s == 7: s = 6
            hit.add(s)
            if len(hit) == 7: break
        if len(hit) == 7:
            tot += (b - a)
    return tot

# ---------------- EXACT engine (confirm) ----------------
def frac(x):
    r = x - int(x)
    return r + 1 if r < 0 else r

def cell_iter(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(e):
            for j in range(7):
                bps.add(F(7 * m + j, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    for a, b in zip(bps, bps[1:]):
        if b <= a: continue
        xm = (a + b) / 2
        hit = set()
        for e in E:
            s = int(frac(e * xm) * 7)
            hit.add(6 if s == 7 else s)
        yield (b - a, hit)

def meas_S7_exact(E):
    return sum((L for L, hit in cell_iter(E) if len(hit) == 7), F(0))

def P_N6_exact(E):
    tot = F(0)
    for L, hit in cell_iter(E):
        if set(s for s in hit if s != 0) == set():
            tot += L
    return tot

def n_components_GE(E):
    cells = [(L, set(s for s in hit if s != 0) == set()) for L, hit in cell_iter(E)]
    runs = 0; prev = False
    for L, cur in cells:
        if cur and not prev: runs += 1
        prev = cur
    if cells and cells[0][1] and cells[-1][1] and runs >= 2:
        runs -= 1
    return runs

def is_primitive(E):
    g = 0
    for e in E:
        if e: g = gcd(g, e)
    return g == 1

cap = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91), 11: F(66, 91), 12: F(6, 7)}
capf = {k: float(v) for k, v in cap.items()}

# calibrate float vs exact on consec
print("[calibration float vs exact meas_S7 on consec]")
for k in [8, 9, 10]:
    fe = meas_S7_float(list(range(k))); ex = float(meas_S7_exact(list(range(k))))
    print(f"  k={k}: float={fe:.9f} exact={ex:.9f} diff={abs(fe-ex):.2e}")

EPS = 1e-9
print("\n" + "=" * 72)
print("B2 — wide exhaustive bounded primitive sweep (FLOAT screen, EXACT confirm beaters)")
print("=" * 72)
def sweep(k, B):
    cvf = meas_S7_float(list(range(k)))
    cve = meas_S7_exact(list(range(k)))
    bestf = cvf; argf = tuple(range(k))
    beaters = []; overcaps = []
    cnt = 0
    for rest in combinations(range(1, B + 1), k - 1):
        E = (0,) + rest
        if not is_primitive(E): continue
        cnt += 1
        v = meas_S7_float(list(E))
        if v > bestf + EPS:
            bestf = v; argf = E
        if v > cvf + EPS:
            beaters.append(E)
        if v > capf[k] + EPS:
            overcaps.append(E)
    # exact-confirm the float-flagged candidates
    real_beaters = [E for E in beaters if meas_S7_exact(list(E)) > cve]
    real_over = [E for E in overcaps if meas_S7_exact(list(E)) > cap[k]]
    return cnt, bestf, argf, len(beaters), len(real_beaters), real_beaters[:5], len(real_over), real_over[:5]

for k, B in [(8, 20), (9, 18), (10, 16), (11, 15)]:
    cnt, bestf, argf, nbf, nbr, rb, nor, ro = sweep(k, B)
    print(f" k={k} B={B} prim#={cnt}: consec={float(meas_S7_exact(list(range(k)))):.6f} "
          f"SUPfloat={bestf:.6f} arg={argf} float-beaters={nbf} EXACT-beaters={nbr} "
          f"EXACT-overcap={nor}")
    if rb: print("    CONFIRMED BEATERS:", rb)
    if ro: print("    CONFIRMED OVERCAP:", ro)

print("\n" + "=" * 72)
print("B3 — wide / resonant(w=0 mod7) / short-relation RANDOM hunt")
print("=" * 72)
random.seed(20260619)

def gen_wide(k):
    span = random.randint(20, 300)
    s = {0}
    while len(s) < k: s.add(random.randint(1, span))
    return sorted(s)

def gen_resonant(k):
    # apex-prime-7 resonance: many multiples of 7 + a few strangers, then made primitive by a stranger
    nres = random.randint(k - 3, k - 1)
    s = {0}
    while len([x for x in s if x % 7 == 0]) < nres:
        s.add(7 * random.randint(1, 25))
    while len(s) < k:
        x = random.randint(1, 250)
        if x % 7 != 0: s.add(x)
    return sorted(s)[:k]

def gen_shortrel(k):
    N = random.randint(5, 150); s = {0}; j = 0
    while len(s) < k:
        s.add(j * N); s.add(j * N + 1); j += 1
    return sorted(s)[:k]

def hunt(name, gen, k, trials):
    cve = meas_S7_exact(list(range(k))); cvf = float(cve)
    bestf = cvf; argf = tuple(range(k)); flagged = []; over_flag = []
    seen = set()
    for _ in range(trials):
        E = tuple(gen(k))
        if E in seen: continue
        seen.add(E)
        if not is_primitive(E): continue
        v = meas_S7_float(list(E))
        if v > bestf + EPS: bestf = v; argf = E
        if v > cvf + EPS: flagged.append(E)
        if v > capf[k] + EPS: over_flag.append(E)
    rb = [E for E in flagged if meas_S7_exact(list(E)) > cve]
    ro = [E for E in over_flag if meas_S7_exact(list(E)) > cap[k]]
    print(f"  [{name}] k={k} tried={len(seen)}: consec={cvf:.6f} SUPfloat={bestf:.6f} cap={capf[k]:.6f} "
          f"float-beat={len(flagged)} EXACT-beat={len(rb)} EXACT-overcap={len(ro)}")
    if rb: print("      CONFIRMED BEATER:", rb[:3], "vals:", [float(meas_S7_exact(list(E))) for E in rb[:3]])
    if ro: print("      CONFIRMED OVERCAP:", ro[:3])

for k in [8, 9, 10]:
    hunt("wide", gen_wide, k, 40000)
    hunt("resonant7", gen_resonant, k, 25000)
    hunt("short-rel", gen_shortrel, k, 25000)

print("\n" + "=" * 72)
print("C — LEMMA A: P(N=6) <= 1/(7(k-1)) for primitive E + component-count bound")
print("=" * 72)
def lemmaA(k, B):
    bound = F(1, 7 * (k - 1)); maxv = F(0); arg = None
    nviol = 0; cc_viol = 0; cc_ex = []; cnt = 0
    for rest in combinations(range(1, B + 1), k - 1):
        E = (0,) + rest
        if not is_primitive(E): continue
        cnt += 1
        pn6 = P_N6_exact(list(E))
        if pn6 > maxv: maxv = pn6; arg = E
        if pn6 > bound: nviol += 1
        nc = n_components_GE(list(E))
        if F(nc) > F(max(E), k - 1):
            cc_viol += 1; cc_ex.append((nc, max(E), E))
    return bound, maxv, arg, nviol, cc_viol, cc_ex, cnt

for k, B in [(8, 22), (9, 18), (10, 16), (11, 15)]:
    bound, maxv, arg, nviol, ccv, cce, cnt = lemmaA(k, B)
    print(f" k={k} B={B} prim#={cnt}: bound={bound}={float(bound):.6f} maxP(N=6)={maxv}={float(maxv):.6f} "
          f"arg={arg} #viol_len={nviol} #viol_compcount={ccv}")
    if cce: print("    COMP-COUNT VIOLATIONS:", cce[:5])

# wide Lemma A
print("\n [Lemma A wide/resonant hunt]")
random.seed(99)
for k in [8, 9, 10, 11]:
    bound = F(1, 7 * (k - 1)); maxv = F(0); arg = None; nviol = 0; ccv = 0; cce = []
    for _ in range(20000):
        E = gen_wide(k)
        if not is_primitive(E): continue
        pn6 = P_N6_exact(E)
        if pn6 > maxv: maxv = pn6; arg = tuple(E)
        if pn6 > bound: nviol += 1
        nc = n_components_GE(E)
        if F(nc) > F(max(E), k - 1): ccv += 1; cce.append((nc, max(E), tuple(E)))
    print(f"  k={k}: bound={float(bound):.6f} maxP(N=6)={float(maxv):.6f} arg={arg} #viol_len={nviol} #viol_cc={ccv}")
    if cce: print("     CC VIOL:", cce[:3])

print("\n" + "=" * 72)
print("E — k=12 Lemma-B failure witness")
print("=" * 72)
E12 = list(range(11)) + [12]
v12 = meas_S7_exact(E12); c12 = meas_S7_exact(list(range(12)))
print(f"  E=[0..10,12]: meas_S7={v12}={float(v12):.6f} consec_12={float(c12):.6f} "
      f"beats_consec={v12>c12} over_cap12={v12>cap[12]}")

print("\nDONE")
