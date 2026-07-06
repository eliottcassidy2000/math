#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S5 -- HYP-4272: THE COUPLED 2-TORUS BRACKET CENSUS for (A).

(A) [the S4 J-K reduction]: no coupled proper 2-subtorus of (R/Z)^12 has
M(U) in (1/13, 2/25].  For a coupling (u, v) (12 coordinate pairs):
    M(U) = max_{(t,s) in T^2} min_i || u_i t + v_i s ||.

METHOD (rigorous brackets, no exactness needed to clear the window):
  f(t,s) = min_i ||u_i t + v_i s|| is Lipschitz with constant
  L1 = max|u_i| in t and L2 = max|v_i| in s (sup-norm: |f(x)-f(y)| <=
  L1|dt| + L2|ds|).  On an N1 x N2 grid with steps h1 = 1/N1, h2 = 1/N2:
    grid_max <= M(U) <= grid_max + (L1 h1 + L2 h2)/2.
  Choose N's so the slack (L1h1 + L2h2)/2 < the distance from the bracket
  to the window edges; refine locally where the bracket straddles.
  VERDICT per U: SAFE-LOW (upper < 1/13 + eps? impossible: M >= 1/13 only
  for 12-dim... actually 2-torus values can be anything in [0, 1/2]),
  SAFE-ABOVE (lower > 2/25), IN-WINDOW-RISK (bracket intersects
  (1/13, 2/25]) -> refined; final report.

DIRECTIONS CENSUSED (u = (1..12) the AP base throughout -- the lift-limit
tori; the S4 reduction's dangerous families are exactly these):
  - v = e_i (single-coordinate rays, 12): the THM-621 ladder limits;
  - v = e_i + e_j (66): the double-lift limits (incl. the {4,6} attainer's
    direction -- S4 measured its 1-dim species at ~2/17);
  - v = consecutive blocks (len 3..6);
  - v = the 14r-ladder direction (v_i = i on {7..12}-supported patterns);
  - 200 random sparse v (entries in {0,1,2}, support 1..6).
Coupled = v not identically 0 and (u,v) rank 2 (v not a rational multiple
of u -- always true for sparse v).
"""
from fractions import Fraction as F
from math import gcd
import random, time

T0 = time.time()
def log(m=""):
    print(m, flush=True)
random.seed(5)

U = list(range(1, 13))
LO, HI = F(1, 13), F(2, 25)

def dist1(x):
    x = x - int(x)
    if x < 0:
        x += 1
    return min(x, 1 - x)

def bracket(u, v, N1=400, N2=400, depth=0):
    """returns (lower, upper) rigorous bracket for M(U)."""
    L1, L2 = max(abs(a) for a in u), max(abs(b) for b in v)
    best = 0.0
    bt = bs = 0.0
    for i1 in range(N1):
        t = i1 / N1
        # precompute u_i t fractional parts
        ut = [a * t for a in u]
        for i2 in range(N2):
            s = i2 / N2
            m = 1.0
            for k in range(12):
                d = dist1(ut[k] + v[k] * s)
                if d < m:
                    m = d
                    if m <= best:
                        break
            if m > best:
                best = m
                bt, bs = t, s
    slack = (L1 / N1 + L2 / N2) / 2
    lower, upper = best, best + slack
    # refine locally if the bracket straddles the window
    if depth < 2 and not (lower > float(HI) or upper < float(LO)):
        # local refine around (bt, bs)
        span1, span2 = 2.0 / N1, 2.0 / N2
        M1 = M2 = 200
        bb = lower
        for i1 in range(M1 + 1):
            t = bt - span1 / 2 + span1 * i1 / M1
            ut = [a * t for a in u]
            for i2 in range(M2 + 1):
                s = bs - span2 / 2 + span2 * i2 / M2
                m = 1.0
                for k in range(12):
                    d = dist1(ut[k] + v[k] * s)
                    if d < m:
                        m = d
                        if m <= bb:
                            break
                if m > bb:
                    bb = m
        lower = max(lower, bb)
        slack2 = (L1 * span1 / M1 + L2 * span2 / M2) / 2
        # the true max is either near (bt,bs) (bracketed by bb + slack2) or
        # elsewhere on the coarse grid (bracketed by best + slack); take max:
        upper = max(bb + slack2, best + slack)  # conservative
        # tighten the global upper by a finer full pass if still straddling
        if not (lower > float(HI) or upper < float(LO)):
            return bracket(u, v, N1 * 2, N2 * 2, depth + 1)
    return lower, upper

def verdict(lo, up):
    if lo > float(HI):
        return "SAFE-ABOVE"
    if up <= float(LO):
        return "BELOW-13TH"       # <= 1/13: outside the open window's interior... (G) cares about (1/13, 2/25]
    if up <= float(HI) and lo > float(LO):
        return "IN-WINDOW !!"
    return "UNRESOLVED"

log("(A) coupled 2-torus bracket census; window (1/13, 2/25] = (0.076923, 0.080000]")
cases = []
for i in range(12):
    v = [0] * 12
    v[i] = 1
    cases.append((f"e_{i+1}", v))
for i in range(12):
    for j in range(i + 1, 12):
        v = [0] * 12
        v[i] = 1
        v[j] = 1
        cases.append((f"e_{i+1}+e_{j+1}", v))
for ln in range(3, 7):
    for a in range(0, 13 - ln):
        v = [0] * 12
        for x in range(a, a + ln):
            v[x] = 1
        cases.append((f"block[{a+1}..{a+ln}]", v))
v = [0] * 12
for r in range(7, 13):
    v[r - 1] = r
cases.append(("ladder14 (v_r = r, r=7..12)", v))
for _ in range(200):
    v = [0] * 12
    for i in random.sample(range(12), random.randint(1, 6)):
        v[i] = random.choice([1, 1, 2])
    cases.append(("rand", v))

stats = {"SAFE-ABOVE": 0, "BELOW-13TH": 0, "IN-WINDOW !!": 0, "UNRESOLVED": 0}
inwin = []
unres = []
for name, v in cases:
    lo, up = bracket(U, v)
    vd = verdict(lo, up)
    stats[vd] += 1
    if vd == "IN-WINDOW !!":
        inwin.append((name, v, lo, up))
        log(f"  !! IN-WINDOW: {name} v={v} bracket=({lo:.6f},{up:.6f})")
    elif vd == "UNRESOLVED":
        unres.append((name, v, lo, up))
        log(f"  ?? UNRESOLVED: {name} v={v} bracket=({lo:.6f},{up:.6f})")
    elif name.startswith(("e_", "block", "ladder")) or (name == "e_4+e_6"):
        log(f"  {name:>22}: bracket=({lo:.6f},{up:.6f}) {vd}")
log(f"\ncensus: {len(cases)} directions; verdicts: {stats}")
log("(A) VERDICT on the structured landscape: " +
    ("NO coupled torus in the window -- (A) holds on every censused direction"
     if stats["IN-WINDOW !!"] == 0 and stats["UNRESOLVED"] == 0 else "see flagged cases"))
log(f"[t = {time.time()-T0:.0f}s]")
