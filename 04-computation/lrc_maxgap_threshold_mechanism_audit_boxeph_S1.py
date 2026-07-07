"""
LRC(14) density-floor MECHANISM AUDIT  (boxeph-2026-07-07-S1, HYP-4760)
Pure-Python (no numpy). Resolve the 1/7-vs-2/7 threshold question and pin down
which object is load-bearing for the density->reach reduction.

Objects:
  (A) THM-527 rho*: co-offset config {frac(e_i x)} WITH 0 (e=Vmax-u),
      circular-maxgap > 2/7, intersect G_P. rho*>0 => M(S)>=1/14 via Vmax-ruler.
  (B) recent mu_{1/7}(E): SPEED config {frac(e x)} (no forced 0), maxgap > 1/7.

Questions Q1..Q4 (see comments).
"""
from fractions import Fraction as F
import math, random

def frac(y):  # fractional part
    return y - math.floor(y)

def dist_int(y):
    f = frac(y)
    return min(f, 1.0 - f)

# ---- loneliness M(S)=max_t min_i ||v_i t|| by adaptive sampling ----
def M_of(S, N=None):
    S = sorted(set(int(v) for v in S))
    Vmax = max(abs(v) for v in S)
    if N is None:
        N = max(120000, 40 * Vmax)
    best = -1.0; bt = 0.0
    for a in range(N):
        t = (a + 0.5) / N
        m = 1.0
        for v in S:
            d = dist_int(v * t)
            if d < m: m = d
            if d < best: break
        if m > best:
            best = m; bt = t
    lo, hi = bt - 2.0/N, bt + 2.0/N
    for _ in range(4):
        for a in range(2000):
            t = lo + (hi-lo)*(a+0.5)/2000
            m = 1.0
            for v in S:
                d = dist_int(v*t)
                if d < m: m = d
            if m > best:
                best = m; bt = t
        span = (hi-lo)/2000
        lo, hi = bt-2*span, bt+2*span
    return best

# ---- circular max-gap of config {frac(e x): e in E} ----
def maxgap_at(E, x):
    pts = sorted(frac(e * x) for e in E)
    if len(pts) == 1:
        return 1.0
    mg = 0.0
    for i in range(1, len(pts)):
        g = pts[i]-pts[i-1]
        if g > mg: mg = g
    wrap = 1.0 - pts[-1] + pts[0]
    if wrap > mg: mg = wrap
    return mg

def mu_threshold(E, thr, G=40000):
    c = 0
    for a in range(G):
        x = (a+0.5)/G
        if maxgap_at(E, x) > thr: c += 1
    return c / G

def E_maxgap(E, G=40000):
    s = 0.0
    for a in range(G):
        x = (a+0.5)/G
        s += maxgap_at(E, x)
    return s / G

# ---- EXACT mu_{thr} for AP-like sets (piecewise linear, rational breakpoints) ----
def mu_exact(E, thr):
    """EXACT meas{x in (0,1): circular max-gap {frac(e x): e in E} > thr}."""
    E = list(E)
    thr = F(thr)
    denoms = set()
    for i in range(len(E)):
        for j in range(len(E)):
            if E[i] != E[j]:
                denoms.add(abs(E[i]-E[j]))
        denoms.add(abs(E[i]) if E[i] != 0 else 1)
    bps = {F(0), F(1)}
    for d in denoms:
        for m in range(0, d+1):
            bps.add(F(m, d))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        mid = (a+b)/2
        # phases p_e(x)=e*x-floor(e*mid); order fixed on (a,b)
        fl = {e: (e*mid).__floor__() for e in E}
        order = sorted(E, key=lambda e: e*mid - fl[e])
        m = len(order)
        # gaps g_s(x) between order[s], order[s+1]
        subs = []
        for s in range(m):
            e1 = order[s]; e2 = order[(s+1) % m]
            if s < m-1:
                c = F(e2-e1); b0 = F(-(fl[e2]-fl[e1]))
            else:
                c = F(order[0]-order[-1]); b0 = F(-(fl[order[0]]-fl[order[-1]]) + 1)
            # region c*x+b0 > thr on (a,b)
            if c == 0:
                if b0 > thr: subs.append((a, b))
            elif c > 0:
                x0 = (thr-b0)/c; lo = max(a, x0)
                if lo < b: subs.append((lo, b))
            else:
                x0 = (thr-b0)/c; hi = min(b, x0)
                if a < hi: subs.append((a, hi))
        # union
        subs.sort()
        cur = None
        for lo, hi in subs:
            if lo >= hi: continue
            if cur is None: cur = [lo, hi]
            elif lo <= cur[1]: cur[1] = max(cur[1], hi)
            else:
                total += cur[1]-cur[0]; cur = [lo, hi]
        if cur is not None:
            total += cur[1]-cur[0]
    return total

print("="*70)
print("Q1/Q2: sharp threshold of the Vmax-ruler (co-offset config WITH 0)")
print("  A ruler-period at slow x is GOOD iff exists phi in (1/14,13/14) with")
print("  ||phi - e x||>=1/14 for all e.  Config seen = {frac(e x): e in E} incl 0.")
print("  Witness phi EXISTS <=> some gap > 1/7 (each endpoint eats 1/14).")
print("  => predict SHARP threshold = 1/7 on co-offset-with-0 config.")
print("="*70)

def actual_good_fraction(S, Vmax, npts=400):
    good = 0
    for j in range(Vmax):
        lo = (14*j+1)/(14*Vmax); hi = (14*j+13)/(14*Vmax)
        mx = 0.0
        for a in range(npts):
            t = lo + (hi-lo)*(a+0.5)/npts
            m = 1.0
            for v in S:
                d = dist_int(v*t)
                if d < m: m = d
            if m > mx: mx = m
        if mx >= 1/14 - 1e-9: good += 1
    return good / Vmax

test_offsets = [
    ("consec5",   [0,1,2,3,4]),
    ("perf7",     [0,2,3,4,5,6,8]),
    ("APcoff12",  list(range(13))),
]
for name, E in test_offsets:
    g17 = float(mu_exact(E, F(1,7)))
    g27 = float(mu_exact(E, F(2,7)))
    print(f"{name:9s} |E|={len(E):2d}  limit good-density (EXACT): thr=1/7 -> {g17:.4f}   thr=2/7 -> {g27:.4f}")
    for Vmax in [200, 1000]:
        S = [Vmax] + [Vmax - e for e in E if e > 0]
        af = actual_good_fraction(S, Vmax)
        print(f"          Vmax={Vmax:5d}: ACTUAL good-period frac = {af:.4f}")
    print()

print("="*70)
print("Q3: SPEED-config {1..k} vs CO-OFFSET-config {0..k-1}: mu_{1/7} & mu_{2/7} (EXACT)")
print("="*70)
print(f"{'k':>3} | {'mu17(speed)':>11} {'mu17(coff0)':>11} | {'mu27(speed)':>11} {'mu27(coff0)':>11}")
for k in range(8, 14):
    speed = list(range(1, k+1)); coff = list(range(0, k))
    print(f"{k:>3} | {float(mu_exact(speed,F(1,7))):>11.4f} {float(mu_exact(coff,F(1,7))):>11.4f} "
          f"| {float(mu_exact(speed,F(2,7))):>11.4f} {float(mu_exact(coff,F(2,7))):>11.4f}")

print("="*70)
print("Q4: mu_{1/7}(AP) exact vs known, E[maxgap](AP), adversarial AP-minimality")
print("="*70)
known = {8:F(691,735),9:F(247,294),10:F(38,49),11:F(1381,2205),12:F(13823,24255),13:F(477,1078)}
print(f"{'k':>3} | {'mu17 exact':>18} {'= known?':>9} | {'E[maxgap]~':>10}")
for k in range(8, 14):
    ap = list(range(1, k+1))
    me = mu_exact(ap, F(1,7))
    ok = "YES" if me == known[k] else f"NO({float(known[k]):.4f})"
    print(f"{k:>3} | {str(me):>18} {ok:>9} | {E_maxgap(ap,60000):>10.4f}")

print("\nAdversarial AP-minimality (k=13): beat mu17(AP)=477/1078=%.4f or E[maxgap](AP)?" % (477/1078))
ap13 = list(range(1,14))
base_mu = float(mu_exact(ap13, F(1,7)))
base_emg = E_maxgap(ap13, 120000)
print(f"  AP13: mu17={base_mu:.4f}  E[maxgap]={base_emg:.4f}")
trials = [
    ("2*AP",              [2*j for j in range(1,14)]),
    ("primes<=41",        [2,3,5,7,11,13,17,19,23,29,31,37,41]),
    ("GW {1..11,13,24}",  [1,2,3,4,5,6,7,8,9,10,11,13,24]),
    ("split 1..6,20..26", [1,2,3,4,5,6,20,21,22,23,24,25,26]),
    ("fib-ish",           [1,2,3,5,8,13,21,34,55,89,144,233,377]),
]
for nm, E in trials:
    m = float(mu_exact(E, F(1,7))); e = E_maxgap(E, 60000)
    fm = "  <<BELOW AP mu!" if m < base_mu-1e-4 else ""
    fe = "  <<BELOW AP E!"  if e < base_emg-1e-4 else ""
    print(f"  {nm:20s} mu17={m:.4f}{fm}   E[maxgap]={e:.4f}{fe}")

random.seed(20260707)
worst_mu, wm = base_mu, ap13
worst_e, we = base_emg, ap13
bm = be = 0
for _ in range(120):
    E = sorted(random.sample(range(1, 60), 13))
    m = mu_threshold(E, 1/7, 20000); e = E_maxgap(E, 20000)
    if m < worst_mu: worst_mu, wm = m, E
    if e < worst_e:  worst_e, we = e, E
    if m < base_mu-3e-3: bm += 1
    if e < base_emg-3e-3: be += 1
print(f"  random(120,coarse): {bm} beat AP mu17, {be} beat AP E[maxgap]")
print(f"  worst mu17={worst_mu:.4f} at {wm}")
print(f"  worst E[maxgap]={worst_e:.4f} at {we}")
