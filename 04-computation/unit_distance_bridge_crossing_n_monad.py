"""
monad-explorer-2026-06-13  --  GATE-3 PROBE for OPEN-Q-057 (extends THM-494/HYP-2461).

QUESTION: is t=3 (sqrt(-11)) UNIQUE for the 3N crossing, or just FIRST?
HYPOTHESIS (minimal-transverse-distance principle): t=3 crosses at n=28 because
sqrt(3) is the MINIMAL Loeschian transverse distance (>1), hence the only one
realized in the SMALL optimal factors landing at n~28.  Larger Loeschian t have
transverse distance sqrt(t) that only appears in LARGER factors, so their crossing
(if any) should be at LARGER n -- mirroring HYP-2301's single-lattice sqrt7@32.

So: for each transverse-bearing L_t, find the SMALLEST n where a free densest patch
beats 3n.  Prediction: n_cross(t=3) = 28 < n_cross(t=13) < ...  If a larger t also
crosses (at larger n), t=3 is FIRST not UNIQUE; the gate-3 selector is "smallest n",
governed by minimal sqrt(t).

PART A (cheap, exact): the minimal number of points of a unit-distance graph (triangular
patch) that contains a norm-t Eisenstein pair, for small Loeschian t.  Confirms sqrt(3)
is realized in the smallest factors.

PART B (annealing, lower bound): extend the bridge-family free search to larger n for a
few transverse-bearing t, locating each one's first crossing of 3n.

Reuses the engine of unit_distance_bridge_lattice_family_monad.py verbatim.
"""
import sys, random, importlib.util as IU
from itertools import product as iproduct

random.seed(20260613)
spec = IU.spec_from_file_location("blf",
        "04-computation/unit_distance_bridge_lattice_family_monad.py")
blf = IU.module_from_spec(spec)
# prevent its __main__ from running:
import builtins; _name = blf.__name__
spec.loader.exec_module(blf)

# ---------- PART A: minimal factor carrying a norm-t Eisenstein pair ----------
# Eisenstein points a + b*zeta6 ; squared length (norm) = a^2 + a*b + b^2.
def eis_norm(a, b):
    return a*a + a*b + b*b

def loeschian(t):
    for a in range(-t, t+1):
        for b in range(-t, t+1):
            if eis_norm(a, b) == t:
                return True
    return False

def min_udg_with_normt_pair(t, maxpts=8):
    """Smallest k such that SOME k-point connected triangular patch contains a pair
    at Eisenstein norm t.  We look for a short unit-edge path realizing the displacement."""
    # An Eisenstein vector of norm t = sum of unit Eisenstein steps; the minimal number
    # of unit steps to realize a displacement of norm t is its "Eisenstein word length".
    units = [(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]  # the 6 unit Eisenstein vectors
    # BFS over displacements by #unit steps
    from collections import deque
    target_disps = {(a,b) for a in range(-t,t+1) for b in range(-t,t+1) if eis_norm(a,b)==t}
    seen = {(0,0):0}; dq = deque([(0,0)])
    while dq:
        x = dq.popleft(); d = seen[x]
        if x in target_disps:
            return d+1, d  # k points = d edges + 1 ; (#points, #edges = word length)
        if d >= maxpts: continue
        for u in units:
            y = (x[0]+u[0], x[1]+u[1])
            if y not in seen:
                seen[y] = d+1; dq.append(y)
    return None, None

def part_A():
    print("PART A — minimal unit-distance path realizing a norm-t (=dist sqrt(t)) pair", flush=True)
    print("  (Eisenstein word length = min # unit steps; #points = length+1)", flush=True)
    print(f"  {'t':>3} {'sqrt(t)':>8} {'Loeschian':>10} {'min#pts':>8} {'min#edges':>9}", flush=True)
    for t in range(1, 32):
        if not loeschian(t): continue
        pts, edges = min_udg_with_normt_pair(t)
        print(f"  {t:>3} {t**0.5:>8.3f} {'yes':>10} {str(pts):>8} {str(edges):>9}", flush=True)
    print("  => the MINIMAL transverse distance >1 is sqrt(3) (t=3), realized by a 3-point", flush=True)
    print("     2-edge path; t=4 needs a colinear length-2 (3 pts collinear); larger t need", flush=True)
    print("     progressively larger factors. The optimal 3-pt UDG is K3 (sqrt3-FREE), which", flush=True)
    print("     is why 27=3^3 gets zero bonus but 28=4*7 (rhombus has a sqrt3 pair) crosses.", flush=True)

# ---------- PART B: first crossing n for each transverse-bearing L_t ----------
def first_crossing(t, box, n_lo, n_hi, iters, restarts):
    blf.UNITS = blf.enumerate_units(t, box); blf.UNITSET = set(blf.UNITS)
    assert len(blf.UNITS) == blf.predicted_units(t), "raise box"
    print(f"\nL_{t}: #units={len(blf.UNITS)} (sqrt(-{4*t-1}))  scanning n={n_lo}..{n_hi}", flush=True)
    print(f"   {'n':>3} {'best':>5} {'3n':>4} {'U-3n':>5}  {'cross?':>6}", flush=True)
    cross_n = None
    for n in range(n_lo, n_hi+1):
        b, S = blf.densest(n, iters=iters, restarts=restarts)
        rc = blf.recount(S, t); assert rc == b
        tn = 3*n; tag = ''
        if b > tn:
            tag = '  *** CROSSES ***'
            if cross_n is None: cross_n = n
        elif b == tn: tag = '  tie'
        print(f"   {n:>3} {b:>5} {tn:>4} {b-tn:>5} {tag}", flush=True)
    print(f"   -> first crossing for t={t}: n = {cross_n}", flush=True)
    return cross_n

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("all","A"):
        part_A()
    if mode in ("all","B"):
        print("\nPART B — first 3N crossing per transverse-bearing L_t (annealing lower bound)", flush=True)
        # t=3 (sqrt11): confirm 28.   t=13 (sqrt51, 24 units): where does it first cross?
        results = {}
        results[3]  = first_crossing(3,  box=4, n_lo=26, n_hi=32, iters=90000,  restarts=10)
        results[13] = first_crossing(13, box=5, n_lo=28, n_hi=40, iters=120000, restarts=10)
        print("\nSUMMARY first-crossing n:", results, flush=True)
        print("If n_cross(13) > 28 = n_cross(3): t=3 is FIRST (minimal sqrt(t)), and larger-t", flush=True)
        print("crossings exist at larger n -> gate 3 = 'smallest n', minimal-distance governed.", flush=True)
    print("\nDONE.", flush=True)
