"""
monad-explorer-2026-06-13  --  OPEN-Q-057 handoff (HYP-2461 next-explorer probe).

THE EXACT-30-DEGREE BISECTOR LATTICE  Z[zeta_12]  --  DOES IT CROSS 3N AT 28?

Handoff (HYP-2461 / reflection 'the-unit-distance-tie-is-carrier-robust...'):
the L_t bridge family rotates the second triangular copy by the MOSER angle
arccos((2t-1)/2t) -- t=3 -> 33.6 deg (CROSSES at 28, u=85), t=4 -> 29.0 deg
(TIES 81 at 27 but caps 83<84 at 28).  The "perfect bisector" 30 deg sits BETWEEN
them and is NOT an L_t (cos 30 = sqrt3/2 is never (2t-1)/2t).  Z[zeta_12] is the
natural exact-30deg carrier: 12 unit vectors = the 12th roots of unity = two
triangular hexagons offset by exactly 30 deg.  Does it cross at 28 (=> "good-angle
band") or fail (=> sqrt(-11) arithmetically singular)?

KEY ARITHMETIC FACT (this session, proved):  Z[zeta_12] has EXACTLY 12 unit
vectors and ZERO transverse units.
  Proof: if z in Z[zeta_12], |z|=1, then z*conj(z) in Z[sqrt3] (the real subfield),
  and equals 1 under the chosen real embedding; since sqrt3 is irrational the only
  element of Z[sqrt3] mapping to 1 is 1 itself, so z*conj(z)=1 EXACTLY.  Then ALL
  archimedean conjugates of z have absolute value 1, so by Kronecker z is a root of
  unity, hence z in mu_12.  The 12 split as mu_6 (even powers, one hexagon) and
  zeta_12*mu_6 (odd powers, the 30deg-rotated hexagon).  NO unit vector mixes the
  two hexagons => transverse-free, like the t=2,5 (sqrt7,sqrt19) bridge lattices.

PREDICTION (from the transverse-free analogy): Z[zeta_12] should behave like t=2,5
-- cap near 78 at n=27, NOT reach the 81 tie, NOT cross at 28.  This REFRAMES the
handoff: Z[zeta_12] fails not because of a "bad angle band" but because the exact
30deg bisector has NO transverse unit vectors at all (cross-hexagon distances are
2*sin(odd*15deg), never 1).  The transverse engine REQUIRES an irrational Moser
angle arccos((2t-1)/2t); the rational-cosine 30deg cannot carry it.

EXACT INTEGER ADJACENCY (basis 1,zeta,zeta^2,zeta^3 ; zeta=e^{i pi/6}; zeta^4=zeta^2-1):
  z=(a,b,c,d).  X=2*Re = (2a+c)+b*sqrt3,  Y=2*Im = (b+2d)+c*sqrt3.
  Let P1=2a+c, Q1=b, P2=b+2d, Q2=c.  Then 4|z|^2 = (P1^2+3Q1^2+P2^2+3Q2^2)
     + 2(P1Q1+P2Q2)*sqrt3.   So  |z|^2 = 1  <=>  P1Q1+P2Q2 = 0  AND
     P1^2+3Q1^2+P2^2+3Q2^2 = 4.   (Pure integer test; sqrt3 only on the cross term.)

Engine (greedy seed + simulated annealing + independent exact pairwise recount) is
COPIED VERBATIM from unit_distance_bridge_lattice_family_monad.py so the numbers are
directly comparable to the peer's L_t table.
"""
import sys, random
from itertools import product

random.seed(20260613)

# ---- exact integer adjacency for Z[zeta_12] -------------------------------
def is_unit_z12(v):
    a, b, c, d = v
    P1 = 2*a + c; Q1 = b
    P2 = b + 2*d; Q2 = c
    if P1*Q1 + P2*Q2 != 0:           # sqrt3 coefficient must vanish
        return False
    return P1*P1 + 3*Q1*Q1 + P2*P2 + 3*Q2*Q2 == 4   # rational part == 4|z|^2 = 4

# the 12 unit offsets (zeta^k, k=0..11), as integer 4-tuples in basis 1,zeta,zeta^2,zeta^3
Z12_UNITS = [
    ( 1, 0, 0, 0),  # zeta^0  = 1
    ( 0, 1, 0, 0),  # zeta^1
    ( 0, 0, 1, 0),  # zeta^2  (=zeta_6)
    ( 0, 0, 0, 1),  # zeta^3  (=i)
    (-1, 0, 1, 0),  # zeta^4  = zeta^2 - 1
    ( 0,-1, 0, 1),  # zeta^5
    (-1, 0, 0, 0),  # zeta^6  = -1
    ( 0,-1, 0, 0),  # zeta^7
    ( 0, 0,-1, 0),  # zeta^8
    ( 0, 0, 0,-1),  # zeta^9
    ( 1, 0,-1, 0),  # zeta^10
    ( 0, 1, 0,-1),  # zeta^11
]
# the 6 even powers = triangular sub-hexagon (for the Harborth calibration)
TRI_UNITS = [Z12_UNITS[k] for k in (0,2,4,6,8,10)]

def brute_units(box):
    return [v for v in product(range(-box, box+1), repeat=4)
            if v != (0,0,0,0) and is_unit_z12(v)]

# ---- search machinery (identical to the bridge-lattice engine) ------------
UNITS = []; UNITSET = set()

def deg_in(p, S):
    return sum((p[0]+u[0],p[1]+u[1],p[2]+u[2],p[3]+u[3]) in S for u in UNITS)

def U_exact(S):
    s = S if isinstance(S, set) else set(S)
    return sum(deg_in(p, s) for p in s) // 2

def cands(S):
    c = set()
    for p in S:
        for u in UNITS:
            q = (p[0]+u[0],p[1]+u[1],p[2]+u[2],p[3]+u[3])
            if q not in S: c.add(q)
    return c

def greedy_grow_from(S, N):
    S = set(S)
    while len(S) < N:
        cs = cands(S)
        if not cs: break
        S.add(max(cs, key=lambda q: deg_in(q, S)))
    return S

def greedy_grow(N):
    return greedy_grow_from({(0,0,0,0)}, N)

def anneal(S, iters):
    S = set(S); E = U_exact(S); best = E; bestS = set(S)
    for it in range(iters):
        T = max(0.05, 3.0*(1 - it/iters))
        u = random.choice(tuple(S))
        p = random.choice(tuple(S)); v = random.choice(UNITS)
        w = (p[0]+v[0],p[1]+v[1],p[2]+v[2],p[3]+v[3])
        if w in S or w == u: continue
        du = deg_in(u, S)
        diff = (w[0]-u[0],w[1]-u[1],w[2]-u[2],w[3]-u[3])
        dw = deg_in(w, S) - (1 if diff in UNITSET else 0)
        delta = dw - du
        if delta >= 0 or random.random() < pow(2.718281828, delta/max(T,1e-9)):
            S.remove(u); S.add(w); E += delta
            if E > best: best = E; bestS = set(S)
    return best, bestS

def densest(N, iters, restarts):
    best = -1; bestS = None
    base = greedy_grow(N + 8)
    for r in range(restarts):
        if r == 0:
            seed = greedy_grow(N)
        else:
            b = list(base); random.shuffle(b)
            seed = greedy_grow_from(set(b[:1]), N)
        e, S = anneal(seed, iters)
        if e > best: best, bestS = e, S
    return best, bestS

def recount(S):
    L = list(S); n = len(L); c = 0
    for i in range(n):
        for j in range(i+1, n):
            d = (L[i][0]-L[j][0], L[i][1]-L[j][1], L[i][2]-L[j][2], L[i][3]-L[j][3])
            if is_unit_z12(d): c += 1
    return c

ENGEL  = {21:57,22:60,23:64,24:68,25:72,26:76,27:81,28:85,29:89,30:93}
HARBORTH = {n: 3*n - __import__('math').isqrt(12*n-3) for n in range(21,31)}  # floor(3n-sqrt(12n-3))

def run(units, label, ns, iters, restarts, ref):
    global UNITS, UNITSET
    UNITS = list(units); UNITSET = set(UNITS)
    print(f"\n{'='*70}\n{label}: #unit vectors = {len(UNITS)}", flush=True)
    print(f"   {'n':>3} {'best':>5} {'recount':>8} {'ref':>5} {'3n':>4} {'U-3n':>5}  cmp",
          flush=True)
    res = {}
    for n in ns:
        b, S = densest(n, iters=iters, restarts=restarts)
        rc = recount(S)
        assert rc == b, f"recount mismatch at n={n}: {rc} vs {b}"
        res[n] = b
        r = ref.get(n); tn = 3*n
        if r is None: cmp=''
        elif b>r: cmp=f'  <<< +{b-r} over ref'
        elif b==r: cmp='  = ref'
        else: cmp=f'  ref-{r-b}'
        tag = '  *** CROSSES 3N ***' if b>tn else ('  =3N tie' if b==tn else '')
        print(f"   {n:>3} {b:>5} {rc:>8} {str(r):>5} {tn:>4} {b-tn:>5} {cmp}{tag}",
              flush=True)
    return res

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    ns = list(range(21, 31))

    # 0) confirm the arithmetic: brute-force unit enumeration must be EXACTLY the 12 roots.
    bu = brute_units(box=3)
    print(f"Z[zeta_12] brute-force units in box 3: found {len(bu)} (claim: 12)", flush=True)
    assert sorted(bu) == sorted(Z12_UNITS), f"unit set mismatch: {sorted(bu)}"
    print("  CONFIRMED: exactly the 12 roots of unity; transverse-free.", flush=True)

    if mode in ("all","calib"):
        # 1) CALIBRATION: the 6 even units = pure triangular lattice -> must match Harborth.
        run(TRI_UNITS, "CALIBRATION triangular sub-hexagon (expect Harborth floor(3n-sqrt(12n-3)))",
            ns, iters=60000, restarts=8, ref=HARBORTH)

    if mode in ("all","z12"):
        # 2) THE TEST: full Z[zeta_12], 12 unit vectors at exact 30deg spacing.
        run(Z12_UNITS, "Z[zeta_12] exact-30deg bisector (DOES IT CROSS 3N AT 28?)",
            ns, iters=110000, restarts=12, ref=ENGEL)

    print("\nDONE.", flush=True)
