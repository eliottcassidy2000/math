"""
monad-explorer-2026-06-13  —  HYP-2461  (companion to peer HYP-2460/THM-493).

DENSER BRIDGE LATTICES FOR SMALL UNIT-DISTANCE MAXIMA.

Peer THM-493 (same day) showed L_t is the Minkowski product of the triangular
lattice with an w_t-rotated copy, and that the t=3 (Moser) "resonant angle" is what
enables its u(28)=85 crossing via a transverse bonus -- using CURATED 2-FACTOR product
configs.  Here I confirm the same conclusion from a NON-product method: free-patch
simulated annealing over the WHOLE lattice, across the unit-count-increasing family.

THM-434:  the Moser-ladder bridge lattice
    L_t = Z[zeta6] (+) Z[zeta6]*w_t,   w_t = ((2t-1)+i*sqrt(4t-1))/(2t), |w_t|=1
(rank-4 in C iff 4t-1 is NOT 3*square) has EXACTLY  #units(L_t) = 12 + 6*B(t)
unit vectors, B(t) = sum_{d|t} chi(d), chi=(-3/.).  So:
    t=3 -> 18 units (the MOSER lattice M_L; Engel et al. use ONLY this one)
    t=13,21,31,43 -> 24 units;   t=49 -> 30 units;   t=133 -> 36 units.

QUESTION.  The Engel/Moser dense-patch construction reproducing
u(n)=...,81(tie@27),85(>3N @28),... was run ONLY for t=3.  Does a bridge lattice
with MORE unit vectors give DENSER small patches -- u(27)>=82 (=> N* <= 27, refuting
HYP-2299) or u(28)>=86 -- or does t=3 stay best (corroborating N*=28)?

METHOD.  EXACT, IDENTICAL search on every L_t: greedy seed + simulated annealing
(the s4 search that provably reproduces Engel's table on L_3).  Adjacency is the
PURE INTEGER test below; the lattice is infinite (neighbours generated on the fly via
the unit-vector offsets), so no disk/cap heuristic distorts the comparison.  A patch's
unit-distance count = #pairs whose integer difference is a unit vector (exact, since
{1,zeta6,w_t,zeta6 w_t} are Q-independent => distinct 4-tuples are distinct points).

EXACT integer adjacency.  z=(a,b,c,d) => |z|^2=1  iff  bc-ad=0  AND
  2t(a^2+b^2+c^2+d^2)+2t(ab+cd)+2(2t-1)(ac+bd)+(2t-1)(ad+bc) = 2t.
(Derived from the Gram matrix; the only irrationality sqrt(3(4t-1)) sits on bc-ad.)
"""
import sys, random
from math import isqrt
from itertools import product

random.seed(20260613)

def is_three_times_square(d):
    if d % 3: return False
    k = d // 3; r = isqrt(k); return r*r == k

def chi(d):
    d %= 3
    return 0 if d == 0 else (1 if d == 1 else -1)

def B_of(t):
    res = 0; d = 1
    while d*d <= t:
        if t % d == 0:
            res += chi(d)
            if d != t // d: res += chi(t // d)
        d += 1
    return res

def predicted_units(t):
    return 12 + 6 * B_of(t)

def is_unit(v, t):
    a, b, c, d = v
    if b*c - a*d != 0: return False
    return (2*t*(a*a+b*b+c*c+d*d) + 2*t*(a*b+c*d)
            + 2*(2*t-1)*(a*c+b*d) + (2*t-1)*(a*d+b*c)) == 2*t

def enumerate_units(t, box):
    return [v for v in product(range(-box, box+1), repeat=4)
            if v != (0,0,0,0) and is_unit(v, t)]

# ---- search (globals UNITS/UNITSET set per lattice) -----------------------
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

def densest(N, iters=110000, restarts=12):
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

def recount(S, t):
    """fully independent exact pairwise recount (integer)."""
    L = list(S); n = len(L); c = 0
    for i in range(n):
        for j in range(i+1, n):
            d = (L[i][0]-L[j][0], L[i][1]-L[j][1], L[i][2]-L[j][2], L[i][3]-L[j][3])
            if is_unit(d, t): c += 1
    return c

ENGEL = {21:57, 22:60, 23:64, 24:68, 25:72, 26:76, 27:81, 28:85, 29:89, 30:93}

def run_lattice(t, ns, box, iters, restarts):
    global UNITS, UNITSET
    assert not is_three_times_square(4*t-1)
    UNITS = enumerate_units(t, box); UNITSET = set(UNITS)
    pred = predicted_units(t); ok = (len(UNITS) == pred)
    print(f"\n{'='*66}\nL_{t}: 4t-1={4*t-1}, #units={len(UNITS)} "
          f"(THM-434:{pred}) {'OK' if ok else '*** raise box ***'}", flush=True)
    assert ok, "incomplete unit enumeration -- raise box"
    print(f"   {'n':>3} {'best':>5} {'recount':>8} {'Engel':>6} {'3n':>4} "
          f"{'U-3n':>5}  cmp", flush=True)
    res = {}
    for n in ns:
        b, S = densest(n, iters=iters, restarts=restarts)
        rc = recount(S, t)
        assert rc == b, f"recount mismatch at n={n}: {rc} vs {b}"
        res[n] = b
        e = ENGEL.get(n); tn = 3*n
        if e is None: cmp = ''
        elif b > e:   cmp = f'  <<< BEATS ENGEL by {b-e}'
        elif b == e:  cmp = '  = Engel'
        else:         cmp = f'  Engel-{e-b}'
        print(f"   {n:>3} {b:>5} {rc:>8} {str(e):>6} {tn:>4} {b-tn:>5} {cmp}",
              flush=True)
    return res

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    ns = list(range(21, 31))
    print("UNIT-DISTANCE DENSE PATCHES OVER MOSER-LADDER BRIDGE LATTICES L_t",
          flush=True)
    print("Exact integer adjacency; identical annealing search on each lattice.",
          flush=True)
    # (t, box) -- box big enough to capture all 12+6B(t) unit vectors
    plan = [(3,4), (13,5), (21,6), (31,7), (49,8)]
    if which == "validate":
        plan = [(3,4)]
    for t, box in plan:
        run_lattice(t, ns, box=box, iters=110000, restarts=12)
    print("\nDONE.", flush=True)
