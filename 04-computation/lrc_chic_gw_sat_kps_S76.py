#!/usr/bin/env python3
r"""
lrc_chic_gw_sat_kps_S76.py  (kind-pasteur-2026-07-07-S76)

RESOLVE chi_c(G_GW) in (13,14].  G_GW = Cay(Z, +-GW), GW = {1..11,13,24}.
  chi_f(G_GW) = 13 (= 1/mu, mu=1/13, opus-S144)
  chi(G_GW)   = 14 (= 1/M, opus-S145 THM-652, odd-cycle matching obstruction)
  chi_c       in (13,14]  -- the OPEN rung question.

KEY REDUCTION (kps-S76): rotation/linear colorings c(x)=floor(p*frac(x*t)) give exactly
chi_c <= 1/M = 14 (a lonely-runner witness t is a rotation (p,q)-coloring with q/p=M).  So
    chi_c(G_GW) < 14  <=>  G_GW admits a NON-ROTATION circular (p,q)-coloring with p/q<14.
That is the LINEARIZATION GAP (opus-S141 homomorphism ladder) at the GW tight instance.

A (p,q)-coloring of C(Z_N,GW) for N>48 (clean quotient: GW mod N and -GW mod N are 26
distinct nonzero residues) LIFTS to a periodic (p,q)-coloring of G(Z,GW).  So
    exists N>48, (p,q)-coloring of C(Z_N,GW) with p/q<14  =>  chi_c(G_GW) < 14.

SAT ENCODING: vars b[x][a] = "vertex x has color a" (x in Z_N, a in Z_p).  Exactly-one
color per vertex.  For each edge (x, x+s mod N), s in GW, forbid color pairs (a,b) with
circular distance |a-b|_p < q.  Solve; a model = a circular coloring beating 14.

We sweep p/q in (13,14) with small denominators q and N a range of periods > 48.
"""
import sys
from math import gcd
from pysat.solvers import Cadical153 as SatSolver
from pysat.formula import CNF

GW = [1,2,3,4,5,6,7,8,9,10,11,13,24]

def circ_dist(a, b, p):
    d = abs(a - b) % p
    return min(d, p - d)

def build_cnf(N, p, q):
    """CNF for a (p,q)-coloring of the circulant C(Z_N, GW)."""
    cnf = CNF()
    def V(x, a):  # 1-based var id
        return x * p + a + 1
    # exactly-one color per vertex
    for x in range(N):
        cnf.append([V(x, a) for a in range(p)])           # at least one
        for a in range(p):
            for b in range(a + 1, p):
                cnf.append([-V(x, a), -V(x, b)])          # at most one
    # edge constraints: forbid close color pairs
    seen = set()
    for x in range(N):
        for s in GW:
            y = (x + s) % N
            e = (min(x, y), max(x, y))
            if e in seen:
                continue
            seen.add(e)
            for a in range(p):
                for b in range(p):
                    if circ_dist(a, b, p) < q:
                        cnf.append([-V(x, a), -V(y, b)])
    return cnf, V

def try_coloring(N, p, q, conf_budget=200000):
    cnf, V = build_cnf(N, p, q)
    s = SatSolver(bootstrap_with=cnf.clauses)
    if conf_budget:
        s.conf_budget(conf_budget)
        res = s.solve_limited(expect_interrupt=True)
    else:
        res = s.solve()
    model = None
    if res is True:
        m = s.get_model()
        mset = set(v for v in m if v > 0)
        model = []
        for x in range(N):
            for a in range(p):
                if V(x, a) in mset:
                    model.append(a); break
    s.delete()
    return res, model  # True / False / None(=budget exceeded)

def is_rotation(coloring, N, p):
    """check whether the found coloring is (essentially) a linear/rotation coloring
    c(x) = round(a*x/N * p + b) type -- i.e. color increments are near-constant."""
    if coloring is None:
        return None
    diffs = [(coloring[(x+1) % N] - coloring[x]) % p for x in range(N)]
    return len(set(diffs)) <= 2  # near-constant increment => rotation-like

print("=" * 84)
print("chi_c(G_GW): SAT search for a (p,q)-coloring of C(Z_N,GW) with 13 < p/q < 14")
print("  (any hit => chi_c < 14 => linearization FAILS at GW; the rotation bound is 14)")
print("=" * 84)

# --- sanity: validate the encoding against opus THM-652 (chi=14) ---
print("\n  [sanity] integer colorings of C(Z_78,GW): expect (13,1) UNSAT, (14,1) SAT")
for p in (13, 14):
    res, col = try_coloring(78, p, 1, conf_budget=500000)
    print(f"    ({p},1)-coloring: {'SAT' if res is True else 'UNSAT' if res is False else 'budget'}")

# candidate (p,q) with 13 < p/q < 14; ORDER BY p ASCENDING (smallest CNF / fastest first)
CANDS = []
for q in range(2, 9):
    for p in range(13 * q + 1, 14 * q):        # strictly between 13 and 14
        if gcd(p, q) == 1:
            CANDS.append((p, q, p / q))
CANDS.sort(key=lambda t: (t[0], t[1]))          # smallest p first (fast SAT)
print(f"\n  candidates by ascending p: {[f'{p}/{q}={r:.3f}' for p,q,r in CANDS[:14]]} ...")

# periods > 48; multiples of 13/26 (the pattern scale) + a few others for freedom
NS = [52, 65, 78, 91, 104, 130]
found = None
for (p, q, r) in CANDS:
    if r >= 14:
        continue
    row = []
    for N in NS:
        res, coloring = try_coloring(N, p, q)
        if res is True:
            rot = is_rotation(coloring, N, p)
            print(f"  p/q={p}/{q}={r:.4f}  N={N}: *** SAT ***  rotation-like={rot}")
            found = (p, q, r, N, coloring)
            break
        row.append(f"N{N}:{'unsat' if res is False else '?'}")
    if found:
        break
    print(f"  p/q={p}/{q}={r:.4f}: " + " ".join(row))

print()
if found:
    p, q, r, N, col = found
    print(f"=> chi_c(G_GW) <= {p}/{q} = {r:.4f} < 14: LINEARIZATION FAILS at GW.")
    print(f"   coloring (period {N}, {p} colors, gap {q}): {col}")
    # verify
    ok = all(circ_dist(col[x], col[(x+s) % N], p) >= q for x in range(N) for s in GW)
    print(f"   verification (all edges gap>=q): {ok}")
else:
    print("=> no sub-14 coloring found at these (p,q,N): evidence for chi_c(G_GW) = 14")
    print("   (rotation bound tight; the odd-cycle obstruction may extend to chi_c)")
