#!/usr/bin/env python3
"""procgen_kuratowski_20260925, part K: Kohl's Collatz group G_C as a Tait-coloured Schreier graph.

Kohl (J. Group Theory 20 (2017), Prop. 1.2 = Prop. 2.1; primary text read in the atlas lane, file
scratch/procgen_atlas/lit_rcwa/dl/kohl_collatzgroups.txt): G_C = <a,b,c>,
    a = tau_{1(2),4(6)}  (n <-> 3n+1, n odd),
    b = tau_{1(3),2(6)}  (m <-> 2m, m = 1 mod 3),
    c = tau_{2(3),4(6)}  (m <-> 2m, m = 2 mod 3),
acts transitively on N \ 0(6) iff the Collatz conjecture holds.

Sections (every check raises on failure):
  K1  generators: involutions, disjoint classes, fixed sets; the mod-6 vertex-type table (degree,
      semi-edges); the Schreier graph (non-loop edges) equals the undirected C-graph induced on
      Z \\ 0(6) on the window |n| <= X.
  K2  cycles: on N \\ 0(6) the only cycle met by orbits n <= X1 is the rainbow triangle {1,4,2};
      on Z_<0 \\ 0(6): the (a,c) digon {-1,-2}, the 5-cycle through -5, the 18-cycle through -17,
      with their colour words.
  K3  Kempe chains: <b,c> = doubling rays / singletons; <a,b> = P_n = {2n?, n, 3n+1, 6n+2};
      <a,c> = odd runs y -> T^k(y), k = v2(y+1); the ONLY closed Kempe chain on Z \\ 0(6) is the
      <a,c>-digon {-1,-2}; <b,c>-chains are D-orbits (D(x)=2x) and <a,c>-chains are E-orbits
      (E(x)=(2x-1)/3), the two Banach contractions of THM-4471 with fixed points 0 (deleted by
      Kohl) and -1 (the digon); the Mersenne start 2^k-1 opens the chain ending at 3^k-1.
  K4  Kempe quotients: K/<b,c> = undirected Syracuse graph; K/<a,b> = subcubic "sibling" graph on
      odd numbers with edges {n,4n+1}, {n,T(n)} (n=3 mod 4), {n,T^2(n)} (n=1 mod 8), degrees 3/2;
      K/<a,c> = run graph with degree 1 + [y=1 mod 6] + v2(y+1).
  K5  Z2xZ2 boundary of the Tait colouring: delta = c,a,a,0,b on 1,2,3,4,5 mod 6.
  K6  switches: a Kempe switch on an rcwa union of <a,b>-chains keeps the graph (new rcwa
      generators); nu-conjugation gives (a,b,c) -> (a^h0, c, b) with h0 = tau_{2(6),4(6)}, the
      3n-1 group; single-colour class-transposition census (moduli 6,12 in full, 6..24 counted).
  K7  Kohl's mod-3 shape needs 3 | q (explicit q = 5 cut); the mod-p generalisation exists iff p | q and
      -1 is a power of 2 mod p; the mod-5 Tait group of 5n+1 (a DRIFT control) is built and checked.
  K8  Kohl's G_T = <tau_{0(2),1(2)}, tau_{1(2),2(4)}, tau_{1(4),2(6)}>: Schreier graph on [0,X2]
      has cyclomatic number 1 (the b/c digon {1,2}); colour a is the consecutive pairing {2k,2k+1}.
Run: python3 04-computation/experiments/procgen_kuratowski_20260925_kohl_tait.py  (~1 min, < 300 MB)
"""
import itertools
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction as Fr

X = 60000        # window |n| <= X for structural checks
X1 = 200000      # Collatz orbit check on N
XN = 200000      # negative orbit check


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)


def ct(r1, m1, r2, m2):
    """class transposition tau_{r1(m1), r2(m2)} on Z (classes assumed disjoint)."""
    def f(n):
        if n % m1 == r1:
            return r2 + (n - r1) // m1 * m2
        if n % m2 == r2:
            return r1 + (n - r2) // m2 * m1
        return n
    return f


a = ct(1, 2, 4, 6)
b = ct(1, 3, 2, 6)
c = ct(2, 3, 4, 6)
GEN = {"a": a, "b": b, "c": c}


def C(n):
    return 3 * n + 1 if n % 2 else n // 2


def T(n):
    return (3 * n + 1) // 2 if n % 2 else n // 2


def v2(n):
    n = abs(n)
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def oddpart(n):
    while n % 2 == 0:
        n //= 2
    return n


def colour_of_edge(u, v):
    """colours (set) of the generator(s) joining u and v."""
    return {g for g, f in GEN.items() if f(u) == v and u != v}


print("=" * 100)
print("K1  generators, vertex types, Schreier graph = undirected C-graph on Z \\ 0(6)")
print("=" * 100)
W = [n for n in range(-X, X + 1)]
for g, f in GEN.items():
    for n in W:
        check(f(f(n)) == n, f"{g} involution at {n}")
print("  a, b, c are involutions on |n| <=", X)
# disjointness of the class pairs (gcd test)
for (r1, m1, r2, m2) in ((1, 2, 4, 6), (1, 3, 2, 6), (2, 3, 4, 6)):
    check((r1 - r2) % math.gcd(m1, m2) != 0, "classes disjoint")
print("  class pairs 1(2)|4(6), 1(3)|2(6), 2(3)|4(6) are disjoint")
typetab = {}
for n in W:
    r = n % 6
    moved = tuple(sorted(g for g, f in GEN.items() if f(n) != n))
    typetab.setdefault(r, set()).add(moved)
for r in range(6):
    check(len(typetab[r]) == 1, f"type of residue {r} mod 6 is constant")
TYPE = {r: next(iter(typetab[r])) for r in range(6)}
print("  residue mod 6 -> moving generators (degree), semi-edges:")
for r in range(6):
    mv = TYPE[r]
    print(f"    {r}: moved by {''.join(mv) or '-':3s} degree {len(mv)}   semi-edges {''.join(sorted(set('abc') - set(mv))) or '-'}")
check(TYPE == {0: (), 1: ("a", "b"), 2: ("b", "c"), 3: ("a",), 4: ("a", "b", "c"), 5: ("a", "c")},
      "vertex type table")
# Schreier edges (non-loop) vs undirected C-graph induced on Z \ 0(6), both endpoints in window
Wset = set(W)
V6 = [n for n in W if n % 6 != 0]
sch = set()
for n in V6:
    for g, f in GEN.items():
        m = f(n)
        if m != n and m in Wset:
            sch.add((min(n, m), max(n, m), g))
sch_edges = Counter((u, v) for (u, v, g) in sch)
cg = Counter()
for n in V6:
    m = C(n)
    if m in Wset and m % 6 != 0 and m != n:
        cg[(min(n, m), max(n, m))] += 1
check(sch_edges == cg, "Schreier multigraph == undirected C multigraph on Z \\ 0(6) (window)")
print(f"  Schreier multigraph == undirected C-multigraph on (Z \\ 0(6)) cap [-{X},{X}]: {sum(cg.values())} edges")
print("  the only multi-edge:", [e for e, k in cg.items() if k > 1], "(the C 2-cycle -1 -> -2 -> -1)")
check([e for e, k in cg.items() if k > 1] == [(-2, -1)], "unique double edge")
# C maps N \ 0(6) into itself and C^{-1}(0(3)) = 0(6)
for n in range(1, X + 1):
    if n % 6:
        check(C(n) % 6 != 0 and C(n) >= 1, "C maps N\\0(6) into itself")
    if C(n) % 3 == 0:
        check(n % 6 == 0, "C^{-1}(0(3)) = 0(6)")
print("  C maps N \\ 0(6) into itself; C^{-1}(0(3)) = 0(6): multiples of 3 are transient (0(6) rays hang off the leaves 3(6))")

print()
print("=" * 100)
print("K2  cycles and their colours")
print("=" * 100)


def find_cycle_from(n, f, bound=10 ** 60):
    seen = {}
    path = []
    x = n
    while x not in seen:
        if abs(x) > bound:
            return None
        seen[x] = len(path)
        path.append(x)
        x = f(x)
    return path[seen[x]:]


def canon_cycle(cyc):
    k = cyc.index(min(cyc, key=abs))
    return tuple(cyc[k:] + cyc[:k])


cycles_pos = set()
reach1 = 0
known = {1, 2, 4}
for n in range(1, X1 + 1):
    x = n
    path = []
    while x not in known:
        path.append(x)
        x = C(x)
    known.update(path)
    reach1 += 1
check(reach1 == X1, "every n <= X1 reaches {1,2,4}")
print(f"  every n in [1,{X1}] reaches the cycle (1 4 2); on N \\ 0(6) the only cycle met is the triangle")
tri = [(1, 4), (4, 2), (2, 1)]
tri_cols = [colour_of_edge(u, v) for u, v in tri]
print("  triangle colours:", [(u, v, ''.join(sorted(s))) for (u, v), s in zip(tri, tri_cols)])
check(tri_cols == [{"a"}, {"c"}, {"b"}], "rainbow triangle a,c,b")
neg_cycles = set()
cyc_of = {}
for n in range(1, XN + 1):
    x = -n
    path = []
    pos = {}
    while x not in cyc_of and x not in pos:
        pos[x] = len(path)
        path.append(x)
        x = C(x)
    if x in pos:
        cz = canon_cycle(path[pos[x]:])
        neg_cycles.add(cz)
        cid = cz
    else:
        cid = cyc_of[x]
    for y in path:
        cyc_of[y] = cid
neg_cycles = sorted(neg_cycles, key=len)
print(f"  C-cycles met from [-{XN},-1]: {len(neg_cycles)}")
for cyc in neg_cycles:
    cols = []
    L = len(cyc)
    for i in range(L):
        u, v = cyc[i], cyc[(i + 1) % L]
        s = colour_of_edge(u, v)
        if len(s) == 2:   # digon: both colours, take the one belonging to u -> v in the C-direction
            s = {g for g in s if (g == "a") == (u % 2 != 0)}
        cols.append(''.join(sorted(s)))
    print(f"    length {L:2d}: {cyc}")
    print(f"               colour word (C-direction): {''.join(cols)}   counts {dict(Counter(cols))}")
check([len(z) for z in neg_cycles] == [2, 5, 18], "negative cycles have lengths 2, 5, 18")
check(neg_cycles[0] == (-1, -2), "digon")
print("  (0 is the C-fixed point; 0 in 0(6) is removed by Kohl)")

print()
print("=" * 100)
print("K3  Kempe chains (orbits of the dihedral groups <b,c>, <a,b>, <a,c>)")
print("=" * 100)


def orbit2(n, f, g, cap=10 ** 7, maxlen=4000):
    """orbit of the dihedral group <f,g> containing n (finite chains); returns (set, closed?) or None if too long."""
    seen = {n}
    frontier = [n]
    while frontier:
        x = frontier.pop()
        for h in (f, g):
            y = h(x)
            if y not in seen:
                if len(seen) > maxlen or abs(y) > cap ** 3:
                    return None
                seen.add(y)
                frontier.append(y)
    # closed chain: every vertex has two distinct neighbours inside and no semi-edge
    closed = all(f(x) != x and g(x) != x for x in seen)
    return seen, closed


def P_index(v):
    """odd n with v in the <a,b>-chain P_n."""
    if v % 2:
        return v
    if v % 6 == 4:
        return (v - 1) // 3
    if v % 12 == 2:
        return v // 2
    if v % 12 == 8:
        return (v - 2) // 6
    raise ValueError(v)


def P_chain(n):
    s = {n, 3 * n + 1, 6 * n + 2}
    if n % 6 == 1:
        s.add(2 * n)
    return s


def run_start(v):
    """start y (odd, y = 1 or 3 mod 6) of the <a,c>-chain through v (v != -1,-2)."""
    if v % 6 == 2:          # end of a run: walk back via c then a
        v = (2 * v - 1) // 3  # c: v -> 2v (in 4(6)); a: 2v -> (2v-1)/3
    elif v % 6 == 4:
        v = (v - 1) // 3
    while v % 6 == 5:
        v = (2 * v - 1) // 3
    return v


def run_chain(y):
    k = v2(y + 1)
    s = set()
    x = y
    for _ in range(k):
        s.add(x)
        s.add(3 * x + 1)
        x = (3 * x + 1) // 2
    s.add(x)
    return s, k, x


nbc_ray = nab = nac = 0
closed_found = []
for n in V6:
    # <b,c>: doubling ray through the odd part
    o = oddpart(n)
    if o % 3 == 0:
        check(n == o and b(n) == n and c(n) == n, "<b,c> singleton at odd multiples of 3")
    else:
        # walk down to o by the unique b/c edge towards the odd base, check alternation
        x = n
        while x % 2 == 0:
            y = b(x) if b(x) == x // 2 else c(x)
            check(y == x // 2, "b/c down-step is halving")
            x = y
        check(x == o, "ray base")
        col0 = "b" if o % 3 == 1 else "c"
        col1 = "c" if col0 == "b" else "b"
        check(GEN[col0](o) == 2 * o and GEN[col1](o) == o, "ray colour at base")
        nbc_ray += 1
    # <a,b>
    res = orbit2(n, a, b)
    check(res is not None and not res[1], "<a,b> chains are finite paths")
    check(res[0] == P_chain(P_index(n)), f"<a,b> chain formula at {n}")
    nab += 1
    # <a,c>
    if n in (-1, -2):
        res = orbit2(n, a, c)
        check(res[0] == {-1, -2} and res[1], "the digon is a closed <a,c> chain")
        closed_found.append(("ac", tuple(sorted(res[0]))))
        continue
    res = orbit2(n, a, c)
    check(res is not None and not res[1], f"<a,c> chain through {n} is a finite path")
    y = run_start(n)
    ch, k, end = run_chain(y)
    check(y % 2 == 1 and y % 6 in (1, 3), "run start type")
    check(res[0] == ch, f"<a,c> chain formula at {n}")
    check(len(ch) == 2 * k + 1 and end % 6 == 2, "run length 2k and end in 2(6)")
    nac += 1
print(f"  checked every vertex of (Z \\ 0(6)) cap [-{X},{X}] ({len(V6)} vertices):")
print("   <b,c>: the doubling ray {o 2^j} of the odd part o (b/c alternate: 2 = -1 mod 3), singleton if 3 | o")
print("   <a,b>: P_n = {2n (if n = 1 mod 6), n, 3n+1, 6n+2}, a path with 3 or 2 edges and exactly one a-edge")
print("   <a,c>: the odd run y, 3y+1, T(y), ..., T^k(y) with k = v2(y+1), y = 1,3 mod 6; 2k edges")
print("  closed Kempe chains found:", closed_found)
check(closed_found == [("ac", (-2, -1)), ("ac", (-2, -1))], "unique closed chain (seen from both its vertices)")
# Banach reading: <a,c> backwards = E, <b,c> = D
for y in range(-X, X):
    if y % 6 == 5 and y != -1:
        # c then a from y gives E(y)
        check(a(c(y)) == (2 * y - 1) // 3 and (2 * y - 1) % 3 == 0, "c then a = E on 5(6)")
print("  two Kempe steps c,a from x = 5 mod 6 give E(x) = (2x-1)/3; b/c steps along a ray are D(x) = 2x:")
print("  <a,c>-chains are maximal E-orbit segments, <b,c>-chains are D-orbits; D fixes 0 (in 0(6), deleted),")
print("  E fixes -1 (the digon). v2(y+1) = -log2|y+1|_2: the <a,c> chain through y has 2k edges with 2^-k = |y+1|_2.")
for k in range(1, 12):
    ch, kk, end = run_chain(2 ** k - 1)
    check(kk == k and end == 3 ** k - 1, "Mersenne chain ends at 3^k - 1")
print("  Mersenne start y = 2^k - 1 opens the <a,c> chain of length 2k ending at 3^k - 1 (k = 1..11 checked)")
# distribution of <a,c> run parameter k over starts y in [1, X]
kd = Counter(v2(y + 1) for y in range(1, X + 1) if y % 6 in (1, 3))
print("  run-parameter k = v2(y+1) over starts y <= X:", dict(sorted(kd.items())[:8]))

print()
print("=" * 100)
print("K4  Kempe quotients")
print("=" * 100)
# K/<b,c>: a-edges between rays
syr_ok = all(oddpart(3 * n + 1) == oddpart(a(n)) for n in range(-X + 1, X, 2))
check(syr_ok, "a-edge joins ray(n) to ray(S(n))")
print("  K/<b,c>: vertices = odd numbers (rays), edges = a-edges {n, S(n)}, S = Syracuse map: the undirected Syracuse graph")
# K/<a,b>: c-edges between P-chains
qab = Counter()
for m in range(-3 * X, 3 * X + 1):
    if m % 3 == 2 and abs(2 * m) <= 6 * X:
        u, v = P_index(m), P_index(2 * m)
        qab[(u, v)] += 1


def Tn(n):
    return (3 * n + 1) // 2 if n % 2 else n // 2


def qab_formula(u, v):
    return (v == 4 * u + 1 or u == 4 * v + 1
            or (u % 4 == 3 and v == Tn(u)) or (v % 4 == 3 and u == Tn(v))
            or (u % 8 == 1 and v == Tn(Tn(u))) or (v % 8 == 1 and u == Tn(Tn(v))))


bad = [(u, v) for (u, v) in qab if not qab_formula(u, v)]
check(not bad, f"K/<a,b> edge formula {bad[:5]}")
deg = Counter()
for (u, v), k in qab.items():
    deg[u] += k
    deg[v] += k
for n in range(-X // 8, X // 8):
    if n % 2:
        check(deg[n] == (3 if n % 3 else 2), f"K/<a,b> degree at {n}: {deg[n]}")
print("  K/<a,b>: vertices = odd n (chains P_n), edges = c-edges; every edge is {n,4n+1}, {n,T(n)} (n=3 mod 4)")
print("           or {n,T^2(n)} (n=1 mod 8); degree 3 if 3 does not divide n, 2 if n = 3 mod 6 (checked |n| <= X/8)")
print("           P_1 carries a loop (the c-edge {2,4} of the triangle) and the edge {1,5}")
check(qab[(1, 1)] == 1 and deg[1] == 3, "loop at P_1")
# K/<a,c>: b-edges between runs, counted on the independently computed chain (BFS orbit of <a,c>)
for y in range(1, X // 4):
    if y % 6 in (1, 3):
        chain, closed = orbit2(y, a, c)
        nb = sum(1 for v in chain if b(v) != v)
        pred = 1 + (1 if y % 6 == 1 else 0) + v2(y + 1)
        check(nb == pred, f"K/<a,c> degree at {y}: {nb} vs {pred}")
print("  K/<a,c>: vertices = odd runs (start y = 1,3 mod 6), edges = b-edges; degree 1 + [y = 1 mod 6] + v2(y+1)")
print("  connectivity of K <=> connectivity of each quotient (chains are connected); the loop/triangle survives in each")

print()
print("=" * 100)
print("K5  the Z2 x Z2 boundary of the Tait colouring (a=(1,0), b=(0,1), c=(1,1))")
print("=" * 100)
VEC = {"a": (1, 0), "b": (0, 1), "c": (1, 1)}
NAME = {(0, 0): "0", (1, 0): "a", (0, 1): "b", (1, 1): "c"}
dtab = {}
for n in V6:
    s = (0, 0)
    for g, f in GEN.items():
        if f(n) != n:
            s = (s[0] ^ VEC[g][0], s[1] ^ VEC[g][1])
    dtab.setdefault(n % 6, set()).add(NAME[s])
print("  boundary delta(v) = sum of the colours at v = sum of the missing colours:")
for r in range(1, 6):
    print(f"    v = {r} mod 6: delta = {sorted(dtab[r])}")
check({r: dtab[r] for r in range(1, 6)} == {1: {"c"}, 2: {"a"}, 3: {"a"}, 4: {"0"}, 5: {"b"}}, "boundary table")
print("  the colouring is a nowhere-zero Z2^2-flow exactly at the branch points 4(6) (two C-preimages);")
print("  every other vertex is a source: 1,2,5 mod 6 (degree 2) and the leaves 3(6) (degree 1)")

print()
print("=" * 100)
print("K6  switches: Kempe switches keep the graph; nu and single-colour conjugations change it")
print("=" * 100)
# (i) Kempe switch on the <a,b>-chains P_n with n = 1 mod 4 (an rcwa-definable union)
inS = lambda v: v % 6 != 0 and P_index(v) % 4 == 1
a2 = lambda v: b(v) if inS(v) else a(v)
b2 = lambda v: a(v) if inS(v) else b(v)
for n in V6:
    check(a2(a2(n)) == n and b2(b2(n)) == n, "switched generators are involutions")
E1 = {(min(n, f(n)), max(n, f(n))) for n in V6 for f in (a, b, c) if f(n) != n}
E2 = {(min(n, f(n)), max(n, f(n))) for n in V6 for f in (a2, b2, c) if f(n) != n}
check(E1 == E2, "Kempe switch preserves the Schreier graph")
diff = sum(1 for n in V6 if a2(n) != a(n))
# modulus of a2: find the least M such that a2 is affine on each class mod M
def rcwa_modulus(f, cands=(2, 4, 6, 8, 12, 24, 48, 72, 96, 144, 288)):
    for M in cands:
        ok = True
        for r in range(M):
            pts = [r + M * k for k in range(1, 40) if (r + M * k) % 6 != 0]
            if len(pts) < 3:
                continue
            x0, x1, x2 = pts[0], pts[1], pts[2]
            if (f(x1) - f(x0)) * (x2 - x1) != (f(x2) - f(x1)) * (x1 - x0):
                ok = False
                break
        if ok:
            return M
    return None
print(f"  Kempe switch a<->b on the chains P_n, n = 1 mod 4: new involutions a', b' (rcwa, modulus of a' = {rcwa_modulus(a2)}),")
print(f"  a' != a at {diff} window points, Schreier graph identical ({len(E1)} edges): transitivity is unchanged.")
# (ii) nu-conjugation and the sheet switch h0
nu = lambda n: -n
h0 = ct(2, 6, 4, 6)
a_m = ct(1, 2, 2, 6)   # n <-> 3n-1, n odd
for n in W:
    check(nu(a(nu(n))) == a_m(n), "nu a nu = tau_{1(2),2(6)}")
    check(nu(b(nu(n))) == c(n) and nu(c(nu(n))) == b(n), "nu b nu = c, nu c nu = b")
    check(h0(a(h0(n))) == a_m(n), "a^h0 = tau_{1(2),2(6)}")
print("  nu a nu = a^h0 = tau_{1(2),2(6)} (n <-> 3n-1), nu b nu = c, nu c nu = b, with h0 = tau_{2(6),4(6)}:")
print("  nu G_C nu = <a^h0, c, b> = G_C^- (the 3n-1 group); the sheet swap = one-colour conjugation of a by the")
print("  '+-2' class transposition h0, composed with the global Kempe relabelling b <-> c (which keeps the graph).")
# G_C^- on N: cycles of 3n-1
Cm = lambda n: 3 * n - 1 if n % 2 else n // 2
cyc_m = sorted({canon_cycle(find_cycle_from(n, Cm)) for n in range(1, 20001)}, key=len)
print("  3n-1 C-cycles met from [1,20000]:", [cz[:6] + (('...',) if len(cz) > 6 else ()) for cz in cyc_m])
check([len(z) for z in cyc_m] == [2, 5, 18], "three 3n-1 cycles")
for cz in cyc_m:
    check(all(x % 6 != 0 for x in cz), "cycles avoid 0(6)")
print("  => G_C^- = <tau_{1(2),2(6)}, b, c> has >= 3 orbits on N \\ 0(6) (PROVED: three distinct cycles)")
# (iii) census of single-colour switches a -> a^h, h = tau of two disjoint classes inside 2(6) u 4(6)


def classes_inside(m):
    return [(r, m) for r in range(m) if r % 6 in (2, 4)]


def disjoint(c1, c2):
    return (c1[0] - c2[0]) % math.gcd(c1[1], c2[1]) != 0


def meets46(cl):
    r, m = cl
    return any((r + m * k) % 6 == 4 for k in range(6))


def orbit_stats(F, starts, cap=10 ** 30, maxsteps=30000):
    known = {}
    cycles = {}
    esc = 0
    for s in starts:
        x = s
        pos = {}
        path = []
        st = 0
        res = None
        while True:
            if x in known:
                res = known[x]
                if res == "ESC":
                    esc += 1
                break
            if x in pos:
                cyc = path[pos[x]:]
                res = min(cyc)
                cycles[res] = len(cyc)
                break
            if abs(x) > cap or st > maxsteps:
                esc += 1
                res = "ESC"
                break
            pos[x] = len(path)
            path.append(x)
            x = F(x)
            st += 1
        for y in path:
            known[y] = res
    return cycles, esc


def Fraction_str(p, q):
    g = math.gcd(p, q)
    p, q = p // g, q // g
    return f"{p}" if q == 1 else f"{p}/{q}"


def haar_drift(F, L, depth=18):
    """exact Haar mean over odd n of log(F(n)/n) - v2(F(n)) log 2, i.e. the log-drift of the odd-to-odd map,
    computed on the odd classes mod L*2^depth (v2 capped at depth, error <= 2^-depth)."""
    M = L * 2 ** depth
    s = 0.0
    cnt = 0
    for r in range(1, M, 2):
        n1, n2 = r + M, r + 2 * M
        slope = (F(n2) - F(n1)) / M
        k = min(v2(F(n1)), depth)
        s += math.log(slope) - k * math.log(2)
        cnt += 1
    return s / cnt


def descent_certificate(F, L, depth=8, nsmall=5000):
    """PROVED-transitivity test for F (odd n -> F(n), even n -> n/2): on every odd class n = r mod M
    (M = L 2^depth) F(n) = s n + t and v2(F(n)) >= k := min(v2(F(r)), depth); if s < 2^k on every class,
    the next odd iterate is < n for n > n0 = max t/(2^k - s); then every orbit enters [1, n0], and a direct
    check that [1, max(n0, nsmall)] flows into the cycle through 1 proves a single grand orbit on N."""
    M = L * 2 ** depth
    n0 = 0
    for r in range(1, M, 2):
        s_ = Fr(F(r + M) - F(r), M)
        t_ = F(r) - s_ * r
        k = min(v2(F(r)) if F(r) else depth, depth)
        if s_ >= 2 ** k:
            return False, None
        n0 = max(n0, math.ceil(t_ / (2 ** k - s_)) if t_ > 0 else 0)
    top = max(n0, nsmall)
    good = {1}
    for n in range(1, top + 1):
        x = n
        path = []
        while x not in good:
            if x <= 0 or len(path) > 10 ** 5:
                return False, None
            path.append(x)
            x = F(x)
        good.update(path)
    return True, n0


def ascent_certificate(F, L, depth=8):
    """PROVED-divergence test: on every odd class n = r mod M (M = L 2^depth) the valuation k = v2(F(n)) is
    constant (v2(F(r)) < depth) and s > 2^k; then the next odd iterate exceeds n for n > n1, so every odd
    start above n1 has a strictly increasing odd subsequence (a divergent orbit).  Returns (ok, n1)."""
    M = L * 2 ** depth
    n1 = 0
    for r in range(1, M, 2):
        s_ = Fr(F(r + M) - F(r), M)
        t_ = F(r) - s_ * r
        k = v2(F(r)) if F(r) else depth
        if k >= depth or s_ <= 2 ** k:
            return False, None
        n1 = max(n1, math.ceil(-t_ / (s_ - 2 ** k)) if t_ < 0 else 0)
    return True, n1


print("  single-colour switch census: h = tau_{c1,c2}, c1,c2 disjoint classes inside 2(6) u 4(6), support meeting 4(6);")
print("  F_h(n) = h(3n+1) (n odd), n/2 (n even) is the map whose functional graph is the Schreier graph of <a^h,b,c>.")
print("  columns: up-slopes | Haar log-drift of the odd->odd map | cycles met from [1,3000] | escapes past 1e30")
rows = []
cls12 = [cl for m in (6, 12) for cl in classes_inside(m)]
pairs12 = [(c1, c2) for c1, c2 in itertools.combinations(cls12, 2) if disjoint(c1, c2) and (meets46(c1) or meets46(c2))]
for c1, c2 in pairs12:
    h = ct(c1[0], c1[1], c2[0], c2[1])
    F = lambda n, h=h: h(3 * n + 1) if n % 2 else n // 2
    L = 2 * math.lcm(c1[1], c2[1])
    upmap = []
    for r in range(1, L, 2):
        s_ = Fr(F(r + L) - F(r), L)
        t_ = F(r) - s_ * r
        upmap.append((r, s_, t_))
    pieces = sorted({(s_, t_) for (r, s_, t_) in upmap})
    upstr = []
    for (s_, t_) in pieces:
        cls_ = sorted(r for (r, s2, t2) in upmap if (s2, t2) == (s_, t_))
        # express the residues of n (odd) mod L that use this piece, reduced to the smallest modulus
        for mod in (2, 4, 8, 12, 24, 48):
            if L % mod == 0 and all(((r % mod) in {x % mod for x in cls_}) for r in cls_) and \
                    len({x % mod for x in cls_}) * (L // mod) == len(cls_):
                cl_txt = ",".join(str(x) for x in sorted({x % mod for x in cls_})) + f" mod {mod}"
                break
        num = s_.numerator
        den = s_.denominator
        tt = t_ * den
        expr = f"{num}n{int(tt):+d}" if den == 1 else f"({num}n{int(tt):+d})/{den}"
        upstr.append(f"{expr} [n={cl_txt}]")
    upmap = upstr
    dr = haar_drift(F, L, depth=12)
    cycles, esc = orbit_stats(F, range(1, 3001))
    okd, n0 = descent_certificate(F, L) if (esc == 0 and len(cycles) == 1) else (False, None)
    oka, n1 = ascent_certificate(F, L) if esc > 0 else (False, None)
    if oka:
        kind = f"DRIFT-type, PROVED divergent beyond n1 = {n1} (one-step ascent) => PROVED intransitive"
    elif esc > 0 and dr > 0:
        kind = "DRIFT-type (escapes; heuristic)"
    elif esc == 0 and len(cycles) >= 2:
        kind = "SHEET-type (>=2 cycles, PROVED intransitive)"
    elif okd:
        kind = f"PROVED transitive (odd-step descent beyond n0 = {n0})"
    elif esc == 0 and len(cycles) == 1:
        kind = "Collatz-like (1 cycle; OPEN)"
    else:
        kind = "mixed"
    rows.append((c1, c2, upmap, dr, cycles, esc, kind))
    cyc_show = sorted(cycles.items())[:4]
    print(f"   tau_{{{c1[0]}({c1[1]}),{c2[0]}({c2[1]})}}: odd n -> {'; '.join(upmap)}")
    print(f"        Haar drift {dr:+.3f}  cycles {len(cycles)} (min, length) {cyc_show}  escapes {esc:4d}  => {kind}")
check(rows[0][0] == (2, 6) and rows[0][1] == (4, 6) and rows[0][6].startswith("SHEET"), "the h0 switch is the SHEET control")
check(not descent_certificate(lambda n: 3 * n + 1 if n % 2 else n // 2, 2)[0], "Collatz itself has no descent certificate")
print("  (control: Collatz itself fails the descent test -- n = 3 mod 4 has s/2^k = 3/2)")
# the four mod-4 sign strategies of the 3n+-1 game as one-colour switches
sq = {}
for (c1, c2, upmap, dr_, cyc_, esc_, kind_) in rows:
    sq[(c1, c2)] = kind_
print("  the four mod-4 sign strategies of Althofer's 3n+-1 game are one-colour switches of G_C:")
print("    (+,+) identity           : 3n+1                       -> Collatz (OPEN)")
print("    (-,-) tau_{2(6),4(6)}    : 3n-1                       ->", sq[((2, 6), (4, 6))])
print("    (-,+) tau_{2(12),4(12)}  : 3n-1 (n=1 mod 4), 3n+1 (n=3 mod 4): one halving each ->", sq[((2, 12), (4, 12))])
print("    (+,-) tau_{8(12),10(12)} : 3n+1 (n=1 mod 4), 3n-1 (n=3 mod 4): 4 | 3n+-1        ->", sq[((8, 12), (10, 12))])
check(sq[((2, 12), (4, 12))].startswith("DRIFT-type, PROVED") and sq[((8, 12), (10, 12))].startswith("PROVED transitive"),
      "the strategy square")
print("  => Collatz is one corner of a square whose three other corners are PROVED (two intransitive, one transitive).")
# larger census counts
cls24 = [cl for m in (6, 12, 18, 24) for cl in classes_inside(m)]
pairs24 = [(c1, c2) for c1, c2 in itertools.combinations(cls24, 2) if disjoint(c1, c2) and (meets46(c1) or meets46(c2))]
kinds = Counter()
for c1, c2 in pairs24:
    h = ct(c1[0], c1[1], c2[0], c2[1])
    F = lambda n, h=h: h(3 * n + 1) if n % 2 else n // 2
    L = 2 * math.lcm(c1[1], c2[1])
    dr = haar_drift(F, L, depth=8)
    cycles, esc = orbit_stats(F, range(1, 1501))
    if esc > 0 and ascent_certificate(F, L, depth=6)[0]:
        kinds["DRIFT-type, PROVED divergent (one-step ascent)"] += 1
    elif esc > 0 and dr > 0:
        kinds["DRIFT-type (positive drift, escapes; heuristic)"] += 1
    elif esc == 0 and len(cycles) >= 2:
        kinds["SHEET-type (>=2 cycles: PROVED intransitive)"] += 1
    elif esc == 0 and len(cycles) == 1 and descent_certificate(F, L, depth=6, nsmall=1500)[0]:
        kinds["PROVED transitive (odd-step descent)"] += 1
    elif esc == 0 and len(cycles) == 1:
        kinds["Collatz-like (1 cycle, no escapes; OPEN)"] += 1
    else:
        kinds["mixed (escapes with negative drift, or none)"] += 1
print(f"  moduli 6..24: {len(pairs24)} switches:", dict(kinds))
print("  no class-transposition switch is a DEFECT control: a nonidentity rcwa change has positive-density support.")

print()
print("=" * 100)
print("K7  Kohl's mod-3 shape <a_q, b, c> tracks C_q (qn+1) only when 3 | q; the mod-p generalisation")
print("=" * 100)


for q in (3, 5, 9):
    aq = ct(1, 2, (q + 1) % (2 * q), 2 * q)
    # which up-images land in 0(3)?
    hits = [n for n in range(1, 200, 2) if (q * n + 1) % 3 == 0][:4]
    print(f"  q = {q}: a_q = tau_{{1(2),{(q + 1) % (2 * q)}({2 * q})}}; odd n with qn+1 = 0 mod 3: {hits or 'none'}")
# q = 5 explicit cut: the C_5 path 7 -> 36 -> 18 -> 9 -> 46
aq5 = ct(1, 2, 6, 10)
print("  q = 5: C_5 path 7 -> 36 -> 18 -> 9 -> 46; in <a_5,b,c>: 36 has neighbours",
      sorted({f(36) for f in (aq5, b, c)} - {36}), ", 18 has", sorted({f(18) for f in (aq5, b, c)} - {18}) or "none (isolated)")
check(sorted({f(18) for f in (aq5, b, c)} - {18}) == [], "18 isolated in the q=5 Kohl shape")
check(sorted({f(36) for f in (aq5, b, c)} - {36}) == [7], "36 a leaf")
print("  so the q = 5 Kohl-shape (mod-3) group cuts C_5-orbits (18 is a fixed point of all three generators):")
print("  Kohl's mod-3 shape needs 3 | q (multiples of 3 transient) as well as 2 = -1 mod 3 (doubling flips the sign).")
# the mod-p generalisation: colour {m,2m} by the parity of the discrete log of m in (Z/p)^x / <2>-cosets
def modp_presentation(q, p):
    """doubling colours b_p, c_p from residues mod p (needs -1 in <2> mod p), up colour a_q = tau_{1(2),(q+1)(2q)}."""
    orb = {}
    for r0 in range(1, p):
        if r0 in orb:
            continue
        x, k = r0, 0
        while x not in orb:
            orb[x] = k % 2
            x = (2 * x) % p
            k += 1
    ok = all(orb[(2 * r) % p] != orb[r] for r in range(1, p))
    return ok, orb


tab = []
for pp in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
    ok, _ = modp_presentation(pp, pp)
    minus1 = any(pow(2, k, pp) == pp - 1 for k in range(1, pp))
    check(ok == minus1, "proper residue 2-colouring of doubling mod p <=> -1 in <2> mod p")
    tab.append((pp, ok))
print("  mod-p generalisation: {m,2m} (p not dividing m) can be 2-coloured properly by m mod p iff -1 is a power of 2")
print("  mod p (the doubling cycles on (Z/p)^x have even length):", ", ".join(f"p={a}:{'yes' if b else 'no'}" for a, b in tab))
ok5, col5 = modp_presentation(5, 5)
b5 = lambda m: (2 * m if (m % 5 and col5[m % 5] == 0) else (m // 2 if (m % 10 and m % 2 == 0 and m % 5 and col5[(m // 2) % 5] == 0) else m))
c5 = lambda m: (2 * m if (m % 5 and col5[m % 5] == 1) else (m // 2 if (m % 10 and m % 2 == 0 and m % 5 and col5[(m // 2) % 5] == 1) else m))
# b5, c5 as involutions: m <-> 2m for m of colour 0 (resp. 1); a point is either the small or the large end
def invol(colour):
    def f(m):
        if m % 5 == 0:
            return m
        up = 2 * m if col5[m % 5] == colour else None
        dn = m // 2 if (m % 2 == 0 and col5[(m // 2) % 5] == colour) else None
        check(not (up is not None and dn is not None), "involution well defined")
        return up if up is not None else (dn if dn is not None else m)
    return f
b5, c5 = invol(0), invol(1)
a5 = ct(1, 2, 6, 10)
C5 = lambda n: 5 * n + 1 if n % 2 else n // 2
W5 = set(range(-20000, 20001))
sch5 = Counter()
for n in W5:
    if n % 10 == 0:
        continue
    for f in (a5, b5, c5):
        check(f(f(n)) == n, "q=5 generators are involutions")
        m = f(n)
        if m != n and m in W5 and n < m:
            sch5[(n, m)] += 1
cg5 = Counter()
for n in W5:
    if n % 10 == 0:
        continue
    m = C5(n)
    if m in W5 and m % 10 != 0 and m != n:
        cg5[(min(n, m), max(n, m))] += 1
check(sch5 == cg5, "mod-5 Tait group of 5n+1 has Schreier graph = C_5-graph on Z \\ 0(10)")
for n in range(1, 20001):
    if C5(n) % 5 == 0:
        check(n % 10 == 0, "C_5^{-1}(0(5)) = 0(10)")
print("  q = 5 with p = 5: <tau_{1(2),6(10)}, b_5, c_5> (b_5, c_5 = doubling split by the parity of log_2 m mod 5)")
print("  has Schreier graph = the C_5-graph on Z \\ 0(10) (window |n| <= 20000), multiples of 5 transient;")
print("  it is intransitive on N \\ 0(10) (cycles through 1, 13, 17): the Tait framework is DRIFT-blind as well.")

print()
print("=" * 100)
print("K8  Kohl's second group G_T")
print("=" * 100)
aT, bT, cT = ct(0, 2, 1, 2), ct(1, 2, 2, 4), ct(1, 4, 2, 6)
X2 = 20000
edgesT = Counter()
for n in range(0, X2 + 1):
    for g, f in (("a", aT), ("b", bT), ("c", cT)):
        m = f(n)
        if m > n and m <= X2:
            edgesT[(n, m, g)] += 1
parent = list(range(X2 + 1))


def findT(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


multi = Counter((u, v) for (u, v, g) in edgesT)
for (u, v, g) in edgesT:
    ra, rb = findT(u), findT(v)
    if ra != rb:
        parent[ra] = rb
ncomp = len({findT(x) for x in range(0, X2 + 1)})
beta = len(edgesT) - (X2 + 1) + ncomp
print(f"  Schreier graph of G_T on [0,{X2}]: {len(edgesT)} edges, {ncomp} window components, cyclomatic number {beta}")
print("  multi-edges:", [(e, [g for (u, v, g) in edgesT if (u, v) == e]) for e, k in multi.items() if k > 1])
check(beta == 1 and [e for e, k in multi.items() if k > 1] == [(1, 2)], "G_T: unique cycle = b/c digon {1,2}")
check(all(aT(2 * k) == 2 * k + 1 for k in range(1000)), "a = consecutive pairing")
print("  colour a of G_T is the consecutive pairing {2k, 2k+1} (THM-4470's 3n-1 pairing); the graph is a Tait-")
print("  coloured tree plus one digon on the window (FINITE-EXACT).")
print()
print("ALL CHECKS PASSED (kohl_tait)")
