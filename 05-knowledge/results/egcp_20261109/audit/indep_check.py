"""Independent audit checks for THM-5040 (audit code only; nothing imported from egcp_lib.py).

usage: python indep_check.py [NBUILD]   (NBUILD = largest order for the end-to-end builder check, default 11)
Sections:
  A  validation of the audit machinery (generator-based evenness vs all n! permutations; two d algorithms)
  B  Lemma 0 (complement formula, union additivity, C_4-free formula)
  C  Lemma 5' (multi-join sign rule) on random multi-joins with repeated blocks
  D  Lemma 5 (E_2-multi-join), including G with a universal vertex and G a join
  E  Lemma 7 (SR), including G with isolated vertices and K_2 components
  F  Lemma 8 (jump family) value formula and criterion (iii), including H with universal vertices other than z
  G  Lemma 9 frontier families F-a, F-b, F-c (and the negative control rho = 4 for F-c)
  H  SR chains (Theorem T old j, Lemma 10) from the census L graphs up to order 15
  I  end-to-end: every x needing exactly n vertices, 10 <= n <= NBUILD, built by the Section 4 recursion
     (re-implemented from the text) and checked for order, d and evenness
"""
import sys, os, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_lib import *

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, '..', '..', '..', '..', '04-computation', 'experiments', 'egcp_20261109', 'base_witnesses.txt')
NBUILD = int(sys.argv[1]) if len(sys.argv) > 1 else 11
rng = random.Random(5040)
FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print("  FAIL:", msg, flush=True)


def rand_graph(n, p):
    return graph(n, [(u, v) for u in range(n) for v in range(u + 1, n) if rng.random() < p])


def leaves(g):
    return [v for v, d in enumerate(degrees(g)) if d == 1]


t0 = time.time()
# ---------------- A ----------------
print("A. machinery validation")
cnt = 0
for _ in range(400):
    n = rng.randint(1, 7)
    g = rand_graph(n, rng.random())
    check(is_even(g) == is_even_bruteforce(g), "evenness mismatch %s" % (g,))
    check(dom_count(g) == dom_count_ie(g), "d mismatch %s" % (g,))
    cnt += 1
# structured symmetric graphs too
for g in [Kn(5), En(6), Ck(6), Pk(6), cocktail(6), join(En(3), En(3)), join(En(2), En(2), En(2)),
          disjoint_union(Pk(3), Pk(3)), complement(Ck(7)), disjoint_union(Ck(4), Kn(2))]:
    check(is_even(g) == is_even_bruteforce(g), "evenness mismatch structured %s" % (g,))
    cnt += 1
print("  %d graphs: generator evenness == brute-force evenness, subset-DP d == inclusion-exclusion d" % cnt)
print("  known parities: K_2 odd %s, K_3 odd %s, K_4-e odd %s, P_3 even %s, K_33 odd %s, K_22 even %s" % (
    not is_even(Kn(2)), not is_even(Kn(3)), not is_even(complement(graph(4, [(0, 1)]))), is_even(Pk(3)),
    not is_even(join(En(3), En(3))), is_even(join(En(2), En(2)))))
print("  d(K_2), d(K_3), d(K_4 - e) =", dom_count(Kn(2)), dom_count(Kn(3)), dom_count(complement(graph(4, [(0, 1)]))))

# ---------------- B ----------------
print("B. Lemma 0")
bad = 0
for _ in range(300):
    n = rng.randint(1, 10)
    g = rand_graph(n, rng.random())
    bad += dom_count(g) != 2 ** n - D_count(complement(g))
check(bad == 0, "complement formula")
bad = 0
for _ in range(200):
    h1 = rand_graph(rng.randint(1, 6), rng.random())
    h2 = rand_graph(rng.randint(1, 6), rng.random())
    bad += D_count(disjoint_union(h1, h2)) != D_count(h1) + D_count(h2) - 1
check(bad == 0, "union additivity")


def c4free(g):
    n, adj = g
    for u in range(n):
        for v in range(u + 1, n):
            if popcount(adj[u] & adj[v]) >= 2:
                return False
    return True


bad = tested = 0
for _ in range(3000):
    n = rng.randint(1, 11)
    g = rand_graph(n, rng.uniform(0.05, 0.45))
    if not c4free(g):
        continue
    tested += 1
    deg = degrees(g)
    formula = 1 + sum(1 for d in deg if d > 0) + sum(2 ** d - 1 - d for d in deg)
    bad += formula != D_count(g)
check(bad == 0, "C4-free formula")
print("  complement formula 300, union 200, C_4-free formula %d random C_4-free graphs: failures %d" % (tested, len(FAIL)))

# ---------------- C ----------------
print("C. Lemma 5' multi-join sign rule")


def canon_small(g):
    n = g[0]
    import itertools
    best = None
    E = edges_of(g)
    for p in itertools.permutations(range(n)):
        key = tuple(sorted(tuple(sorted((p[u], p[v]))) for u, v in E))
        if best is None or key < best:
            best = key
    return (n, best)


def coconnected(g):
    h = complement(g)
    n, adj = h
    if n == 0:
        return False
    seen = 1
    frontier = 1
    while frontier:
        nxt = 0
        x = frontier
        while x:
            low = x & -x
            nxt |= adj[low.bit_length() - 1]
            x ^= low
        frontier = nxt & ~seen
        seen |= nxt
    return seen == (1 << n) - 1


pool = []
while len(pool) < 40:
    b = rand_graph(rng.randint(1, 4), rng.random())
    if coconnected(b):
        pool.append(b)
pool += [b for b in [graph(1), En(2), En(3), Pk(3), Pk(4), disjoint_union(Pk(2), graph(1))] if coconnected(b)]
assert all(coconnected(b) for b in pool)
bad = tested = 0
for _ in range(400):
    r = rng.randint(2, 4)
    blocks = []
    for _ in range(r):
        if blocks and rng.random() < 0.4:
            blocks.append(rng.choice(blocks))  # repeat an isomorphic block
        else:
            blocks.append(rng.choice(pool))
    rng.shuffle(blocks)
    J = join(*blocks)
    if J[0] > 10:
        continue
    tested += 1
    keys = [canon_small(b) for b in blocks]
    pred = all(is_even(b) for b in blocks)
    for i in range(r):
        for j in range(i + 1, r):
            if keys[i] == keys[j] and blocks[i][0] % 2 == 1:
                pred = False
    bad += pred != is_even(J)
check(bad == 0, "Lemma 5' rule")
print("  %d random multi-joins (with repeated blocks, blocks in arbitrary positions): mismatches %d" % (tested, bad))

# ---------------- D ----------------
print("D. Lemma 5 (E_2-multi-join)")
bad = tested = special = 0
for _ in range(1500):
    g0 = rand_graph(rng.randint(1, 6), rng.random())
    kind = rng.random()
    if kind < 0.3:
        g0 = cone(g0)  # G with a universal vertex
    elif kind < 0.5:
        g0 = join(g0, rand_graph(rng.randint(1, 3), rng.random()))  # G a join
    if g0[0] > 7 or not is_even(g0):
        continue
    t = rng.randint(1, 2)
    J = e2_multijoin(g0, t)
    n = J[0]
    tested += 1
    special += kind < 0.5
    ok = dom_count(J) == 2 ** n - 2 ** g0[0] - 2 * t + dom_count(g0) and is_even(J)
    bad += not ok
    if tested >= 400:
        break
check(bad == 0, "Lemma 5")
print("  %d even G (%d with a universal vertex or a join), t in {1,2}: failures %d" % (tested, special, bad))

# ---------------- E ----------------
print("E. Lemma 7 (SR)")
bad = tested = withiso = withk2 = 0
for _ in range(20000):
    n = rng.randint(3, 8)
    g = rand_graph(n, rng.uniform(0.1, 0.6))
    if rng.random() < 0.4:
        g = pad(g, n + rng.randint(1, 2))  # isolated vertices become leaves at the hub
    if not has_leaf(g) or has_universal(g) or not is_even(g):
        continue
    for q in leaves(g):
        h = SR(g, q)
        tested += 1
        withiso += 0 in degrees(g)
        nq = [v for v in range(g[0]) if (g[1][q] >> v) & 1][0]
        withk2 += degrees(g)[nq] == 1
        ok = (dom_count(h) == 3 * 2 ** (g[0] - 1) + dom_count(g) and is_even(h) and has_leaf(h)
              and not has_universal(h))
        bad += not ok
    if tested >= 500:
        break
check(bad == 0, "Lemma 7")
print("  %d (G, q) pairs (%d with isolated vertices in G, %d with q in a K_2 component): failures %d" % (
    tested, withiso, withk2, bad))

# ---------------- F ----------------
print("F. Lemma 8 (jump family)")
badv = badc = tested = unique = withW = kone = nev = 0
for _ in range(6000):
    rho = rng.randint(2, 6)
    H = rand_graph(rho, rng.uniform(0.1, 0.9))
    if rng.random() < 0.3:  # force a universal vertex of H (a degree-1 neighbour of the hub in J)
        H = cone(rand_graph(rho - 1, rng.random()))
    z = rng.randrange(rho)
    J = jump(H, z)
    h = popcount(H[1][z])
    m = J[0]
    tested += 1
    badv += m != 2 * rho + 1 - h or dom_count(J) != 2 ** (m - 1) - D_count(H) or not has_leaf(J) or has_universal(J)
    dg = degrees(J)
    if dg[rho] == m - 2 and dg.count(m - 2) == 1:  # hypothesis of (iii)
        unique += 1
        W = [v for v in range(rho) if v != z and popcount(H[1][v]) == rho - 1]
        withW += bool(W)
        kone += rho - h == 1
        crit = all(edge_sign(H, p) == perm_sign(p) for p in stabiliser_generators(H, z))
        ev = is_even(J)
        nev += ev
        badc += crit != ev
    if tested >= 2500:
        break
check(badv == 0, "Lemma 8 value")
check(badc == 0, "Lemma 8 (iii)")
print("  %d (H, z): value/leaf/no-universal failures %d; criterion (iii) on %d cases with a unique hub "
      "(%d with universal vertices of H other than z, %d with k = 1; %d even): failures %d" % (
          tested, badv, unique, withW, kone, nev, badc))
# the hypothesis rho > h + 1 implies uniqueness
bad = 0
for _ in range(500):
    rho = rng.randint(2, 7)
    H = rand_graph(rho, rng.random())
    z = rng.randrange(rho)
    h = popcount(H[1][z])
    if rho > h + 1:
        dg = degrees(jump(H, z))
        bad += not (dg[rho] == max(dg) and dg.count(max(dg)) == 1)
check(bad == 0, "rho > h + 1 => unique hub")

# ---------------- G ----------------
print("G. Lemma 9 frontier families")
rows = []
for rho in range(3, 10):
    H = path_plus_K1(rho)  # path 0-1-...-(rho-2), isolated rho-1
    fam = [('F-a', H, 0, 2 * rho, 2 ** (2 * rho - 1) - (2 * rho - 3), 2 * rho - 3)]
    if rho >= 5:
        S, ends = spider([rho - 4, 1, 2])
        fam.append(('F-b', S, ends[0], 2 * rho, 2 ** (2 * rho - 1) - (2 * rho + 1), 2 * rho + 1))
    if rho >= 4:
        fam.append(('F-c', H, 1, 2 * rho - 1, 2 ** (2 * rho - 2) - (2 * rho - 3), 2 * rho - 3))
    for name, Hh, z, m, dpred, Dpred in fam:
        if m > 17:
            continue
        J = jump(Hh, z)
        h = popcount(Hh[1][z])
        stab = stabiliser_generators(Hh, z)
        dJ = dom_count(J)
        ev = is_even(J)
        rows.append((name, rho, m, J[0], D_count(Hh), Dpred, dJ, dpred, h, len(stab), ev, has_leaf(J), has_universal(J)))
        expected_even = not (name == 'F-c' and rho == 4)
        check(J[0] == m and D_count(Hh) == Dpred and dJ == dpred and ev == expected_even and has_leaf(J)
              and not has_universal(J) and rho > h + 1 and (len(stab) == 0) == expected_even,
              "Lemma 9 %s rho=%d" % (name, rho))
for r in rows:
    print("  %-3s rho=%d m=%2d |J|=%2d D(H)=%3d (pred %3d) d=%6d (pred %6d) h=%d |stab gens|=%d even=%s leaf=%s univ=%s" % r)
# F-b as X(m-1), F-a as j = m - 3, F-c as j = m - 2
check(all(2 ** (m - 1) - (m + 1) == X(m - 1) for m in range(10, 40, 2)), "F-b = X(m-1)")

# ---------------- H ----------------
print("H. SR chains from the census L graphs")
Wt = {}
Lb = {}
for line in open(BASE):
    if line.startswith('#') or not line.strip():
        continue
    parts = line.split()
    tag, n, d = parts[0], int(parts[1]), int(parts[2])
    g = from_hex(n, parts[3] if len(parts) > 3 else '')
    if tag == 'W':
        Wt[d] = g
    else:
        Lb[(n, d)] = g
# re-verify base witnesses with audit code
badW = 0
for d, g in Wt.items():
    n = g[0]
    least = min(k for k in range(1, 10) if d in O(k))
    badW += not (dom_count(g) == d and is_even(g) and n == least)
badL = 0
for (n, d), g in Lb.items():
    badL += not (g[0] == n and dom_count(g) == d and is_even(g) and has_leaf(g) and not has_universal(g))
covered = all((m, 2 ** (m - 1) - j) in Lb for m in (8, 9) for j in Tj(m)) and \
    all((m, 2 ** (m - 1) + c) in Lb for m in (8, 9) for c in (1, 3, 5, 7, 11, 13))
check(badW == 0 and badL == 0 and covered and len(Wt) == 248 and set(Wt) == O(9), "base witnesses")
print("  base: %d W witnesses (set == O(9): %s; order == least order of O: bad %d), %d L graphs (bad %d), "
      "L covers T(8), T(9), offsets: %s" % (len(Wt), set(Wt) == O(9), badW, len(Lb), badL, covered))


def leaf_top(m, j):
    if m <= 9:
        return Lb[(m, 2 ** (m - 1) - j)]
    if j in Tj(m - 2):
        g = leaf_top(m - 2, j)
        return SR(g, leaves(g)[0])
    if m % 2 == 1 and j == m - 2:
        return jump(path_plus_K1((m + 1) // 2), 1)
    if m % 2 == 0 and j == m - 3:
        return jump(path_plus_K1(m // 2), 0)
    if m % 2 == 0 and j == m + 1:
        S, ends = spider([m // 2 - 4, 1, 2])
        return jump(S, ends[0])
    raise ValueError((m, j))


def leaf_up(m, c):
    if m <= 9:
        return Lb[(m, 2 ** (m - 1) + c)]
    g = leaf_up(m - 2, c)
    return SR(g, leaves(g)[0])


# Tj bookkeeping (Theorem T step)
for m in range(10, 200):
    new = Tj(m) - Tj(m - 2)
    check(new == ({m - 2} if m % 2 else {m - 3, m + 1}), "Tj difference m=%d" % m)
cntH = 0
for m in range(10, 16):
    for j in sorted(Tj(m)):
        g = leaf_top(m, j)
        cntH += 1
        check(g[0] == m and dom_count(g) == 2 ** (m - 1) - j and is_even(g) and has_leaf(g) and not has_universal(g),
              "leaf_top m=%d j=%d" % (m, j))
    for c in (1, 3, 5, 7, 11, 13):
        g = leaf_up(m, c)
        cntH += 1
        check(g[0] == m and dom_count(g) == 2 ** (m - 1) + c and is_even(g) and has_leaf(g) and not has_universal(g),
              "leaf_up m=%d c=%d" % (m, c))
print("  %d top-band and offset L-graphs on 10..15 vertices verified (d, even, leaf, no universal vertex)" % cntH)
print("  time so far %.1fs" % (time.time() - t0), flush=True)

# ---------------- I ----------------
print("I. end-to-end builder (Section 4 re-implemented) for 10 <= n <= %d" % NBUILD)
memo = {}


def build(n, x):
    key = (n, x)
    if key in memo:
        return memo[key]
    assert x in O(n), (n, x)
    if n <= 9:
        g = pad(Wt[x], n)
    elif x == M(n):
        g = cocktail(n) if n % 2 == 0 else cone(cocktail(n - 1))  # Lemma 6
    elif x <= 2 ** (n - 1) - n - 2:
        g = pad(build(n - 1, x), n)  # Claim_3(n-1)
    elif x < 2 ** (n - 1):  # Case 1
        j = 2 ** (n - 1) - x
        if j in Tj(n):
            g = leaf_top(n, j)
        else:
            assert x == M(n - 1)
            g = pad(build(n - 1, x), n)
    else:
        D = 2 ** n - x
        assert D >= n + 3
        if D <= 2 ** (n - 2) + 1:  # Case 2a
            g0 = 1 if n % 2 else 2
            while 2 ** g0 + n - g0 - 1 < D:
                g0 += 2
            assert g0 <= n - 2
            dG = 2 ** g0 - D + n - g0
            assert dG % 2 == 1 and 1 <= dG <= 2 ** g0 - g0 - 3
            if dG not in (3, 7, 13):
                g = e2_multijoin(build(g0, dG), (n - g0) // 2)
            elif g0 + 2 <= n - 2:
                d2 = dG + 3 * 2 ** g0 - 2
                assert d2 > 13 and d2 <= 2 ** (g0 + 2) - (g0 + 2) - 3 and d2 == 2 ** (g0 + 2) - D + n - g0 - 2
                g = e2_multijoin(build(g0 + 2, d2), (n - g0 - 2) // 2)
            else:
                assert g0 == n - 2
                c = x - 3 * 2 ** (n - 2)
                assert c in (1, 5, 11)
                g = cone(leaf_up(n - 1, c))
        else:  # Case 2b
            y = 2 ** (n - 1) - D
            assert 1 <= y <= 2 ** (n - 2) - 3
            if y <= 2 ** (n - 2) - n - 1 and y not in (3, 7, 13):
                g = cone(pad(build(n - 2, y), n - 1))
            elif y in (3, 7, 13):
                g = leaf_up(n, y)
            else:
                j = 2 ** (n - 2) - y
                assert 3 <= j <= n
                if j in Tj(n - 1):
                    g = cone(leaf_top(n - 1, j))
                else:
                    assert y == M(n - 2)
                    mg = cocktail(n - 2) if n % 2 == 0 else cone(cocktail(n - 3))
                    g = cone(pad(mg, n - 1))
    memo[key] = g
    return g


for n in range(10, NBUILD + 1):
    vals = sorted(O(n) - O(n - 1))
    bad = 0
    for x in vals:
        g = build(n, x)
        ok = g[0] == n and dom_count(g) == x and is_even(g)
        if not ok:
            bad += 1
            check(False, "build n=%d x=%d" % (n, x))
    print("  n=%d: all %d values of O(n) - O(n-1) built and verified: failures %d  (%.1fs)" % (
        n, len(vals), bad, time.time() - t0), flush=True)

print("TOTAL FAILURES:", len(FAIL))
for f in FAIL[:20]:
    print("  ", f)
