#!/usr/bin/env python3
"""
ADVERSARIAL RE-CHECK of THM-478 (digit-one degree law), mac-mini 2026-06-11.

Fresh code, written from definitions with its OWN conventions (deliberately
different tile enumeration order and bit convention from the explorer /
rm_digit_ladder script). Checks:

  P1. Pairing identity:  prod z_e + prod(1+z_e) = sum over PROPER subsets
      (mod 2), degree <= L-1.  Symbolic + brute-force evaluation, L = 3, 5.
  P2. Directed-cycle pairing: all directed odd cycles group into reversal
      pairs {C->, C<-} with NO fixed points; for L=3 the two orientations of
      a triple are distinct directed cycles; per-pair indicator identity.
  P3. Reed-Muller flat annihilation at m=4, D=1:
      (a) deg f <= 1  <=>  XOR of f over every 2-flat is 0   (all 2^16 f)
      (b) min-weight codewords of RM(2,4) are exactly the 140 2-flat
          indicators, and they span RM(2,4) (rank 11).
  P4. Grinberg-Stanley Thm 7.1 (arXiv:2307.05569):
      H = 1 + 2*c_odd  (mod 4)  over ALL tournaments at n = 4, 5, 6
      (not just the tiling cube).  c_odd = # directed cycles of odd length>1.
  P5. Tiling-cube digit ladder, fresh conventions, n = 4, 5, 6, 7:
      digit_0 == 1; GF(2) ANF degrees of digit_k; real multilinear (Walsh)
      degree of integer H; pointwise digit_1 == c_odd mod 2 on the cube at
      n = 4, 5, 6.

Conventions (MINE, on purpose):
  - vertices 0..n-1; base Hamiltonian path arcs (i+1) -> i, i = 0..n-2.
  - tiles = pairs (x,y), x - y >= 2, enumerated lexicographically by (y,x)
    ASCENDING (explorer uses a different order).
  - tile bit b = 1 means arc x -> y (high beats low); b = 0 means y -> x.
    (rm_digit script / explorer may use the opposite; ANF and real degrees
    are invariant under complementing variables, so results must agree.)
"""
import itertools
import math
import numpy as np

ok_all = True


def check(label, cond):
    global ok_all
    print(("PASS" if cond else "FAIL") + "  " + label)
    if not cond:
        ok_all = False


# ---------------------------------------------------------------- P1
print("=" * 72)
print("P1: pairing identity  prod z + prod(1+z) = sum_{S proper} prod_S z mod 2")
for L in (3, 5):
    # symbolic GF(2) polynomials: poly = set of frozensets (monomials)
    def pmul(p, q):
        out = {}
        for a in p:
            for b in q:
                mset = a | b
                out[mset] = out.get(mset, 0) ^ 1
        return {mn for mn, c in out.items() if c}

    edges = tuple(range(L))
    P_fwd = {frozenset(edges)}                 # prod z_e   (C-> indicator)
    P_bwd = {frozenset()}
    for e in edges:
        P_bwd = pmul(P_bwd, {frozenset(), frozenset([e])})  # prod (1+z_e)
    lhs = P_fwd ^ P_bwd                        # GF(2) sum = symmetric diff
    proper = {frozenset(s) for k in range(L)
              for s in itertools.combinations(edges, k)}
    check(f"L={L}: symbolic LHS == sum over proper subsets", lhs == proper)
    check(f"L={L}: degree == L-1 = {L-1}",
          max(len(mn) for mn in lhs) == L - 1)
    # brute force evaluation on all 2^L points
    brute = True
    for bits in itertools.product((0, 1), repeat=L):
        fwd = all(bits)                # all arcs forward
        bwd = not any(bits)            # all arcs backward
        val = sum(all(bits[e] for e in mn) for mn in lhs) % 2
        if val != (fwd + bwd) % 2:
            brute = False
    check(f"L={L}: brute-force evaluation matches on all 2^{L} points", brute)

# ---------------------------------------------------------------- P2
print("=" * 72)
print("P2: directed odd cycles pair up under reversal, no fixed points")


def directed_cycles(subset):
    """All directed cycles (rotation classes) on a vertex subset."""
    vs = sorted(subset)
    return [(vs[0],) + p for p in itertools.permutations(vs[1:])]


def canon(cyc):
    """Canonical form of a cyclic sequence: rotate so min vertex first."""
    i = cyc.index(min(cyc))
    return cyc[i:] + cyc[:i]


def rev(cyc):
    return canon((cyc[0],) + tuple(reversed(cyc[1:])))


# L=3 explicit: two orientations of a triple are distinct directed cycles
tri = directed_cycles((0, 1, 2))
check("L=3: a vertex triple carries exactly (3-1)! = 2 directed cycles",
      len(tri) == 2)
check("L=3: the two orientations are distinct and mutually reversed",
      tri[0] != tri[1] and rev(tri[0]) == tri[1] and rev(tri[1]) == tri[0])

for n, L in ((6, 3), (7, 5)):
    allc = []
    for sub in itertools.combinations(range(n), L):
        allc.extend(directed_cycles(sub))
    allc = [canon(c) for c in allc]
    check(f"n={n},L={L}: count = C(n,L)*(L-1)! = "
          f"{len(list(itertools.combinations(range(n), L))) * math.factorial(L-1)}",
          len(allc) == len(set(allc)) and
          len(allc) == len(list(itertools.combinations(range(n), L)))
          * math.factorial(L - 1))
    fixed = [c for c in allc if rev(c) == c]
    check(f"n={n},L={L}: reversal has NO fixed points", len(fixed) == 0)
    paired = all(rev(rev(c)) == c and rev(c) in set(allc) for c in allc)
    check(f"n={n},L={L}: reversal is an involution onto the same set, "
          f"so cycles group into {len(allc)//2} disjoint pairs", paired)

# per-pair indicator identity, brute force, L=3 and L=5:
for L in (3, 5):
    good = True
    for bits in itertools.product((0, 1), repeat=L):
        # z_e = 1 iff arc e oriented along C->.  C-> present iff all 1;
        # C<- present iff all 0.  Pair sum mod 2 must equal the proper-subset
        # polynomial, which has degree <= L-1 (checked in P1); here just
        # confirm pair sum = [all] + [none] and never 2 (orientations exclusive)
        s = all(bits) + (not any(bits))
        if s == 2:
            good = False
    check(f"L={L}: C-> and C<- never simultaneously present "
          f"(distinct directed cycles, exclusive indicators)", good)

# ---------------------------------------------------------------- P3
print("=" * 72)
print("P3: RM flat-annihilation equivalence at m=4, D=1")
m4 = 4
pts = list(range(16))
pc16 = np.array([bin(x).count("1") for x in range(1 << 16)], dtype=np.uint8)

# all 2-dim linear subspaces of F_2^4
subspaces = set()
for d1 in range(1, 16):
    for d2 in range(d1 + 1, 16):
        if d2 == d1:
            continue
        sp = frozenset({0, d1, d2, d1 ^ d2})
        if len(sp) == 4:
            subspaces.add(sp)
check("number of 2-dim subspaces of F_2^4 == Gaussian[4,2]_2 == 35",
      len(subspaces) == 35)
flats = set()
for sp in subspaces:
    for a in range(16):
        flats.add(frozenset(x ^ a for x in sp))
check("number of 2-dim affine flats == 35 * 4 == 140", len(flats) == 140)
flat_masks = []
for fl in flats:
    msk = 0
    for x in fl:
        msk |= 1 << x
    flat_masks.append(msk)

# (a) equivalence over ALL 2^16 boolean functions on F_2^4
F = np.arange(1 << 16, dtype=np.uint32)
ann = np.ones(1 << 16, dtype=bool)          # annihilated by every 2-flat?
for msk in flat_masks:
    ann &= (pc16[F & np.uint32(msk)] & 1) == 0
# degree <= 1 functions: f(x) = c + sum a_i x_i  -> 32 of them
deg_le_1 = np.zeros(1 << 16, dtype=bool)
for c in (0, 1):
    for a in range(16):  # a = linear form selector over 4 vars
        fv = 0
        for x in range(16):
            val = c ^ (bin(a & x).count("1") & 1)
            fv |= val << x
        deg_le_1[fv] = True
check("for ALL 65536 functions:  (deg <= 1)  <=>  (XOR over every 2-flat = 0)",
      bool(np.all(ann == deg_le_1)))

# (b) min-weight codewords of RM(2,4) = 2-flat indicators; they span RM(2,4)
basis = []
for k in range(3):                       # monomials of degree <= 2
    for S in itertools.combinations(range(4), k):
        fv = 0
        for x in range(16):
            if all((x >> i) & 1 for i in S):
                fv |= 1 << x
        basis.append(fv)
check("dim RM(2,4) == 11", len(basis) == 11)
codewords = {0}
for b in basis:
    codewords |= {c ^ b for c in codewords}
check("|RM(2,4)| == 2^11", len(codewords) == 2048)
w4 = [c for c in codewords if bin(c).count("1") == 4]
nonzero_w = [bin(c).count("1") for c in codewords if c]
check("min weight of RM(2,4) == 4", min(nonzero_w) == 4)
flat_mask_set = set(flat_masks)
check("ALL weight-4 codewords are 2-flat indicators and vice versa "
      f"({len(w4)} == 140)", set(w4) == flat_mask_set and len(w4) == 140)
# rank of the 140 flat indicators over GF(2) (pivot on highest set bit)
pivots = {}
rank = 0
for r in w4:
    cur = r
    while cur:
        hb = cur.bit_length() - 1
        if hb in pivots:
            cur ^= pivots[hb]
        else:
            pivots[hb] = cur
            rank += 1
            break
check("2-flat indicators span RM(2,4)  (GF(2) rank == 11)", rank == 11)

# ---------------------------------------------------------------- machinery
def ham_counts(n, arc, T):
    """Vectorised #Hamiltonian-path count for T tournaments at once.
    arc[u][v] = uint8 array (T,) with 1 iff arc u->v present."""
    full = (1 << n) - 1
    layers = {}
    for v in range(n):
        layers[(1 << v, v)] = np.ones(T, dtype=np.int32)
    for size in range(2, n + 1):
        for S in itertools.combinations(range(n), size):
            Sm = 0
            for v in S:
                Sm |= 1 << v
            for v in S:
                acc = np.zeros(T, dtype=np.int32)
                prev = Sm ^ (1 << v)
                for u in S:
                    if u == v:
                        continue
                    acc += layers[(prev, u)] * arc[u][v]
                layers[(Sm, v)] = acc
    H = np.zeros(T, dtype=np.int64)
    for v in range(n):
        H += layers[(full, v)]
    return H


def anf_degree(bits, m):
    a = np.array(bits, dtype=np.uint8).copy()
    for i in range(m):
        step = 1 << i
        a = a.reshape(-1, 2 * step)
        a[:, step:] ^= a[:, :step]
        a = a.reshape(-1)
    nz = np.nonzero(a)[0]
    if len(nz) == 0:
        return -1
    return max(bin(int(x)).count("1") for x in nz)


def real_degree(vals, m):
    c = np.array(vals, dtype=np.int64).copy()
    for i in range(m):
        step = 1 << i
        c = c.reshape(-1, 2 * step)
        c[:, step:] -= c[:, :step]
        c = c.reshape(-1)
    nz = np.nonzero(c)[0]
    if len(nz) == 0:
        return -1
    return max(bin(int(x)).count("1") for x in nz)


def count_odd_directed_cycles(n, arc, T):
    tot = np.zeros(T, dtype=np.int64)
    for L in range(3, n + 1, 2):
        for sub in itertools.combinations(range(n), L):
            for cyc in directed_cycles(sub):
                ind = np.ones(T, dtype=np.uint8)
                for a, b in zip(cyc, cyc[1:] + (cyc[0],)):
                    ind = ind * arc[a][b]
                tot += ind
    return tot


# ---------------------------------------------------------------- P4
print("=" * 72)
print("P4: Grinberg-Stanley Thm 7.1 over ALL tournaments (independent check)")
for n in (4, 5, 6):
    pairs = [(u, v) for u in range(n) for v in range(u)]   # u > v
    mfull = len(pairs)
    T = 1 << mfull
    idx = np.arange(T)
    arc = [[None] * n for _ in range(n)]
    for i, (u, v) in enumerate(pairs):
        b = ((idx >> i) & 1).astype(np.uint8)
        arc[u][v] = b          # bit 1: u -> v
        arc[v][u] = (1 - b).astype(np.uint8)
    H = ham_counts(n, arc, T)
    codd = count_odd_directed_cycles(n, arc, T)
    check(f"n={n}: Redei (H odd) for all {T} tournaments",
          bool(np.all(H & 1 == 1)))
    check(f"n={n}: H mod 4 == (1 + 2*c_odd) mod 4 for ALL {T} tournaments",
          bool(np.all(H % 4 == (1 + 2 * codd) % 4)))

# ---------------------------------------------------------------- P5
print("=" * 72)
print("P5: tiling-cube digit ladder (fresh conventions)")
results = {}
for n in (4, 5, 6, 7):
    tiles = [(x, y) for y in range(n) for x in range(y + 2, n)]
    m = len(tiles)
    assert m == (n - 1) * (n - 2) // 2
    T = 1 << m
    idx = np.arange(T)
    arc = [[None] * n for _ in range(n)]
    zero = np.zeros(T, dtype=np.uint8)
    one = np.ones(T, dtype=np.uint8)
    for u in range(n):
        for v in range(n):
            if u == v:
                continue
            arc[u][v] = zero
    for i in range(n - 1):                       # base path arcs (i+1)->i
        arc[i + 1][i] = one
        arc[i][i + 1] = zero
    for i, (x, y) in enumerate(tiles):           # tile bit 1: x->y
        b = ((idx >> i) & 1).astype(np.uint8)
        arc[x][y] = b
        arc[y][x] = (1 - b).astype(np.uint8)
    H = ham_counts(n, arc, T)
    check(f"n={n}: H odd on entire tiling cube (digit_0 == 1, Redei)",
          bool(np.all(H & 1 == 1)))
    ndig = int(H.max()).bit_length()
    degs = []
    for k in range(ndig):
        degs.append(anf_degree((H >> k) & 1, m))
    rdeg = real_degree(H, m)
    D = 2 * ((n - 1) // 2)
    results[n] = (m, degs, rdeg)
    print(f"  n={n}  m={m}  D=2*floor((n-1)/2)={D}   max H={H.max()}")
    print(f"      ANF degrees of digits 0..{ndig-1}: {degs}")
    print(f"      real multilinear (Walsh) degree of H: {rdeg}")
    check(f"n={n}: digit_0 ANF degree == 0", degs[0] == 0)
    check(f"n={n}: digit_1 ANF degree == D == {D}  (claimed equality)",
          degs[1] == D)
    check(f"n={n}: digit_1 ANF degree <= D (the PROVED bound)", degs[1] <= D)
    check(f"n={n}: real Walsh degree of H == D == {D}", rdeg == D)
    if ndig > 3:
        check(f"n={n}: digits k>=3 all have ANF degree >= m-1 = {m-1}",
              all(d >= m - 1 for d in degs[3:]))
    # pointwise digit_1 == c_odd mod 2 on the cube (n <= 6 cheap; n=7 too big
    # for the cycle enumeration to stay cheap -> do n=4,5,6)
    if n <= 6:
        codd = count_odd_directed_cycles(n, arc, T)
        check(f"n={n}: pointwise digit_1 == c_odd mod 2 on all {T} tilings",
              bool(np.all(((H >> 1) & 1) == (codd & 1))))
        check(f"n={n}: pointwise H mod 4 == 1 + 2*c_odd mod 4 on the cube",
              bool(np.all(H % 4 == (1 + 2 * codd) % 4)))
    # convention robustness: flip ALL tile bits (complement variables) and
    # re-derive digit_1 degree -- must be identical
    if n <= 6:
        Hf = H[::-1].copy()      # t -> complement(t) reverses index order
        check(f"n={n}: digit_1 degree invariant under bit-convention flip",
              anf_degree((Hf >> 1) & 1, m) == degs[1])

print("=" * 72)
expected = {4: (3, [0, 2, 3], 2),
            5: (6, [0, 4, 6, 5], 4),
            6: (10, [0, 4, 7], 4),     # digits 0..2 claimed; k>=3 near-full
            7: (15, [0, 6, 11], 6)}    # digits 0..2 claimed
for n in (4, 5, 6, 7):
    m, degs, rdeg = results[n]
    em, edegs, erdeg = expected[n]
    check(f"n={n}: claimed digit degrees {edegs} match computed "
          f"{degs[:len(edegs)]}, claimed Walsh {erdeg} == {rdeg}",
          m == em and degs[:len(edegs)] == edegs and rdeg == erdeg)

print()
print("OVERALL:", "ALL CHECKS PASSED" if ok_all else "SOME CHECKS FAILED")
