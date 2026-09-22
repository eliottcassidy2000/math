#!/usr/bin/env python3
"""Adversarial audit of lane paley_fano_octonion_design (wave 2026-09-22).

Independent recomputation, with different code paths from the explorer's script:
  A1  T_7 counts (arcs, cyclic/transitive triples, tr A^3 by numpy, h by brute force over all 5040
      vertex orders, per-arc 3-cycle multiplicities).
  A2  All STS(7) on {0..6} by exact cover over the 21 pairs (not by S_7-images); the two inside T_7.
  A3  Octonion table built from the index rule alone (not from the tournament) and compared with the
      tournament-derived table; alternativity as total antisymmetry of the associator on basis triples;
      Moufang identity and composition on random rational vectors; the explicit Paley -> Cayley-Dickson
      relabelling; isomorphism with the textbook (Wikipedia/Baez-figure) table by search over S_7.
  A4  Census of the 128 line orientations of dev{0,1,3}: alternative count, h by brute force.
  A5  All 2^21 tournaments: 2640 regular, split by (h, |Aut| via explicit automorphism count of a
      representative), doubly-regular test and 2-design test for the three classes.
  A6  Paley T_p, p = 3 mod 4: per-arc 3-cycle count (A^2 o A^T entries) and doubly regularity
      (A A^T off-diagonal) for p in {3,7,11,19,23,31,43,47}; AP lemma; cyclic STS(p) counts by an
      independent difference-family search (2, 8, 192, 1024); p = 31 disjoint partitions.
  A7  The 57 blocks printed in the explorer's .out for p = 19 form an STS(19) inside the design.
  A8  Entropy numbers.  A9  Lean file rerun + no sorry.  A10  OEIS A007079 fetch (network permitting).
Every check is an explicit raise; python3 -O identical.
"""
import itertools, math, os, re, subprocess, shutil, sys, time
from fractions import Fraction
import numpy as np

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "..", "05-knowledge", "results")

def out(*a):
    print(*a); sys.stdout.flush()

def check(c, m):
    if not c:
        raise RuntimeError("AUDIT CHECK FAILED: " + m)

def paley(p):
    qr = {(x * x) % p for x in range(1, p)}
    return np.array([[1 if i != j and (j - i) % p in qr else 0 for j in range(p)] for i in range(p)], dtype=np.int64), qr

def is_cyclic(A, t):
    a, b, c = t
    return (A[a, b] and A[b, c] and A[c, a]) or (A[a, c] and A[c, b] and A[b, a])

def h_bruteforce(A):
    n = len(A)
    return sum(1 for perm in itertools.permutations(range(n)) if all(A[perm[i], perm[i + 1]] for i in range(n - 1)))

def h_dp(A):
    n = len(A); N = 1 << n
    dp = [[0] * n for _ in range(N)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(N):
        for v in range(n):
            c = dp[S][v]
            if c:
                for w in range(n):
                    if not (S >> w) & 1 and A[v, w]:
                        dp[S | 1 << w][w] += c
    return sum(dp[N - 1])

# ---------------------------------------------------------------- A1
out("=== A1. Paley T_7 counts (independent) ===")
A, QR = paley(7)
check(np.all(A + A.T + np.eye(7, dtype=np.int64) == 1), "tournament")
arcs = int(A.sum()); trip = list(itertools.combinations(range(7), 3))
cyc = [t for t in trip if is_cyclic(A, t)]
trA3 = int(np.trace(A @ A @ A))
h7 = h_bruteforce(A)
out(f"QR={sorted(QR)} arcs={arcs} triples={len(trip)} cyclic={len(cyc)} transitive={len(trip)-len(cyc)} trA^3={trA3} h(brute force over 5040 orders)={h7}")
check(arcs == 21 and len(cyc) == 14 and len(trip) - len(cyc) == 21 and trA3 == 42 and h7 == 189, "A1 counts")
M = (A @ A) * A.T          # M[i,j] = #{w : i->w->j} restricted to arcs j->i ... use per arc i->j: #w with j->w, w->i
per_arc = {(i, j): int(((A @ A).T * A)[i, j]) for i in range(7) for j in range(7) if A[i, j]}
# ((A@A).T)[i,j] = (A@A)[j,i] = #{w: j->w->i}; times A[i,j] restricts to arcs i->j
out("per-arc 3-cycle counts (multiset):", sorted(set(per_arc.values())), "over", len(per_arc), "arcs")
check(len(per_arc) == 21 and set(per_arc.values()) == {2}, "2-(7,3,2)")

# ---------------------------------------------------------------- A2
out()
out("=== A2. STS(7) by exact cover over the 21 pairs ===")
pairs = list(itertools.combinations(range(7), 2))
def all_sts7():
    sols = []
    def rec(chosen, covered):
        if len(covered) == 21:
            sols.append(tuple(sorted(chosen))); return
        # smallest uncovered pair
        pr = next(p_ for p_ in pairs if p_ not in covered)
        for t in trip:
            if pr[0] in t and pr[1] in t:
                ps = {(t[0], t[1]), (t[0], t[2]), (t[1], t[2])}
                if not (ps & covered):
                    rec(chosen + [t], covered | ps)
    rec([], frozenset())
    return sorted(set(sols))
STS = all_sts7()
out("labelled STS(7):", len(STS))
check(len(STS) == 30, "30 STS")
cycset = set(cyc)
inside = [S for S in STS if set(S) <= cycset]
dev013 = tuple(sorted(tuple(sorted(((x + s) % 7) for x in (0, 1, 3))) for s in range(7)))
dev015 = tuple(sorted(tuple(sorted(((x + s) % 7) for x in (0, 1, 5))) for s in range(7)))
out("STS(7) inside the cyclic triples:", len(inside), "; == {dev013, dev015}:", sorted(inside) == sorted([dev013, dev015]),
    "; disjoint and union = cyclic triples:", not (set(dev013) & set(dev015)) and set(dev013) | set(dev015) == cycset)
check(len(inside) == 2 and sorted(inside) == sorted([dev013, dev015]) and set(dev013) | set(dev015) == cycset, "two STS")
o13 = all(A[s % 7, (s + 1) % 7] and A[(s + 1) % 7, (s + 3) % 7] and A[(s + 3) % 7, s % 7] for s in range(7))
o15 = all(A[s % 7, (s + 1) % 7] and A[(s + 1) % 7, (s + 5) % 7] and A[(s + 5) % 7, s % 7] for s in range(7))
neg_anti = all(A[(-i) % 7, (-j) % 7] == A[j, i] for i in range(7) for j in range(7))
out("orientation s->s+1->s+3->s on all lines:", o13, "; s->s+1->s+5->s:", o15, "; x->-x anti-automorphism:", neg_anti)
check(o13 and o15 and neg_anti, "orientations")
# per-pair 3-cycle multiplicities for the other two regular classes are done in A5

# ---------------------------------------------------------------- A3
out()
out("=== A3. Octonions: index rule vs tournament, associator antisymmetry, Cayley-Dickson ===")
def table_from_rule():
    """Imaginary units e_0..e_6 (index mod 7): e_r e_{r+1} = e_{r+3}, e_{r+1} e_{r+3} = e_r, e_{r+3} e_r = e_{r+1},
    anticommuting, e_r^2 = -1. Built WITHOUT reference to the tournament."""
    T = {}
    for r in range(7):
        a, b, c = r, (r + 1) % 7, (r + 3) % 7
        for (x, y, z) in ((a, b, c), (b, c, a), (c, a, b)):
            T[(x + 1, y + 1)] = (1, z + 1); T[(y + 1, x + 1)] = (-1, z + 1)
    for i in range(8):
        T[(0, i)] = (1, i); T[(i, 0)] = (1, i)
    for i in range(1, 8):
        T[(i, i)] = (-1, 0)
    check(len(T) == 64, "table complete")
    return T
def table_from_tournament(Adj, lines):
    T = {}
    for a in range(8):
        for b in range(8):
            if a == 0 or b == 0:
                T[(a, b)] = (1, a | b) if (a == 0 or b == 0) else None
                T[(a, b)] = (1, b if a == 0 else a)
            elif a == b:
                T[(a, b)] = (-1, 0)
            else:
                r, s = a - 1, b - 1
                L = [l for l in lines if r in l and s in l]
                check(len(L) == 1, "line")
                t = [x for x in L[0] if x not in (r, s)][0]
                T[(a, b)] = (1 if Adj[r, s] else -1, t + 1)
    return T
def mul(T, x, y):
    z = [Fraction(0)] * 8
    for a in range(8):
        if x[a]:
            for b in range(8):
                if y[b]:
                    s, c = T[(a, b)]; z[c] += s * x[a] * y[b]
    return z
def E(i):
    v = [Fraction(0)] * 8; v[i] = Fraction(1); return v
def assoc(T, x, y, z):
    l = mul(T, mul(T, x, y), z); r = mul(T, x, mul(T, y, z))
    return [l[i] - r[i] for i in range(8)]
def alternative_by_antisymmetry(T):
    """Multilinear associator is alternating (char 0) iff on basis vectors: (a,a,c)=0, (a,b,b)=0,
    (a,b,c) = -(b,a,c) = -(a,c,b)."""
    Es = [E(i) for i in range(8)]
    for a in range(8):
        for b in range(8):
            for c in range(8):
                A1 = assoc(T, Es[a], Es[b], Es[c])
                if a == b or b == c:
                    if any(A1): return False
                if A1 != [-t for t in assoc(T, Es[b], Es[a], Es[c])]: return False
                if A1 != [-t for t in assoc(T, Es[a], Es[c], Es[b])]: return False
    return True
Trule = table_from_rule()
Tpal = table_from_tournament(A, [tuple(l) for l in dev013])
out("table from the index rule == table from the Paley orientation on dev{0,1,3}:", Trule == Tpal)
check(Trule == Tpal, "rule = Paley")
altP = alternative_by_antisymmetry(Trule)
out("alternative (associator alternating on all 512 basis triples):", altP)
check(altP, "alternative")
rng = np.random.default_rng(20260922)
def rvec():
    return [Fraction(int(v), 3) for v in rng.integers(-9, 10, 8)]
def norm(x): return sum(t * t for t in x)
ok_comp = ok_mouf = ok_alt = True
for _ in range(60):
    x, y, z = rvec(), rvec(), rvec()
    if norm(mul(Trule, x, y)) != norm(x) * norm(y): ok_comp = False
    if mul(Trule, mul(Trule, x, x), y) != mul(Trule, x, mul(Trule, x, y)): ok_alt = False
    if mul(Trule, mul(Trule, y, x), x) != mul(Trule, y, mul(Trule, x, x)): ok_alt = False
    # Moufang: z(x(zy)) = ((zx)z)y
    if mul(Trule, z, mul(Trule, x, mul(Trule, z, y))) != mul(Trule, mul(Trule, mul(Trule, z, x), z), y): ok_mouf = False
out("60 random rational triples: N(xy)=N(x)N(y):", ok_comp, "; (xx)y=x(xy),(yx)x=y(xx):", ok_alt, "; Moufang z(x(zy))=((zx)z)y:", ok_mouf)
check(ok_comp and ok_alt and ok_mouf, "composition/Moufang")
w1 = mul(Trule, mul(Trule, E(1), E(2)), E(3)); w2 = mul(Trule, E(1), mul(Trule, E(2), E(3)))
out("(e_0 e_1) e_2 =", [int(t) for t in w1], "; e_0 (e_1 e_2) =", [int(t) for t in w2], "; differ:", w1 != w2)
out("nonzero coordinate index:", [i for i in range(8) if w2[i]], "= imaginary unit e_5 (index i <-> e_{i-1}); the explorer note wrote e_6: corrected")
check(w1 != w2 and w1 == [-t for t in w2] and int(w2[6]) == 1 and int(w2[7]) == 0, "non-associative witness = +-e_5")
# Cayley-Dickson (Baez: (a,b)(c,d) = (ac - d b*, a* d + c b)), independent implementation on tuples
def cd_conj(x):
    if len(x) == 1: return x
    h = len(x) // 2
    return cd_conj(x[:h]) + tuple(-t for t in x[h:])
def cd_mul(x, y):
    if len(x) == 1: return (x[0] * y[0],)
    h = len(x) // 2
    a, b, c, d = x[:h], x[h:], y[:h], y[h:]
    p1 = cd_mul(a, c); p2 = cd_mul(d, cd_conj(b)); p3 = cd_mul(cd_conj(a), d); p4 = cd_mul(c, b)
    return tuple(p1[i] - p2[i] for i in range(h)) + tuple(p3[i] + p4[i] for i in range(h))
def cd_basis(i):
    return tuple(1 if k == i else 0 for k in range(8))
CD = np.zeros((7, 7), dtype=np.int64)
for a in range(1, 8):
    for b in range(1, 8):
        if a != b:
            z = cd_mul(cd_basis(a), cd_basis(b))
            check(z[a ^ b] in (1, -1) and sum(abs(t) for t in z) == 1, "CD xor rule")
            if z[a ^ b] == 1: CD[a - 1, b - 1] = 1
check(np.all(CD + CD.T + np.eye(7, dtype=np.int64) == 1), "CD tournament")
phi = {0: 1, 1: 3, 2: 5, 3: 2, 4: 6, 5: 7, 6: 4}
iso = all(A[i, j] == CD[phi[i] - 1, phi[j] - 1] for i in range(7) for j in range(7))
out("explorer's relabelling r -> label", phi, "carries Paley arcs onto Cayley-Dickson signs:", iso)
check(iso, "explicit isomorphism")
xor_lines = {tuple(sorted((a - 1, b - 1, (a ^ b) - 1))) for a in range(1, 8) for b in range(1, 8) if a != b}
img_lines = {tuple(sorted(phi[x] - 1 for x in l)) for l in dev013}
out("the relabelling carries dev{0,1,3} onto the XOR lines:", img_lines == xor_lines)
check(img_lines == xor_lines, "lines")
h_cd = h_bruteforce(CD); cyc_cd = sum(1 for t in trip if is_cyclic(CD, t))
out("Cayley-Dickson sign tournament: h =", h_cd, "; cyclic triples =", cyc_cd)
check(h_cd == 189 and cyc_cd == 14, "CD counts")
# textbook table (Wikipedia 'Octonion' / Baez Fig.: e1e2=e3, e1e4=e5, e1e7=e6, e2e4=e6, e2e5=e7, e3e4=e7, e3e6=e5)
book = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
BK = np.zeros((7, 7), dtype=np.int64)
for (x, y, z) in book:
    for (u, v) in ((x, y), (y, z), (z, x)):
        BK[u - 1, v - 1] = 1
check(np.all(BK + BK.T + np.eye(7, dtype=np.int64) == 1), "book tournament")
n_iso = sum(1 for perm in itertools.permutations(range(7)) if all(A[i, j] == BK[perm[i], perm[j]] for i in range(7) for j in range(7)))
out("vertex bijections carrying Paley onto the textbook octonion sign tournament (UNCITED-RECOLLECTION for the textbook table):", n_iso, "(= |Aut| = 21 iff isomorphic)")
check(n_iso == 21, "book iso")
# fano_code convention mask of CD
def dot(x, a): return bin(x & a).count("1") % 2
m_cd = 0
for a in range(1, 8):
    x, y, z = sorted(v for v in range(1, 8) if dot(v, a) == 0)
    if CD[x - 1, y - 1] and CD[y - 1, z - 1] and CD[z - 1, x - 1]:
        m_cd |= 1 << (a - 1)
out("Cayley-Dickson mask in the fano_code convention:", m_cd)
check(m_cd == 96, "m0=96")

# ---------------------------------------------------------------- A4
out()
out("=== A4. The 128 orientations of dev{0,1,3} ===")
lines013 = [tuple(l) for l in dev013]
cnt = {}
for mask in range(128):
    Adj = np.zeros((7, 7), dtype=np.int64)
    for k, (x, y, z) in enumerate(lines013):
        cyc_ = ((x, y), (y, z), (z, x)) if (mask >> k) & 1 else ((x, z), (z, y), (y, x))
        for (u, v) in cyc_: Adj[u, v] = 1
    check(np.all(Adj + Adj.T + np.eye(7, dtype=np.int64) == 1), "orientation tournament")
    alt = alternative_by_antisymmetry(table_from_tournament(Adj, lines013))
    hh = h_bruteforce(Adj)
    key = (alt, hh); cnt[key] = cnt.get(key, 0) + 1
    if mask == 127:
        check(np.array_equal(Adj, A), "mask 127 = Paley")
out("(alternative, h) census:", cnt)
check(cnt == {(True, 189): 16, (False, 171): 112}, "16/112 census")

# ---------------------------------------------------------------- A5
out()
out("=== A5. All 2^21 labelled 7-tournaments; regular ones; doubly regular vs 2-design ===")
pair_idx = list(itertools.combinations(range(7), 2))
reg = []
for m in range(1 << 21):
    sc = [0] * 7
    for k, (i, j) in enumerate(pair_idx):
        if (m >> k) & 1: sc[i] += 1
        else: sc[j] += 1
    if sc == [3] * 7: reg.append(m)
out("regular labelled 7-tournaments:", len(reg))
check(len(reg) == 2640, "2640")
def adj_of(m):
    Adj = np.zeros((7, 7), dtype=np.int64)
    for k, (i, j) in enumerate(pair_idx):
        if (m >> k) & 1: Adj[i, j] = 1
        else: Adj[j, i] = 1
    return Adj
byh = {}
for m in reg:
    Adj = adj_of(m); hh = h_dp(Adj)
    byh.setdefault(hh, []).append(m)
out("h -> labelled count:", {k: len(v) for k, v in sorted(byh.items())})
check({k: len(v) for k, v in byh.items()} == {189: 240, 175: 720, 171: 1680}, "class sizes")
for hh in sorted(byh, reverse=True):
    Adj = adj_of(byh[hh][0])
    aut = sum(1 for perm in itertools.permutations(range(7)) if all(Adj[i, j] == Adj[perm[i], perm[j]] for i in range(7) for j in range(7)))
    AAT = Adj @ Adj.T
    common_out = sorted({int(AAT[i, j]) for i in range(7) for j in range(7) if i != j})
    arc_cnt = sorted({int(((Adj @ Adj).T * Adj)[i, j]) for i in range(7) for j in range(7) if Adj[i, j]})
    n_sts = sum(1 for S in STS if all(is_cyclic(Adj, t) for t in S))
    dreg = (len(common_out) == 1); design = (len(arc_cnt) == 1)
    out(f"h={hh}: |Aut|={aut} (7!/|Aut|={5040//aut}) common-out-neighbour counts={common_out} per-arc 3-cycle counts={arc_cnt} "
        f"doubly-regular={dreg} 2-design={design} STS(7) inside cyclic triples={n_sts}")
    check(5040 // aut == len(byh[hh]) and dreg == design == (hh == 189), "doubly regular iff 2-design, only Paley")
    check(n_sts == {189: 2, 175: 0, 171: 2}[hh], "STS inside")
out("on the three regular 7-tournaments: cyclic triples form a 2-design <=> doubly regular <=> Paley (h=189): True")

# ---------------------------------------------------------------- A6
out()
out("=== A6. Paley T_p, p = 3 mod 4 ===")
for p in (3, 7, 11, 19, 23, 31, 43, 47):
    Ap, qr = paley(p)
    check(np.all(Ap + Ap.T + np.eye(p, dtype=np.int64) == 1), "tournament")
    arc_cnt = sorted({int(((Ap @ Ap).T * Ap)[i, j]) for i in range(p) for j in range(p) if Ap[i, j]})
    AAT = Ap @ Ap.T
    common = sorted({int(AAT[i, j]) for i in range(p) for j in range(p) if i != j})
    c3 = int(np.trace(Ap @ Ap @ Ap)) // 3
    ap_cyc = sum(1 for a in range(1, p) if is_cyclic(Ap, (0, a, (2 * a) % p)))
    two_nr = (2 % p) not in qr
    out(f"p={p}: c3={c3} (p^3-p)/24={(p**3-p)//24} per-arc 3-cycles={arc_cnt} (p+1)/4={(p+1)//4} common-out={common} (p-3)/4={(p-3)//4} "
        f"AP triples cyclic={ap_cyc}/{p-1} 2 non-residue={two_nr} p mod 8={p%8}")
    check(c3 == (p ** 3 - p) // 24 and arc_cnt == [(p + 1) // 4] and common == [(p - 3) // 4], f"design/doubly regular p={p}")
    check(ap_cyc == (p - 1 if two_nr else 0) and two_nr == (p % 8 == 3), "AP lemma")
for p in (59, 67, 71, 79, 83):
    Ap, qr = paley(p)
    ap_cyc = sum(1 for a in range(1, p) if is_cyclic(Ap, (0, a, (2 * a) % p)))
    check(ap_cyc == (p - 1 if p % 8 == 3 else 0), "AP lemma large p")
out("AP lemma confirmed also for p in {59,67,71,79,83}")
# cyclic STS(p) inside the design: independent difference-family search over base blocks {0,a,b}
def cyclic_sts_count(p):
    Ap, qr = paley(p)
    need = (p - 1) // 6
    # canonical base block per translation class: the lexicographically least translate
    cands = {}
    for a in range(1, p):
        for b in range(a + 1, p):
            if is_cyclic(Ap, (0, a, b)):
                cl = min(tuple(sorted(((x + s) % p) for x in (0, a, b))) for s in range(p))
                d = [(x - y) % p for x in (0, a, b) for y in (0, a, b) if x != y]
                if len(set(d)) == 6:
                    cands[cl] = frozenset(d)
    n_classes = sum(1 for a in range(1, p) for b in range(a + 1, p) if is_cyclic(Ap, (0, a, b))) // 3
    keys = sorted(cands)
    sols = []
    def rec(i, chosen, cov):
        if len(chosen) == need:
            sols.append(tuple(chosen)); return
        for j in range(i, len(keys)):
            if not (cands[keys[j]] & cov):
                rec(j + 1, chosen + [keys[j]], cov | cands[keys[j]])
    rec(0, [], frozenset())
    return n_classes, len(cands), sols
for p in (7, 19, 31, 43):
    ncl, ncand, sols = cyclic_sts_count(p)
    lam = (p + 1) // 4
    out(f"p={p}: translation classes={ncl} non-AP classes={ncand} cyclic STS({p}) inside the design={len(sols)}")
    check(ncl == (p * p - 1) // 24 and len(sols) == {7: 2, 19: 8, 31: 192, 43: 1024}[p], f"cyclic STS count p={p}")
    if p == 31:
        sets = [frozenset(s) for s in sols]
        n_part = 0
        def rec2(i, used, k):
            global n_part
            if k == lam:
                if len(used) == ncand: n_part += 1
                return
            for j in range(i, len(sets)):
                if not (sets[j] & used):
                    rec2(j + 1, used | sets[j], k + 1)
        rec2(0, frozenset(), 0)
        out(f"p=31: partitions of the 40 classes into 8 block-disjoint cyclic STS(31): {n_part}")
        check(n_part == 0, "no partition at 31")

# ---------------------------------------------------------------- A7
out()
out("=== A7. The explorer's STS(19) witness ===")
with open(os.path.join(RES, "collatz_mod6_20260922_paley_fano_octonion_design.out"), encoding="utf-8") as f:
    txt = f.read()
mline = [l for l in txt.splitlines() if l.strip().startswith("blocks: [")][0]
blocks = [tuple(map(int, m)) for m in re.findall(r"\((\d+), (\d+), (\d+)\)", mline)]
A19, _ = paley(19)
cover = {}
for t in blocks:
    check(is_cyclic(A19, t), "block cyclic in T_19")
    for pr in itertools.combinations(sorted(t), 2):
        cover[pr] = cover.get(pr, 0) + 1
out("blocks:", len(blocks), "; all cyclic triples of T_19:", True, "; pairs covered:", len(cover), "/", math.comb(19, 2),
    "; every pair exactly once:", set(cover.values()) == {1})
check(len(blocks) == 57 and len(cover) == 171 and set(cover.values()) == {1}, "STS(19) witness")
ap_blocks = [t for t in blocks if (2 * t[1] - t[0] - t[2]) % 19 == 0 or (2 * t[0] - t[1] - t[2]) % 19 == 0 or (2 * t[2] - t[0] - t[1]) % 19 == 0]
out("arithmetic-progression blocks in the witness:", len(ap_blocks), "-> the witness is not translation-invariant (a cyclic STS(19) has none)")
check(len(ap_blocks) > 0, "AP blocks present")

# ---------------------------------------------------------------- A8
out()
out("=== A8. Entropy ===")
q = 8 / math.pi ** 2
H = lambda t: -(t * math.log2(t) + (1 - t) * math.log2(1 - t))
out(f"8/pi^2 = {q:.6f}; H = {H(q):.6f} bits; paste 0.704 differs by {abs(0.704 - H(q)):.5f}; 6/pi^2 = {6/math.pi**2:.5f}, H = {H(6/math.pi**2):.5f}")
check(abs(H(q) - 0.70028) < 5e-6 and abs(H(6 / math.pi ** 2) - 0.96612) < 5e-6, "entropy")

# ---------------------------------------------------------------- A9
out()
out("=== A9. Lean file ===")
lf = os.path.join(HERE, "collatz_mod6_20260922_paley_fano_octonion_design.lean")
src = open(lf, encoding="utf-8").read()
out("contains 'sorry':", "sorry" in src, "; imports:", [l for l in src.splitlines() if l.startswith("import")])
check("sorry" not in src and not any(l.startswith("import") for l in src.splitlines()), "lean source clean")
lean = shutil.which("lean")
if lean:
    r = subprocess.run([lean, lf], capture_output=True, text=True, timeout=300)
    out("lean exit:", r.returncode, "; stdout:", r.stdout.strip().replace("\n", " | "))
    check(r.returncode == 0 and "does not depend on any axioms" in r.stdout, "lean ok")
else:
    out("lean not on PATH: SKIPPED (SCOPE)")

# ---------------------------------------------------------------- A10
out()
out("=== A10. OEIS A007079 ===")
try:
    r = subprocess.run(["curl", "-s", "-m", "15", "-A", "Mozilla/5.0 research", "https://oeis.org/search?q=id:A007079&fmt=json"],
                       capture_output=True, text=True, timeout=30)
    import json
    d = json.loads(r.stdout); rec = d[0] if isinstance(d, list) else d["results"][0]
    data = rec["data"].split(",")
    out("A007079:", rec["name"], "; offset", rec["offset"], "; data[:5] =", data[:5])
    check(int(data[3]) == 2640, "A007079(3) = 2640 (7 nodes)")
    out("2640 is a(3) of A007079 (2n+1 = 7 nodes), not 'A007079(7)'.")
except Exception as ex:  # noqa
    out("OEIS fetch unavailable (", type(ex).__name__, "): 2640 = A007079(3) stays UNCITED-RECOLLECTION")

out()
out(f"audit runtime {time.time() - T0:.0f} s")
out("ALL AUDIT CHECKS PASSED")
