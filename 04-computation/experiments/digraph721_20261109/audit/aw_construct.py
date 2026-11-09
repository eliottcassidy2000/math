"""AUDIT D1 (THM-5055): build explicit graphs EXACTLY as PROOF.md prescribes and emit them for brute-force checking
by aw_verify (d by inclusion-exclusion, evenness from a stabiliser-chain generating set, flags).
  realize(n, k):   section 5 (ranges 1-4) for n >= 13; base n <= 12 from the session's witness file (already
                   verified independently by aw_verify);
  good_pos(c, n):  Lemma 7 + Lemma 6 chains (R: +2, PH(3,2): +3) from the n in {8, 9} witnesses;
  good_neg(j, n):  Lemma 8 - base L(j) <= 12 from the seed witnesses + chains; L(j) >= 13 by cases A1/A2/A3/B1/B3
                   (seed A + R, padding of a Lemma 7 graph, double padding, seed B + padding + PH(3,2), padding).
Every recursive call asserts the hypotheses the proof states (ranges, bounds on i, j, x, y, z; orders).
usage: python3 aw_construct.py witnesses.txt out.txt NMAX NRAND"""
import sys, random
WFILE, OUT, NMAX, NRAND = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
random.seed(1729)
def M(n): return 2**n - n if n % 2 else 2**n - n - 1
def Xv(n): return 2**n - n - 2 if (n % 2 == 1 and n >= 5) else None
def inF(n, k):
    if k < 1 or k % 2 == 0 or k > M(n) or k in (3, 7, 13): return False
    return k != Xv(n)
def L(k): return k.bit_length()
Ew = {}; Gw = {}
for line in open(WFILE):
    left, _, right = line.partition('|')
    tag, st, n, k = left.split()[:4]; n = int(n); k = int(k)
    edges = [tuple(map(int, e.split('-'))) for e in right.split()]
    if st == 'E': Ew[k] = (n, edges)
    if st == 'G': Gw[(n, k)] = edges
def cone(g, m): return g + [(v, m) for v in range(m)]
def e2join(g, m): return g + [(v, m) for v in range(m)] + [(v, m + 1) for v in range(m)]
def ph(g, m, u, l): return g + [(v, m) for v in range(m) if v != u] + [(m, m + 1 + i) for i in range(l)]
def degs(g, m):
    dg = [0] * m
    for a, b in g: dg[a] += 1; dg[b] += 1
    return dg
def nbr(g, v): return [b if a == v else a for a, b in g if v in (a, b)]
def R(g, m):          # PH(G, p, 1), p a pendant with non-universal neighbour
    dg = degs(g, m)
    for v in range(m):
        if dg[v] == 1 and dg[nbr(g, v)[0]] != m - 1: return ph(g, m, v, 1)
    raise AssertionError("R: no pendant with non-universal neighbour")
def PH32(g, m):
    dg = degs(g, m)
    for v in range(m):
        if dg[v] == 2: return ph(g, m, v, 2)
    raise AssertionError("PH32: no degree-2 vertex")
stats = {}
def note(s): stats[s] = stats.get(s, 0) + 1
def realize(n, k):
    assert inF(n, k), (n, k)
    if n <= 12:
        m, g = Ew[k]; assert m <= n; note('base<=12'); return g
    if k <= M(n - 1) and k != Xv(n - 1):
        note('r1'); return realize(n - 1, k)
    if (M(n - 1) < k < 2**(n - 1)) or k == Xv(n - 1):
        j = 2**(n - 1) - k
        assert j <= n - 1 or (k == Xv(n - 1) and j == n + 1)
        assert n >= L(j) + 8
        note('r2'); return good_neg(j, n)
    if 2**(n - 1) < k <= 3 * 2**(n - 2) - 3:
        x = k - 2**(n - 1); assert 1 <= x <= 2**(n - 2) - 3
        if inF(n - 2, x): note('r3a'); return cone(realize(n - 2, x), n - 1)        # padded to n - 1, then cone
        if x in (3, 7, 13): note('r3c'); return good_pos(x, n)
        assert x > M(n - 2) or x == Xv(n - 2)
        jp = 2**(n - 2) - x; assert jp <= n and n - 1 >= L(jp) + 8
        note('r3b'); return cone(good_neg(jp, n - 1), n - 1)
    assert 3 * 2**(n - 2) - 1 <= k <= M(n)
    y = k - 3 * 2**(n - 2) + 2; assert 1 <= y <= M(n - 2)
    if inF(n - 2, y): note('r4a'); return e2join(realize(n - 2, y), n - 2)
    assert y in (3, 7, 13), (n, k, y)
    c = y - 2; note('r4b'); return cone(good_pos(c, n - 1), n - 1)
def chain(g, n1, n, src):    # Lemma 6: from a good graph on n1 vertices to n, with n - n1 in {0, 2, 3, 4, ...}
    d = n - n1; assert d == 0 or d >= 2, (n1, n)
    while d > 0:
        if d == 3 or (d >= 5 and random.random() < 0.3): g = PH32(g, n1); n1 += 3; d -= 3; note(src + ':PH32')
        else: g = R(g, n1); n1 += 2; d -= 2; note(src + ':R')
    return g
def good_pos(c, n):
    assert c in (1, 3, 5, 7, 11, 13) and n >= 8
    if n in (8, 9): note('L7base'); return Gw[(n, 2**(n - 1) + c)]
    if n - 3 >= 8 and random.random() < 0.4: note('L7:PH32'); return PH32(good_pos(c, n - 3), n - 3)
    note('L7:R'); return R(good_pos(c, n - 2), n - 2)
def good_neg(j, n, top=True):
    J = L(j); assert j % 2 == 1
    if top: assert n >= J + 8, (j, n)
    if J <= 12:
        n0 = min(m for (m, k) in Gw if k == 2**(m - 1) - j and (m + 1, 2**m - j) in Gw)
        assert n0 <= J + 8 and n >= n0
        if n in (n0, n0 + 1): note('L8base'); return Gw[(n, 2**(n - 1) - j)]
        if n - 3 >= n0 and random.random() < 0.4: note('L8b:PH32'); return PH32(good_neg(j, n - 3, False), n - 3)
        note('L8b:R'); return R(good_neg(j, n - 2, False), n - 2)
    assert 2**(J - 1) < j < 2**J
    if j > 3 * 2**(J - 2):
        x = 2**J - j; assert 1 <= x < 2**(J - 2)
        if inF(J - 2, x):                                  # A1: seed A (s = J), R -> good on J + 4
            X = realize(J - 2, x)
            S = cone(X, J)                                 # K_1 v (X u E_2) on J + 1 vertices; plus K_1 -> J + 2
            g, n1 = R(S, J + 2), J + 4; note('A1')
        elif x in (3, 7, 13):                              # A2: Lemma 7 graph on J + 1, padded -> J + 2
            g, n1 = good_pos(x, J + 1), J + 2; note('A2')
        else:                                              # A3
            assert x > M(J - 2) or x == Xv(J - 2)
            i = 2**(J - 2) - x; assert 1 <= i <= J and i % 2 == 1 and J - 1 >= L(i) + 8 and j == 3 * 2**(J - 2) + i
            g, n1 = good_neg(i, J - 1), J + 1; note('A3')
    else:
        y = 2**J - j; z = y + 2; assert 2**(J - 2) + 2 < z < 2**(J - 1) + 2
        if inF(J - 1, z):                                  # B1: seed B (s = J), pad, PH(3,2) -> J + 6
            X = realize(J - 1, z)
            S = e2join(X, J)                               # E_2 v (X u K_1) on J + 2; padded -> J + 3
            g, n1 = PH32(S, J + 3), J + 6; note('B1')
        else:                                              # B3
            assert z > 13 and (z > M(J - 1) or z == Xv(J - 1))
            i = j - 2**(J - 1); assert i % 2 == 1 and (i <= J or (z == Xv(J - 1) and i == J + 3)) and J >= L(i) + 8
            g, n1 = good_neg(i, J), J + 1; note('B3')
    return chain(g, n1, n, 'L8')
claims = []
def emit(st, n, k, g): claims.append((st, n, k, g))
# (1) section 5 for 13 <= n <= NMAX: every delicate value (ranges 2, 3b, 3c, 4b) plus NRAND random values per order
for n in range(13, NMAX + 1):
    ks = set()
    for k in range(M(n - 1) + 2, 2**(n - 1), 2): ks.add(k)
    if Xv(n - 1): ks.add(Xv(n - 1))
    for x in list(range(M(n - 2) + 2, 2**(n - 2), 2)) + ([Xv(n - 2)] if Xv(n - 2) else []) + [3, 7, 13]:
        ks.add(2**(n - 1) + x)
    for c in (1, 5, 11): ks.add(3 * 2**(n - 2) + c)
    ks.add(M(n))
    vals = [k for k in range(1, M(n) + 1, 2) if inF(n, k)]
    ks |= set(random.sample(vals, NRAND))
    for k in sorted(ks):
        if inF(n, k): emit('E', n, k, realize(n, k))
# (2) Lemma 8 inductive cases at J = 13 (j in (2^12, 2^13)): all A2/A3/B3 j's and samples of A1/B1, at n = 21
J = 13; picks = []
for j in range(2**(J - 1) + 1, 2**J, 2):
    if j > 3 * 2**(J - 2):
        x = 2**J - j; case = 'A1' if inF(J - 2, x) else ('A2' if x in (3, 7, 13) else 'A3')
    else:
        z = 2**J - j + 2; case = 'B1' if inF(J - 1, z) else 'B3'
    picks.append((case, j))
sel = [j for c, j in picks if c in ('A2', 'A3', 'B3')]
sel += random.sample([j for c, j in picks if c == 'A1'], 6) + random.sample([j for c, j in picks if c == 'B1'], 6)
for j in sel:
    for n in (J + 8, J + 9):
        if n <= NMAX + 9: emit('G', n, 2**(n - 1) - j, good_neg(j, n))
with open(OUT, 'w') as f:
    for st, n, k, g in claims:
        es = sorted(set(tuple(sorted(e)) for e in g))
        assert all(0 <= a < n and 0 <= b < n and a != b for a, b in es), (st, n, k)
        f.write(f"cons-{st}-{n}-{k} {st} {n} {k} | " + " ".join(f"{a}-{b}" for a, b in es) + "\n")
print("graphs:", len(claims), " max order:", max(n for _, n, _, _ in claims))
print("construction steps used:", dict(sorted(stats.items())))
