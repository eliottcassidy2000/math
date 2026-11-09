"""AUDIT D1 (THM-5055): independent value-level closure, written from scratch (Python standard library only; big-int
bitsets), started from the AUDITOR's exhaustive base (aw_enum ... rep: every flag tuple of an even graph on <= NB
vertices with a representative graph), with explicit witness reconstruction.

Classes per order n (sets of odd values d):
  E  even graph on n vertices          W  ... no universal vertex       P  ... has a pendant vertex
  S1 ... a pendant whose neighbour is not universal                    S2 ... two such pendants (Seed2)
  S3 ... a vertex of degree 2 (Seed3)  G  ... S2 and S3 ("good")
Moves (PROOF.md Lemmas 1-4, Corollary 5, seeds A/B; each re-proved by the auditor):
  padding (all classes; E -> W; P -> S1), cone over W, E_2-join (E -> E, W -> W), cone over (G u K_1) -> P,
  seed A: E_{n-4} + 2^(n-2) -> S2, seed B: E_{n-3} + 3*2^(n-2) - 2 -> S3,
  R: S1_{n-2} + 3*2^(n-3) -> S1, S2_{n-2} + 3*2^(n-3) -> G;  PH(3,2): S3_{n-3} + 7*2^(n-4) -> G.
Outputs: E(n) vs F(n); Gamma(c) for c in {1,3,5,7,11,13}; least n0 with n0, n0+1 in Gamma(-j) for odd j < 2^JM;
and a witness file (aw_verify format) for a sample of claims.
usage: python3 aw_closure.py base_rep_file NB NM JM witness_out NSAMPLE"""
import sys, random
basefile, NB, NM, JM, wout, NS = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], int(sys.argv[6])
CL = ['E', 'W', 'P', 'S1', 'S2', 'S3', 'G']
B = {c: [0] * (NM + 1) for c in CL}
rep = {}
for line in open(basefile):
    left, _, right = line.partition('|')
    p = left.split()
    if len(p) < 7: continue
    n, d, U, pend, pnu, iso, r3 = map(int, p[:7])
    if n > NB or n > NM: continue
    edges = [tuple(map(int, e.split('-'))) for e in right.split()]
    cls = ['E'] + (['W'] if not U else []) + (['P'] if pend >= 1 else []) + (['S1'] if pnu >= 1 else []) \
        + (['S2'] if pnu >= 2 else []) + (['S3'] if r3 else []) + (['G'] if pnu >= 2 and r3 else [])
    for c in cls:
        B[c][n] |= 1 << d
        rep.setdefault((c, n, d), edges)
def lim(n): return (1 << (1 << n)) - 1          # values < 2^n
def sh(x, off, n): return (x << off) & lim(n)
for n in range(2, NM + 1):
    for c in CL: B[c][n] |= B[c][n - 1]
    B['W'][n] |= B['E'][n - 1]
    B['S1'][n] |= B['P'][n - 1]
    if n - 1 >= 2: B['E'][n] |= sh(B['W'][n - 1], 1 << (n - 1), n)
    if n - 2 >= 1:
        B['E'][n] |= sh(B['E'][n - 2], 3 * (1 << (n - 2)) - 2, n)
        B['W'][n] |= sh(B['W'][n - 2], 3 * (1 << (n - 2)) - 2, n)
    if n - 1 >= 2: B['P'][n] |= sh(B['E'][n - 2], 1 << (n - 1), n)
    if n - 4 >= 1: B['S2'][n] |= sh(B['E'][n - 4], 1 << (n - 2), n)
    if n - 3 >= 1: B['S3'][n] |= sh(B['E'][n - 3], 3 * (1 << (n - 2)) - 2, n)
    if n - 2 >= 3:
        B['S1'][n] |= sh(B['S1'][n - 2], 3 * (1 << (n - 3)), n)
        B['G'][n] |= sh(B['S2'][n - 2], 3 * (1 << (n - 3)), n)
    if n - 3 >= 3: B['G'][n] |= sh(B['S3'][n - 3], 7 * (1 << (n - 4)), n)
    B['S2'][n] |= B['G'][n]; B['S3'][n] |= B['G'][n]; B['S1'][n] |= B['S2'][n]; B['P'][n] |= B['S1'][n]
    if n >= 3: B['W'][n] |= B['S1'][n]
    B['E'][n] |= B['W'][n] | B['P'][n] | B['S3'][n]
def has(c, n, d): return 1 <= n <= NM and d >= 1 and (B[c][n] >> d) & 1 == 1
def M(n): return 2**n - n if n % 2 else 2**n - n - 1
def inF(n, k):
    if k < 1 or k % 2 == 0 or k > M(n) or k in (3, 7, 13): return False
    return not (n % 2 == 1 and n >= 5 and k == 2**n - n - 2)
tot_miss = tot_extra = 0
for n in range(1, NM + 1):
    Fm = int.from_bytes(b'\xaa' * ((M(n) >> 3) + 1), 'little') & ((1 << (M(n) + 1)) - 1)   # odd k <= M(n)
    for k in (3, 7, 13): Fm &= ~(1 << k)
    if n % 2 == 1 and n >= 5: Fm &= ~(1 << (2**n - n - 2))
    Em = B['E'][n]
    miss = bin(Fm & ~Em).count('1'); extra = bin(Em & ~Fm).count('1')
    tot_miss += miss; tot_extra += extra
    print(f"E({n}): missing {miss}, extra {extra}")
print(f"TOTAL n <= {NM}: missing {tot_miss}, extra {tot_extra}")
for c in (1, 3, 5, 7, 11, 13):
    print(f"Gamma(+{c}):", [n for n in range(3, NM + 1) if has('G', n, 2**(n - 1) + c)])
n0s = {}; bad = 0; worst = -99; undecided = 0
for j in range(1, 1 << JM, 2):
    n0 = None
    for n in range(3, NM):
        if 2**(n - 1) > j and has('G', n, 2**(n - 1) - j) and has('G', n + 1, 2**n - j): n0 = n; break
    if n0 is None:
        if j.bit_length() + 9 <= NM: bad += 1; print("seed FAIL j =", j)
        else: undecided += 1
    else:
        n0s[j] = n0; worst = max(worst, n0 - j.bit_length())
        if n0 > j.bit_length() + 8: bad += 1; print("seed TOO LATE j =", j, n0)
print(f"seed base, odd j < 2^{JM}: failures {bad}, undecided (L(j)+9 > NM) {undecided}, max n0 - L(j) = {worst}")
# ---------------- witness reconstruction ----------------
sys.setrecursionlimit(100000)
def pad(g, n): return g  # edge lists on vertex set range(n); padding adds an isolated vertex = same edges
def cone(g, m):          # g on m vertices -> m + 1 vertices, apex m
    return g + [(v, m) for v in range(m)]
def e2join(g, m):        # g on m vertices -> m + 2, new a = m, b = m + 1 (non-adjacent)
    return g + [(v, m) for v in range(m)] + [(v, m + 1) for v in range(m)]
def ph(g, m, u, l):      # hub h = m adjacent to all of range(m) except u; leaves m+1..m+l
    return g + [(v, m) for v in range(m) if v != u] + [(m, m + 1 + i) for i in range(l)]
def degs(g, m):
    dg = [0] * m
    for a, b in g: dg[a] += 1; dg[b] += 1
    return dg
def pick_pnu(g, m):
    dg = degs(g, m)
    for v in range(m):
        if dg[v] == 1:
            q = [b if a == v else a for a, b in g if v in (a, b)][0]
            if dg[q] != m - 1: return v
    raise RuntimeError("no pendant with non-universal neighbour")
def pick_r3(g, m):
    dg = degs(g, m)
    for v in range(m):
        if dg[v] == 2: return v
    raise RuntimeError("no degree-2 vertex")
memo = {}
def wit(c, n, d):
    key = (c, n, d)
    if key in memo: return memo[key]
    assert has(c, n, d), key
    g = None
    if (c, n, d) in rep: g = rep[(c, n, d)]
    elif has(c, n - 1, d): g = wit(c, n - 1, d)
    elif c == 'W' and has('E', n - 1, d): g = wit('E', n - 1, d)
    elif c == 'S1' and has('P', n - 1, d): g = wit('P', n - 1, d)
    else:
        tries = []
        if c == 'E': tries = [('W', n), ('P', n), ('S3', n), ('cone', 'W', n - 1, 1 << (n - 1)), ('e2', 'E', n - 2, 3 * (1 << (n - 2)) - 2)]
        if c == 'W': tries = [('S1', n), ('e2', 'W', n - 2, 3 * (1 << (n - 2)) - 2)]
        if c == 'P': tries = [('S1', n), ('coneK1', 'E', n - 2, 1 << (n - 1))]
        if c == 'S1': tries = [('S2', n), ('R', 'S1', n - 2, 3 * (1 << (n - 3)))]
        if c == 'S2': tries = [('G', n), ('A', 'E', n - 4, 1 << (n - 2))]
        if c == 'S3': tries = [('G', n), ('Bs', 'E', n - 3, 3 * (1 << (n - 2)) - 2)]
        if c == 'G': tries = [('R', 'S2', n - 2, 3 * (1 << (n - 3))), ('PH32', 'S3', n - 3, 7 * (1 << (n - 4)))]
        for t in tries:
            if len(t) == 2:
                if has(t[0], n, d) and t[0] != c: g = wit(t[0], n, d); break
                continue
            op, src, m, off = t
            if m < 1 or not has(src, m, d - off): continue
            h = wit(src, m, d - off)
            if op == 'cone': g = cone(h, m)
            elif op == 'e2': g = e2join(h, m)
            elif op == 'coneK1': g = cone(h, m + 1)                      # G u K_1 on m + 1 vertices
            elif op == 'R': g = ph(h, m, pick_pnu(h, m), 1)
            elif op == 'PH32': g = ph(h, m, pick_r3(h, m), 2)
            elif op == 'A': g = cone(h, m + 2)                           # K_1 v (X u E_2) on m + 3, then u K_1
            elif op == 'Bs': g = e2join(h, m + 1)                        # E_2 v (Z u K_1) on m + 3
            break
    if g is None: raise RuntimeError(f"no rule for {key}")
    assert all(0 <= a < n and 0 <= b < n for a, b in g), key
    memo[key] = g
    return g
random.seed(5055)
claims = []
for n in range(1, NM + 1):                                   # sample of E(n) values plus every top-band value
    vals = [k for k in range(1, 1 << n, 2) if inF(n, k)]
    Xm = 2**(n - 1) - (n - 1) - 2 if (n - 1 >= 5 and (n - 1) % 2 == 1) else -1
    top = [k for k in vals if M(n - 1) < k < 2**(n - 1) or k == Xm]                      # top band + X(n-1)
    top += [k for k in vals if 2**(n - 1) < k < 3 * 2**(n - 2) and (k - 2**(n - 1) > M(n - 2) or k - 2**(n - 1) in (3, 7, 13))]
    top += [k for k in vals if k - 3 * 2**(n - 2) in (1, 5, 11)]
    for k in set(top + random.sample(vals, min(len(vals), NS))): claims.append(('E', n, k))
for c in (1, 3, 5, 7, 11, 13):
    for n in range(8, NM + 1): claims.append(('G', n, 2**(n - 1) + c))
for j, n0 in n0s.items():
    if j < 128 or random.random() < 0.05: claims += [('G', n0, 2**(n0 - 1) - j), ('G', n0 + 1, 2**n0 - j)]
with open(wout, 'w') as f:
    for c, n, d in claims:
        g = wit(c, n, d)
        f.write(f"aw-{c}-{n}-{d} {'S' if c == 'S1' else c} {n} {d} | " + " ".join(f"{a}-{b}" for a, b in sorted(set(tuple(sorted(e)) for e in g))) + "\n")
print("witnesses written:", len(claims))
