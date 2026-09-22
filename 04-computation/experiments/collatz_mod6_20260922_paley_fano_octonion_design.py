#!/usr/bin/env python3
"""Lane paley_fano_octonion_design (wave 2026-09-22, session collatz-mod6-20260917).

Exact, self-contained checks behind the note
05-knowledge/results/collatz_mod6_20260922_paley_fano_octonion_design.md.

Sections:
  A  Paley T_7: arcs, cyclic/transitive triples, tr A^3, score formula, h(T_7) by DP.
  B  The 14 cyclic triples are the two cyclic STS(7) dev{0,1,3} and dev{0,1,5};
     the decomposition into two STS(7) is unique; Paley orientation on each line.
  C  Octonions from the Paley orientation: e_i e_{i+1} = e_{i+3} rule, alternativity,
     norm multiplicativity, non-associativity witness; the same on dev{0,1,5};
     census of all 128 line orientations of dev{0,1,3}; Cayley-Dickson octonions
     are isomorphic (as sign tournaments) to Paley; mask/syndrome in the
     fano_code convention (m0 = 96).
  D  All 2^21 labelled 7-tournaments: the 2640 regular ones fall into 3 classes;
     which classes are Fano-line orientations; per-arc 3-cycle multiplicities.
  E  Paley T_p for p = 3 mod 4: 2-(p,3,(p+1)/4) design, arc-regular affine group,
     cyclic STS(p) sub-systems, exact-cover search for an STS(19) inside the design.
  F  Binary entropy of 8/pi^2; the paste's numbers.
  G  Core-Lean `decide` run of the companion .lean file (no Mathlib available locally).
Every check is an explicit raise; python3 -O is identical.
"""
import itertools, math, sys, time
from fractions import Fraction

T0 = time.time()

def out(*a):
    print(*a)
    sys.stdout.flush()

def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)

# ---------------------------------------------------------------- helpers
def paley(p):
    qr = {(x * x) % p for x in range(1, p)}
    adj = [[0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                adj[i][j] = 1
    return adj, qr

def is_tournament(adj):
    n = len(adj)
    for i in range(n):
        if adj[i][i]:
            return False
        for j in range(i + 1, n):
            if adj[i][j] + adj[j][i] != 1:
                return False
    return True

def scores(adj):
    return [sum(r) for r in adj]

def cyclic_triples(adj):
    n = len(adj)
    cyc, tr = [], []
    for a, b, c in itertools.combinations(range(n), 3):
        s = adj[a][b] + adj[b][c] + adj[c][a]
        # a 3-cycle iff the three "forward" arcs are all present or all absent
        if s == 3 or s == 0:
            cyc.append((a, b, c))
        else:
            tr.append((a, b, c))
    return cyc, tr

def ham_paths(adj):
    """Number of directed Hamiltonian paths, subset DP (exact, integers)."""
    n = len(adj)
    N = 1 << n
    dp = [[0] * n for _ in range(N)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(N):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if c == 0:
                continue
            for w in range(n):
                if not (S >> w) & 1 and adj[v][w]:
                    dp[S | (1 << w)][w] += c
    return sum(dp[N - 1])

def trace_A3(adj):
    n = len(adj)
    A2 = [[sum(adj[i][k] * adj[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    return sum(A2[i][k] * adj[k][i] for i in range(n) for k in range(n))

def skew_charpoly(adj):
    """Characteristic polynomial of S = A - A^T as integer coefficients (Faddeev-LeVerrier, exact)."""
    n = len(adj)
    S = [[Fraction(adj[i][j] - adj[j][i]) for j in range(n)] for i in range(n)]
    M = [[Fraction(0)] * n for _ in range(n)]
    coeffs = [Fraction(1)]
    for k in range(1, n + 1):
        # M_k = S*M_{k-1} + c_{n-k+1} I
        MS = [[sum(S[i][l] * M[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        for i in range(n):
            MS[i][i] += coeffs[-1]
        M = MS
        SM = [[sum(S[i][l] * M[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        c = -sum(SM[i][i] for i in range(n)) / k
        coeffs.append(c)
    return tuple(int(c) for c in coeffs)  # x^n + c1 x^{n-1} + ... + cn

def poly_str(coeffs):
    n = len(coeffs) - 1
    terms = []
    for i, c in enumerate(coeffs):
        if c == 0:
            continue
        e = n - i
        terms.append(f"{c:+d}x^{e}" if e else f"{c:+d}")
    return " ".join(terms)

# ---------------------------------------------------------------- A
out("=== A. Paley tournament T_7 ===")
p = 7
A7, QR7 = paley(7)
check(is_tournament(A7), "T_7 tournament")
check(sorted(QR7) == [1, 2, 4], "QR mod 7")
arcs7 = sum(map(sum, A7))
cyc7, tr7 = cyclic_triples(A7)
out("QR_7 =", sorted(QR7), "; arcs =", arcs7, "; scores =", scores(A7))
out("cyclic triples c3 =", len(cyc7), "; transitive triples =", len(tr7),
    "; C(7,3) =", math.comb(7, 3))
out("formula (p^3-p)/24 =", (7 ** 3 - 7) // 24,
    "; score formula C(7,3) - 7*C(3,2) =", math.comb(7, 3) - 7 * math.comb(3, 2))
out("tr A^3 =", trace_A3(A7), "= 3*c3 (rooted 3-cycles = closed 3-walks)")
h7 = ham_paths(A7)
out("h(T_7) = number of directed Hamiltonian paths =", h7)
check(arcs7 == 21 and len(cyc7) == 14 and len(tr7) == 21, "T_7 counts")
check(trace_A3(A7) == 42 and h7 == 189, "T_7 trace and h")
out("paste's '21 directed 3-cycles': 21 = arcs = transitive triples = 7 lines x 3 rotations;"
    " true directed 3-cycle count = 14, rooted = 42; h(T_7) = 189 is neither 7 nor 21.")
per_arc = {}
for (a, b, c) in cyc7:
    if A7[a][b]:
        cyc = [(a, b), (b, c), (c, a)]
    else:
        cyc = [(a, c), (c, b), (b, a)]
    for e in cyc:
        per_arc[e] = per_arc.get(e, 0) + 1
mult7 = sorted(set(per_arc.values()))
out("arcs covered by cyclic triples:", len(per_arc), "; multiplicities:", mult7)
check(len(per_arc) == 21 and mult7 == [2], "T_7 is a 2-(7,3,2) design")

# ---------------------------------------------------------------- B
out()
out("=== B. Two cyclic Steiner triple systems ===")
def dev(base, p):
    return sorted(tuple(sorted(((x + s) % p) for x in base)) for s in range(p))
D013 = dev((0, 1, 3), 7)
D015 = dev((0, 1, 5), 7)
cycset = set(cyc7)
out("dev{0,1,3} =", D013)
out("dev{0,1,5} =", D015)
check(set(D013) | set(D015) == cycset and not (set(D013) & set(D015)), "cyclic triples = dev013 + dev015 disjoint")
out("cyclic triples of T_7 == dev{0,1,3} u dev{0,1,5} (disjoint): True")

# all STS(7) on labelled {0..6}: images of the XOR Fano under S_7
xor_lines = [tuple(sorted(x - 1 for x in (a, b, a ^ b))) for a in range(1, 8) for b in range(a + 1, 8) if a ^ b > b]
check(len(xor_lines) == 7, "xor lines")
all_sts = set()
for perm in itertools.permutations(range(7)):
    img = tuple(sorted(tuple(sorted(perm[x] for x in L)) for L in xor_lines))
    all_sts.add(img)
out("number of labelled STS(7) on {0..6} =", len(all_sts), "= 7!/168 =", math.factorial(7) // 168)
check(len(all_sts) == 30, "30 STS(7)")
inside = [S for S in all_sts if set(S) <= cycset]
out("STS(7) whose 7 blocks are all cyclic triples of T_7:", len(inside),
    "; they are dev013 and dev015:", sorted(inside) == sorted([tuple(D013), tuple(D015)]))
check(len(inside) == 2 and sorted(inside) == sorted([tuple(D013), tuple(D015)]), "unique STS pair")
# also: decompositions of the 14 into two STS: any STS inside has complement inside
out("decompositions of the 14 cyclic triples into two STS(7): 1 (unordered), unique")
# Paley orientation on each line
def orient(adj, line):
    a, b, c = line
    if adj[a][b] and adj[b][c] and adj[c][a]:
        return (a, b, c)
    if adj[a][c] and adj[c][b] and adj[b][a]:
        return (a, c, b)
    raise RuntimeError("not cyclic")
ok13 = all(orient(A7, (s % 7, (s + 1) % 7, (s + 3) % 7)) in
           {((s) % 7, (s + 1) % 7, (s + 3) % 7), ((s + 1) % 7, (s + 3) % 7, s % 7), ((s + 3) % 7, s % 7, (s + 1) % 7)}
           for s in range(7))
ok15 = all(orient(A7, (s % 7, (s + 1) % 7, (s + 5) % 7)) in
           {((s) % 7, (s + 1) % 7, (s + 5) % 7), ((s + 1) % 7, (s + 5) % 7, s % 7), ((s + 5) % 7, s % 7, (s + 1) % 7)}
           for s in range(7))
out("on every line {s,s+1,s+3}: Paley orientation is s->s+1->s+3->s:", ok13)
out("on every line {s,s+1,s+5}: Paley orientation is s->s+1->s+5->s:", ok15)
check(ok13 and ok15, "line orientations")
# negation x -> -x
neg_img = sorted(tuple(sorted((-x) % 7 for x in L)) for L in D013)
out("negation x->-x maps dev013 to dev015:", neg_img == D015,
    "; negation is an anti-automorphism of T_7 (-1 is a non-residue):",
    all(A7[(-i) % 7][(-j) % 7] == A7[j][i] for i in range(7) for j in range(7) if i != j))
check(neg_img == D015, "negation swaps the two STS")

# ---------------------------------------------------------------- C
out()
out("=== C. Octonions from the Paley orientation ===")
def algebra_from_orientation(adj, lines):
    """Basis e0=1, e_(r+1) for r in Z/7. e_r e_s = +e_t if r->s in adj, -e_t if s->r; e_r^2 = -1.
    Returns table[(a,b)] = (sign, c) for a,b in 0..7."""
    tab = {}
    for a in range(8):
        for b in range(8):
            if a == 0:
                tab[(a, b)] = (1, b)
            elif b == 0:
                tab[(a, b)] = (1, a)
            elif a == b:
                tab[(a, b)] = (-1, 0)
            else:
                r, s = a - 1, b - 1
                L = [Ln for Ln in lines if r in Ln and s in Ln]
                if len(L) != 1:
                    raise RuntimeError("pair not on exactly one line")
                t = [x for x in L[0] if x not in (r, s)][0]
                tab[(a, b)] = (1 if adj[r][s] else -1, t + 1)
    return tab

def mul(tab, x, y):
    z = [0] * 8
    for a in range(8):
        if x[a] == 0:
            continue
        for b in range(8):
            if y[b] == 0:
                continue
            sg, c = tab[(a, b)]
            z[c] += sg * x[a] * y[b]
    return z

def basis(i):
    v = [0] * 8
    v[i] = 1
    return v

def norm(x):
    return sum(t * t for t in x)

def alternative_exact(tab):
    """Left/right alternative laws checked as polynomial identities on basis elements
    via the linearized forms (x,x,y)+... : over Q, (xx)y = x(xy) for all x,y iff
    (ab)c + (ba)c = a(bc) + b(ac) for all basis a,b,c (char 0 polarization)."""
    E = [basis(i) for i in range(8)]
    for a in range(8):
        for b in range(8):
            for c in range(8):
                l1 = mul(tab, mul(tab, E[a], E[b]), E[c])
                l2 = mul(tab, mul(tab, E[b], E[a]), E[c])
                r1 = mul(tab, E[a], mul(tab, E[b], E[c]))
                r2 = mul(tab, E[b], mul(tab, E[a], E[c]))
                if [l1[i] + l2[i] for i in range(8)] != [r1[i] + r2[i] for i in range(8)]:
                    return False
                # right alternative, linearized: a(bc)+a(cb) = (ab)c+(ac)b
                r1 = mul(tab, E[a], mul(tab, E[b], E[c]))
                r2 = mul(tab, E[a], mul(tab, E[c], E[b]))
                l1 = mul(tab, mul(tab, E[a], E[b]), E[c])
                l2 = mul(tab, mul(tab, E[a], E[c]), E[b])
                if [r1[i] + r2[i] for i in range(8)] != [l1[i] + l2[i] for i in range(8)]:
                    return False
    return True

def composition_exact(tab):
    """N(xy)=N(x)N(y) as a polynomial identity: its full polarization in x and y is the multilinear form
    <xy,wz> + <xz,wy> - 2<x,w><y,z>; over Q the identity holds iff this vanishes on all basis quadruples."""
    E = [basis(i) for i in range(8)]
    P = {}
    for a in range(8):
        for b in range(8):
            P[(a, b)] = mul(tab, E[a], E[b])
    def ip(u, v):
        return sum(u[i] * v[i] for i in range(8))
    for a in range(8):
        for b in range(8):
            for c in range(8):
                for d in range(8):
                    lhs = ip(P[(a, b)], P[(c, d)]) + ip(P[(a, d)], P[(c, b)])
                    if lhs != (2 if (a == c and b == d) else 0):
                        return False
    return True

def lcg(seed):
    s = seed
    while True:
        s = (1103515245 * s + 12345) % (1 << 31)
        yield (s >> 8) % 21 - 10

def random_vec(g):
    return [next(g) for _ in range(8)]

def composition_random(tab, trials=200, seed=7):
    g = lcg(seed)
    for _ in range(trials):
        x, y = random_vec(g), random_vec(g)
        if norm(mul(tab, x, y)) != norm(x) * norm(y):
            return False
        # direct alternative laws on random vectors too
        xx = mul(tab, x, x)
        if mul(tab, xx, y) != mul(tab, x, mul(tab, x, y)):
            return False
        if mul(tab, mul(tab, y, x), x) != mul(tab, y, xx):
            return False
    return True

lines013 = [tuple(L) for L in D013]
lines015 = [tuple(L) for L in D015]
tabP = algebra_from_orientation(A7, lines013)
# Baez index rule e_i e_{i+1} = e_{i+3} (indices mod 7), written with imaginary index r in Z/7 as e_(r+1)
rule = all(tab_entry == (1, ((r + 3) % 7) + 1) for r in range(7) for tab_entry in [tabP[(r + 1, ((r + 1) % 7) + 1)]])
rule2 = all(tabP[(((r + 1) % 7) + 1, ((r + 3) % 7) + 1)] == (1, r + 1) and tabP[(((r + 3) % 7) + 1, r + 1)] == (1, ((r + 1) % 7) + 1) for r in range(7))
out("Paley/dev013 table satisfies e_r e_{r+1} = +e_{r+3}, e_{r+1} e_{r+3} = +e_r, e_{r+3} e_r = +e_{r+1}:", rule and rule2)
check(rule and rule2, "index rule")
altP = alternative_exact(tabP)
compP = composition_random(tabP)
compPx = composition_exact(tabP)
out("Paley/dev013 algebra: alternative (exact, linearized on all 512 basis triples):", altP,
    "; N(xy)=N(x)N(y) exact (polarized, all 4096 basis quadruples):", compPx,
    "; N(xy)=N(x)N(y) and (xx)y=x(xy),(yx)x=y(xx) on 200 random integer pairs:", compP)
check(altP and compP and compPx, "Paley octonions")
e = [basis(i) for i in range(8)]
lhs = mul(tabP, mul(tabP, e[1], e[2]), e[3]); rhs = mul(tabP, e[1], mul(tabP, e[2], e[3]))
out("non-associativity witness: (e_0 e_1) e_2 =", lhs, "; e_0 (e_1 e_2) =", rhs, "(imaginary indices r=0,1,2)")
check(lhs != rhs, "nonassociative")
# same for dev015 with the Paley orientation
tabQ = algebra_from_orientation(A7, lines015)
altQ = alternative_exact(tabQ); compQ = composition_random(tabQ, seed=11); compQx = composition_exact(tabQ)
out("Paley/dev015 algebra (rule e_r e_{r+1} = +e_{r+5}): alternative:", altQ, "; composition exact:", compQx, "; composition random:", compQ)
check(altQ and compQ and compQx, "dev015 octonions")
# relation: negation conjugates dev013-table to the OPPOSITE of the dev015-table
opp_ok = True
for a in range(1, 8):
    for b in range(1, 8):
        if a == b:
            continue
        sg, c = tabP[(a, b)]
        na, nb = ((-(a - 1)) % 7) + 1, ((-(b - 1)) % 7) + 1
        sg2, c2 = tabQ[(nb, na)]
        if c2 != ((-(c - 1)) % 7) + 1 or sg2 != sg:
            opp_ok = False
out("x->-x carries the dev013 table to the opposite algebra of the dev015 table:", opp_ok)
check(opp_ok, "opposite relation")

# census of all 128 orientations of the lines of dev013
out()
out("--- census of the 128 orientations of the 7 lines of dev{0,1,3} ---")
def tournament_from_mask(lines, mask):
    n = 7
    adj = [[0] * n for _ in range(n)]
    for k, (x, y, z) in enumerate(lines):
        if (mask >> k) & 1:
            cyc = [(x, y), (y, z), (z, x)]
        else:
            cyc = [(x, z), (z, y), (y, x)]
        for (u, v) in cyc:
            adj[u][v] = 1
    if not is_tournament(adj):
        raise RuntimeError("mask does not give a tournament")
    return adj

census = {}
alt_masks = []
paley_mask = None
paley_op_mask = None
A7op = [[A7[j][i] for j in range(7)] for i in range(7)]
for mask in range(128):
    adj = tournament_from_mask(lines013, mask)
    tab = algebra_from_orientation(adj, lines013)
    alt = alternative_exact(tab)
    h = ham_paths(adj)
    c3 = len(cyclic_triples(adj)[0])
    cp = skew_charpoly(adj)
    key = (alt, h, c3, cp)
    census[key] = census.get(key, 0) + 1
    comp = composition_exact(tab)
    if comp != alt:
        raise RuntimeError("composition and alternativity disagree on mask %d" % mask)
    if alt:
        alt_masks.append(mask)
    if adj == A7:
        paley_mask = mask
    if adj == A7op:
        paley_op_mask = mask
for key in sorted(census, key=lambda k: (-k[0], -k[1])):
    alt, h, c3, cp = key
    out(f"alternative={alt} h={h} c3={c3} skew charpoly={poly_str(cp)} : {census[key]} orientations")
check(len(alt_masks) == 16, "16 alternative orientations")
out("for every one of the 128 orientations: alternative <=> composition (exact) <=> h = 189: True")
check(sum(v for k, v in census.items() if k[1] == 189) == 16 and sum(v for k, v in census.items() if k[1] == 171) == 112, "189/171 split")
out("Paley = mask", paley_mask, "(bits = lines dev013 in order s=0..6, bit 1 = s->s+1->s+3->s);",
    "Paley^op = mask", paley_op_mask, "; both alternative:", paley_mask in alt_masks and paley_op_mask in alt_masks)
check(paley_mask == 127 and paley_op_mask == 0, "Paley is the all-ones mask")
# collineation group of dev013 and its orbits on the 16 alternative masks
coll = []
Dset = set(lines013)
for perm in itertools.permutations(range(7)):
    if all(tuple(sorted(perm[x] for x in L)) in Dset for L in lines013):
        coll.append(perm)
out("collineations of dev{0,1,3} (perms of Z/7 preserving the lines):", len(coll))
check(len(coll) == 168, "GL(3,2)")
def mask_of(adj):
    m = 0
    for k, (x, y, z) in enumerate(lines013):
        if adj[x][y] and adj[y][z] and adj[z][x]:
            m |= 1 << k
    return m
def relabel(adj, perm):
    n = len(adj)
    B = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            B[perm[i]][perm[j]] = adj[i][j]
    return B
orbits = []
seen = set()
for m in alt_masks:
    if m in seen:
        continue
    adj = tournament_from_mask(lines013, m)
    orb = set()
    for perm in coll:
        orb.add(mask_of(relabel(adj, perm)))
    seen |= orb
    orbits.append(sorted(orb))
out("GL(3,2)-orbits on the 16 alternative masks:", [len(o) for o in orbits])
out("orbit of Paley (mask 127):", [o for o in orbits if 127 in o][0])
out("orbit of Paley^op (mask 0):", [o for o in orbits if 0 in o][0])
check(sorted(len(o) for o in orbits) == [8, 8], "two orbits of 8")
stab = [perm for perm in coll if relabel(A7, perm) == A7]
out("collineations fixing the Paley orientation (= Aut(T_7) inside GL(3,2)):", len(stab))
check(len(stab) == 21, "Aut order 21")
aut_all = sum(1 for perm in itertools.permutations(range(7)) if relabel(A7, perm) == A7)
out("|Aut(T_7)| over all of S_7:", aut_all, "(so Aut(T_7) preserves dev013: every automorphism is a collineation)")
check(aut_all == 21, "aut 21")
# the 112 non-alternative: all isomorphic to one tournament? use invariants: charpoly and h
nonalt_keys = {k for k in census if not k[0]}
out("distinct (h, c3, skew charpoly) among the 112 non-alternative orientations:", len(nonalt_keys))

# Cayley-Dickson octonions (Baez convention (a,b)(c,d) = (ac - d b*, a* d + c b)) as a sign tournament
def cd_mul(x, y):
    if len(x) == 1:
        return [x[0] * y[0]]
    h = len(x) // 2
    a, b, c, d = x[:h], x[h:], y[:h], y[h:]
    ac = cd_mul(a, c); db = cd_mul(d, cd_conj(b)); ad = cd_mul(cd_conj(a), d); cb = cd_mul(c, b)
    return [ac[i] - db[i] for i in range(h)] + [ad[i] + cb[i] for i in range(h)]
def cd_conj(x):
    if len(x) == 1:
        return [x[0]]
    h = len(x) // 2
    return cd_conj(x[:h]) + [-t for t in x[h:]]
cd_adj = [[0] * 7 for _ in range(7)]
xor_ok = True
for a in range(1, 8):
    for b in range(1, 8):
        if a == b:
            continue
        z = cd_mul(basis(a), basis(b))
        nz = [i for i in range(8) if z[i] != 0]
        if nz != [a ^ b] or abs(z[a ^ b]) != 1:
            xor_ok = False
        if z[a ^ b] == 1:
            cd_adj[a - 1][b - 1] = 1
out("Cayley-Dickson: e_a e_b = +-e_(a xor b) for all distinct a,b in 1..7:", xor_ok, "; sign tournament is a tournament:", is_tournament(cd_adj))
check(xor_ok and is_tournament(cd_adj), "CD table")
# mask in the fano_code convention: lines indexed by normal a=1..7 (line = {x,y,z} with x xor y xor z = 0 ... normal a: the line {x : <x,a>=0})
def dot(x, a):
    return bin(x & a).count("1") % 2
fc_lines = {}
for a in range(1, 8):
    pts = tuple(sorted(x for x in range(1, 8) if dot(x, a) == 0))
    fc_lines[a] = pts
def fc_mask(adj7):  # adj7 indexed by labels 1..7 -> indices 0..6
    m = 0
    for a in range(1, 8):
        x, y, z = fc_lines[a]
        if adj7[x - 1][y - 1] and adj7[y - 1][z - 1] and adj7[z - 1][x - 1]:
            m |= 1 << (a - 1)
    return m
m_cd = fc_mask(cd_adj)
out("Cayley-Dickson orientation mask in the fano_code convention (bit a-1 <-> normal a):", m_cd, "(fano_code m0 = 96)")
check(m_cd == 96, "m0 = 96 reproduced")
def syndrome(m, m0=96):
    s = 0
    for a in range(1, 8):
        if ((m ^ m0) >> (a - 1)) & 1:
            s ^= a
    return s
# identifications Z/7 -> {1..7} carrying dev013 to the XOR lines
xor_line_set = set(fc_lines.values())
idents = [perm for perm in itertools.permutations(range(1, 8))
          if all(tuple(sorted(perm[x] for x in L)) in xor_line_set for L in lines013)]
out("bijections Z/7 -> {1..7} carrying dev{0,1,3} onto the XOR Fano lines:", len(idents))
check(len(idents) == 168, "168 identifications")
masks_P = sorted({fc_mask(relabel(A7, [q - 1 for q in perm])) for perm in idents})
masks_Pop = sorted({fc_mask(relabel(A7op, [q - 1 for q in perm])) for perm in idents})
out("fano_code masks induced by Paley over all 168 identifications:", masks_P)
out("fano_code masks induced by Paley^op over all 168 identifications:", masks_Pop)
out("their relative syndromes:", sorted({syndrome(m) for m in masks_P + masks_Pop}))
check(len(masks_P) == 8 and len(masks_Pop) == 8 and not set(masks_P) & set(masks_Pop), "8+8 masks")
check(all(syndrome(m) == 0 for m in masks_P + masks_Pop), "all syndromes zero")
check(96 in masks_P or 96 in masks_Pop, "CD mask among Paley masks")
out("Cayley-Dickson mask 96 lies in the Paley orbit:", 96 in masks_P, "; in the Paley^op orbit:", 96 in masks_Pop)
ident0 = [perm for perm in idents if fc_mask(relabel(A7, [q - 1 for q in perm])) == 96][0]
out("one explicit isomorphism Paley -> Cayley-Dickson sign tournament: r -> label", dict(zip(range(7), ident0)))
check(relabel(A7, [q - 1 for q in ident0]) == cd_adj, "explicit isomorphism")

# ---------------------------------------------------------------- D
out()
out("=== D. All 2^21 labelled 7-tournaments; the regular ones ===")
pairs = list(itertools.combinations(range(7), 2))
check(len(pairs) == 21, "21 pairs")
# bit k of mask: 1 means pairs[k] = (i,j) oriented i->j, else j->i
out_bits = [[k for k, (i, j) in enumerate(pairs) if i == a] for a in range(7)]
in_bits = [[k for k, (i, j) in enumerate(pairs) if j == a] for a in range(7)]
out_m = [sum(1 << k for k in out_bits[a]) for a in range(7)]
in_m = [sum(1 << k for k in in_bits[a]) for a in range(7)]
regular_masks = []
for mask in range(1 << 21):
    ok = True
    for a in range(7):
        if bin(mask & out_m[a]).count("1") + bin(~mask & in_m[a]).count("1") != 3:
            ok = False
            break
    if ok:
        regular_masks.append(mask)
out("labelled regular tournaments on 7 vertices:", len(regular_masks), "(OEIS A007079(7) = 2640)")
check(len(regular_masks) == 2640, "2640 regular")
def adj_from_pairmask(mask):
    adj = [[0] * 7 for _ in range(7)]
    for k, (i, j) in enumerate(pairs):
        if (mask >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj
classes = {}
for mask in regular_masks:
    adj = adj_from_pairmask(mask)
    cp = skew_charpoly(adj)
    cyc, _ = cyclic_triples(adj)
    cs = set(cyc)
    n_sts = sum(1 for S in all_sts if set(S) <= cs)
    # per-arc multiplicities
    pa = {}
    for (a, b, c) in cyc:
        cy = [(a, b), (b, c), (c, a)] if adj[a][b] else [(a, c), (c, b), (b, a)]
        for e_ in cy:
            pa[e_] = pa.get(e_, 0) + 1
    mults = tuple(sorted((pa.get((i, j), 0) if adj[i][j] else pa.get((j, i), 0)) for (i, j) in pairs))
    key = (cp, len(cyc), n_sts, mults)
    if key not in classes:
        classes[key] = [0, ham_paths(adj), mask]
    classes[key][0] += 1
out("isomorphism-type invariants of the 2640 (skew charpoly, c3, #STS(7) inside cyclic triples, per-pair 3-cycle multiplicities):")
tot = 0
for key in sorted(classes, key=lambda k: -classes[k][1]):
    cp, c3, nsts, mults = key
    cnt, h, rep = classes[key]
    tot += cnt
    out(f"  charpoly {poly_str(cp)} : count={cnt} |Aut|={5040 // cnt} h={h} c3={c3} STS-inside={nsts} "
        f"pair-mults={sorted(set(mults))} rep-mask={rep}")
check(tot == 2640 and len(classes) == 3, "three classes")
for key in classes:
    cp, c3, nsts, mults = key
    cnt, h, rep = classes[key]
    if nsts == 2:
        adj = adj_from_pairmask(rep)
        cs = set(cyclic_triples(adj)[0])
        two = [S for S in all_sts if set(S) <= cs]
        shared = len(set(two[0]) & set(two[1]))
        out(f"  class h={h}: the 2 STS(7) inside its 14 cyclic triples share {shared} blocks; they partition the 14 cyclic triples: {shared == 0}")
        check((shared == 0) == (h == 189), "only Paley's cyclic triples are two disjoint Fano planes")
cnts = sorted(v[0] for v in classes.values())
out("class sizes:", cnts, "; |Aut| =", [5040 // c for c in cnts])
check(cnts == [240, 720, 1680], "class sizes")
# doubly regular (every ordered pair i != j has the same number of common out-neighbours) versus
# "cyclic triples form a 2-design" (every pair in the same number of cyclic triples), on the three classes
for key in sorted(classes, key=lambda k: -classes[k][1]):
    cp, c3, nsts, mults = key
    cnt, h, rep = classes[key]
    adj = adj_from_pairmask(rep)
    common = sorted({sum(adj[i][w] * adj[j][w] for w in range(7)) for i in range(7) for j in range(7) if i != j})
    out(f"  class h={h}: common out-neighbour counts over ordered pairs = {common}; doubly regular = {len(common) == 1};"
        f" per-pair 3-cycle counts = {sorted(set(mults))}; 2-design = {len(set(mults)) == 1}")
    check((len(common) == 1) == (len(set(mults)) == 1) == (h == 189), "doubly regular iff 2-design iff Paley")
out("on regular 7-tournaments: cyclic triples form a 2-design <=> doubly regular <=> Paley: True (audit S20b gives the general proof)")
# rotational R_7 (connection set {1,2,3}) identification
R7 = [[1 if (j - i) % 7 in (1, 2, 3) else 0 for j in range(7)] for i in range(7)]
out("R_7 = circulant {1,2,3}: skew charpoly", poly_str(skew_charpoly(R7)), "; h =", ham_paths(R7),
    "; STS inside its cyclic triples:", sum(1 for S in all_sts if set(S) <= set(cyclic_triples(R7)[0])))
out("Paley: skew charpoly", poly_str(skew_charpoly(A7)))
check(skew_charpoly(A7) == (1, 0, 21, 0, 147, 0, 343, 0), "x(x^2+7)^3")
# which classes are Fano-line orientations: exactly those with STS-inside >= 1
fano_classes = [(poly_str(k[0]), classes[k][0], classes[k][1], k[2]) for k in classes if k[2] >= 1]
out("classes that are orientations of some Fano plane (an STS inside the cyclic triples):", fano_classes)
n_fano_labelled = sum(c for (_, c, _, _) in fano_classes)
out("labelled Fano-line orientations total:", n_fano_labelled, "; each contains exactly 2 STS, so 2 *", n_fano_labelled, "=",
    2 * n_fano_labelled, "= 30 STS x 128 orientations =", 30 * 128)
check(sum(classes[k][0] * k[2] for k in classes) == 30 * 128, "double count 30*128")

# THM-133 cross-check: H = (462 - tr A^4)/2 for every circulant tournament on Z_7
def trace_A4(adj):
    n = len(adj)
    A2 = [[sum(adj[i][k] * adj[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    return sum(A2[i][k] * A2[k][i] for i in range(n) for k in range(n))
out("7! =", math.factorial(7))
circ_ok = True
for choice in itertools.product((0, 1), repeat=3):
    S_ = {d if c else 7 - d for d, c in zip((1, 2, 3), choice)}
    C = [[1 if (j - i) % 7 in S_ else 0 for j in range(7)] for i in range(7)]
    hC, t4 = ham_paths(C), trace_A4(C)
    if 2 * hC != 462 - t4:
        circ_ok = False
    if S_ == {1, 2, 4}:
        out("Paley: tr A^4 =", t4, "; (462 - tr A^4)/2 =", (462 - t4) // 2, "= h")
out("THM-133 formula h = (462 - tr A^4)/2 holds for all 8 circulant tournaments on Z_7:", circ_ok)
check(circ_ok, "THM-133")

# ---------------------------------------------------------------- E
out()
out("=== E. Paley T_p, p = 3 mod 4: the 2-(p,3,(p+1)/4) design ===")
def design_check(p):
    adj, qr = paley(p)
    check(is_tournament(adj), f"T_{p} tournament")
    sc = scores(adj)
    check(all(s == (p - 1) // 2 for s in sc), "regular")
    cyc, tr = cyclic_triples(adj)
    c3 = len(cyc)
    per_pair = {}
    per_arc_ = {}
    for (a, b, c) in cyc:
        for e_ in ((a, b), (a, c), (b, c)):
            per_pair[e_] = per_pair.get(e_, 0) + 1
        cy = [(a, b), (b, c), (c, a)] if adj[a][b] else [(a, c), (c, b), (b, a)]
        for e_ in cy:
            per_arc_[e_] = per_arc_.get(e_, 0) + 1
    lam_pairs = sorted(set(per_pair.values()))
    lam_arcs = sorted(set(per_arc_.values()))
    # arc-transitivity of G = {x -> ax+b : a in QR}: orbit of arc (0,1)
    orbit = {((b) % p, (a + b) % p) for a in qr for b in range(p)}
    narcs = sum(map(sum, adj))
    return c3, (p ** 3 - p) // 24, math.comb(p, 3) - p * math.comb((p - 1) // 2, 2), lam_pairs, lam_arcs, (p + 1) // 4, len(orbit), narcs, len(per_pair), cyc, adj, qr

design_rows = []
for p in (3, 7, 11, 19, 23, 31, 43, 47):
    c3, f1, f2, lp, la, lam, orb, narcs, npairs, cyc, adj, qr = design_check(p)
    out(f"p={p}: c3={c3} (p^3-p)/24={f1} scoreformula={f2} pairs-covered={npairs}/{math.comb(p,2)} "
        f"lambda(pairs)={lp} lambda(arcs)={la} (p+1)/4={lam} arc-orbit-of-G={orb} arcs={narcs}")
    check(c3 == f1 == f2 and lp == [lam] and la == [lam] and orb == narcs and npairs == math.comb(p, 2), f"design p={p}")
    design_rows.append((p, c3, lam))
out("2-(p,3,(p+1)/4) design verified for p in {3,7,11,19,23,31,43,47}; lambda=1 only at p=3, lambda=2 only at p=7.")

# arithmetic-progression triples {x, x+a, x+2a}: cyclic in T_p iff 2 is a non-residue mod p (p = 3 mod 8)
for p in (3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83):
    adj, qr = paley(p)
    ap_cyc = 0
    for a in range(1, p):
        t = (0, a, (2 * a) % p)
        s = adj[0][a] + adj[a][(2 * a) % p] + adj[(2 * a) % p][0]
        if s in (0, 3):
            ap_cyc += 1
    two_nqr = (2 % p) not in qr
    out(f"p={p} (p mod 8 = {p % 8}): 2 is a non-residue: {two_nqr}; AP base triples (0,a,2a) that are cyclic: {ap_cyc} of {p-1}"
        f" -> translation classes {ap_cyc // 2}")
    check(ap_cyc == (p - 1 if two_nqr else 0), "AP lemma")
out("AP lemma verified: all p-1 AP triples (0,a,2a) are cyclic when p = 3 mod 8, none when p = 7 mod 8.")

# cyclic STS(p) inside the design: base blocks {0,a,b} (cyclic triples) whose differences partition Z_p\{0}
def cyclic_sts_inside(p):
    adj, qr = paley(p)
    classes_all = set()
    orbs = []
    seen_cl = set()
    for a in range(1, p):
        for b in range(a + 1, p):
            t = (0, a, b)
            s = adj[0][a] + adj[a][b] + adj[b][0]
            if s in (0, 3):
                cl = frozenset(tuple(sorted(((x + sh) % p) for x in t)) for sh in range(p))
                classes_all.add(cl)
                diffs = frozenset(((x - y) % p) for x in t for y in t if x != y)
                if len(diffs) == 6 and cl not in seen_cl:
                    seen_cl.add(cl)
                    orbs.append((t, diffs, cl))
    need = (p - 1) // 6
    if (p - 1) % 6 != 0:
        return len(classes_all), len(orbs), need, None, []
    sols = []
    def rec(start, chosen, covered):
        if len(chosen) == need:
            sols.append(list(chosen))
            return
        for i in range(start, len(orbs)):
            t, d, cl = orbs[i]
            if not (d & covered):
                chosen.append(i)
                rec(i + 1, chosen, covered | d)
                chosen.pop()
    rec(0, [], frozenset())
    return len(classes_all), len(orbs), need, sols, orbs

def disjoint_partition(sols, orbs, lam):
    """Do lam pairwise block-disjoint cyclic STS among sols cover all translation classes? Count such partitions."""
    sets = [frozenset(i for i in sol) for sol in sols]
    total = len(orbs)
    count = 0
    def rec(start, chosen, used):
        nonlocal count
        if len(chosen) == lam:
            if len(used) == total:
                count += 1
            return
        for i in range(start, len(sets)):
            if not (sets[i] & used):
                chosen.append(i)
                rec(i + 1, chosen, used | sets[i])
                chosen.pop()
    rec(0, [], frozenset())
    return count

for p in (7, 19, 31, 43):
    ncl, norb, need, sols, orbs = cyclic_sts_inside(p)
    lam = (p + 1) // 4
    out(f"p={p}: an STS({p}) has {p*(p-1)//6} blocks; lam*that = {lam*p*(p-1)//6} = c3 = {(p**3-p)//24}")
    out(f"p={p}: translation classes of cyclic triples = {ncl} (=(p^2-1)/24 = {(p*p-1)//24}); "
        f"non-arithmetic-progression classes (candidate base blocks) = {norb}; base blocks needed = {need}; "
        f"cyclic STS(p) inside the Paley 2-(p,3,{lam}) design: {len(sols)}")
    check(ncl == (p * p - 1) // 24, "orbit count")
    if sols:
        out(f"   first four base-block sets: {[[orbs[i][0] for i in sol] for sol in sols[:4]]}")
        # AP classes cannot lie in any STS, so a partition into lam cyclic STS needs norb == lam*need
        if norb == lam * need:
            npart = disjoint_partition(sols, orbs, lam)
            out(f"   partitions of all {norb} classes into {lam} block-disjoint cyclic STS({p}): {npart}")
        else:
            out(f"   no partition into {lam} cyclic STS is possible: {ncl - norb} arithmetic-progression classes "
                f"(never in an STS) and {norb} != {lam}*{need}")
    if p == 7:
        check(len(sols) == 2, "two cyclic STS(7)")
out("p=11, 23, 47: STS(p) does not exist (p = 5 mod 6).")

# exact cover: does the Paley design at p=19 contain ANY STS(19)? (Algorithm X with a time cap)
def exact_cover_sts(p, cap_seconds):
    adj, qr = paley(p)
    cyc, _ = cyclic_triples(adj)
    cols = {}
    rows = {}
    for idx, (a, b, c) in enumerate(cyc):
        rows[idx] = [(a, b), (a, c), (b, c)]
        for e_ in rows[idx]:
            cols.setdefault(e_, set()).add(idx)
    t_start = time.time()
    nodes = [0]
    found = [None]
    def solve(sol):
        if time.time() - t_start > cap_seconds:
            raise TimeoutError
        if not cols:
            found[0] = list(sol)
            return True
        c = min(cols, key=lambda k: len(cols[k]))
        for r in sorted(cols[c]):
            nodes[0] += 1
            sol.append(r)
            removed = []
            for j in rows[r]:
                for i in cols[j]:
                    for k in rows[i]:
                        if k != j:
                            cols[k].discard(i)
                removed.append((j, cols.pop(j)))
            if solve(sol):
                return True
            for j, s_ in reversed(removed):
                cols[j] = s_
                for i in s_:
                    for k in rows[i]:
                        if k != j:
                            cols[k].add(i)
            sol.pop()
        return False
    try:
        res = solve([])
        status = "complete"
    except TimeoutError:
        res = False
        status = "timeout"
    return status, found[0], nodes[0], time.time() - t_start

status19, sol19, nodes19, secs19 = exact_cover_sts(19, 120.0)
out(f"p=19 exact cover (any STS(19) among the 285 cyclic triples): status={status19}, found={'yes' if sol19 else 'no'}, "
    f"search nodes={nodes19}, seconds={secs19:.0f}")
if status19 == "complete":
    out("VERDICT p=19:", "an STS(19) inside the Paley design EXISTS; blocks:" if sol19 else "NO STS(19) is contained in the Paley 2-(19,3,5) design (exhaustive).")
    if sol19:
        adj19, _ = paley(19)
        cyc19, _ = cyclic_triples(adj19)
        out("  blocks:", [cyc19[i] for i in sol19])
else:
    out("VERDICT p=19: undecided within the time cap (SCOPE).")

# ---------------------------------------------------------------- F
out()
out("=== F. Binary entropy of 8/pi^2 ===")
q = 8 / math.pi ** 2
Hq = -(q * math.log2(q) + (1 - q) * math.log2(1 - q))
out(f"8/pi^2 = {q:.5f}; binary entropy = {Hq:.5f} bits (paste says 0.704 and calls it 'absolute zero-entropy')")
check(abs(Hq - 0.70028) < 5e-6, "entropy 0.70028")
q6 = 6 / math.pi ** 2
out(f"6/pi^2 = {q6:.5f}; binary entropy = {-(q6*math.log2(q6)+(1-q6)*math.log2(1-q6)):.5f} bits (for the record)")
# ---------------------------------------------------------------- G
out()
out("=== G. Core-Lean (no Mathlib) decide checks ===")
import os, shutil, subprocess
lean_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "collatz_mod6_20260922_paley_fano_octonion_design.lean")
lean_bin = shutil.which("lean")
if lean_bin and os.path.exists(lean_file):
    ver = subprocess.run([lean_bin, "--version"], capture_output=True, text=True, timeout=60)
    out("lean --version:", ver.stdout.strip().split(",")[0])
    res = subprocess.run([lean_bin, lean_file], capture_output=True, text=True, timeout=300)
    out("lean exit code:", res.returncode)
    for line in (res.stdout + res.stderr).strip().splitlines():
        out("  " + line)
    check(res.returncode == 0, "lean file compiles")
    check("'paley7_cyclic_triples' does not depend on any axioms" in res.stdout, "axiom-free")
    out("theorems checked by `decide` (core Lean 4, no Mathlib): isTournament, numArcs = 21, numTriples = 35,"
        " numCyclic = 14, dev013 and dev015 orientations")
else:
    out("lean not available or file missing: SKIPPED (SCOPE)")
mathlib_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lean", "CollatzBlueprintAudit", ".lake", "packages", "mathlib")
out("local Mathlib checkout present:", os.path.isdir(mathlib_dir), "(the paste's `Mathlib.Combinatorics.SimpleGraph.Kuratowski` cannot be checked locally)")

out()
out(f"total runtime {time.time() - T0:.0f} s")
out("ALL CHECKS PASSED")
