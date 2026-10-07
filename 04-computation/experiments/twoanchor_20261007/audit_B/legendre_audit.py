#!/usr/bin/env python3
"""Audit B (independent) of THM-4602 (3) and its 'Reading' remarks: Legendre chirotopes of P^1(F_p), p = 3 mod 4.
Written from scratch (no code shared with tournament_chirotope.py):
  * quadratic residues from the set of squares (not Euler's criterion);
  * Singer cycle from a companion matrix in GL2(F_p) with irreducible characteristic polynomial (not F_p[s]/(s^2-r));
  * a generic tournament isomorphism test (joint colour refinement + individualisation), applied to the
    sink-normalised switching representatives (two tournaments are switching-isomorphic iff for a fixed v there is
    a w with T1^(v) - v ~ T2^(w) - w, where T^(v) is the unique switching of T in which v is a sink).
"""
import itertools, math, random, sys
from collections import Counter

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)

def primes_3mod4(lo, hi):
    return [q for q in range(lo, hi + 1) if q % 4 == 3 and all(q % d for d in range(2, int(q**0.5) + 1))]

class Fp:
    def __init__(self, p):
        self.p = p
        self.sq = {(x*x) % p for x in range(1, p)}
    def chi(self, a):
        a %= self.p
        return 0 if a == 0 else (1 if a in self.sq else -1)
    def inv(self, a): return pow(a % self.p, self.p - 2, self.p)

def det(u, v, p): return (u[0]*v[1] - u[1]*v[0]) % p

def pattern(F, V):
    n = len(V); B = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                c = F.chi(det(V[i], V[j], F.p)); assert c != 0; B[i][j] = c
    return B

def paley_plus_sink(F):
    p = F.p; n = p + 1
    B = [[0]*n for _ in range(n)]
    for x in range(p):
        for y in range(p):
            if x != y: B[x][y] = F.chi(y - x)          # x -> y iff y - x is a residue
        B[x][p], B[p][x] = 1, -1                       # infinity (index p) is a sink
    return B

def sink_normalise(B, v):
    n = len(B)
    sw = [(-1 if (j != v and B[v][j] == 1) else 1) for j in range(n)]   # switch every j with v -> j
    C = [[B[i][j]*sw[i]*sw[j] if i != j else 0 for j in range(n)] for i in range(n)]
    assert all(C[j][v] == 1 for j in range(n) if j != v)
    keep = [j for j in range(n) if j != v]
    return [[C[i][j] for j in keep] for i in keep], keep, sw

def find_iso(A, B):
    """Isomorphism phi: V(A) -> V(B) with A[i][j] = B[phi i][phi j], or None. Joint refinement + backtracking."""
    n = len(A)
    if n != len(B): return None
    outA = [[u for u in range(n) if A[v][u] == 1] for v in range(n)]
    inA = [[u for u in range(n) if A[v][u] == -1] for v in range(n)]
    outB = [[u for u in range(n) if B[v][u] == 1] for v in range(n)]
    inB = [[u for u in range(n) if B[v][u] == -1] for v in range(n)]
    def refine(cA, cB):
        k = len(set(cA) | set(cB))
        while True:
            sA = [(cA[v], tuple(sorted(cA[u] for u in outA[v])), tuple(sorted(cA[u] for u in inA[v]))) for v in range(n)]
            sB = [(cB[v], tuple(sorted(cB[u] for u in outB[v])), tuple(sorted(cB[u] for u in inB[v]))) for v in range(n)]
            uniq = sorted(set(sA) | set(sB)); idx = {s: i for i, s in enumerate(uniq)}
            cA, cB = [idx[s] for s in sA], [idx[s] for s in sB]
            if len(uniq) == k: return cA, cB
            k = len(uniq)
    nodes = [0]
    def search(cA, cB):
        nodes[0] += 1
        cA, cB = refine(cA, cB)
        ca, cb = Counter(cA), Counter(cB)
        if ca != cb: return None
        non = [c for c in ca if ca[c] > 1]
        if not non:
            pos = {c: w for w, c in enumerate(cB)}
            phi = [pos[cA[v]] for v in range(n)]
            ok = all(A[i][j] == B[phi[i]][phi[j]] for i in range(n) for j in range(n) if i != j)
            return phi if ok else None
        c = min(non, key=lambda c: (ca[c], c))
        v = cA.index(c); new = max(ca) + 1
        for w in [j for j in range(n) if cB[j] == c]:
            a2, b2 = cA[:], cB[:]; a2[v] = new; b2[w] = new
            r = search(a2, b2)
            if r is not None: return r
        return None
    phi = search([0]*n, [0]*n)
    return phi, nodes[0]

def switching_iso(T1, T2):
    """Certificate that T1, T2 are switching-isomorphic: returns (phi, eps) with eps_i eps_j T1[i][j] = T2[phi i][phi j]."""
    n = len(T1)
    A, keepA, swA = sink_normalise(T1, 0)
    order = [n - 1] + list(range(n - 1))              # try the sink of T2 first, then the rest
    for w in order:
        Bw, keepB, swB = sink_normalise(T2, w)
        phi0, nodes = find_iso(A, Bw)
        if phi0 is not None:
            phi = [None]*n; phi[0] = w
            for a_idx, a in enumerate(keepA): phi[a] = keepB[phi0[a_idx]]
            # explicit verification of the full switching isomorphism
            eps = [swA[i]*swB[phi[i]] for i in range(n)]
            assert sorted(phi) == list(range(n))
            assert all(eps[i]*eps[j]*T1[i][j] == T2[phi[i]][phi[j]] for i in range(n) for j in range(n) if i != j)
            return phi, nodes, w
    return None, None, None

def normalise_strip(points, p):
    """consecutive determinant 1: v_0 = points[0], v_(k+1) = mu_k * points[k+1] with det(v_k, v_(k+1)) = 1."""
    V = [points[0]]
    for P in points[1:]:
        d = det(V[-1], P, p); assert d != 0
        mu = pow(d, p - 2, p)
        V.append(((mu*P[0]) % p, (mu*P[1]) % p))
    return V

def closure_constant(V, p):
    """v_(n) := normalised lift of the first point after v_(n-1); returns c with v_n = c v_0."""
    d = det(V[-1], V[0], p); mu = pow(d, p - 2, p)
    return mu % p          # v_n = mu * v_0

def proj(v, p):
    return ('inf',) if v[0] % p == 0 else ((v[1]*pow(v[0], p - 2, p)) % p,)

rnd = random.Random(424242)
PR = primes_3mod4(3, 131)

# ---------------------------------------------------------------- (a) and (b)
say("== (a) Legendre pattern of the standard chart = Paley + sink; (b) rescaling switches ==")
for p in PR:
    F = Fp(p)
    S = [(1, x) for x in range(p)] + [(0, 1)]
    assert pattern(F, S) == paley_plus_sink(F)
    for _ in range(20):
        V = S[:]; i = rnd.randrange(p + 1); lam = rnd.randrange(1, p)
        V[i] = ((lam*V[i][0]) % p, (lam*V[i][1]) % p)
        B0, B1 = pattern(F, S), pattern(F, V)
        flipped = {(a, b) for a in range(p + 1) for b in range(p + 1) if a != b and B0[a][b] != B1[a][b]}
        expect = {(a, b) for a in range(p + 1) for b in range(p + 1) if a != b and (a == i or b == i)} if F.chi(lam) == -1 else set()
        assert flipped == expect
say(f"  CONFIRMED for all {len(PR)} primes p = 3 mod 4, p <= 131")

# ---------------------------------------------------------------- (c) random frieze strips
say("== (c) random SL2 frieze strips (distinct points) = switched induced subtournaments ==")
cnt = 0
for p in [7, 11, 19, 23, 31, 43]:
    F = Fp(p); PS = paley_plus_sink(F)
    for _ in range(60):
        m = rnd.randint(3, min(p + 1, 12))
        idx = rnd.sample(range(p + 1), m)
        pts = [((1, x) if x < p else (0, 1)) for x in idx]
        # random representatives, then normalise
        pts = [((l*a) % p, (l*b) % p) for (a, b), l in zip(pts, [rnd.randrange(1, p) for _ in range(m)])]
        V = normalise_strip(pts, p)
        assert all(det(V[k], V[k + 1], p) == 1 for k in range(m - 1))
        T = pattern(F, V)
        # V[k] = lam_k * standard rep; the pattern must be the switching by chi(lam_k) of the induced Paley+sink
        lam = []
        for k in range(m):
            s = (1, idx[k]) if idx[k] < p else (0, 1)
            l = (V[k][0]*pow(s[0], p - 2, p)) % p if s[0] else (V[k][1]*pow(s[1], p - 2, p)) % p
            assert ((l*s[0]) % p, (l*s[1]) % p) == V[k]
            lam.append(F.chi(l))
        assert all(T[a][b] == lam[a]*lam[b]*PS[idx[a]][idx[b]] for a in range(m) for b in range(m) if a != b)
        cnt += 1
say(f"  CONFIRMED on {cnt} random strips (p = 7..43, 3..12 points)")

# ---------------------------------------------------------------- (d) Singer strips, own isomorphism test
say("== (d) Singer strips through all p+1 points; independent switching-isomorphism test ==")
def singer_matrix(p):
    F = Fp(p)
    for t in range(1, p):
        for d in range(1, p):
            if F.chi(t*t - 4*d) != -1: continue          # need irreducible x^2 - t x + d
            g = ((0, (-d) % p), (1, t))                  # companion matrix, columns (0,1)->..., rows listed
            # projective order of g
            M = ((1, 0), (0, 1)); k = 0
            while True:
                M = (((M[0][0]*g[0][0] + M[0][1]*g[1][0]) % p, (M[0][0]*g[0][1] + M[0][1]*g[1][1]) % p),
                     ((M[1][0]*g[0][0] + M[1][1]*g[1][0]) % p, (M[1][0]*g[0][1] + M[1][1]*g[1][1]) % p))
                k += 1
                if M[0][1] == 0 and M[1][0] == 0 and M[0][0] == M[1][1]: break
            if k == p + 1: return g, t, d
    return None
def apply(g, v, p): return ((g[0][0]*v[0] + g[0][1]*v[1]) % p, (g[1][0]*v[0] + g[1][1]*v[1]) % p)

rows = []
for p in PR:
    F = Fp(p)
    g, t, d = singer_matrix(p)
    pts = [(1, 0)]
    for _ in range(p): pts.append(apply(g, pts[-1], p))
    assert len({proj(v, p) for v in pts}) == p + 1                    # all of P^1(F_p)
    V = normalise_strip(pts, p)
    c = closure_constant(V, p)
    quid = [det(V[k - 1], V[k + 1], p) for k in range(1, p)]
    T = pattern(F, V)
    phi, nodes, w = switching_iso(T, paley_plus_sink(F))
    assert phi is not None
    # random orderings of all p+1 points (no Singer structure at all)
    rand_ok, closes = 0, 0
    for _ in range(5 if p <= 47 else 2):
        order = list(range(p + 1)); rnd.shuffle(order)
        P2 = [((1, x) if x < p else (0, 1)) for x in order]
        V2 = normalise_strip(P2, p)
        phi2, _, _ = switching_iso(pattern(F, V2), paley_plus_sink(F))
        assert phi2 is not None; rand_ok += 1
    for _ in range(400):
        order = list(range(p + 1)); rnd.shuffle(order)
        V2 = normalise_strip([((1, x) if x < p else (0, 1)) for x in order], p)
        closes += (closure_constant(V2, p) == p - 1)
    rows.append((p, (t, d), c == p - 1, sorted(set(quid)), (quid[0]*quid[1]) % p == (t*t*pow(d, p - 2, p)) % p,
                 nodes, rand_ok, closes))
    say(f"  p={p:3d}: Singer g = companion(x^2-{t}x+{d}); strip covers P^1; closes antiperiodically (v_(p+1) = -v_0): {c == p - 1};"
        f" quiddity values {sorted(set(quid))} (2-periodic: {len(set(quid[0::2])) == 1 and len(set(quid[1::2])) == 1},"
        f" a*b = t^2/d: {rows[-1][4]}); switching-iso to Paley+sink: yes ({nodes} search nodes);"
        f" random orderings also switching-iso: {rand_ok}/{rand_ok}; random orderings closing as frieze: {closes}/400")

# ---------------------------------------------------------------- (e) constant quiddity
say("== (e) constant-quiddity strips: orbit of e1 under g = [[a,-1],[1,0]] on P^1(F_p) ==")
def mat_order_mod_pm(g, p):
    M = ((1, 0), (0, 1)); k = 0
    ordS, ordP = None, None
    while True:
        M = (((M[0][0]*g[0][0] + M[0][1]*g[1][0]) % p, (M[0][0]*g[0][1] + M[0][1]*g[1][1]) % p),
             ((M[1][0]*g[0][0] + M[1][1]*g[1][0]) % p, (M[1][0]*g[0][1] + M[1][1]*g[1][1]) % p))
        k += 1
        if ordP is None and M[0][1] == 0 and M[1][0] == 0 and M[0][0] == M[1][1] and M[0][0] in (1, p - 1): ordP = k
        if M == ((1, 0), (0, 1)): ordS = k; break
    return ordS, ordP
for p in [7, 11, 19, 23, 31, 43, 47]:
    F = Fp(p)
    summ = {'elliptic': [], 'parabolic': [], 'hyperbolic': []}
    for a in range(p):
        disc = F.chi(a*a - 4)
        kind = 'parabolic' if disc == 0 else ('elliptic' if disc == -1 else 'hyperbolic')
        g = ((a, p - 1), (1, 0))
        orb, v = [], (1, 0)
        seen = set()
        while proj(v, p) not in seen:
            seen.add(proj(v, p)); v = apply(g, v, p)
        # Chebyshev strip v_(k+1) = a v_k - v_(k-1) from (1,0),(0,1): distinct points before repetition
        u0, u1 = (1, 0), (0, 1); pts = [proj(u0, p)]
        while proj(u1, p) not in pts:
            pts.append(proj(u1, p)); u0, u1 = u1, ((a*u1[0] - u0[0]) % p, (a*u1[1] - u0[1]) % p)
        oS, oP = mat_order_mod_pm(g, p)
        assert len(seen) == len(pts)
        if kind == 'elliptic': assert len(seen) == oP and (p + 1) % oS == 0
        summ[kind].append((a, len(seen), oS, oP))
    ell = summ['elliptic']
    say(f"  p={p}: elliptic a: orbit sizes {sorted(set(x[1] for x in ell))} (max {max(x[1] for x in ell)} = (p+1)/2 = {(p+1)//2}),"
        f" SL2 orders {sorted(set(x[2] for x in ell))}, PSL2 orders {sorted(set(x[3] for x in ell))};"
        f" parabolic a=+-2: orbit {sorted(set(x[1] for x in summ['parabolic']))} (= p, not (p+1)/2);"
        f" hyperbolic max orbit {max(x[1] for x in summ['hyperbolic'])} (<= (p-1)/2 = {(p-1)//2})")
    assert max(x[1] for x in ell) == (p + 1)//2
    assert all(x[1] == p for x in summ['parabolic'])
    assert max(x[1] for x in summ['hyperbolic']) <= (p - 1)//2
# parabolic a = -2: g = -(unipotent) has g^p = -I, so the constant strip CLOSES as an SL2 frieze through p points
for p in PR[1:]:
    F = Fp(p); a = p - 2
    V = [(1, 0), (0, 1)]
    while len(V) < p + 2: V.append(((a*V[-1][0] - V[-2][0]) % p, (a*V[-1][1] - V[-2][1]) % p))
    assert V[p] == (p - 1, 0) and V[p + 1] == (0, p - 1)            # v_p = -v_0, v_(p+1) = -v_1
    assert len({proj(v, p) for v in V[:p]}) == p                    # p distinct points (all but [1:-1])
    assert all(det(V[k], V[k + 1], p) == 1 for k in range(p + 1))
    QR = [[F.chi(y - x) if x != y else 0 for y in range(p)] for x in range(p)]
    phi, nodes, w = switching_iso(pattern(F, V[:p]), QR)
    assert phi is not None
say("  CONFIRMED: quiddity a = -2 (parabolic) gives a CLOSED SL2 frieze of width p-3 through p of the p+1 points,"
    " whose Legendre pattern is switching-isomorphic to QR_p itself (all p = 3 mod 4, 7 <= p <= 131)")
# elliptic torus of SL2(F_p): cyclic of order p+1, image of order (p+1)/2 in PSL2
for p in [7, 11, 19, 23]:
    F = Fp(p)
    best = max(mat_order_mod_pm(((a, p - 1), (1, 0)), p) for a in range(p) if F.chi(a*a - 4) == -1)
    assert best == (p + 1, (p + 1)//2)
say("  CONFIRMED: elliptic torus of SL2(F_p) has order p+1, image (p+1)/2 in PSL2 (p = 7, 11, 19, 23)")

# ---------------------------------------------------------------- chi-positive friezes and transitive subtournaments
say("== Reading: chi-positive strips (all pairwise determinants residues) vs transitive subtournaments ==")
def tt_paley(F, reduce=True):
    """largest transitive subtournament of QR_p (x -> y iff y - x residue)."""
    p = F.p
    out = [frozenset((x + s) % p for s in F.sq) for x in range(p)]
    best = [0]
    def rec(depth, cand):
        if depth > best[0]: best[0] = depth
        if depth + len(cand) <= best[0]: return
        for v in sorted(cand):
            rec(depth + 1, cand & out[v])
            if depth + len(cand) <= best[0]: return
    if reduce:
        # arc-transitivity of x -> a x + b (a residue): WLOG the chain starts 0 -> 1
        rec(2, out[0] & out[1])
    else:
        rec(0, frozenset(range(p)))
    return best[0]
def max_chi_positive_strip(F):
    """brute force over the 2(p+1) 'half-points' {+-s}: longest chain with chi(det(v_i,v_j)) = +1 for all i<j."""
    p = F.p
    H = [(1, x) for x in range(p)] + [(0, 1)]
    H = H + [((-a) % p, (-b) % p) for (a, b) in H]
    N = len(H)
    out = [frozenset(j for j in range(N) if j != i and det(H[i], H[j], p) != 0 and F.chi(det(H[i], H[j], p)) == 1) for i in range(N)]
    best = [0]
    def rec(depth, cand):
        if depth > best[0]: best[0] = depth
        if depth + len(cand) <= best[0]: return
        for v in sorted(cand):
            rec(depth + 1, cand & out[v])
            if depth + len(cand) <= best[0]: return
    rec(0, frozenset(range(N)))
    return best[0]
tab = []
for p in primes_3mod4(3, 199):
    F = Fp(p)
    tt = tt_paley(F)
    if p <= 43: assert tt == tt_paley(F, reduce=False)
    L = max_chi_positive_strip(F) if p <= 43 else None
    if L is not None: assert L == tt + 1
    tab.append((p, tt, L, round(math.log2(p), 2), round(2*math.sqrt(p) + 1, 1)))
say("  (p, tt(QR_p), max chi-positive strip [brute force, p<=43], log2 p, 2 sqrt(p)+1):")
for r in tab: say("   ", r)
say("  CONFIRMED: max chi-positive strip = tt(QR_p) + 1 = tt(Paley+sink) for p <= 43 (brute force).")
say("  NOTE: proved upper bound is tt(QR_p) <= 2 sqrt(p) + 1 (bilinear character sum); a logarithmic bound is OPEN.")
say("ALL LEGENDRE CHECKS PASSED")
with open(__file__.replace('.py', '.out'), 'w') as f: f.write("\n".join(OUT) + "\n")
