# Verify residue colourings on concrete unit-distance graphs with EXACT edge detection.
import sys, time, random, itertools, math
from fractions import Fraction as Fr
from mq import MQ
random.seed(1)

def hensel_sqrt_2adic(r, K):
    """s with s^2 = r mod 2^K, r odd int with r = 1 mod 8; s = 1 mod 4."""
    assert r % 8 == 1
    s = 1
    for k in range(3, K):
        if (s*s - r) % (1 << (k+1)) != 0:
            s += 1 << (k-1)
    assert (s*s - r) % (1 << K) == 0
    return s

def hensel_sqrt_padic(r, p, K, s0):
    """odd p: s^2 = r mod p^K with s = s0 mod p."""
    mod = p
    s = s0 % p
    assert (s*s - r) % p == 0
    for k in range(1, K):
        mod *= p
        # Newton step
        s = (s - (s*s - r) * pow(2*s, -1, mod)) % mod
    return s

def rat_mod(a, p, K):
    """a Fraction with denominator prime to p -> int mod p^K"""
    M = p**K
    assert a.denominator % p != 0, a
    return (a.numerator * pow(a.denominator, -1, M)) % M

def two_adic_scaled(a, M, K):
    """2^M * a as int mod 2^K (requires v_2(denominator) <= M)."""
    d = a.denominator; m = 0
    while d % 2 == 0: d //= 2; m += 1
    assert m <= M
    return (a.numerator * (1 << (M - m)) * pow(d, -1, 1 << K)) % (1 << K)

# ---------------- residue maps ----------------
def make_res_2adic(F, dlist, K=80, M=40):
    """L = Q(sqrt(-d_i)), all d_i = 3 mod 8 (d_1 may be 3).  P | 2 with L_P = Q_2(w), sqrt(-d_i) -> (2w+1)*s_i, s_i^2 = d_i/3.
       Returns x -> residue in F_4 encoded as (bit0, bit1) for  bit0 + bit1*w ."""
    MOD = 1 << K
    s = []
    for d in dlist:
        r = (d * pow(3, -1, MOD)) % MOD     # d/3 in Z_2
        s.append(hensel_sqrt_2adic(r, K))
    k = len(dlist)
    def res(x):
        A = 0; B = 0
        for S, a in enumerate(x):
            if a == 0: continue
            bits = [i for i in range(k) if S >> i & 1]
            prod = 1
            for i in bits: prod = prod * s[i] % MOD
            h = len(bits)
            fac = pow(-3, h//2, MOD) % MOD
            val = two_adic_scaled(Fr(a), M, K) * prod * fac % MOD
            if h % 2 == 0: A = (A + val) % MOD
            else: B = (B + val) % MOD
        AB = (A + B) % MOD; B2 = (2*B) % MOD
        # integrality: 2^M * (A+B) and 2^M*(2B) must be divisible by 2^M
        assert AB % (1 << M) == 0 and B2 % (1 << M) == 0, "not P-integral"
        return ((AB >> M) & 1, (B2 >> M) & 1)
    return res

# ---------------- patches ----------------
def build_patch(F, gens, depth, cap, extra_units=()):
    units = list(gens) + list(extra_units)
    allu = units + [tuple(-c for c in u) for u in units]
    zero = tuple([Fr(0)]*F.n)
    pts = {zero: 0}; frontier = [zero]
    for _ in range(depth):
        new = []
        for x in frontier:
            for u in allu:
                y = F.add(x, u)
                if y not in pts:
                    pts[y] = len(pts); new.append(y)
                    if len(pts) >= cap: break
            if len(pts) >= cap: break
        frontier = new
        if len(pts) >= cap: break
    return list(pts)

def exact_edges(F, pts, tol=1e-7):
    vals = [F.numval(p) for p in pts]
    grid = {}
    for i, z in enumerate(vals):
        grid.setdefault((math.floor(z.real), math.floor(z.imag)), []).append(i)
    cand = 0; edges = []; diffs = set()
    for i, z in enumerate(vals):
        gx, gy = math.floor(z.real), math.floor(z.imag)
        for dx in (-1, 0, 1, 2):
            for dy in (-1, 0, 1, 2):
                for j in grid.get((gx+dx-1+1-1, gy+dy-1+1-1), []):
                    pass
        for dx in range(-2, 3):
            for dy in range(-2, 3):
                for j in grid.get((gx+dx, gy+dy), []):
                    if j <= i: continue
                    if abs(abs(vals[j]-z) - 1) < tol:
                        cand += 1
                        d = F.sub(pts[j], pts[i])
                        if F.is_one(F.norm2(d)):
                            edges.append((i, j)); diffs.add(d if d > tuple(-c for c in d) else tuple(-c for c in d))
    return edges, cand, diffs

def check(name, F, pts, color, edges, diffs, ncolors_claim):
    cols = [color(p) for p in pts]
    bad = sum(1 for i, j in edges if cols[i] == cols[j])
    used = len(set(cols))
    print(f"[{name}] points={len(pts)} exact unit edges={len(edges)} distinct unit directions(+-)={len(diffs)} "
          f"colours used={used} (claim <= {ncolors_claim}) monochromatic edges={bad}", flush=True)
    return bad == 0

def sat_colorable(n, edges, k):
    from pysat.solvers import Cadical153
    v = lambda i, c: i*k + c + 1
    s = Cadical153()
    for i in range(n):
        s.add_clause([v(i, c) for c in range(k)])
    for a, b in edges:
        for c in range(k): s.add_clause([-v(a, c), -v(b, c)])
    r = s.solve(); s.delete(); return r

def hilbert90_units(F, count, rng=3):
    out = []
    for _ in range(count):
        z = F.elt({S: random.randint(-rng, rng) for S in range(F.n)})
        if all(c == 0 for c in z): continue
        u = F.mul(z, F.inv(F.conj(z)))
        assert F.is_one(F.norm2(u))
        out.append(u)
    return out

if __name__ == "__main__":
    allok = True
    # ---- A. Moser field Q(sqrt-3, sqrt-11): F_4 colouring at P|2 (Gibbs 2018 for the ring; here the field) ----
    F = MQ([-3, -11])
    w1 = F.elt({0: Fr(1, 2), 1: Fr(1, 2)}); w3 = F.elt({0: Fr(5, 6), 2: Fr(1, 6)})
    w3b = F.conj(w3)
    gens = []
    z = F.one()
    for j in range(6):
        gens += [z, F.mul(z, w3), F.mul(z, w3b)]
        z = F.mul(z, w1)
    res = make_res_2adic(F, [3, 11])
    for u in gens: assert res(u) != (0, 0)
    pts = build_patch(F, gens, 3, 2500)
    edges, cand, diffs = exact_edges(F, pts)
    ok = check("Moser ring Z[w1,w3] patch", F, pts, res, edges, diffs, 4); allok &= ok
    # sub-patch 3-colourability (should be UNSAT: contains spindles)
    sub = pts[:400]; idx = {p: i for i, p in enumerate(sub)}
    e2 = [(i, j) for i, j in edges if i < 400 and j < 400]
    print("   3-colourable (first 400 pts)?", sat_colorable(400, e2, 3), " 4-colourable?", sat_colorable(400, e2, 4))
    # field points: Hilbert-90 unit vectors (not in the ring), mixed patch
    hu = hilbert90_units(F, 25)
    for u in hu:
        r = res(u); assert r != (0, 0)
    pts2 = build_patch(F, gens[:6] + hu[:10], 2, 2500)
    edges2, _, diffs2 = exact_edges(F, pts2)
    ok = check("Moser FIELD patch (incl. Hilbert-90 units)", F, pts2, res, edges2, diffs2, 4); allok &= ok
    print("ALL OK so far:", allok, flush=True)
