#!/usr/bin/env python3
"""
THE TOURNAMENT SIMPLEX: simplicial geometry x odd maps
macmini-2026-06-10-S2  (exploratory; NO canon claims)

For a tournament T with adjacency A, S = A - A^T, M = I + S.
Rows of M are n vertices of the cube {±1}^n.  The TOURNAMENT SIMPLEX is
conv{0, row_1, ..., row_n} (the cone simplex; 0 = cube center), with
vol(T) = det(M)/n!   (det(M) > 0 always: eigenvalues of I+S are 1+i*mu).

Part 1: volume spectrum per iso class, n=4..7; ratios to max-volume simplices;
        switching = coordinate reflection + selective antipodal moves (verified).
Part 2: Gram geometry. Gram(M) = M M^T = I + S S^T exactly.
        Distance dictionary: |row_i - row_j|^2 = 4*(1 + disagree(i,j))
        where disagree(i,j) = #{k != i,j : exactly one of i,j beats k}.
        Regular / orthocentric / equifacetal census; shape spectrum.
Part 3: odd maps. x -> Sx is an odd tangent field. Pfaffian vector w
        (adj(S) = w w^T, odd n) = the hairy-ball singularity. Sign/score laws.

All determinants/adjugates in EXACT integer arithmetic (fraction-free Bareiss).
"""

import subprocess, sys
from math import comb, factorial, isqrt, gcd, sqrt
from fractions import Fraction
from itertools import permutations, combinations
from collections import Counter, defaultdict

sys.stdout.reconfigure(line_buffering=True)

# ----------------------------------------------------------------------
# exact linear algebra
# ----------------------------------------------------------------------

def det_int(mat):
    """Fraction-free Bareiss determinant on Python ints. Exact."""
    a = [row[:] for row in mat]
    n = len(a)
    if n == 0: return 1
    sign = 1
    prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            for r in range(k + 1, n):
                if a[r][k] != 0:
                    a[k], a[r] = a[r], a[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) // prev
        prev = a[k][k]
    return sign * a[n - 1][n - 1]

def minor(mat, r, c):
    return [[mat[i][j] for j in range(len(mat)) if j != c]
            for i in range(len(mat)) if i != r]

def adjugate(mat):
    """adj[i][j] = (-1)^(i+j) det(mat with row j, col i deleted). Exact."""
    n = len(mat)
    return [[(-1) ** (i + j) * det_int(minor(mat, j, i)) for j in range(n)]
            for i in range(n)]

def matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)]
            for i in range(n)]

# ----------------------------------------------------------------------
# tournament enumeration (gentourng, one rep per iso class)
# ----------------------------------------------------------------------

def tournaments(n):
    m = comb(n, 2)
    P = [(i, j) for i in range(n) for j in range(i + 1, n)]
    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    out = []
    for line in (r.stdout or '').split('\n'):
        line = line.strip()
        if len(line) == m and set(line) <= {'0', '1'}:
            A = [[0] * n for _ in range(n)]
            for k, (i, j) in enumerate(P):
                if line[k] == '1':
                    A[i][j] = 1
                else:
                    A[j][i] = 1
            out.append(A)
    return out

def skew(A):
    n = len(A)
    return [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]

def Mmat(A):
    n = len(A)
    S = skew(A)
    return [[(1 if i == j else S[i][j]) for j in range(n)] for i in range(n)]

def scores(A):
    return [sum(r) for r in A]

# ----------------------------------------------------------------------
# shape machinery
# ----------------------------------------------------------------------

MAXDET_PM1 = {1: 1, 2: 2, 3: 4, 4: 16, 5: 48, 6: 160, 7: 576, 8: 4096}
# Hadamard maximal determinant for +-1 matrices of order n (classical values)

def gram(M):
    n = len(M)
    return [[sum(M[i][k] * M[j][k] for k in range(n)) for j in range(n)]
            for i in range(n)]

def dist2_matrix(G):
    """squared distances between rows: d2_ij = G_ii + G_jj - 2 G_ij"""
    n = len(G)
    return [[G[i][i] + G[j][j] - 2 * G[i][j] for j in range(n)] for i in range(n)]

def disagree_count(A, i, j):
    n = len(A)
    c = 0
    for k in range(n):
        if k == i or k == j:
            continue
        # exactly one of i,j beats k  <=>  S_ik != S_jk
        if (A[i][k] - A[k][i]) != (A[j][k] - A[k][j]):
            c += 1
    return c

def codegrees(A, i, j):
    n = len(A)
    cp = sum(1 for k in range(n) if k != i and k != j and A[i][k] and A[j][k])
    cm = sum(1 for k in range(n) if k != i and k != j and A[k][i] and A[k][j])
    return cp, cm

def is_orthocentric(d2):
    """n-point simplex with squared distances d2: exists t_i with
    d2_ij = t_i + t_j for all i<j?  Exact (Fractions)."""
    n = len(d2)
    if n < 3:
        return True
    t0 = Fraction(d2[0][1] + d2[0][2] - d2[1][2], 2)
    t = [t0] + [Fraction(d2[0][i]) - t0 for i in range(1, n)]
    for i in range(1, n):
        for j in range(i + 1, n):
            if t[i] + t[j] != d2[i][j]:
                return False
    return True

_PERM_CACHE = {}
def facet_canon(d2sub):
    """canonical form of a k-point distance matrix under vertex permutation"""
    k = len(d2sub)
    if k not in _PERM_CACHE:
        _PERM_CACHE[k] = list(permutations(range(k)))
    best = None
    for p in _PERM_CACHE[k]:
        f = tuple(d2sub[p[i]][p[j]] for i in range(k) for j in range(k))
        if best is None or f < best:
            best = f
    return best

def is_equifacetal(d2):
    """(n-1)-simplex on n rows: all n facets congruent?"""
    n = len(d2)
    facets = []
    for drop in range(n):
        idx = [i for i in range(n) if i != drop]
        sub = [[d2[a][b] for b in idx] for a in idx]
        facets.append(sub)
    # cheap invariant: multiset of all distances per facet
    inv = [tuple(sorted(sub[a][b] for a in range(n - 1) for b in range(a + 1, n - 1)))
           for sub in facets]
    if len(set(inv)) > 1:
        return False
    # cheap invariant 2: per-vertex distance multisets
    inv2 = [tuple(sorted(tuple(sorted(row[b] for b in range(n - 1) if b != a))
                         for a, row in enumerate(sub)))
            for sub in facets]
    if len(set(inv2)) > 1:
        return False
    canons = {facet_canon(sub) for sub in facets}
    return len(canons) == 1

def shape_spectrum(d2):
    """multiset of squared distances /4 = Hamming distances (sorted tuple)"""
    n = len(d2)
    return tuple(sorted(d2[i][j] // 4 for i in range(n) for j in range(i + 1, n)))

# ----------------------------------------------------------------------
# Pfaffian vector (odd n) and harmonic direction (even n)
# ----------------------------------------------------------------------

def pfaffian_vector(S):
    """odd n: adj(S) = w w^T; return w with sum(w) > 0. Exact."""
    n = len(S)
    adj = adjugate(S)
    w0sq = adj[0][0]
    assert w0sq > 0, "adj_00 = w_0^2 should be positive"
    w0 = isqrt(w0sq)
    assert w0 * w0 == w0sq, "adj_00 not a perfect square?!"
    w = [w0]
    for i in range(1, n):
        assert adj[i][0] % w0 == 0
        w.append(adj[i][0] // w0)
    # verify rank-1 factorization
    for i in range(n):
        for j in range(n):
            assert adj[i][j] == w[i] * w[j], "adj != w w^T"
    # verify kernel
    for i in range(n):
        assert sum(S[i][j] * w[j] for j in range(n)) == 0, "S w != 0"
    if sum(w) < 0:
        w = [-x for x in w]
    return w

def harmonic_direction(S):
    """even n: u = adj(S) . 1  (so S u = det(S) 1).  Exact."""
    n = len(S)
    adj = adjugate(S)
    u = [sum(adj[i][j] for j in range(n)) for i in range(n)]
    d = det_int(S)
    for i in range(n):
        assert sum(S[i][j] * u[j] for j in range(n)) == d, "S u != det * 1"
    assert sum(u) == 0, "1^T adj(S) 1 != 0"
    return u, d

# ----------------------------------------------------------------------
# main analysis per n
# ----------------------------------------------------------------------

def pearson(xs, ys):
    n = len(xs)
    if n < 2: return float('nan')
    mx = sum(xs) / n; my = sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0: return float('nan')
    return cov / sqrt(vx * vy)

def banner(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)

def analyze(n, full_table):
    banner(f"n = {n}")
    Ts = tournaments(n)
    print(f"iso classes: {len(Ts)}")
    nf = factorial(n)
    recs = []
    for cid, A in enumerate(Ts):
        S = skew(A)
        M = Mmat(A)
        dM = det_int(M)
        assert dM > 0
        assert dM % (1 << (n - 1)) == 0, "det(I+S) not divisible by 2^(n-1)?!"
        dT = dM >> (n - 1)
        G = gram(M)
        # Gram identity: M M^T = I + S S^T  (since S + S^T = 0)
        SSt = matmul(S, [[S[j][i] for j in range(n)] for i in range(n)])
        for i in range(n):
            for j in range(n):
                expect = (1 if i == j else 0) + SSt[i][j]
                assert G[i][j] == expect, "Gram != I + S S^T"
        assert det_int(G) == dM * dM, "det(M M^T) != det(M)^2"
        d2 = dist2_matrix(G)
        # distance dictionary check
        for i in range(n):
            for j in range(i + 1, n):
                dis = disagree_count(A, i, j)
                cp, cm = codegrees(A, i, j)
                assert d2[i][j] == 4 * (1 + dis), "d2 != 4(1+disagree)"
                assert dis == (n - 2) - (cp + cm)
                assert d2[i][j] == 4 * (n - 1) - 4 * (cp + cm), "co-degree dict fails"
                # Hamming check
                ham = sum(1 for k in range(n) if M[i][k] != M[j][k])
                assert d2[i][j] == 4 * ham, "d2 != 4*Hamming"
        spec = shape_spectrum(d2)
        regular = (len(set(spec)) == 1)
        ortho_rows = is_orthocentric(d2)
        # cone simplex {0, rows}: d2 to origin = n for every row
        d2c = [[0] * (n + 1) for _ in range(n + 1)]
        for i in range(n):
            d2c[0][i + 1] = d2c[i + 1][0] = n
            for j in range(n):
                d2c[i + 1][j + 1] = d2[i][j]
        ortho_cone = is_orthocentric(d2c)
        equif = is_equifacetal(d2)
        rec = dict(cid=cid, A=A, S=S, M=M, det=dM, dT=dT, sc=tuple(scores(A)),
                   spec=spec, regular=regular, ortho=ortho_rows,
                   ortho_cone=ortho_cone, equif=equif, G=G, d2=d2)
        recs.append(rec)

    # ---------------- Part 1: volume spectrum ----------------
    md_n = MAXDET_PM1[n]
    md_n1 = MAXDET_PM1[n + 1]
    if n % 2 == 1:
        ceil_t = (n + 1) ** ((n - 1) // 2)   # DRT ceiling
        ceil_lbl = f"(n+1)^((n-1)/2) = {ceil_t} [DRT ceiling]"
    else:
        ceil_t = isqrt(n ** n)
        ceil_lbl = f"n^(n/2) = {ceil_t} [skew-conference ceiling]"
    dets = sorted({r['det'] for r in recs})
    print(f"\nPART 1 — VOLUME SPECTRUM  (vol = det(M)/n!, det(M) = 2^(n-1) d(T))")
    print(f"  distinct det values ({len(dets)}): {dets}")
    print(f"  distinct d(T) values: {sorted({r['dT'] for r in recs})}")
    cnt = Counter(r['det'] for r in recs)
    print(f"  det multiplicities: {dict(sorted(cnt.items()))}")
    mx = max(dets)
    print(f"  tournament max det = {mx};  ceiling {ceil_lbl}; attained: {mx == ceil_t}")
    print(f"  max-vol simplex {{center + n cube vertices}}: maxdet_pm1({n}) = {md_n}"
          f"  -> tournament/maxdet ratio = {mx}/{md_n} = {mx/md_n:.4f}")
    print(f"  true max-vol inscribed simplex (n+1 cube vertices): maxdet_pm1({n+1})/n!"
          f" -> ratio = {mx}/{md_n1} = {mx/md_n1:.4f}")
    if n % 2 == 0:
        print(f"  regular-simplex ratio (b): det/n^(n/2): max class = {mx}/{ceil_t} = {mx/ceil_t:.4f}"
              f"  (edge sqrt(2n)={sqrt(2*n):.4f} regular simplex fits iff skew-conference exists)")
    else:
        print(f"  NOTE: at odd n the regular tournament simplex has edge sqrt(2n+2)={sqrt(2*n+2):.4f},"
              f" NOT sqrt(2n) — Gram = (n+1)I - J (DRT), d2 = 2n+2.")

    # switching experiment
    print(f"\n  SWITCHING (M -> D M D, D = diag(+-1)):")
    n_spec_changed = 0
    example = None
    for r in recs:
        changed = False
        for mask in range(1, 1 << (n - 1)):
            d = [1] + [(-1 if (mask >> (i - 1)) & 1 else 1) for i in range(1, n)]
            Msw = [[d[i] * r['M'][i][j] * d[j] for j in range(n)] for i in range(n)]
            assert det_int(Msw) == r['det'], "switching changed det!"
            # geometric claim: row_i' = d_i * (coordinate-reflected row_i)
            for i in range(n):
                for j in range(n):
                    assert Msw[i][j] == d[i] * (r['M'][i][j] * d[j])
            Gsw = gram(Msw)
            spec_sw = shape_spectrum(dist2_matrix(Gsw))
            if spec_sw != r['spec']:
                changed = True
                if example is None:
                    example = (r['cid'], mask, r['spec'], spec_sw)
        if changed:
            n_spec_changed += 1
    print(f"  det(M) invariant under ALL 2^(n-1) switchings of ALL classes: VERIFIED")
    print(f"  classes where some switching changes the distance spectrum: "
          f"{n_spec_changed}/{len(recs)}")
    if example:
        cid, mask, s_old, s_new = example
        print(f"  example: class {cid}, switch-set mask {mask}: spectrum {s_old} -> {s_new}")

    # ---------------- Part 2: shape census ----------------
    print(f"\nPART 2 — GRAM / SHAPE GEOMETRY  (Gram = I + S S^T; d2_ij = 4(1+disagree))")
    n_reg = sum(r['regular'] for r in recs)
    n_ortho = sum(r['ortho'] for r in recs)
    n_oc = sum(r['ortho_cone'] for r in recs)
    n_eq = sum(r['equif'] for r in recs)
    print(f"  regular row-simplex: {n_reg}/{len(recs)} classes")
    print(f"  orthocentric row-simplex (n rows): {n_ortho}/{len(recs)}")
    print(f"  orthocentric cone simplex (0 + rows): {n_oc}/{len(recs)}"
          f"   [derived: cone orthocentric <=> rows regular: {'CONSISTENT' if n_oc == n_reg and all(r['ortho_cone'] == r['regular'] for r in recs) else 'VIOLATED'}]")
    print(f"  equifacetal row-simplex: {n_eq}/{len(recs)}")
    nonreg_eq = [r['cid'] for r in recs if r['equif'] and not r['regular']]
    nonreg_or = [r['cid'] for r in recs if r['ortho'] and not r['regular']]
    print(f"  equifacetal but NOT regular: {nonreg_eq if nonreg_eq else 'none'}")
    print(f"  orthocentric rows but NOT regular: {nonreg_or if nonreg_or else 'none'}")
    specs = Counter(r['spec'] for r in recs)
    print(f"  distinct shape spectra: {len(specs)} over {len(recs)} classes"
          f"  (shape spectrum {'SEPARATES' if len(specs) == len(recs) else 'does NOT separate'} iso classes)")
    coll = {s: c for s, c in specs.items() if c > 1}
    if coll:
        worst = max(coll.values())
        print(f"  largest spectrum collision: {worst} classes share one spectrum"
              f"  ({len(coll)} colliding spectra)")
        # does (spectrum, det) separate?
        sd = Counter((r['spec'], r['det']) for r in recs)
        print(f"  (spectrum, det) pairs: {len(sd)} distinct "
              f"({'separates' if len(sd) == len(recs) else 'still collides'})")
    if full_table:
        print(f"\n  per-class table (Ham spectrum = d2/4 multiset):")
        print(f"  {'cid':>3} {'scores':<16} {'det':>5} {'d(T)':>4} {'reg':>3} {'orth':>4} {'equi':>4}  spectrum")
        for r in recs:
            print(f"  {r['cid']:>3} {str(r['sc']):<16} {r['det']:>5} {r['dT']:>4}"
                  f" {'Y' if r['regular'] else '.':>3} {'Y' if r['ortho'] else '.':>4}"
                  f" {'Y' if r['equif'] else '.':>4}  {r['spec']}")
    else:
        # condensed: histogram of spectra with counts, plus special classes
        print(f"\n  spectrum histogram (top 12 by multiplicity):")
        for s, c in specs.most_common(12):
            print(f"    x{c:>3}  {s}")
        for r in recs:
            if r['regular'] or r['equif'] or r['ortho']:
                print(f"  special class {r['cid']}: scores {r['sc']} det {r['det']}"
                      f" reg={r['regular']} ortho={r['ortho']} equif={r['equif']} spec {r['spec']}")

    # ---------------- Part 3: odd maps ----------------
    print(f"\nPART 3 — ODD TANGENT FIELD x -> Sx")
    if n % 2 == 0:
        print(f"  even n: det(S) = Pf(S)^2 odd => ker S = 0 for ALL classes:")
        allinv = all(det_int(r['S']) != 0 for r in recs)
        print(f"    verified S nonsingular on all {len(recs)} classes: {allinv}")
        print(f"  projected field P(Sx) on sum-zero sphere S^{n-2} (EVEN-dim) must vanish (hairy ball).")
        print(f"  unique antipodal zero pair: x* ∝ S^(-1) 1  (in 1^perp automatically since 1^T S^-1 1 = 0).")
        rows = []
        for r in recs:
            u, dS = harmonic_direction(r['S'])
            g = 0
            for x in u: g = gcd(g, abs(x))
            r['u'] = u; r['ugcd'] = g; r['detS'] = dS
            rows.append(r)
        # correlations with score
        xs = []; ys = []
        sign_score = Counter()
        for r in rows:
            for i in range(n):
                xs.append(r['sc'][i]); ys.append(r['u'][i])
                if r['u'][i] != 0:
                    sign_score[(r['sc'][i], 1 if r['u'][i] > 0 else -1)] += 1
        print(f"  harmonic direction u = adj(S)·1: sum(u)=0 verified; S u = det(S)·1 verified.")
        print(f"  Pearson corr(u_i, s_i) over all (class,vertex): {pearson(xs, ys):+.4f}")
        nz = sum(1 for r in rows for i in range(n) if r['u'][i] == 0)
        print(f"  zero entries of u: {nz} of {len(rows)*n}")
        print(f"  gcd(u) distribution: {dict(Counter(r['ugcd'] for r in rows))}")
        par = Counter(tuple(x % 2 for x in r['u']) for r in rows)
        print(f"  parity patterns of u (per class): {dict(par)}")
        # test linear law u_i ~ a*(s_i) + b per class?
        lin = 0
        for r in rows:
            # check if u_i is an affine function of s_i within the class
            pairs = {}
            ok = True
            for i in range(n):
                if r['sc'][i] in pairs and pairs[r['sc'][i]] != r['u'][i]:
                    ok = False; break
                pairs[r['sc'][i]] = r['u'][i]
            if ok and len(pairs) >= 2:
                ss = sorted(pairs)
                a = Fraction(pairs[ss[1]] - pairs[ss[0]], ss[1] - ss[0])
                b = pairs[ss[0]] - a * ss[0]
                ok = all(pairs[s] == a * s + b for s in ss)
            lin += ok
        print(f"  classes where u_i is a well-defined affine function of score s_i: {lin}/{len(rows)}")
        if full_table:
            print(f"  per-class u:")
            for r in rows:
                print(f"    cid {r['cid']:>2} scores {r['sc']} det(S)={r['detS']:>4} u={r['u']}")
        # numeric index of the zero (sample): sign of det of projected linearization
        import random
        random.seed(7)
        idxsigns = Counter()
        for r in rows:
            u = r['u']
            nu = sqrt(sum(x * x for x in u))
            xstar = [x / nu for x in u]
            # tangent space T = 1^perp ∩ xstar^perp; build orthonormal basis numerically
            basis = []
            cand = [[1.0 if k == i else 0.0 for k in range(n)] for i in range(n)]
            ones = [1 / sqrt(n)] * n
            for v in cand:
                w2 = v[:]
                for b in ([ones, xstar] + basis):
                    dp = sum(a * c for a, c in zip(w2, b))
                    w2 = [a - dp * c for a, c in zip(w2, b)]
                nv = sqrt(sum(a * a for a in w2))
                if nv > 1e-9:
                    basis.append([a / nv for a in w2])
                if len(basis) == n - 2:
                    break
            # field F(x) = P Sx ; linearization on T: L[a][b] = e_a^T P S e_b (P = proj 1^perp)
            L = [[0.0] * (n - 2) for _ in range(n - 2)]
            for bidx, bv in enumerate(basis):
                Sv = [sum(r['S'][i][j] * bv[j] for j in range(n)) for i in range(n)]
                mean = sum(Sv) / n
                PSv = [x - mean for x in Sv]
                # subtract radial part along xstar then project to T
                dpx = sum(a * c for a, c in zip(PSv, xstar))
                PSv = [a - dpx * c for a, c in zip(PSv, xstar)]
                for aidx, av in enumerate(basis):
                    L[aidx][bidx] = sum(p * q for p, q in zip(PSv, av))
            # numeric determinant sign (n-2 small)
            import copy
            mat = [row[:] for row in L]
            sgn = 1.0
            for k in range(n - 2):
                piv = max(range(k, n - 2), key=lambda r2: abs(mat[r2][k]))
                if abs(mat[piv][k]) < 1e-9:
                    sgn = 0.0; break
                if piv != k:
                    mat[k], mat[piv] = mat[piv], mat[k]; sgn = -sgn
                for i2 in range(k + 1, n - 2):
                    f = mat[i2][k] / mat[k][k]
                    for j2 in range(k, n - 2):
                        mat[i2][j2] -= f * mat[k][j2]
                sgn *= (1.0 if mat[k][k] > 0 else -1.0)
            idxsigns[int(sgn)] += 1
        print(f"  index of zero x* (sign det of linearized field on tangent space): {dict(idxsigns)}"
              f"  [each antipodal zero same index; sum must be chi(S^{n-2}) = 2]")
    else:
        print(f"  odd n: rank(S) = n-1 always; ker S = span(w), adj(S) = w w^T (Pfaffian vector).")
        print(f"  full sphere S^{n-1} is EVEN-dim: hairy ball forces a zero; zeros are exactly ±w/|w|.")
        print(f"  sum-zero sphere S^{n-2} is ODD-dim: field x->Sx NONVANISHING there because")
        print(f"  sum(w) = sqrt(det(I+2A)) is ODD (never 0) => w not in 1^perp. Verify below.")
        all_odd = True
        sign_law = Counter()       # epsilon_i = sign(w_i) * (-1)^{s_i} constant per class?
        mod4_law = Counter()
        sumw_vals = []
        gcds = Counter()
        xs_abs = []; ys_score = []; ys_dev = []
        rows = []
        for r in recs:
            w = pfaffian_vector(r['S'])
            r['w'] = w
            # check det(I+2A) = (sum w)^2
            I2A = [[(1 if i == j else 2 * r['A'][i][j]) for j in range(n)] for i in range(n)]
            assert det_int(I2A) == sum(w) ** 2, "det(I+2A) != (1^T w)^2"
            assert sum(w) % 2 == 1
            sumw_vals.append(sum(w))
            if any(x % 2 == 0 for x in w):
                all_odd = False
            g = 0
            for x in w: g = gcd(g, abs(x))
            gcds[g] += 1
            eps = {(1 if w[i] > 0 else -1) * (-1) ** (r['sc'][i]) for i in range(n)}
            sign_law['constant' if len(eps) == 1 else 'mixed'] += 1
            m4 = {(w[i] - (-1) ** (r['sc'][i])) % 4 for i in range(n)}
            m4b = {(w[i] + (-1) ** (r['sc'][i])) % 4 for i in range(n)}
            mod4_law['w_i ≡ ±(-1)^s_i mod 4 (global sign)' if (m4 == {0} or m4b == {0})
                      else 'fails'] += 1
            for i in range(n):
                xs_abs.append(abs(w[i]))
                ys_score.append(r['sc'][i])
                ys_dev.append((r['sc'][i] - (n - 1) / 2) ** 2)
            rows.append(r)
        print(f"  all w_i odd (never zero): {all_odd}  [w_i = ±Pf(S minus i), Pf odd]")
        print(f"  sum(w) values (odd, = sqrt(det(I+2A))): {sorted(set(sumw_vals))}")
        print(f"  gcd(w) distribution: {dict(gcds)}")
        print(f"  SIGN-SCORE LAW  sign(w_i)·(-1)^(s_i) constant per class: {dict(sign_law)}")
        print(f"  MOD-4 LAW  w_i ≡ ±(-1)^(s_i) (mod 4) with one global sign: {dict(mod4_law)}")
        print(f"  Pearson corr(|w_i|, s_i): {pearson(xs_abs, ys_score):+.4f}")
        print(f"  Pearson corr(|w_i|, (s_i-mean)^2): {pearson(xs_abs, ys_dev):+.4f}")
        # flow-balance interpretation: sum_{j: i->j} w_j = sum_{j: j->i} w_j
        print(f"  flow balance (S w = 0): for every vertex i, sum of w over out-neighbors")
        print(f"    equals sum of w over in-neighbors — w is the unique 'parity-balanced' weighting.")
        if full_table:
            print(f"  per-class w (sum(w)>0 normalization):")
            for r in rows:
                print(f"    cid {r['cid']:>2} scores {r['sc']} w={r['w']} sum={sum(r['w'])} det(I+2A)={sum(r['w'])**2}")
        else:
            # most interesting: regular classes, transitive, extremes of |w|
            mx = max(rows, key=lambda r: max(abs(x) for x in r['w']))
            print(f"  largest |w_i| at n={n}: class {mx['cid']} scores {mx['sc']} w={mx['w']}")
            for r in rows:
                if len(set(r['sc'])) == 1:  # regular tournament
                    print(f"  regular tournament class {r['cid']}: w={r['w']}"
                          f"  (constant iff arc-transitive-ish)")
    return recs

# ----------------------------------------------------------------------
# QR_7 spot check
# ----------------------------------------------------------------------

def qr7_check():
    banner("QR_7 (Paley) SPOT CHECK")
    n = 7
    QR = {1, 2, 4}
    A = [[1 if ((j - i) % 7) in QR else 0 for j in range(7)] for i in range(7)]
    S = skew(A)
    M = Mmat(A)
    dM = det_int(M)
    G = gram(M)
    d2 = dist2_matrix(G)
    offd = {G[i][j] for i in range(7) for j in range(7) if i != j}
    dd = {d2[i][j] for i in range(7) for j in range(7) if i != j}
    print(f"det(M) = {dM} (DRT ceiling (n+1)^3 = {8**3}); vol = {dM}/{factorial(7)}")
    print(f"Gram off-diagonal values: {offd}  => Gram = 8I - J: {offd == {-1}}")
    print(f"pairwise squared distances: {dd}  => all = 2n+2 = 16, edge = 4 = sqrt(2n+2)")
    print(f"REGULAR simplex: {len(dd) == 1}; edge sqrt(2n)=sqrt(14) does NOT fit at odd n.")
    w = pfaffian_vector(S)
    print(f"Pfaffian vector w = {w} (vertex-transitivity => w ∝ 1)")
    # SS^T = nI - J check
    SSt = matmul(S, [[S[j][i] for j in range(7)] for i in range(7)])
    ok = all(SSt[i][j] == (7 * (i == j) - 1) for i in range(7) for j in range(7))
    print(f"SS^T = 7I - J (DRT): {ok}")

# ----------------------------------------------------------------------
# addendum: law verification and special classes
# ----------------------------------------------------------------------

def c3_through(A, i):
    n = len(A)
    c = 0
    for j in range(n):
        for k in range(j + 1, n):
            if j == i or k == i:
                continue
            trio = (i, j, k)
            # cyclic iff every vertex has out-degree 1 within the triple
            if all(sum(A[a][b] for b in trio) == 1 for a in trio):
                c += 1
    return c

def canon_form(A):
    n = len(A)
    best = None
    for p in permutations(range(n)):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or f < best:
            best = f
    return best

def addendum():
    banner("ADDENDUM — LAW VERIFICATION & SPECIAL CLASSES")

    # A. sign-free mod-4 law (odd n): adj(S)_ij ≡ (-1)^(s_i+s_j)  (mod 4)
    for n in (5, 7):
        ok = True
        for A in tournaments(n):
            S = skew(A); sc = scores(A); adj = adjugate(S)
            for i in range(n):
                for j in range(n):
                    if (adj[i][j] - (-1) ** (sc[i] + sc[j])) % 4 != 0:
                        ok = False
        print(f"A. odd n={n}: adj(S)_ij ≡ (-1)^(s_i+s_j) (mod 4), ALL entries incl diagonal,"
              f" ALL classes: {ok}")

    # B. even-n analog on the harmonic vector u = adj(S)·1
    for n in (4, 6):
        cp = cm = cx = 0
        adjm4 = Counter()
        for A in tournaments(n):
            S = skew(A); sc = scores(A)
            u, _ = harmonic_direction(S)
            rp = all((u[i] - (-1) ** sc[i]) % 4 == 0 for i in range(n))
            rm = all((u[i] + (-1) ** sc[i]) % 4 == 0 for i in range(n))
            if rp: cp += 1
            elif rm: cm += 1
            else: cx += 1
            adj = adjugate(S)
            for i in range(n):
                for j in range(i + 1, n):
                    adjm4[(adj[i][j] % 4, (sc[i] + sc[j]) % 2)] += 1
        print(f"B. even n={n}: u_i ≡ +(-1)^(s_i) mod 4 globally: {cp};"
              f" ≡ -(-1)^(s_i): {cm}; neither: {cx}")
        print(f"   (adj is SKEW at even n; off-diag adj_ij mod 4 vs parity(s_i+s_j):"
              f" {dict(sorted(adjm4.items()))})")

    # C. mod-8 refinement exploration (odd n): after the mod-4 law, the residue
    #    r_i = eps * w_i * (-1)^(s_i) mod 8 lies in {1,5}. What controls it?
    for n in (5, 7):
        tab_c3 = Counter(); tab_s2 = Counter(); tab_tri = Counter()
        for A in tournaments(n):
            S = skew(A); sc = scores(A)
            w = pfaffian_vector(S)
            eps = 1 if all((w[i] - (-1) ** sc[i]) % 4 == 0 for i in range(n)) else -1
            for i in range(n):
                r = (eps * w[i] * (-1) ** sc[i]) % 8
                assert r in (1, 5)
                tab_c3[(r, c3_through(A, i) % 2)] += 1
                tab_s2[(r, sc[i] % 4)] += 1
                tab_tri[(r, (sc[i] * (sc[i] - 1) // 2) % 2)] += 1
        print(f"C. odd n={n}: residue r_i = eps*w_i*(-1)^s_i mod 8 contingency:")
        print(f"   vs c3(i) mod 2:        {dict(sorted(tab_c3.items()))}")
        print(f"   vs s_i mod 4:          {dict(sorted(tab_s2.items()))}")
        print(f"   vs C(s_i,2) mod 2:     {dict(sorted(tab_tri.items()))}")

    # D. orthocenter parameters for orthocentric-but-not-regular classes
    print("D. orthocentric-not-regular classes: exact t_i (d2_ij = t_i + t_j):")
    for n in (5, 7):
        for cid, A in enumerate(tournaments(n)):
            M = Mmat(A)
            d2 = dist2_matrix(gram(M))
            spec = shape_spectrum(d2)
            if len(set(spec)) == 1:
                continue
            if not is_orthocentric(d2):
                continue
            t0 = Fraction(d2[0][1] + d2[0][2] - d2[1][2], 2)
            t = [t0] + [Fraction(d2[0][i]) - t0 for i in range(1, n)]
            print(f"   n={n} class {cid}: scores {tuple(scores(A))} det {det_int(M)}"
                  f" t = {[str(x) for x in t]}")

    # E. are the det-512 classes at n=7 exactly the switching class of QR_7?
    print("E. switching class of QR_7 vs det=512 classes at n=7:")
    QR = {1, 2, 4}
    Aq = [[1 if ((j - i) % 7) in QR else 0 for j in range(7)] for i in range(7)]
    Sq = skew(Aq)
    sw_canons = set()
    for mask in range(64):
        d = [1] + [(-1 if (mask >> (i - 1)) & 1 else 1) for i in range(1, 7)]
        Ssw = [[d[i] * Sq[i][j] * d[j] for j in range(7)] for i in range(7)]
        Asw = [[1 if Ssw[i][j] == 1 else 0 for j in range(7)] for i in range(7)]
        sw_canons.add(canon_form(Asw))
    det512 = []
    for cid, A in enumerate(tournaments(7)):
        if det_int(Mmat(A)) == 512:
            det512.append((cid, canon_form(A)))
    match = {c for _, c in det512} == sw_canons
    print(f"   switching orbit of QR_7 hits {len(sw_canons)} iso classes;"
          f" det-512 classes: {[c for c, _ in det512]};"
          f" sets EQUAL: {match}")

    # F. equifacetal-but-not-regular classes: distance matrices
    print("F. equifacetal-not-regular classes (Hamming distance matrices d2/4):")
    for n, targets in ((5, None), (6, None), (7, None)):
        for cid, A in enumerate(tournaments(n)):
            M = Mmat(A)
            d2 = dist2_matrix(gram(M))
            spec = shape_spectrum(d2)
            if len(set(spec)) == 1:
                continue
            if not is_equifacetal(d2):
                continue
            print(f"   n={n} class {cid} scores {tuple(scores(A))} det {det_int(M)}:")
            for row in d2:
                print(f"      {[x // 4 for x in row]}")

    # G. harden the mod-4 laws at larger n via random sampling
    import random
    random.seed(20260610)
    print("G. mod-4 laws on RANDOM tournaments at larger n (200 samples each):")
    for n in (8, 9, 10, 11):
        bad = 0
        for _ in range(200):
            A = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    if random.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
            S = skew(A); sc = scores(A)
            if n % 2 == 1:
                adj = adjugate(S)
                for i in range(n):
                    for j in range(n):
                        if (adj[i][j] - (-1) ** (sc[i] + sc[j])) % 4 != 0:
                            bad += 1
            else:
                adj = adjugate(S)
                u = [sum(adj[i][j] for j in range(n)) for i in range(n)]
                for i in range(n):
                    if (u[i] - (-1) ** sc[i]) % 4 != 0:
                        bad += 1
        law = ("adj_ij ≡ (-1)^(s_i+s_j) mod 4" if n % 2 == 1
               else "u_i = (adj(S)·1)_i ≡ (-1)^(s_i) mod 4, sign forced +")
        print(f"   n={n}: {law}: violations = {bad}/200 samples")

def main():
    print("THE TOURNAMENT SIMPLEX — macmini-2026-06-10-S2")
    print("M = I + S rows = cube vertices; vol = det(M)/n!; Gram = I + SS^T")
    print("d2_ij = 4(1 + disagree(i,j)) = 4(n-1) - 4(cod+ + cod-) = 4*Hamming(row_i,row_j)")

    qr7_check()

    all_recs = {}
    for n in (4, 5, 6, 7):
        all_recs[n] = analyze(n, full_table=(n <= 5))

    addendum()

    banner("SUMMARY OF FINDINGS (exploratory, no canon claims)")
    print("""
1. SWITCHING GEOMETRY (verified exhaustively n=4..7):
   Switching at vertex set X maps M -> DMD whose i-th row is d_i * R(row_i),
   where R reflects the cube coordinates in X. I.e. the simplex is reflected
   (a cube isometry) and then each SWITCHED vertex is replaced by its antipode.
   Reflection preserves everything; the per-vertex antipodal moves preserve
   det (volume) but change pairwise distances: VOLUME is the switching-class
   invariant, SHAPE is not.

2. DISTANCE DICTIONARY (verified exhaustively n=4..7):
   |row_i - row_j|^2 = 4 * (1 + #{k : exactly one of i,j beats k})
                     = 4(n-1) - 4*(common out-neighbors + common in-neighbors)
                     = 4 * Hamming(row_i, row_j).
   Simplex geometry IS co-degree combinatorics, scaled by 4.

3. REGULARITY: cone simplex orthocentric <=> row simplex regular <=>
   SS^T has constant off-diagonal <=> (odd n) DRT switching class
   [edge sqrt(2n+2)], (even n) skew-conference [edge sqrt(2n)].

4. ODD MAPS: every tournament's x -> Sx is an odd tangent field.
   Odd n: zeros on S^(n-1) are exactly ±w (Pfaffian vector) — the hairy-ball
   singularity HAS a name. On the sum-zero sphere the field never vanishes
   because sum(w) = sqrt(det(I+2A)) is odd.
   Even n: ker S = 0, but the projected field on the sum-zero sphere S^(n-2)
   (even-dim) vanishes exactly at ±adj(S)·1 — a canonical 'harmonic direction'
   every even tournament carries.
""")

if __name__ == '__main__':
    main()
