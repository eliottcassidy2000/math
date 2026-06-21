"""
conn_verify_hopfieldfubi_kps-Sx-wf.py
=====================================
ADVERSARIAL, INDEPENDENT verification of the 'hopfield-fubini-hebbian' connection.

Strategy: do NOT reuse the original script's THM-554 algebraic quadratic-form
construction. Instead:
  - compute c3 (=alpha_1) by pure brute force over EVERY tiling (engine from prompt)
  - INDEPENDENTLY decide whether c3 is exactly a quadratic (degree<=2) function of
    the tile bits, using MULTILINEAR / finite-difference Mobius coefficients over
    GF rationals. A boolean function f:{0,1}^F -> Z is degree<=2 iff ALL multilinear
    coefficients of order >=3 vanish. This is a basis-free test, no Hopfield assumed.
  - if degree<=2, EXTRACT the unique Hopfield couplings W[e,f] and compare ground
    state (max c3) to the regular tournament (THM-027).
  - INDEPENDENTLY test whether H (full Hamiltonian-path count = OCF) is degree<=2,
    and report exact max residual of best quadratic fit (least squares, exact rationals).
  - Re-run Hopfield greedy ascent basin fractions with fresh code.

All arithmetic EXACT (fractions.Fraction / Python ints).
"""
from itertools import combinations, product
from fractions import Fraction
import random

# ---------- engine (copied from prompt) ----------
def tiles(n):
    return [(a, b) for a in range(3, n + 1) for b in range(1, a - 1)]

def adj(n, bits, T):
    A = [[0] * (n + 1) for _ in range(n + 1)]
    for k in range(n, 1, -1):
        A[k][k - 1] = 1
    for (a, b), bit in zip(T, bits):
        if bit == 0:
            A[a][b] = 1
        else:
            A[b][a] = 1
    return A

def c3(A, n):
    t = 0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            for k in range(j + 1, n + 1):
                if (A[i][j] + A[i][k], A[j][i] + A[j][k], A[k][i] + A[k][j]) == (1, 1, 1):
                    t += 1
    return t

def scores(A, n):
    return tuple(sorted(sum(A[v][u] for u in range(1, n + 1)) for v in range(1, n + 1)))

def ham_paths(A, n):
    # count Hamiltonian paths (directed) = Redei H = OCF. DP over subsets.
    # dp[mask][v] = # directed paths covering 'mask' ending at v.
    full = (1 << n) - 1
    dp = [[0] * (n + 1) for _ in range(1 << n)]
    for v in range(1, n + 1):
        dp[1 << (v - 1)][v] = 1
    for mask in range(1 << n):
        for v in range(1, n + 1):
            if not (mask & (1 << (v - 1))):
                continue
            cur = dp[mask][v]
            if cur == 0:
                continue
            for u in range(1, n + 1):
                if mask & (1 << (u - 1)):
                    continue
                if A[v][u] == 1:  # arc v->u, extend path
                    dp[mask | (1 << (u - 1))][u] += cur
    return sum(dp[full][v] for v in range(1, n + 1))

# ---------- independent multilinear (Mobius) degree test ----------
def multilinear_coeffs(values, F):
    """
    values: dict bitvector(tuple of 0/1 length F) -> int.
    Return dict subset(frozenset of indices) -> coefficient c_S such that
    f(b) = sum_S c_S * prod_{i in S} b_i  (multilinear extension on the cube).
    c_S = sum_{T subseteq S} (-1)^{|S|-|T|} f(1_T) where 1_T has bits set on T.
    Compute via the standard 'subset zeta/Mobius' over the cube. Only feasible
    for small F; we use F<=10 (n<=6) plus c3 up to n=7 with a smarter route.
    """
    coeffs = {}
    # iterate subsets S by size; use inclusion-exclusion restricted to S
    idxs = list(range(F))
    for r in range(F + 1):
        for S in combinations(idxs, r):
            Sset = set(S)
            tot = 0
            for r2 in range(r + 1):
                for Tt in combinations(S, r2):
                    b = [0] * F
                    for i in Tt:
                        b[i] = 1
                    sign = (-1) ** (r - r2)
                    tot += sign * values[tuple(b)]
            if tot != 0:
                coeffs[frozenset(S)] = tot
    return coeffs

def max_degree(coeffs):
    return max((len(S) for S in coeffs), default=0)

# ---------- exact best-quadratic least-squares residual ----------
def exact_quadratic_fit_residual(values, F):
    """
    Fit best quadratic q(b)=c0+sum lin_i b_i + sum_{i<j} q_ij b_i b_j to f over all
    2^F points, minimizing max|f-q| ... we instead report the EXACT max abs residual
    of the least-squares fit. Use exact rational normal equations via Gaussian elim.
    Features: 1, b_i (F), b_i b_j (C(F,2)).
    """
    feats = [()]  # constant
    for i in range(F):
        feats.append((i,))
    for i in range(F):
        for j in range(i + 1, F):
            feats.append((i, j))
    M = len(feats)
    pts = list(values.keys())
    # design rows
    def featval(b, f):
        p = 1
        for i in f:
            p *= b[i]
        return p
    # Normal equations A = X^T X (M x M), rhs = X^T y, exact
    A = [[Fraction(0)] * M for _ in range(M)]
    rhs = [Fraction(0)] * M
    for b in pts:
        row = [featval(b, f) for f in feats]
        y = values[b]
        for a in range(M):
            if row[a] == 0:
                continue
            for c in range(M):
                if row[c]:
                    A[a][c] += row[a] * row[c]
            rhs[a] += row[a] * y
    # solve A x = rhs via Gaussian elimination (A may be singular if features dependent;
    # use least-norm by pivoting). We'll do standard elimination, skipping zero pivots.
    aug = [A[i][:] + [rhs[i]] for i in range(M)]
    piv_cols = []
    r = 0
    for col in range(M):
        # find pivot
        pr = None
        for rr in range(r, M):
            if aug[rr][col] != 0:
                pr = rr
                break
        if pr is None:
            continue
        aug[r], aug[pr] = aug[pr], aug[r]
        pv = aug[r][col]
        aug[r] = [x / pv for x in aug[r]]
        for rr in range(M):
            if rr != r and aug[rr][col] != 0:
                fac = aug[rr][col]
                aug[rr] = [aug[rr][k] - fac * aug[r][k] for k in range(M + 1)]
        piv_cols.append(col)
        r += 1
        if r == M:
            break
    x = [Fraction(0)] * M
    for i, col in enumerate(piv_cols):
        x[col] = aug[i][M]
    # compute residuals
    maxres = Fraction(0)
    for b in pts:
        q = sum(x[a] * featval(b, feats[a]) for a in range(M))
        res = abs(Fraction(values[b]) - q)
        if res > maxres:
            maxres = res
    return maxres

# ---------- Hopfield greedy ascent on c3 ----------
def hopfield_basin(values, F, cmax, starts=300, seed=1):
    rng = random.Random(seed)
    reached = 0
    hist = {}
    for _ in range(starts):
        b = [rng.randint(0, 1) for _ in range(F)]
        improved = True
        while improved:
            improved = False
            cur = values[tuple(b)]
            order = list(range(F))
            rng.shuffle(order)
            for i in order:
                b[i] ^= 1
                nv = values[tuple(b)]
                if nv > cur:
                    cur = nv
                    improved = True
                else:
                    b[i] ^= 1
        fv = values[tuple(b)]
        hist[fv] = hist.get(fv, 0) + 1
        if fv == cmax:
            reached += 1
    return reached, starts, hist

# ---------- main ----------
def build_values(n, what="c3"):
    T = tiles(n)
    F = len(T)
    vals = {}
    for bits in product((0, 1), repeat=F):
        A = adj(n, bits, T)
        if what == "c3":
            vals[bits] = c3(A, n)
        elif what == "H":
            vals[bits] = ham_paths(A, n)
        elif what == "score":
            vals[bits] = scores(A, n)
    return F, vals

OUT = []
def p(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.append(s)

p("=" * 70)
p("INDEPENDENT ADVERSARIAL VERIFICATION: hopfield-fubini-hebbian")
p("=" * 70)

p("\n--- TEST A: is c3 EXACTLY degree<=2 in tile bits? (basis-free Mobius test) ---")
p("(no Hopfield assumed; we just compute multilinear coeffs of brute-force c3)")
for n in range(3, 8):
    T = tiles(n)
    F = len(T)
    if F > 15:
        p(f"  n={n}: F={F} too big to enumerate; skipping"); continue
    # enumerate (feasible up to F=15 => 32768)
    vals = {}
    for bits in product((0, 1), repeat=F):
        vals[bits] = c3(adj(n, bits, T), n)
    if F <= 10:
        coeffs = multilinear_coeffs(vals, F)
        deg = max_degree(coeffs)
        n_order3plus = sum(1 for S in coeffs if len(S) >= 3)
        p(f"  n={n}: F={F}, 2^F={2**F} tilings, max multilinear degree of c3 = {deg}"
          f"  ; #coeffs of order>=3 = {n_order3plus}  => {'QUADRATIC (deg<=2)' if deg<=2 else 'NOT quadratic'}")
    else:
        # F=15: full Mobius over all subsets is too big; instead check degree<=2
        # by verifying all 3rd finite differences vanish on a random sample of triples.
        rng = random.Random(7)
        idxs = list(range(F))
        bad = 0
        trials = 2000
        for _ in range(trials):
            i, j, k = rng.sample(idxs, 3)
            base = [rng.randint(0, 1) for _ in range(F)]
            # third mixed difference over coords i,j,k must be 0 for deg<=2
            s = 0
            for di, dj, dk in product((0, 1), repeat=3):
                b = base[:]
                b[i], b[j], b[k] = di, dj, dk
                sign = (-1) ** (di + dj + dk)
                s += sign * vals[tuple(b)]
            if s != 0:
                bad += 1
        p(f"  n={n}: F={F}, 2^F={2**F}; sampled {trials} 3rd-diff triples, nonzero (deg>=3 evidence)={bad}"
          f"  => {'QUADRATIC (deg<=2)' if bad==0 else 'NOT quadratic'}")

p("\n--- TEST B: c3 ground state (max) = regular tournament? (THM-027) ---")
for n in range(3, 8):
    T = tiles(n)
    F = len(T)
    if F > 15:
        continue
    best = -1
    bestscores = set()
    for bits in product((0, 1), repeat=F):
        A = adj(n, bits, T)
        cv = c3(A, n)
        if cv > best:
            best = cv
            bestscores = {scores(A, n)}
        elif cv == best:
            bestscores.add(scores(A, n))
    reg = tuple(sorted([(n - 1) // 2] * n)) if n % 2 == 1 else None
    p(f"  n={n}: max c3={best}; maximizer score multiset(s)={sorted(bestscores)}"
      + (f"  regular={reg} -> {'MATCH' if n%2==1 and bestscores=={reg} else ('only-regular' if n%2==1 and reg in bestscores else 'see set')}" if n%2==1 else "  (even n: near-regular)"))

p("\n--- TEST C: is H (full OCF Hamiltonian-path count) degree<=2? residual = magic gap ---")
for n in range(3, 7):
    T = tiles(n)
    F = len(T)
    vals = {}
    for bits in product((0, 1), repeat=F):
        vals[bits] = ham_paths(adj(n, bits, T), n)
    # sanity: H = 1 + 2*c3 + ... so for n=3 H should be small
    res = exact_quadratic_fit_residual(vals, F)
    # also exact degree via Mobius if small
    if F <= 10:
        coeffs = multilinear_coeffs(vals, F)
        deg = max_degree(coeffs)
    else:
        deg = "?"
    p(f"  n={n}: F={F}; best-quadratic max|H-Qfit| = {res} (={float(res)})  ; exact multilinear deg(H)={deg}")

p("\n--- TEST D: Hopfield greedy ascent basin fractions on c3 (fresh code) ---")
for n in range(4, 8):
    T = tiles(n)
    F = len(T)
    if F > 15:
        continue
    vals = {}
    cmax = -1
    for bits in product((0, 1), repeat=F):
        cv = c3(adj(n, bits, T), n)
        vals[bits] = cv
        cmax = max(cmax, cv)
    reached, starts, hist = hopfield_basin(vals, F, cmax, starts=300, seed=42)
    p(f"  n={n}: global c3-max={cmax}; greedy reached it {reached}/{starts} = {reached/starts:.3f}"
      f" ; final-c3 hist={dict(sorted(hist.items()))}")

p("\n--- TEST E: sanity check H = 1 + 2*c3 + 4*c5disjoint... at the regular/Paley point ---")
# verify OCF base identity H>=1+2*c3 numerically on a few tilings
for n in [5, 7]:
    T = tiles(n)
    F = len(T)
    rng = random.Random(3)
    ok = True
    for _ in range(200):
        bits = tuple(rng.randint(0,1) for _ in range(F))
        A = adj(n, bits, T)
        H = ham_paths(A, n)
        cc = c3(A, n)
        if H < 1 + 2*cc:  # OCF: H = 1+2alpha1+4alpha2+... all terms >=0
            ok = False
            break
    p(f"  n={n}: H >= 1+2*c3 on 200 random tilings: {'OK' if ok else 'VIOLATED'}")

# write output
with open("C:/Users/Eliott/Documents/GitHub/math/05-knowledge/results/conn_verify_hopfieldfubi_kps-Sx-wf.out", "w") as fh:
    fh.write("\n".join(OUT) + "\n")
p("\n[written to 05-knowledge/results/conn_verify_hopfieldfubi_kps-Sx-wf.out]")
