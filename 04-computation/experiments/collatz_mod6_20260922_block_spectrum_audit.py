#!/usr/bin/env python3
"""
Lane block_spectrum_audit (wave 2026-09-22, session collatz-mod6-20260917).

Audits the pasted "Adjunction-Bipartite Tournament Operator"
    M_B = [[A,0,1,0],[0,B,0,1],[K_{m,n},0,0,1],[0,K_{n,m},0,0]]
and its three spectral claims (nilpotent core, purely imaginary pairs from the
bipartite blocks, Re(lambda) < 0 for every trajectory because of a sink).

Sections of the output:
  S1  well-typing of the literal block layout (row/column count)
  S2  the consistent (m+n+2)-vertex reading and the sweep over all tournament
      iso classes m,n <= 4 and three cross-link orientations
  S3  claim (i): nilpotency of tournament blocks
  S4  claim (ii): spectra of one-directional / two-directional / skew
      bipartite blocks
  S5  claim (iii): Perron-Frobenius, sink/source spectrum shift
  S6  Brauer-Gentry bounds on all labelled tournaments n <= 6, regular
      tournaments, Paley T_7
  S7  K_5 / K_{3,3} in the underlying graph (Euler bound + explicit witness)
  S8  verdict table

Deterministic; no timing output. RAM well under 1 GB; runtime ~1-2 min.
"""
import itertools
import math
import sys

import numpy as np
import sympy as sp

np.set_printoptions(linewidth=200)


def fail(msg):
    raise RuntimeError(msg)


def out(*a):
    print(*a)
    sys.stdout.flush()


# --------------------------------------------------------------------------
# tournaments
# --------------------------------------------------------------------------
def all_labelled_tournaments(n):
    pairs = list(itertools.combinations(range(n), 2))
    for bits in itertools.product((0, 1), repeat=len(pairs)):
        A = np.zeros((n, n), dtype=np.int64)
        for (i, j), b in zip(pairs, bits):
            if b:
                A[i, j] = 1
            else:
                A[j, i] = 1
        yield A


def canon(A):
    n = A.shape[0]
    best = None
    for p in itertools.permutations(range(n)):
        key = tuple(int(A[p[i], p[j]]) for i in range(n) for j in range(n))
        if best is None or key < best:
            best = key
    return best


def iso_classes(n):
    seen = {}
    for A in all_labelled_tournaments(n):
        c = canon(A)
        if c not in seen:
            seen[c] = A.copy()
    return list(seen.values())


def is_transitive(A):
    n = A.shape[0]
    s = sorted(int(x) for x in A.sum(axis=1))
    return s == list(range(n))


def transitive(n):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            A[i, j] = 1
    return A


def charpoly_str(M):
    x = sp.Symbol('x')
    P = sp.Matrix(M.tolist()).charpoly(x)
    return sp.factor(P.as_expr())


def spectrum(M):
    ev = np.linalg.eigvals(M.astype(float))
    ev = np.round(ev, 6) + 0.0
    return sorted(ev, key=lambda z: (-z.real, -abs(z.imag), z.imag))


def fmt(z):
    r, i = float(z.real), float(z.imag)
    if abs(i) < 1e-9:
        return "%.4f" % r
    return "%.4f%+.4fi" % (r, i)


# --------------------------------------------------------------------------
out("=" * 78)
out("S1  Well-typing of the literal block layout")
out("=" * 78)
out("Literal layout: block rows (A | 0 | 1 | 0), (0 | B | 0 | 1), (K_{m,n} | 0 | 0 | 1), (0 | K_{n,m} | 0 | 0)")
out("Block row heights: m (A), n (B), m (K_{m,n} is m x n), n (K_{n,m} is n x m)  -> rows = 2m + 2n")
out("Block column widths: m (A above K_{m,n}), n (B above K_{n,m}), 1, 1          -> cols = m + n + 2")
out("Column 1 carries A (width m) and K_{m,n} (width n): consistent only if n = m.")
out("Column 2 carries B (width n) and K_{n,m} (width m): consistent only if m = n.")
out("Square iff 2m+2n = m+n+2 iff m+n = 2 iff m = n = 1.")
for m in range(1, 5):
    for n in range(1, 5):
        out("  m=%d n=%d : rows=%d cols=%d square=%s col-widths-consistent=%s"
            % (m, n, 2 * m + 2 * n, m + n + 2, 2 * m + 2 * n == m + n + 2, m == n))
out("Only square literal instance: m=n=1 with A=B=[0] (the 1-vertex tournament), K=[k], K'=[k'] with k,k' in {0,1}:")
x = sp.Symbol('x')
for k in (0, 1):
    for kp in (0, 1):
        M = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [k, 0, 0, 1], [0, kp, 0, 0]])
        P = M.charpoly(x).as_expr()
        out("  k=%d k'=%d : charpoly = %s ; eigenvalues = %s ; max real eigenvalue = %s" % (k, kp, sp.factor(P), sp.roots(P), max(r for r in sp.roots(P) if r.is_real)))

# --------------------------------------------------------------------------
out()
out("=" * 78)
out("S2  Consistent (m+n+2)-vertex reading and the exhaustive sweep")
out("=" * 78)
out("Vertex order (V_A, V_B, alpha, beta).  M = [[A, K, 0, 0], [K', B, 0, 1_n], [1_m^T, 0, 0, 0], [0, 0, 0, 0]].")
out("alpha = source (arcs alpha -> every vertex of V_A, no in-arcs); beta = sink (arcs every vertex of V_B -> beta, no out-arcs).")
out("Cross-link orientations: ONE = K = J_{m,n}, K' = 0 ; TWO = K = J, K' = J^T (K + K^T both ways) ; ORI = K_ij = [i+j even], K' = J^T - K^T (an oriented bipartite cross-link, union is a tournament on m+n).")


def build(A, B, mode):
    m, n = A.shape[0], B.shape[0]
    N = m + n + 2
    M = np.zeros((N, N), dtype=np.int64)
    M[:m, :m] = A
    M[m:m + n, m:m + n] = B
    J = np.ones((m, n), dtype=np.int64)
    if mode == "ONE":
        K, Kp = J, np.zeros((n, m), dtype=np.int64)
    elif mode == "TWO":
        K, Kp = J, J.T.copy()
    elif mode == "ORI":
        K = np.array([[1 if (i + j) % 2 == 0 else 0 for j in range(n)] for i in range(m)], dtype=np.int64)
        Kp = (J - K).T.copy()
    else:
        fail("mode")
    M[:m, m:m + n] = K
    M[m:m + n, :m] = Kp
    M[m + n, :m] = 1          # alpha -> V_A
    M[m:m + n, m + n + 1] = 1  # V_B -> beta
    return M


classes = {n: iso_classes(n) for n in range(1, 5)}
for n in range(1, 5):
    out("tournament iso classes on %d vertices: %d" % (n, len(classes[n])))
if [len(classes[n]) for n in range(1, 5)] != [1, 1, 2, 4]:
    fail("iso class counts")

rows = []
count = {"ONE": 0, "TWO": 0, "ORI": 0}
stats = {}
for mode in ("ONE", "TWO", "ORI"):
    st = dict(total=0, all_zero=0, has_nonzero_pure_imag=0, has_negative_real_eig=0,
              all_re_negative=0, perron_positive=0, block_triangular_factor_ok=0,
              zero_mult_ge2=0, max_re_ge0=0)
    for m in range(1, 5):
        for n in range(1, 5):
            for ia, A in enumerate(classes[m]):
                for ib, B in enumerate(classes[n]):
                    M = build(A, B, mode)
                    ev = spectrum(M)
                    st["total"] += 1
                    zero_mult = sum(1 for z in ev if abs(z) < 1e-9)
                    if zero_mult == len(ev):
                        st["all_zero"] += 1
                    if any(abs(z.real) < 1e-9 and abs(z.imag) > 1e-9 for z in ev):
                        st["has_nonzero_pure_imag"] += 1
                    if any(z.real < -1e-9 for z in ev):
                        st["has_negative_real_eig"] += 1
                    if all(z.real < -1e-9 for z in ev):
                        st["all_re_negative"] += 1
                    rho = max(abs(z) for z in ev)
                    if rho > 1e-9:
                        st["perron_positive"] += 1
                    if zero_mult >= 2:
                        st["zero_mult_ge2"] += 1
                    if max(z.real for z in ev) >= -1e-9:
                        st["max_re_ge0"] += 1
                    # exact charpoly and block-triangular factorisation check
                    PM = sp.Matrix(M.tolist()).charpoly(x).as_expr()
                    PA = sp.Matrix(A.tolist()).charpoly(x).as_expr()
                    PB = sp.Matrix(B.tolist()).charpoly(x).as_expr()
                    factor_ok = sp.expand(PM - x ** 2 * PA * PB) == 0
                    if factor_ok:
                        st["block_triangular_factor_ok"] += 1
                    # exact multiplicity of the zero eigenvalue
                    zm_exact = 0
                    Pp = sp.Poly(PM, x)
                    while Pp.eval(0) == 0:
                        Pp = sp.Poly(sp.quo(Pp.as_expr(), x), x)
                        zm_exact += 1
                    rows.append((mode, m, n, ia, ib, is_transitive(A), is_transitive(B), zm_exact,
                                 sp.factor(PM), ev))
    stats[mode] = st

out()
out("Sweep statistics (64 tournament pairs per orientation, matrices of size m+n+2 <= 10):")
keys = ["total", "all_zero", "zero_mult_ge2", "has_nonzero_pure_imag", "has_negative_real_eig",
        "all_re_negative", "max_re_ge0", "perron_positive", "block_triangular_factor_ok"]
out("| orientation | " + " | ".join(keys) + " |")
for mode in ("ONE", "TWO", "ORI"):
    out("| %s | " % mode + " | ".join(str(stats[mode][k]) for k in keys) + " |")

if stats["ONE"]["block_triangular_factor_ok"] != 64:
    fail("ONE must factor as x^2 charpoly(A) charpoly(B)")
if stats["ONE"]["all_re_negative"] or stats["TWO"]["all_re_negative"] or stats["ORI"]["all_re_negative"]:
    fail("no 0/1 matrix can have all eigenvalues with negative real part")
for mode in stats:
    if stats[mode]["max_re_ge0"] != 64:
        fail("Perron root must be >= 0")

out()
out("Exact characteristic polynomials, all 64 pairs, orientation ONE (block triangular):")
out("| m | n | A transitive | B transitive | zero mult | charpoly(M) | numerical spectrum |")
for r in rows:
    if r[0] == "ONE":
        out("| %d | %d | %s | %s | %d | %s | %s |" % (r[1], r[2], r[5], r[6], r[7], r[8], ", ".join(fmt(z) for z in r[9])))
out()
out("Exact characteristic polynomials, all 64 pairs, orientation TWO (K + K^T):")
out("| m | n | A transitive | B transitive | zero mult | charpoly(M) | numerical spectrum |")
for r in rows:
    if r[0] == "TWO":
        out("| %d | %d | %s | %s | %d | %s | %s |" % (r[1], r[2], r[5], r[6], r[7], r[8], ", ".join(fmt(z) for z in r[9])))
out()
out("Exact characteristic polynomials, all 64 pairs, orientation ORI (oriented cross-link):")
out("| m | n | A transitive | B transitive | zero mult | charpoly(M) | numerical spectrum |")
for r in rows:
    if r[0] == "ORI":
        out("| %d | %d | %s | %s | %d | %s | %s |" % (r[1], r[2], r[5], r[6], r[7], r[8], ", ".join(fmt(z) for z in r[9])))

# zero multiplicity in ONE equals 2 + (m if A transitive) + (n if B transitive) + zero-mults of cyclic blocks
out()
out("Zero-eigenvalue multiplicity in orientation ONE versus the block formula 2 + mult_0(A) + mult_0(B):")
bad = 0
for r in rows:
    if r[0] != "ONE":
        continue
    A, B = classes[r[1]][r[3]], classes[r[2]][r[4]]
    zA = sum(1 for z in spectrum(A) if abs(z) < 1e-9)
    zB = sum(1 for z in spectrum(B) if abs(z) < 1e-9)
    if r[7] != 2 + zA + zB:
        bad += 1
out("  mismatches: %d" % bad)
if bad:
    fail("zero multiplicity formula")
max_zero_one = max(r[7] for r in rows if r[0] == "ONE")
min_zero_one = min(r[7] for r in rows if r[0] == "ONE")
out("  ONE: zero multiplicity ranges from %d (both blocks 3-cycle-rich) to %d (both transitive, m=n=4: the whole M is nilpotent)" % (min_zero_one, max_zero_one))
out("  Every ONE matrix with both blocks transitive is nilpotent: %s" %
    all(r[7] == r[1] + r[2] + 2 for r in rows if r[0] == "ONE" and r[5] and r[6]))
out("  Every ONE matrix with a non-transitive block has Perron root > 0: %s" %
    all(max(abs(z) for z in r[9]) > 1e-9 for r in rows if r[0] == "ONE" and not (r[5] and r[6])))

# --------------------------------------------------------------------------
out()
out("=" * 78)
out("S3  Claim (i): nilpotent core")
out("=" * 78)
out("Transitive tournament on n vertices sorted by score: strictly upper triangular, A^n = 0.")
for n in range(1, 7):
    A = transitive(n)
    P = np.linalg.matrix_power(A, n)
    out("  n=%d: A^n == 0: %s ; charpoly = %s" % (n, bool((P == 0).all()), charpoly_str(A)))
out("Non-transitive tournaments: tr A = tr A^2 = 0 (THM-1858) so sum lambda = sum lambda^2 = 0; nilpotent iff acyclic iff transitive.")
out("Perron root of every tournament iso class n<=4:")
for n in range(1, 5):
    for A in classes[n]:
        ev = spectrum(A)
        out("  n=%d scores=%s transitive=%s charpoly=%s rho=%.4f spectrum=%s"
            % (n, sorted(int(s) for s in A.sum(axis=1)), is_transitive(A), charpoly_str(A),
               max(abs(z) for z in ev), ", ".join(fmt(z) for z in ev)))

# --------------------------------------------------------------------------
out()
out("=" * 78)
out("S4  Claim (ii): bipartite blocks")
out("=" * 78)
for (m, n) in [(2, 2), (2, 3), (3, 3), (3, 4), (4, 4)]:
    J = np.ones((m, n), dtype=np.int64)
    Z1 = np.zeros((m, m), dtype=np.int64)
    Z2 = np.zeros((n, n), dtype=np.int64)
    one = np.block([[Z1, J], [np.zeros((n, m), dtype=np.int64), Z2]])
    two = np.block([[Z1, J], [J.T, Z2]])
    skew = np.block([[Z1, J], [-J.T, Z2]])
    sv = np.linalg.svd(J.astype(float), compute_uv=False)
    out("  m=%d n=%d K=J: ONE [[0,K],[0,0]] square is zero: %s, charpoly %s, spectrum %s"
        % (m, n, bool((one @ one == 0).all()), charpoly_str(one), ", ".join(fmt(z) for z in spectrum(one))))
    out("           TWO [[0,K],[K^T,0]] symmetric, charpoly %s, spectrum %s, singular values of K %s, sqrt(mn)=%.4f"
        % (charpoly_str(two), ", ".join(fmt(z) for z in spectrum(two)), np.round(sv, 4).tolist(), math.sqrt(m * n)))
    out("           SKEW [[0,K],[-K^T,0]] (not a 0/1 matrix), charpoly %s, spectrum %s"
        % (charpoly_str(skew), ", ".join(fmt(z) for z in spectrum(skew))))
# general oriented K
out("  Oriented K (ORI pattern) alone, m=n=3: ONE nilpotent %s ; TWO spectrum %s ; SKEW spectrum %s" % (
    bool((np.block([[np.zeros((3, 3), int), np.array([[1 if (i + j) % 2 == 0 else 0 for j in range(3)] for i in range(3)])], [np.zeros((3, 3), int), np.zeros((3, 3), int)]]) ** 2 == 0).all()),
    ", ".join(fmt(z) for z in spectrum(np.block([[np.zeros((3, 3), int), np.array([[1 if (i + j) % 2 == 0 else 0 for j in range(3)] for i in range(3)])], [np.array([[1 if (i + j) % 2 == 0 else 0 for j in range(3)] for i in range(3)]).T, np.zeros((3, 3), int)]]))),
    ", ".join(fmt(z) for z in spectrum(np.block([[np.zeros((3, 3), int), np.array([[1 if (i + j) % 2 == 0 else 0 for j in range(3)] for i in range(3)])], [-np.array([[1 if (i + j) % 2 == 0 else 0 for j in range(3)] for i in range(3)]).T, np.zeros((3, 3), int)]])))))
out("Purely imaginary nonzero eigenvalues in the full sweep (S2 table column has_nonzero_pure_imag): ONE %d, TWO %d, ORI %d of 64 each"
    % (stats["ONE"]["has_nonzero_pure_imag"], stats["TWO"]["has_nonzero_pure_imag"], stats["ORI"]["has_nonzero_pure_imag"]))
pure_im_examples = [r for r in rows if any(abs(z.real) < 1e-9 and abs(z.imag) > 1e-9 for z in r[9])]
out("Examples with a purely imaginary pair (from the tournament blocks, not from K):")
out("  (none)" if not pure_im_examples else "  %d matrices" % len(pure_im_examples))
for r in pure_im_examples[:8]:
    out("  %s m=%d n=%d A-class %d B-class %d charpoly %s spectrum %s" % (r[0], r[1], r[2], r[3], r[4], r[8], ", ".join(fmt(z) for z in r[9])))
# which single tournaments n<=4 have purely imaginary eigenvalues
out("Tournament iso classes n<=4 with a purely imaginary nonzero eigenvalue (none expected: sum lambda^2 = 0 with a Perron root rho>0 needs Re-negative partners, and n<=4 spectra are listed in S3):")
pi_found = 0
for n in range(1, 5):
    for A in classes[n]:
        ev = spectrum(A)
        if any(abs(z.real) < 1e-9 and abs(z.imag) > 1e-9 for z in ev):
            pi_found += 1
            out("  n=%d scores=%s charpoly=%s spectrum=%s" % (n, sorted(int(s) for s in A.sum(axis=1)), charpoly_str(A), ", ".join(fmt(z) for z in ev)))
out("  found: %d" % pi_found)

# --------------------------------------------------------------------------
out()
out("=" * 78)
out("S5  Claim (iii): Perron-Frobenius and the sink/source shift")
out("=" * 78)
out("Every square 0/1 matrix M is nonnegative, so rho(M) is an eigenvalue (Perron-Frobenius); Re(lambda)<0 for all lambda is impossible.")
out("Sweep: matrices whose eigenvalues all have Re<0: ONE %d, TWO %d, ORI %d (of 64 each); matrices with max Re >= 0: %d, %d, %d."
    % (stats["ONE"]["all_re_negative"], stats["TWO"]["all_re_negative"], stats["ORI"]["all_re_negative"],
       stats["ONE"]["max_re_ge0"], stats["TWO"]["max_re_ge0"], stats["ORI"]["max_re_ge0"]))
out("Sink/source shift: appending a sink (zero row) or a source (zero column) multiplies the characteristic polynomial by x.")
shift_ok = 0
shift_tot = 0
for n in range(1, 5):
    for A in classes[n]:
        for extra in ("sink", "source", "both", "both+arc"):
            N = n + (2 if extra.startswith("both") else 1)
            M = np.zeros((N, N), dtype=np.int64)
            M[:n, :n] = A
            if extra == "sink":
                M[:n, n] = 1
            elif extra == "source":
                M[n, :n] = 1
            else:
                M[n, :n] = 1        # source alpha -> all
                M[:n, n + 1] = 1    # all -> sink beta
                if extra == "both+arc":
                    M[n, n + 1] = 1  # alpha -> beta
            PM = sp.Matrix(M.tolist()).charpoly(x).as_expr()
            PA = sp.Matrix(A.tolist()).charpoly(x).as_expr()
            k = 2 if extra.startswith("both") else 1
            shift_tot += 1
            if sp.expand(PM - x ** k * PA) == 0:
                shift_ok += 1
out("  charpoly(M) == x^k charpoly(A) (k = number of appended polar vertices): %d of %d cases" % (shift_ok, shift_tot))
if shift_ok != shift_tot:
    fail("sink/source shift")
out("  Negative real eigenvalues do occur (from cyclic tournament blocks), but never for every eigenvalue: sweep count has_negative_real_eig ONE %d TWO %d ORI %d."
    % (stats["ONE"]["has_negative_real_eig"], stats["TWO"]["has_negative_real_eig"], stats["ORI"]["has_negative_real_eig"]))

out("Integer dichotomy: a nonnegative integer matrix has rho = 0 (nilpotent) or rho >= 1, never 0 < rho < 1 (tr M^k is a nonnegative integer; a closed walk gives tr M^k >= 1 for infinitely many k, and tr M^k <= N rho^k).")
nz = [max(abs(z) for z in r[9]) for r in rows if max(abs(z) for z in r[9]) > 1e-9]
out("  sweep: %d of %d matrices have rho > 0; minimum nonzero rho = %.4f; maximum rho = %.4f (TWO, m=n=4)" % (len(nz), len(rows), min(nz), max(nz)))
if min(nz) < 1 - 1e-9:
    fail("integer dichotomy")
out("  Discrete-time contraction x -> M x needs rho(M) < 1, not Re(lambda) < 0; for a 0/1 matrix that forces nilpotency, i.e. an acyclic digraph.")

# --------------------------------------------------------------------------
out()
out("=" * 78)
out("S6  Brauer-Gentry bounds on all labelled tournaments n<=6; regular tournaments; Paley T_7")
out("=" * 78)
out("Bounds: Re(lambda) >= -1/2 ; |lambda| <= (n-1)/2 ; |Im(lambda)| <= (1/2) cot(pi/(2n)).")
nontrans_total = 0
for n in range(1, 7):
    cnt = 0
    min_re = 10.0
    max_mod = 0.0
    max_im = 0.0
    viol = [0, 0, 0]
    trans = 0
    zero_rho_nontrans = 0
    max_rho = 0.0
    for A in all_labelled_tournaments(n):
        cnt += 1
        ev = np.linalg.eigvals(A.astype(float))
        tr = is_transitive(A)
        trans += tr
        rho = max(abs(ev))
        if not tr and rho < 1e-9:
            zero_rho_nontrans += 1
        max_rho = max(max_rho, rho)
        min_re = min(min_re, ev.real.min())
        max_mod = max(max_mod, rho)
        max_im = max(max_im, abs(ev.imag).max())
        if ev.real.min() < -0.5 - 1e-9:
            viol[0] += 1
        if rho > (n - 1) / 2 + 1e-9:
            viol[1] += 1
        if n >= 2 and abs(ev.imag).max() > 0.5 / math.tan(math.pi / (2 * n)) + 1e-9:
            viol[2] += 1
    nontrans_total += cnt - trans
    cot = 0.5 / math.tan(math.pi / (2 * n)) if n >= 2 else 0.0
    out("  n=%d labelled=%d transitive=%d nontransitive-with-rho=0: %d | min Re=%.4f (bound -0.5) | max |lambda|=%.4f (bound %.4f) | max |Im|=%.4f (bound %.4f) | violations %s"
        % (n, cnt, trans, zero_rho_nontrans, min_re, max_mod, (n - 1) / 2, max_im, cot, viol))
    if any(viol) or zero_rho_nontrans:
        fail("Brauer-Gentry violation or non-transitive nilpotent tournament")
out("Proof of Re(lambda) >= -1/2: A + A^T = J - I; for a unit eigenvector v, 2 Re(lambda) = v*(A+A^T)v = |sum v_i|^2 - 1 >= -1.")
out("Proof of |lambda| <= (n-1)/2: rho is the Perron root with a nonnegative unit Perron vector p, rho = p^T A p = (p^T(A+A^T)p)/2 = ((sum p_i)^2-1)/2 <= (n-1)/2 by Cauchy-Schwarz, and |lambda| <= rho for every lambda.")
out("Regular tournaments (all scores (n-1)/2): A is normal (A commutes with A^T = J - I - A), so every non-Perron eigenvector is orthogonal to 1 and has Re(lambda) = -1/2 exactly.")
for n in (3, 5, 7):
    # circulant regular tournaments: connection set S with S ∪ -S = Z_n^* , S ∩ -S = ∅
    seen = 0
    for S in itertools.combinations(range(1, n), (n - 1) // 2):
        if all((-s) % n not in S for s in S):
            A = np.array([[1 if (j - i) % n in S else 0 for j in range(n)] for i in range(n)], dtype=np.int64)
            ev = spectrum(A)
            non_perron = [z for z in ev if abs(z - (n - 1) / 2) > 1e-9]
            ok = all(abs(z.real + 0.5) < 1e-9 for z in non_perron)
            if not ok:
                fail("regular tournament real part")
            seen += 1
            if seen <= 3:
                out("  n=%d circulant S=%s spectrum %s ; all non-Perron Re = -1/2: %s" % (n, S, ", ".join(fmt(z) for z in ev), ok))
    out("  n=%d circulant regular tournaments checked: %d" % (n, seen))
# Paley T_7
S = {1, 2, 4}
P7 = np.array([[1 if (j - i) % 7 in S else 0 for j in range(7)] for i in range(7)], dtype=np.int64)
arcs = int(P7.sum())
cyc = 0
trans3 = 0
for a, b, c in itertools.combinations(range(7), 3):
    sub = P7[np.ix_([a, b, c], [a, b, c])]
    if sorted(sub.sum(axis=1).tolist()) == [1, 1, 1]:
        cyc += 1
    else:
        trans3 += 1
out("  Paley T_7 (S={1,2,4}): arcs %d, cyclic triples %d = (7^3-7)/24 = %d, transitive triples %d, charpoly %s, spectrum %s"
    % (arcs, cyc, (343 - 7) // 24, trans3, charpoly_str(P7), ", ".join(fmt(z) for z in spectrum(P7))))
if (arcs, cyc, trans3) != (21, 14, 21):
    fail("Paley counts")
# Hamiltonian paths of Paley T_7 (for the paste's '21 forbidden Hamiltonian path count')
def ham_paths(A):
    n = A.shape[0]
    cnt = 0
    for p in itertools.permutations(range(n)):
        if all(A[p[i], p[i + 1]] for i in range(n - 1)):
            cnt += 1
    return cnt
out("  Paley T_7 Hamiltonian paths: %d ; directed 3-cycles: %d (the paste says 21 directed 3-cycles)" % (ham_paths(P7), cyc))
h2 = 0.81057 if False else 8 / math.pi ** 2
H = -(h2 * math.log2(h2) + (1 - h2) * math.log2(1 - h2))
out("  Binary entropy of 8/pi^2 = %.5f is %.5f bits (paste: 0.704)" % (h2, H))

# --------------------------------------------------------------------------
out()
out("=" * 78)
out("S7  K_5 and K_{3,3} in the underlying graph")
out("=" * 78)
out("Underlying simple graph of M (any orientation, all cross pairs linked): complete graph K_{m+n} on V_A u V_B, plus alpha ~ V_A and beta ~ V_B.")
for m in range(1, 5):
    for n in range(1, 5):
        v = m + n + 2
        e = m * (m - 1) // 2 + n * (n - 1) // 2 + m * n + m + n
        euler_fail = e > 3 * v - 6 if v >= 3 else False
        k33 = m >= 3 and n >= 3
        k5 = m + n >= 5
        out("  m=%d n=%d : V=%d E=%d 3V-6=%d Euler-violated=%s K_{3,3}-subgraph=%s K_5-subgraph=%s"
            % (m, n, v, e, 3 * v - 6, euler_fail, k33, k5))
out("  Witness K_{3,3} for m=n=3: parts {a1,a2,a3} and {b1,b2,b3}, all nine cross pairs are arcs of K or K' (one of K_ij, K'_ji is 1 in ONE/TWO/ORI: %s)."
    % all((build(classes[3][0], classes[3][0], mode)[:3, 3:6] + build(classes[3][0], classes[3][0], mode)[3:6, :3].T > 0).all() for mode in ("ONE", "TWO", "ORI")))
out("  Witness K_5 for m+n=5: any five vertices of V_A u V_B are pairwise adjacent (tournament arcs inside a part, cross-link arcs across).")
out("  Non-planarity is therefore a property of m+n >= 5 alone; no vertex of M is an integer and no arc is a Collatz step: no map (SCOPE).")

# --------------------------------------------------------------------------
out()
out("=" * 78)
out("S8  Verdict table")
out("=" * 78)
verdicts = [
    ("literal block layout is a square matrix", "REFUTED unless m=n=1 (rows 2m+2n, cols m+n+2)"),
    ("(i) transitive block is nilpotent", "PROVED (strictly triangular after topological sort)"),
    ("(i) 'maximum algebraic multiplicity' of 0 for a general block", "REFUTED: non-transitive block has Perron root > 0 (all %d nontransitive labelled tournaments n<=6 have rho > 0)" % nontrans_total),
    ("(ii) bipartite blocks give purely imaginary pairs +-i gamma", "REFUTED: one-directional block is nilpotent, two-directional is symmetric (real +-sqrt(mn)); imaginary needs the skew block [[0,K],[-K^T,0]], not a 0/1 matrix"),
    ("(iii) Re(lambda)<0 for every trajectory", "REFUTED: Perron-Frobenius gives a real eigenvalue rho >= 0; 0 of 192 sweep matrices have all Re<0"),
    ("(iii) sink gives a negative eigenvalue", "REFUTED: sink/source multiply charpoly by x (extra zero eigenvalue), %d of %d cases" % (shift_ok, shift_tot)),
    ("Brauer-Gentry Re >= -1/2, |lambda| <= (n-1)/2", "PROVED here (one-line) and verified on all labelled tournaments n<=6"),
    ("|Im| <= (1/2)cot(pi/2n)", "UNCITED-RECOLLECTION, verified n<=6"),
    ("K_{3,3} minor from the cross-link", "PROVED for m,n>=3 (subgraph), content-free for Collatz"),
    ("3N+1 functor replaces odd nodes by 3-cycles", "SCOPE: no map from graph minors to integer orbits is defined"),
]
for a, b in verdicts:
    out("  %-62s | %s" % (a, b))
out()
out("Provenance: wave 2026-09-22, lane block_spectrum_audit; numpy %s, sympy %s." % (np.__version__, sp.__version__))
out("Citations: Brauer-Gentry, Bull. Amer. Math. Soc. 74 (1968) (Re >= -1/2, |lambda| <= (n-1)/2; the cot bound is UNCITED-RECOLLECTION from the same paper); Perron-Frobenius as in Horn-Johnson, Matrix Analysis, Thm 8.3.1; Kuratowski (classical).")
out("Inherited tokens not re-derived here: the 23 mod 256 reset cylinder (guards lane), B^3 = -I (THM-4139), tr A = tr A^2 = 0 (THM-1858), h-spectrum omits 7 and 21 (THM-1370).")
out("DONE")
