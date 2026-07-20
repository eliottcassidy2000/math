#!/usr/bin/env python3
"""
Seidel spectra, the odd-n even-graph projection, and why both are sin/cos
                                                        (mac-mini-2026-07-20-S130)
================================================================================
Owner: "think seidel spectra and even graph bijections that are odd n only, and how
these two interlaced can be seen in sin and cosine."

The claim under test is that the two "odd n only" facts are ONE parity fact wearing
two coats, and that the coat is literally the parity of sin and cos:

  SKEW SIDE (R).  S_ij = +1 if i->j, -1 if j->i.  S^T = -S, so
      p(-x) = (-1)^n p(x)   for p(x) = det(xI - S).
  So the characteristic polynomial is an EVEN function at even n and an ODD function
  at odd n -- exactly cos vs sin.  Hence at odd n it has the factor x: a FORCED ZERO
  eigenvalue, the analogue of sin(0) = 0, which cos does not have.
  And Cauchy interlacing under vertex deletion is the analogue of the fact that the
  zeros of sin and cos INTERLACE.

  GF(2) SIDE.  The repo's projection T_cycle = (I + L(K_n)) T mod 2 lands in the even
  graphs (all degrees even) only for odd n.  Computed here: the degree of the image at
  w is (n-1)*d_w mod 2, so it vanishes identically iff n is ODD.  The companion fact is
  that the BICYCLE space Cut & Cycle of K_n is zero iff n is odd, which is what makes
  Cut (+) Cycle a genuine direct sum there.

  CIRCULANTS make the sine literal: for a circulant tournament with connection set C,
      mu_j = 2i * sum_{k in C} sin(2*pi*j*k/n),
  a pure SINE sum, and the cosine part is pinned: Re(z_j) = -1/2 for every j != 0.
"""
import numpy as np
from itertools import permutations, combinations

np.set_printoptions(precision=4, suppress=True)

def scaffold(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    return pairs, {p: k for k, p in enumerate(pairs)}, len(pairs)

def skew(code, pairs, n):
    """S_ij = +1 if i->j.  bit=1 means the pair (i<j) is REVERSED, i.e. j->i."""
    S = np.zeros((n, n))
    for e, (i, j) in enumerate(pairs):
        if code >> e & 1: S[j, i], S[i, j] = 1.0, -1.0
        else:             S[i, j], S[j, i] = 1.0, -1.0
    return S

def spec(S):
    """eigenvalues of iS (Hermitian, so real), sorted."""
    return np.sort(np.linalg.eigvalsh(1j * S).real)

def canon_codes(n):
    """iso-class representative codes by vertex extension."""
    pairs, idx, E = scaffold(n)
    reps = {0}
    for k in range(2, n + 1):
        pk, ik, Ek = scaffold(k)
        op, _, _ = scaffold(k - 1)
        cand = []
        for r in reps:
            base = 0
            for e, (i, j) in enumerate(op):
                if r >> e & 1: base |= 1 << ik[(i, j)]
            newe = [ik[(i, k - 1)] for i in range(k - 1)]
            for mask in range(1 << (k - 1)):
                v = base
                for b, e in enumerate(newe):
                    if mask >> b & 1: v |= 1 << e
                cand.append(v)
        best = {}
        for v in cand:
            bc = None
            for p in permutations(range(k)):
                c = 0
                for e, (i, j) in enumerate(pk):
                    a, b = p[i], p[j]
                    t = ik[(min(a, b), max(a, b))]
                    bit = (v >> e & 1) ^ (1 if a > b else 0)
                    if bit: c |= 1 << t
                if bc is None or c < bc: bc = c
            best[bc] = True
        reps = set(best)
    return sorted(reps)

# ============================================================ PART A
print("=" * 76)
print("PART A -- the skew-Seidel characteristic polynomial has the PARITY OF sin/cos")
print("=" * 76)
print("  p(-x) = det(-xI - S) = (-1)^n det(xI + S) = (-1)^n det(xI - S^T) = (-1)^n p(x)")
print("  => p EVEN at even n (like cos), ODD at odd n (like sin), so odd n forces a ZERO.")
print()
print(f"{'n':>3} {'classes':>8} {'all eigs imaginary':>19} {'p parity = (-1)^n':>18} "
      f"{'0 in spec (odd n)':>18} {'rank always even':>17}")
for n in range(3, 8):
    pairs, idx, E = scaffold(n)
    reps = canon_codes(n)
    imag_ok = par_ok = zero_ok = rank_ok = True
    for r in reps:
        S = skew(r, pairs, n)
        ev = np.linalg.eigvals(S)
        if np.max(np.abs(ev.real)) > 1e-8: imag_ok = False
        c = np.poly(S).real                      # char poly coefficients, high->low
        # parity: coefficients of the "wrong" parity must vanish
        for k, co in enumerate(c):
            deg = len(c) - 1 - k
            if (deg % 2) != (n % 2) and abs(co) > 1e-7: par_ok = False
        s = spec(S)
        if n % 2 == 1 and np.min(np.abs(s)) > 1e-8: zero_ok = False
        if np.linalg.matrix_rank(S, tol=1e-8) % 2 != 0: rank_ok = False
    print(f"{n:>3} {len(reps):>8} {str(imag_ok):>19} {str(par_ok):>18} "
          f"{str(zero_ok) if n%2 else 'n/a':>18} {str(rank_ok):>17}")

# ============================================================ PART B
print()
print("=" * 76)
print("PART B -- the Seidel spectrum IS switching-invariant (unlike H, FAS, triangles)")
print("=" * 76)
print("  Switching at X sends S -> DSD with D = diag(+-1), and D = D^{-1}, so it is a")
print("  SIMILARITY: the spectrum cannot move.  Contrast THM-1420/THM-1415 addendum,")
print("  where H, min-FAS and the cyclic-triangle count are all NOT switching-invariant.")
print()
A049313 = {1: 1, 2: 1, 3: 1, 4: 2, 5: 2, 6: 6, 7: 12, 8: 79}
print(f"{'n':>3} {'iso classes':>12} {'distinct spectra':>17} {'switching classes':>18} "
      f"{'spectrum complete?':>19}")
for n in range(3, 8):
    pairs, idx, E = scaffold(n)
    reps = canon_codes(n)
    sp = set()
    swok = True
    for r in reps:
        S = skew(r, pairs, n)
        sp.add(tuple(np.round(spec(S), 7)))
        # verify switching invariance directly on a few sign patterns
        for m in range(min(1 << n, 32)):
            D = np.diag([-1.0 if m >> v & 1 else 1.0 for v in range(n)])
            if not np.allclose(np.sort(spec(D @ S @ D)), np.sort(spec(S)), atol=1e-7):
                swok = False
    ns = len(sp); sw = A049313.get(n)
    print(f"{n:>3} {len(reps):>12} {ns:>17} {str(sw):>18} "
          f"{('YES' if ns == sw else 'no -- coarser'):>19}")
    assert swok, "switching moved the spectrum -- impossible"
print("  (switching invariance verified numerically on every class at every n)")

# ============================================================ PART C
print()
print("=" * 76)
print("PART C -- WHY the even-graph projection is odd-n only:  two parity facts")
print("=" * 76)
print("  C1. The repo's map T_cycle = (I + L(K_n)) T mod 2.  Degree of the image at w:")
print("      sum_{e at w} [(I+L)T]_e = (n-1)*d_w + sum_x d_x = (n-1)*d_w  (mod 2),")
print("      since sum_x d_x = sum_f 2 T_f = 0.  So the image is EVEN iff n is ODD.")
print("  C2. The BICYCLE space Cut & Cycle of K_n: a cut delta(S) is even-degree iff")
print("      |S| and |S^c| are both even, impossible when n is odd.  So bicycle = 0 iff")
print("      n odd -- which is exactly when Cut (+) Cycle is a DIRECT sum.")
print()

def bicycle_dim(n):
    pairs, idx, E = scaffold(n)
    # cut space basis: stars
    stars = []
    for v in range(n):
        b = 0
        for e, (i, j) in enumerate(pairs):
            if i == v or j == v: b |= 1 << e
        stars.append(b)
    def rref(vs):
        bas, piv = [], []
        for v in vs:
            for b, p in zip(bas, piv):
                if v >> p & 1: v ^= b
            if v: bas.append(v); piv.append(v.bit_length() - 1)
        return bas, piv
    cb, cp = rref(stars)
    # cycle space = even-degree vectors; intersect by testing each cut for evenness
    inter = []
    for m in range(1, 1 << len(cb)):
        v = 0
        for b in range(len(cb)):
            if m >> b & 1: v ^= cb[b]
        deg = [0] * n
        for e, (i, j) in enumerate(pairs):
            if v >> e & 1: deg[i] ^= 1; deg[j] ^= 1
        if not any(deg): inter.append(v)
    ib, _ = rref(inter)
    return len(cb), len(ib)

def projection_even(n, trials=200):
    """test C1: is (I+L)T mod 2 always an even graph?"""
    pairs, idx, E = scaffold(n)
    rng = np.random.default_rng(130)
    ok = True
    for _ in range(trials):
        T = rng.integers(0, 2, E)
        out = np.zeros(E, dtype=int)
        for e, (i, j) in enumerate(pairs):
            acc = T[e]
            for f, (a, b) in enumerate(pairs):
                if f != e and len({i, j} & {a, b}) == 1: acc += T[f]
            out[e] = acc % 2
        deg = [0] * n
        for e, (i, j) in enumerate(pairs):
            if out[e]: deg[i] ^= 1; deg[j] ^= 1
        if any(deg): ok = False
    return ok

print(f"{'n':>3} {'parity':>7} {'cut dim':>8} {'BICYCLE dim':>12} {'predicted':>10} "
      f"{'(I+L)T even?':>14}")
for n in range(3, 10):
    cd, bd = bicycle_dim(n)
    pred = 0 if n % 2 else n - 2
    ev = projection_even(n) if n <= 8 else None
    print(f"{n:>3} {('odd' if n%2 else 'even'):>7} {cd:>8} {bd:>12} {pred:>10} "
          f"{str(ev):>14}")

# ============================================================ PART D
print()
print("=" * 76)
print("PART D -- circulant tournaments: the eigenvalues ARE pure sine sums")
print("=" * 76)
print("  mu_j = sum_{k in C} (w^{jk} - w^{-jk}) = 2i * sum_{k in C} sin(2 pi j k / n)")
print("  and since C u (-C) = Z_n\\{0} disjointly, z_j = sum_{k in C} w^{jk} obeys")
print("  z_j + conj(z_j) = -1, i.e. Re(z_j) = -1/2 for EVERY j != 0: the cosine part is")
print("  PINNED to a critical line and all the content sits in the sine part.")
print()
for n in (5, 7, 9, 11):
    half = (n - 1) // 2
    worst_ev = worst_re = 0.0
    npass = 0
    for m in range(1 << half):
        C = [(k + 1) if not (m >> k & 1) else (n - k - 1) for k in range(half)]
        S = np.zeros((n, n))
        for i in range(n):
            for k in C: S[i, (i + k) % n] = 1.0; S[(i + k) % n, i] = -1.0
        got = np.sort(np.linalg.eigvals(S).imag)
        pred = np.sort([2 * sum(np.sin(2 * np.pi * j * k / n) for k in C)
                        for j in range(n)])
        worst_ev = max(worst_ev, float(np.max(np.abs(got - pred))))
        for j in range(1, n):
            z = sum(np.exp(2j * np.pi * j * k / n) for k in C)
            worst_re = max(worst_re, abs(z.real + 0.5))
        npass += 1
    print(f"  n={n:>2}: {npass:>4} circulant tournaments   "
          f"max |eig - sine-sum| = {worst_ev:.2e}   max |Re(z_j) + 1/2| = {worst_re:.2e}")

print()
print("  Paley check (n = q = 3 mod 4): spectrum should be {0} u {+-i sqrt(q)}.")
for q in (3, 7, 11, 19, 23):
    QR = {(x * x) % q for x in range(1, q)}
    S = np.zeros((q, q))
    for i in range(q):
        for k in QR: S[i, (i + k) % q] = 1.0; S[(i + k) % q, i] = -1.0
    ev = np.sort(np.abs(np.linalg.eigvals(S).imag))
    print(f"    q={q:>2}: sqrt(q)={np.sqrt(q):.5f}   distinct |eig| = "
          f"{sorted({round(float(v),5) for v in ev})}")

# ============================================================ PART E
print()
print("=" * 76)
print("PART E -- interlacing: vertex deletion reproduces the sin/cos zero interlacing")
print("=" * 76)
print("  iS is Hermitian, so Cauchy interlacing holds for principal submatrices.")
print("  Odd n: spectrum symmetric about 0 and CONTAINS 0   <->  zeros of sin (0, +-pi, ...)")
print("  Even n: spectrum symmetric about 0 and OMITS 0     <->  zeros of cos (+-pi/2, ...)")
print("  Deleting a vertex alternates the two, and the spectra must interlace.")
print()
for n in range(4, 8):
    pairs, idx, E = scaffold(n)
    reps = canon_codes(n)
    ok = True; strict = 0; tot = 0
    for r in reps:
        S = skew(r, pairs, n); a = spec(S)
        for v in range(n):
            keep = [u for u in range(n) if u != v]
            b = np.sort(np.linalg.eigvalsh(1j * S[np.ix_(keep, keep)]).real)
            tot += 1
            for k in range(n - 1):
                if not (a[k] <= b[k] + 1e-7 and b[k] <= a[k + 1] + 1e-7): ok = False
            if all(a[k] < b[k] - 1e-7 and b[k] < a[k+1] - 1e-7 for k in range(n - 1)):
                strict += 1
    print(f"  n={n} -> {n-1}: Cauchy interlacing holds on all {tot} deletions? {ok}"
          f"   (strict in {strict}/{tot})")

print()
print("SUMMARY")
print("  The two 'odd n only' facts are the same parity fact:")
print("    skew side  -- p(-x) = (-1)^n p(x) forces a zero eigenvalue at odd n;")
print("    GF(2) side -- the image degree is (n-1)*d_w, and bicycle = 0, both at odd n.")
print("  sin is odd and vanishes at 0; cos is even and does not.  The tournament")
print("  characteristic polynomial is sin-like at odd n and cos-like at even n, the")
print("  circulant eigenvalues are literal sine sums with the cosine part pinned at")
print("  -1/2, and vertex deletion interlaces the two exactly as sin and cos interlace.")
