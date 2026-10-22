# Paley tournaments -> binary quadratic residue codes -> lattices.
# For prime q = 3 mod 4 with q = -1 mod 8 (2 a QR, eps = (2/q) = +1), the F_2-row space of A (A = adjacency of P_q, i -> j iff
# j - i is a nonzero square; the rows are the translates of QR_q, for q = 7 the 7 lines of the Fano plane) is the binary QR code
# of length q (A + I spans only its even-weight subcode); adding a parity coordinate gives the extended QR code of length
# q + 1, which is doubly even and self-dual (Type II).  q = 7: extended Hamming [8,4,4]; q = 23: extended Golay [24,12,8].
# Construction A (L = {x in Z^n : x mod 2 in C} / sqrt 2) then gives an even unimodular lattice: E_8 for q = 7, the Niemeier
# lattice A_1^24 for q = 23; theta series = sum over codewords via the weight enumerator, compared with E_4 and E_4^3 + (N_2 - 720)Delta.
# mac-mini-2026-10-21-fanopsp.
import itertools
def rank_rows(rows, n):
    piv = []; basis = []
    for r in rows:
        x = r
        for b in basis:
            if x ^ b < x: x ^= b
        if x: basis.append(x); basis.sort(reverse=True)
    return basis
def span_weights(basis, n):
    from collections import Counter
    cnt = Counter()
    k = len(basis)
    # Gray-code enumeration
    x = 0; cnt[0] += 1
    for i in range(1, 1 << k):
        j = (i & -i).bit_length() - 1
        x ^= basis[j]; cnt[bin(x).count('1')] += 1
    return dict(sorted(cnt.items()))
def theta_from_code(wt, n, N):
    # Construction A, scaled: vectors x in Z^n with x mod 2 in C; norm(x)/2.  theta(q) = sum_c prod over coords of
    # (theta_even or theta_odd) where theta_even = sum_{m even} q^{m^2/4}, theta_odd = sum_{m odd} q^{m^2/4} in q^(1/4) units.
    # Work in powers of q^(1/4): norm/2 = sum x_i^2 / 2 -> exponent (sum x_i^2)/2 in q.  Use variable t = q^(1/2): exponent sum x_i^2.
    M = 4 * N + 1
    ev = [0] * M; od = [0] * M
    for m in range(-60, 61):
        e = m * m
        if e < M: (ev if m % 2 == 0 else od)[e] += 1
    def mul(a, b):
        c = [0] * M
        for i, x in enumerate(a):
            if x:
                for j in range(M - i):
                    if b[j]: c[i + j] += x * b[j]
        return c
    def pw(a, k):
        r = [1] + [0] * (M - 1)
        for _ in range(k): r = mul(r, a)
        return r
    tot = [0] * M
    for w, c in wt.items():
        term = mul(pw(od, w), pw(ev, n - w))
        for i in range(M): tot[i] += c * term[i]
    # exponent sum x_i^2 = 2 * (norm in the scaled lattice) ... lattice norm = sum x_i^2 / 2, q-exponent = norm/2 = sum x_i^2 / 4
    return [tot[4 * k] for k in range(N + 1)], [tot[i] for i in range(M) if i % 4]
def eis(k, N):
    from math import comb
    from fractions import Fraction
    B = {4: Fraction(-1, 30), 6: Fraction(1, 42), 12: Fraction(691, 2730)}
    c = -2 * k / B[k] if False else None
    coef = {4: 240, 6: -504}[k]
    return [1] + [coef * sum(d ** (k - 1) for d in range(1, n + 1) if n % d == 0) for n in range(1, N + 1)]
def mulq(a, b, N): return [sum(a[i] * b[n - i] for i in range(n + 1)) for n in range(N + 1)]
N = 4
E4 = eis(4, N); E6 = eis(6, N)
E43 = mulq(mulq(E4, E4, N), E4, N); E62 = mulq(E6, E6, N)
Delta = [(E43[n] - E62[n]) // 1728 for n in range(N + 1)]
print("Delta =", Delta, " (Ramanujan tau: 0, 1, -24, 252, -1472)")
for q in (7, 23, 31, 47):
    QR = {(x * x) % q for x in range(1, q)}
    rows = [sum(1 << ((i + d) % q) for d in QR) for i in range(q)]              # row i of A (out-neighbourhood of i)
    basis = rank_rows(rows, q)
    k = len(basis)
    wt = span_weights(basis, q) if k <= 24 else None
    ext = None
    if wt is not None:
        ebasis = [b | ((bin(b).count('1') % 2) << q) for b in basis]          # extend by parity
        ewt = span_weights(ebasis, q + 1)
        de = min(w for w in ewt if w)
        doubly_even = all(w % 4 == 0 for w in ewt)
        print(f"q={q}: rank of A over F_2 = {k} (QR code [{q},{k}]), weights {wt}")
        print(f"      extended code [{q+1},{k},{de}], doubly even={doubly_even}, self-dual={(2*k == q+1)}; weight enumerator {ewt}")
        if q in (7, 23):
            th, junk = theta_from_code(ewt, q + 1, N)
            print(f"      Construction A theta (q^0..q^{N}): {th}; nonintegral-exponent terms zero: {not any(junk)}")
            if q == 7: print(f"      E_4 = {E4}  -> match: {th == E4}")
            else:
                ref = [E43[n] + (48 - 720) * Delta[n] for n in range(N + 1)]
                print(f"      E_4^3 + (48 - 720) Delta = {ref} -> match: {th == ref}  (Niemeier A_1^24, 48 roots)")
    else:
        print(f"q={q}: rank of A over F_2 = {k} (QR code [{q},{k}]) - weight enumeration skipped")
