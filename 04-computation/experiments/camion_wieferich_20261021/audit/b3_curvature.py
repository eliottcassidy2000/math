# Audit b3 (independent code): THM-4750 statement 3, the Fermat-quotient curvature
#     a_2(m)/p = -q_2(p) X(m-2)  (mod p),  m >= 3,
# and the quadratic law  eps^(-l) D(m + l(p-1)) = a_0 + a_1 l + a_2 C(l,2)  (mod p^2)  for all l >= 0.
# Independent choices:
#  * s(n) mod q from the Hurwitz inverse W = T^(*-1) (numpy, vectorised), not Wright's recurrence;
#  * X is computed through the OTHER identity of the theorem,  X = S'' + 2 S'^2 T  (EGF), i.e.
#        X(n) = s(n+2) + 2 * sum_{i+j+k=n} n!/(i! j! k!) s(i+1) s(j+1) T(k),
#    instead of the session's Delta^2 T * (e - S)^2;
#  * exact integer values of X(0..6) and an independent brute force of the combinatorial meaning
#    (tournaments on [m] that become strong when the pair {0,1} is made two-way) for m <= 6.
# Usage: python3 b3_curvature.py
import numpy as np
from math import comb
from itertools import combinations, product

def s_mod(N, q):
    """s(0..N) mod q via W(n) = -sum_{k=1}^n C(n,k) T(k) W(n-k);  s(n) = -W(n), s(0) = 0."""
    T = np.array([pow(2, n * (n - 1) // 2, q) for n in range(N + 1)], dtype=object)
    Tn = np.array([int(x) for x in T], dtype=np.int64)
    W = np.zeros(N + 1, dtype=np.int64); W[0] = 1
    row = np.zeros(N + 1, dtype=np.int64); row[0] = 1          # C(n, .) mod q
    for n in range(1, N + 1):
        row[1:n + 1] = (row[1:n + 1] + row[0:n]) % q
        # sum_{k=1}^n C(n,k) T(k) W(n-k)
        ct = (row[1:n + 1] * Tn[1:n + 1]) % q
        acc = int(np.dot(ct.astype(object), W[n - 1::-1][:n].astype(object)) % q) if q > 3_000_000 else int(np.dot(ct, W[n - 1::-1][:n]) % q)
        W[n] = (-acc) % q
    s = (-W) % q; s[0] = 0
    return [int(x) for x in s], [int(x) for x in Tn]

def X_mod_p(p, s, T, nmax):
    """X(n) mod p for n <= nmax via X = S'' + 2 S'^2 T  (Hurwitz products, binomials mod p)."""
    sp = [s[n + 1] % p for n in range(nmax + 1)]                # S'(n) = s(n+1)
    # binomial rows mod p
    rows = [[1]]
    for n in range(1, nmax + 1):
        r = rows[-1]; rows.append([1] + [(r[i - 1] + r[i]) % p for i in range(1, n)] + [1])
    sp2 = [sum(rows[n][i] * sp[i] * sp[n - i] for i in range(n + 1)) % p for n in range(nmax + 1)]
    sp2T = [sum(rows[n][i] * sp2[i] * T[n - i] for i in range(n + 1)) % p for n in range(nmax + 1)]
    return [(s[n + 2] + 2 * sp2T[n]) % p for n in range(nmax + 1)]

# ---- exact X(0..6) and brute force of the combinatorial meaning ----
NE = 12
Te = [2 ** comb(n, 2) for n in range(NE + 3)]
We = [1] + [0] * (NE + 2)
for n in range(1, NE + 3): We[n] = -sum(comb(n, k) * Te[k] * We[n - k] for k in range(1, n + 1))
se = [0] + [-We[n] for n in range(1, NE + 3)]
def Xexact(n):
    tot = se[n + 2]
    for i in range(n + 1):
        for j in range(n + 1 - i):
            k = n - i - j
            tot += 2 * (comb(n, i) * comb(n - i, j)) * se[i + 1] * se[j + 1] * Te[k]
    return tot
Xe = [Xexact(n) for n in range(7)]
def brute(m):
    pairs = list(combinations(range(m), 2)); cnt = 0
    for bits in product((0, 1), repeat=len(pairs)):
        out = [set() for _ in range(m)]
        for (i, j), b in zip(pairs, bits):
            if b: out[i].add(j)
            else: out[j].add(i)
        out[0].add(1); out[1].add(0)                       # make the pair {0,1} two-way
        ok = True
        for src in range(m):                               # strong iff every vertex reaches every vertex
            seen = {src}; st = [src]
            while st:
                u = st.pop()
                for v in out[u]:
                    if v not in seen: seen.add(v); st.append(v)
            if len(seen) < m: ok = False; break
        if ok: cnt += 1
    return cnt
bf = [brute(m) for m in range(2, 7)]
print("X(0..6) exact via S''+2S'^2T:", Xe)
print("brute force m=2..6 (tournaments made strong by a two-way pair {0,1}):", bf, "match" if bf == Xe[:5] else "MISMATCH")
# base-derivative identity C(n,2) X(n-2) = 2 d/da s_a(n) at a=2, checked exactly with sympy-free finite differences of the
# polynomial s_a(n) in a (degree C(n,2)): compute s_a symbolically via integer polynomial arithmetic
def s_poly(n):
    # T_a(k) = a^C(k,2) as coefficient lists; W_a = T_a^(-1) (Hurwitz) as integer polynomials in a
    from collections import defaultdict
    Tp = [{comb(k, 2): 1} for k in range(n + 1)]
    Wp = [{0: 1}] + [None] * n
    for N in range(1, n + 1):
        acc = defaultdict(int)
        for k in range(1, N + 1):
            for e1, c1 in Tp[k].items():
                for e2, c2 in Wp[N - k].items(): acc[e1 + e2] -= comb(N, k) * c1 * c2
        Wp[N] = dict(acc)
    return {e: -c for e, c in Wp[n].items() if c}
ok_der = True
for n in range(2, 8):
    sp_ = s_poly(n)
    der = sum(e * c * 2 ** (e - 1) for e, c in sp_.items() if e >= 1)
    if comb(n, 2) * Xe[n - 2] != 2 * der: ok_der = False
print("C(n,2) X(n-2) = 2 (d/da) s_a(n)|_{a=2} for n = 2..7:", ok_der)

# ---- the curvature formula for all classes of all primes 5 <= p <= PMAX, and the direct quadratic law ----
def primes(a, b): return [q for q in range(a, b + 1) if q > 1 and all(q % d for d in range(2, int(q ** 0.5) + 1))]
PMAX = 400; PDIRECT = 61
form_checks = form_fail = div_fail = 0; direct = direct_fail = 0
for p in primes(5, PMAX):
    q = p * p; eps = 1 if p % 8 in (1, 7) else -1
    N = 3 * p + 2 if p > PDIRECT else (p + 1) + (p + 1) * (p - 1) + 2
    s, T = s_mod(N, q)
    def c(n): return pow(2, (n - 1) // 2, q) if n % 2 else 0
    D = [(s[n] - c(n)) % q for n in range(N + 1)]
    X = X_mod_p(p, s, T, p)
    q2 = ((pow(2, p - 1, q) - 1) // p) % p
    for m in range(3, p + 2):
        A = [(eps ** l * D[m + l * (p - 1)]) % q for l in range(3)]
        a1 = (A[1] - A[0]) % q; a2 = (A[2] - 2 * A[1] + A[0]) % q
        if a1 % p or a2 % p: div_fail += 1; continue
        form_checks += 1
        if (a2 // p) % p != (-q2 * X[m - 2]) % p: form_fail += 1
        if p <= PDIRECT:
            for l in range(0, p + 2):
                n = m + l * (p - 1)
                if n > N: break
                Al = (eps ** l * D[n]) % q
                pred = (A[0] + a1 * l + a2 * (l * (l - 1) // 2)) % q
                direct += 1
                if pred != Al: direct_fail += 1
print(f"curvature formula a_2/p = -q_2 X(m-2) mod p: primes 5..{PMAX}, {form_checks} classes, failures {form_fail}; "
      f"p|a_1, p|a_2 failures {div_fail}")
print(f"direct quadratic law mod p^2, l = 0..p+1, all classes, p <= {PDIRECT}: {direct} checks, failures {direct_fail}")

# ---- p = 1091: which classes have v(a_2) >= 2 (claimed: exactly one, with X(m-2) = 0 mod p) ----
for p in (1091, 1097):
    q = p * p; eps = 1 if p % 8 in (1, 7) else -1
    s, T = s_mod(3 * p + 2, q)
    def c(n): return pow(2, (n - 1) // 2, q) if n % 2 else 0
    D = [(s[n] - c(n)) % q for n in range(3 * p + 3)]
    X = X_mod_p(p, s, T, p)
    hi = []; xz = []
    for m in range(3, p + 2):
        A = [(eps ** l * D[m + l * (p - 1)]) % q for l in range(3)]
        a2 = (A[2] - 2 * A[1] + A[0]) % q
        if a2 == 0: hi.append(m)
        if X[m - 2] == 0: xz.append(m)
    print(f"p={p}: classes (base m in [3,p+1]) with p^2 | a_2: {hi};  classes with X(m-2) = 0 mod p: {xz}")
