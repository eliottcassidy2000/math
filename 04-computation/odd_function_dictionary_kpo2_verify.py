# odd_function_dictionary_kpo2_verify.py
# ADVERSARIAL VERIFIER for session kind-pasteur-2026-06-10-S2, thread A
# (A-odd-function-dictionary). Independent re-verification of claims
# A1, A2, A3, A4, A5, A6, A10 (+ the set-identity part of A11).
#
# FRESH code: brute force over all 2^(n-1) subsets / sign vectors (not
# pair-based DP), translation-reduced direct distinct-tuple enumeration for
# cluster integrals (A_L per the reflection's definition:
#   A_L = sum over DISTINCT (x_0..x_L) in Z_n of prod sigma(x_{i+1}-x_i),
# from 07-reflections/why-the-paley-path-ratio-is-e-...md), tr(M^L) via
# exact integer circular convolution cross-checked against direct matrix
# powers, and a generic degree-by-degree FGL solver over Fractions.
# Exact integer / rational arithmetic throughout. verifier: kind-pasteur-o2-verify

import sys, random
from fractions import Fraction
from math import comb

CHECKS = []
def check(name, ok, detail=""):
    CHECKS.append((name, ok))
    print(("PASS " if ok else "FAIL ") + name + ((" | " + detail) if detail else ""))

# ---------------------------------------------------------------- utilities
def is_tournament_set(S, n):
    """S subset of {1..n-1}: exactly one of d, n-d for every d."""
    for d in range(1, n):
        if ((d in S) + (((n - d) % n) in S)) != 1:
            return False
    return True

def sigma_from_set(S, n):
    return {d: (1 if d in S else -1) for d in range(1, n)}

def adj_from_set(S, n):
    return [[1 if ((b - a) % n) in S else 0 for b in range(n)] for a in range(n)]

# ============================================================ PART A1
print("=" * 72)
print("A1: dictionary bijection, counts, T^op, relabeling, even n")
print("=" * 72)
for n in [3, 5, 7, 9, 11, 4, 6, 8]:
    # brute force over ALL subsets of {1..n-1}
    tourn_sets = []
    for mask in range(1 << (n - 1)):
        S = {d for d in range(1, n) if (mask >> (d - 1)) & 1}
        if is_tournament_set(S, n):
            tourn_sets.append(S)
    # brute force over ALL sign vectors sigma: {1..n-1} -> {+-1}
    n_odd_sigma = 0
    for mask in range(1 << (n - 1)):
        sig = {d: (1 if (mask >> (d - 1)) & 1 else -1) for d in range(1, n)}
        if all(sig[(n - d) % n] == -sig[d] for d in range(1, n)):
            n_odd_sigma += 1
    expect = 2 ** ((n - 1) // 2) if n % 2 == 1 else 0
    check(f"A1 count n={n}", len(tourn_sets) == expect and n_odd_sigma == expect,
          f"#circ-tourn={len(tourn_sets)}, #odd-sigma={n_odd_sigma}, expect={expect}")
    if n % 2 == 1:
        # roundtrip + A + A^T = J - I + T^op = -sigma  for EVERY tournament
        ok_rt = ok_j = ok_op = True
        for S in tourn_sets:
            sig = sigma_from_set(S, n)
            if not all(sig[(n - d) % n] == -sig[d] for d in range(1, n)): ok_rt = False
            if {d for d in range(1, n) if sig[d] == 1} != S: ok_rt = False
            A = adj_from_set(S, n)
            for a in range(n):
                for b in range(n):
                    want = 0 if a == b else 1
                    if A[a][b] + A[b][a] != want: ok_j = False
            # T^op: reverse every arc -> adjacency transpose; its connection set
            Sop = {d for d in range(1, n) if A[(0 + d) % n][0] == 1}  # arcs into 0... rebuild honestly:
            Sop = {(n - d) % n for d in S}
            sig_op = sigma_from_set(Sop, n)
            if not is_tournament_set(Sop, n): ok_op = False
            if not all(sig_op[d] == -sig[d] for d in range(1, n)): ok_op = False
            # check Sop really is the arc-reversed tournament
            Aop = adj_from_set(Sop, n)
            for a in range(n):
                for b in range(n):
                    if Aop[a][b] != A[b][a]: ok_op = False
        check(f"A1 roundtrip n={n}", ok_rt)
        check(f"A1 A+A^T=J-I n={n}", ok_j)
        check(f"A1 T^op=-sigma n={n}", ok_op)
if True:
    # relabel action x -> ux : sigma'(d) = sigma(u^{-1} d), ALL units, ALL tournaments, n=7,9
    for n in [7, 9]:
        units = [u for u in range(1, n) if __import__("math").gcd(u, n) == 1]
        tourn_sets = []
        for mask in range(1 << (n - 1)):
            S = {d for d in range(1, n) if (mask >> (d - 1)) & 1}
            if is_tournament_set(S, n):
                tourn_sets.append(S)
        ok = True
        for S in tourn_sets:
            A = adj_from_set(S, n)
            sig = sigma_from_set(S, n)
            for u in units:
                uinv = pow(u, -1, n)
                # permuted adjacency B[a][b] = A[u^{-1}a][u^{-1}b]  (relabel x -> ux)
                B = [[A[(uinv * a) % n][(uinv * b) % n] for b in range(n)] for a in range(n)]
                uS = {(u * d) % n for d in S}
                if B != adj_from_set(uS, n): ok = False
                sig2 = sigma_from_set(uS, n)
                if not all(sig2[d] == sig[(uinv * d) % n] for d in range(1, n)): ok = False
        check(f"A1 relabel sigma'(d)=sigma(u^-1 d) n={n} (all units x all tournaments)", ok)

# ============================================================ PART A2
print("=" * 72)
print("A2: completely multiplicative +-1 functions on F_p^*")
print("=" * 72)
def legendre(a, p):
    r = pow(a % p, (p - 1) // 2, p)
    return 1 if r == 1 else (-1 if r == p - 1 else 0)

PRIMES = [p for p in range(3, 60) if all(p % q for q in range(2, p))]
# exhaustive search p <= 19 (worker did <= 19 as well; independent code)
for p in [3, 5, 7, 11, 13, 17, 19]:
    surv = []
    pairs = [(a, b) for a in range(1, p) for b in range(a, p)]
    for mask in range(1 << (p - 1)):
        f = [0] + [1 if (mask >> (d - 1)) & 1 else -1 for d in range(1, p)]
        good = True
        for a, b in pairs:
            if f[(a * b) % p] != f[a] * f[b]:
                good = False; break
        if good:
            surv.append(tuple(f[1:]))
    triv = tuple(1 for d in range(1, p))
    chi = tuple(legendre(d, p) for d in range(1, p))
    check(f"A2 exhaustive p={p}: CM survivors == {{triv,chi}}",
          sorted(surv) == sorted({triv, chi}), f"#survivors={len(surv)}")
for p in PRIMES:
    # structural: homs C_{p-1} -> C_2; primitive root
    g = next(g for g in range(2, p)
             if all(pow(g, (p - 1) // q, p) != 1
                    for q in {q for q in range(2, p) if (p - 1) % q == 0 and all(q % r for r in range(2, q))}))
    chi = {d: legendre(d, p) for d in range(1, p)}
    triv = {d: 1 for d in range(1, p)}
    mult_ok = all(chi[(a * b) % p] == chi[a] * chi[b] for a in range(1, p) for b in range(1, p))
    triv_even = (triv[p - 1] == 1)                      # f(-1)=+1: even
    chi_odd = all(chi[(p - d) % p] == -chi[d] for d in range(1, p))
    euler = (chi[p - 1] == (-1) ** ((p - 1) // 2))
    QR = {d for d in range(1, p) if chi[d] == 1}
    qr_is_tourn = is_tournament_set(QR, p)
    triv_is_tourn = is_tournament_set(set(range(1, p)), p)
    check(f"A2 p={p:2d} ({p%4} mod 4)",
          mult_ok and triv_even and euler and (not triv_is_tourn)
          and (chi_odd == (p % 4 == 3)) and (qr_is_tourn == (p % 4 == 3)),
          f"chi odd={chi_odd}, QR tourn={qr_is_tourn}, chi(-1)={chi[p-1]}")

# ============================================================ cluster integrals
print("=" * 72)
print("A3/A4/A5: cluster integrals A_1..A_4 (distinct-tuple def, translation-reduced)")
print("=" * 72)

def AL(sig, n, L):
    """A_L = n * sum over distinct (0,x1..xL), exact int. Direct enumeration."""
    others = [x for x in range(1, n)]
    total = 0
    if L == 1:
        for x1 in others:
            total += sig[x1]
    elif L == 2:
        for x1 in others:
            s1 = sig[x1]
            for x2 in others:
                if x2 != x1:
                    total += s1 * sig[(x2 - x1) % n]
    elif L == 3:
        for x1 in others:
            s1 = sig[x1]
            for x2 in others:
                if x2 == x1: continue
                s2 = s1 * sig[(x2 - x1) % n]
                for x3 in others:
                    if x3 == x1 or x3 == x2: continue
                    total += s2 * sig[(x3 - x2) % n]
    elif L == 4:
        for x1 in others:
            s1 = sig[x1]
            for x2 in others:
                if x2 == x1: continue
                s2 = s1 * sig[(x2 - x1) % n]
                for x3 in others:
                    if x3 == x1 or x3 == x2: continue
                    s3 = s2 * sig[(x3 - x2) % n]
                    for x4 in others:
                        if x4 == x1 or x4 == x2 or x4 == x3: continue
                        total += s3 * sig[(x4 - x3) % n]
    return n * total

def trML(sig, n, L):
    """tr(M^L) for circulant M[a,b]=sigma(b-a) via exact circular convolution."""
    c = [0] + [sig[d] for d in range(1, n)]
    def conv(u, v):
        w = [0] * n
        for i in range(n):
            ui = u[i]
            if ui:
                for j in range(n):
                    w[(i + j) % n] += ui * v[j]
        return w
    pw = c
    for _ in range(L - 1):
        pw = conv(pw, c)
    return n * pw[0]

def trML_direct(sig, n, L):
    """independent cross-check: trace of explicit matrix power"""
    M = [[sig[(b - a) % n] if a != b else 0 for b in range(n)] for a in range(n)]
    P = [row[:] for row in M]
    for _ in range(L - 1):
        P = [[sum(P[i][k] * M[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    return sum(P[i][i] for i in range(n))

def paley_sigma(p):
    return {d: legendre(d, p) for d in range(1, p)}

def rot_sigma(n):
    m = (n - 1) // 2
    return {d: (1 if 1 <= d <= m else -1) for d in range(1, n)}

def random_odd_sigma(n, seed):
    rng = random.Random(seed)
    sig = {}
    for d in range(1, (n - 1) // 2 + 1):
        s = rng.choice([1, -1])
        sig[d] = s; sig[n - d] = -s
    return sig

families = []
for p in [7, 11, 19, 23]:
    families.append((f"paley_{p}", p, paley_sigma(p)))
for n in [9, 13, 15, 21, 31]:
    families.append((f"rot_{n}", n, rot_sigma(n)))
families.append(("block9={1,2,5,6}", 9, sigma_from_set({1, 2, 5, 6}, 9)))
for n, sd in [(15, 101), (21, 202), (27, 303)]:
    families.append((f"rand_{n}_s{sd}", n, random_odd_sigma(n, sd)))

# validate convolution trace code against direct matrix powers (small n)
ok = all(trML(sig, n, L) == trML_direct(sig, n, L)
         for (nm, n, sig) in families if n <= 13 for L in (2, 3, 4))
check("conv-trace == direct matrix-power trace (all fams n<=13, L=2,3,4)", ok)

print(f"{'family':>16} | {'A_1':>4} {'A_3':>6} | {'A_2':>6} =n(n-1)? -tr(M^2)? | {'A_4':>9} {'tr(M^4)':>9} {'E/n^3':>7}")
for nm, n, sig in families:
    odd_ok = all(sig[(n - d) % n] == -sig[d] for d in range(1, n))
    a1 = AL(sig, n, 1); a2 = AL(sig, n, 2); a3 = AL(sig, n, 3)
    a4 = AL(sig, n, 4)
    t2 = trML(sig, n, 2); t4 = trML(sig, n, 4)
    E = a4 + t4
    print(f"{nm:>16} | {a1:>4} {a3:>6} | {a2:>6} {a2==n*(n-1)!s:>8} {a2==-t2!s:>9} | {a4:>9} {t4:>9} {E/n**3:>7.3f}")
    check(f"A3 A_1=A_3=0 [{nm}]", a1 == 0 and a3 == 0 and odd_ok)
    check(f"A4 A_2=n(n-1)=-trM^2 [{nm}]", a2 == n * (n - 1) and a2 == -t2)
    check(f"A5 |A_4+trM^4| <= 3 n^3 [{nm}]", abs(E) <= 3 * n ** 3, f"E/n^3={E/n**3:.3f}")

check("A6 exact A_4(rot_31) == -225680 (claimed)", AL(rot_sigma(31), 31, 4) == -225680,
      f"A_4(rot_31)={AL(rot_sigma(31),31,4)}")

# Paley DRT: tr(M^4) = p^2 (p-1)
for p in [7, 11, 19, 23, 31, 43]:
    t4 = trML(paley_sigma(p), p, 4)
    check(f"A5 Paley tr(M^4)=p^2(p-1) p={p}", t4 == p * p * (p - 1), f"tr={t4}")

# Paley A_4 ~ +2p^3 (leading order, from reflection/THM-438)
for p in [19, 23]:
    a4 = AL(paley_sigma(p), p, 4)
    check(f"A6 Paley A_4>0, A_4/2p^3 in (0.5,1.5) p={p}", a4 > 0 and 0.5 < a4 / (2 * p ** 3) < 1.5,
          f"A_4={a4}, A_4/2p^3={a4/(2*p**3):.3f}")

# ============================================================ A6: tanh generators
print("=" * 72)
print("A6: rotation-family generators = tanh Taylor coefficients (exact + numeric)")
print("=" * 72)
# Bernoulli numbers, exact
def bernoulli(N):
    B = [Fraction(0)] * (N + 1)
    for m in range(N + 1):
        # recurrence: sum_{j=0}^{m} C(m+1,j) B_j = 0  (m>=1), B_0=1
        if m == 0:
            B[0] = Fraction(1)
        else:
            s = sum(Fraction(comb(m + 1, j)) * B[j] for j in range(m))
            B[m] = -s / (m + 1)
    return B
B = bernoulli(14)

# tanh coefficients via exact series division sinh/cosh (deg 13)
NDEG = 13
sinh = [Fraction(0)] * (NDEG + 1)
cosh = [Fraction(0)] * (NDEG + 1)
fact = [1] * (NDEG + 1)
for i in range(1, NDEG + 1): fact[i] = fact[i - 1] * i
for k in range(0, NDEG + 1):
    if k % 2 == 1: sinh[k] = Fraction(1, fact[k])
    else: cosh[k] = Fraction(1, fact[k])
tanh = [Fraction(0)] * (NDEG + 1)
acc = sinh[:]
for k in range(NDEG + 1):
    tanh[k] = acc[k]
    if tanh[k]:
        for j in range(k, NDEG + 1):
            acc[j] -= tanh[k] * cosh[j - k]
expected_t = [Fraction(1), Fraction(-1, 3), Fraction(2, 15), Fraction(-17, 315),
              Fraction(62, 2835), Fraction(-1382, 155925)]
check("A6 tanh coeffs (series division) = 1,-1/3,2/15,-17/315,62/2835,-1382/155925",
      [tanh[2 * k - 1] for k in range(1, 7)] == expected_t)

# zeta formula: a_{2k} = (-1)^(k-1) * 2 * (2/pi)^(2k) * (1-2^(-2k)) * zeta(2k)
# with zeta(2k)/pi^(2k) = (-1)^(k+1) B_{2k} 2^(2k-1)/(2k)!  -> exact rational
ok = True
for k in range(1, 7):
    z_over_pi = (-1) ** (k + 1) * B[2 * k] * Fraction(2 ** (2 * k - 1), fact[2 * k])
    f_k = (-1) ** (k - 1) * 2 * Fraction(2 ** (2 * k)) * (1 - Fraction(1, 2 ** (2 * k))) * z_over_pi
    if f_k != tanh[2 * k - 1]: ok = False
    print(f"  k={k}: zeta-formula={f_k} vs [x^{2*k-1}]tanh={tanh[2*k-1]}")
check("A6 zeta/Bernoulli formula == tanh coefficients, k=1..6 (exact Fractions)", ok)

# numeric: -tr(M^{2k})/n^{2k} for rotation, exact integer conv, n = 101..401
def trM_powers(sig, n, Lmax):
    c = [0] + [sig[d] for d in range(1, n)]
    out = {}
    pw = c
    for L in range(2, Lmax + 1):
        nxt = [0] * n
        for i in range(n):
            ui = pw[i]
            if ui:
                for j in range(n):
                    nxt[(i + j) % n] += ui * c[j]
        pw = nxt
        out[L] = n * pw[0]
    return out

vals4 = {}; vals6 = {}
for n in [101, 201, 401]:
    t = trM_powers(rot_sigma(n), n, 6)
    vals4[n] = -t[4] / n ** 4
    vals6[n] = -t[6] / n ** 6
    print(f"  n={n}: -tr(M^4)/n^4 = {vals4[n]:+.6f}   -tr(M^6)/n^6 = {vals6[n]:+.6f}")
rich4 = 2 * vals4[401] - vals4[201]   # error ~ c/n Richardson (step doubling)
rich6 = 2 * vals6[401] - vals6[201]
print(f"  Richardson(1/n): a_4 ~ {rich4:+.6f} (target -1/3 = {-1/3:+.6f}); a_6 ~ {rich6:+.6f} (target 2/15 = {2/15:+.6f})")
check("A6 numeric a_4 -> -1/3 (Richardson within 0.002)", abs(rich4 + 1 / 3) < 0.002, f"{rich4:+.6f}")
check("A6 numeric a_6 -> +2/15 (Richardson within 0.002)", abs(rich6 - 2 / 15) < 0.002, f"{rich6:+.6f}")

# ============================================================ A10: odd formal group
print("=" * 72)
print("A10: F(x,y)=(x+y)/(1+xy) is an odd formal group (exact Fraction series)")
print("=" * 72)
NB = 13  # bivariate truncation (total degree)
def bmul(u, v, N=NB):
    w = {}
    for (i, j), a in u.items():
        if a == 0: continue
        for (k, l), b in v.items():
            if b == 0: continue
            if i + k + j + l <= N:
                w[(i + k, j + l)] = w.get((i + k, j + l), Fraction(0)) + a * b
    return {k: v for k, v in w.items() if v != 0}
def badd(u, v):
    w = dict(u)
    for k, val in v.items():
        w[k] = w.get(k, Fraction(0)) + val
    return {k: v for k, v in w.items() if v != 0}
def bscale(u, c):
    return {k: c * v for k, v in u.items() if c * v != 0}

X = {(1, 0): Fraction(1)}; Y = {(0, 1): Fraction(1)}
XY = bmul(X, Y)
# (1+xy)^{-1} = sum (-xy)^k
inv1xy = {(0, 0): Fraction(1)}
term = {(0, 0): Fraction(1)}
for _ in range(NB // 2 + 1):
    term = bscale(bmul(term, XY), Fraction(-1))
    inv1xy = badd(inv1xy, term)
def cayley_F():
    return bmul(badd(X, Y), inv1xy)
def mult_F():
    return badd(badd(X, Y), XY)
def tangent_F():
    inv = {(0, 0): Fraction(1)}
    t = {(0, 0): Fraction(1)}
    for _ in range(NB // 2 + 1):
        t = bmul(t, XY)
        inv = badd(inv, t)
    return bmul(badd(X, Y), inv)

def biv_from_log(lcoef, N=NB):
    """FGL from odd-or-not log l(x)=x+...: F = l^{-1}(l(x)+l(y)); lcoef[k]=coeff x^k"""
    # compositional inverse of l, univariate, degree N
    inv = [Fraction(0)] * (N + 1); inv[1] = Fraction(1)
    for m in range(2, N + 1):
        # coefficient of x^m in l(inv(x)) must be 0
        c = Fraction(0)
        # compute l(inv) up to x^m with current inv (inv[m]=0 placeholder)
        powv = [Fraction(0)] * (N + 1); powv[0] = Fraction(1)
        tot = [Fraction(0)] * (N + 1)
        cur = [Fraction(0)] * (N + 1)
        cur[0] = Fraction(1)
        for k in range(1, N + 1):
            # cur = inv^k
            new = [Fraction(0)] * (N + 1)
            for i in range(N + 1):
                if cur[i]:
                    for j in range(1, N + 1 - i):
                        new[i + j] += cur[i] * inv[j]
            cur = new
            if k <= len(lcoef) - 1 and lcoef[k]:
                for i in range(N + 1):
                    tot[i] += lcoef[k] * cur[i]
        inv[m] = -tot[m]
    # A(x,y) = l(x)+l(y); F = inv(A)
    A = {}
    for k in range(1, len(lcoef)):
        if lcoef[k]:
            A[(k, 0)] = A.get((k, 0), Fraction(0)) + lcoef[k]
            A[(0, k)] = A.get((0, k), Fraction(0)) + lcoef[k]
    F = {}
    powA = {(0, 0): Fraction(1)}
    for k in range(1, N + 1):
        powA = bmul(powA, A)
        if inv[k]:
            F = badd(F, bscale(powA, inv[k]))
    return F, inv

def neg_series(F, N=NB):
    """solve F(x, i(x)) = 0 degree by degree; returns list i[k]"""
    i = [Fraction(0)] * (N + 1)
    i[1] = Fraction(-1)
    for m in range(2, N + 1):
        # r(x) = F(x, i(x)) truncated to degree m
        # substitute y -> i(x)
        r = [Fraction(0)] * (N + 1)
        # precompute powers of i(x)
        ipow = [[Fraction(0)] * (N + 1) for _ in range(N + 1)]
        ipow[0][0] = Fraction(1)
        for k in range(1, N + 1):
            for a in range(N + 1):
                if ipow[k - 1][a]:
                    for b in range(1, N + 1 - a):
                        ipow[k][a + b] += ipow[k - 1][a] * i[b]
        for (a, b), c in F.items():
            for e in range(N + 1 - a):
                if ipow[b][e]:
                    r[a + e] += c * ipow[b][e]
        i[m] = -r[m]
    return i

def F_props(F, name):
    odd_total = all((a + b) % 2 == 1 for (a, b) in F)
    Fneg = {(a, b): c * ((-1) ** (a + b)) for (a, b), c in F.items()}
    anti = (Fneg == {k: -c for k, c in F.items()})
    i = neg_series(F)
    minus_x = (i[1] == -1 and all(i[k] == 0 for k in range(2, NB + 1)))
    return odd_total, anti, i, minus_x

Fc = cayley_F()
odd_tot, anti, ineg, minus_x = F_props(Fc, "cayley")
check("A10 Cayley F: only odd total degrees", odd_tot)
check("A10 Cayley F: F(-x,-y) = -F(x,y)", anti)
check("A10 Cayley F: [-1](x) = -x exactly (deg<=13, generic solver)", minus_x,
      "i = " + str([str(c) for c in ineg[1:6]]))
check("A10 Cayley F: integer coefficients", all(c.denominator == 1 for c in Fc.values()))

# log = arctanh: check L(F(x,y)) = L(x)+L(y) exactly (bivariate, deg<=13)
def apply_univ(lcoef, F, N=NB):
    out = {}
    powF = {(0, 0): Fraction(1)}
    for k in range(1, N + 1):
        powF = bmul(powF, F)
        if k < len(lcoef) and lcoef[k]:
            out = badd(out, bscale(powF, lcoef[k]))
    return out
arctanh = [Fraction(0)] * (NB + 1)
for k in range(0, NB + 1):
    if k % 2 == 1: arctanh[k] = Fraction(1, k)
lhs = apply_univ(arctanh, Fc)
rhs = {}
for k in range(1, NB + 1):
    if arctanh[k]:
        rhs[(k, 0)] = arctanh[k]; rhs[(0, k)] = arctanh[k]
check("A10 arctanh(F(x,y)) = arctanh x + arctanh y (deg<=13 exact)", lhs == rhs)
check("A10 arctanh has only odd powers (definitionally)", all(arctanh[k] == 0 for k in range(2, NB + 1, 2)))
# round trip: tanh(arctanh(x)) = x
tanh_l = [Fraction(0)] * (NB + 1)
for k in range(NB + 1):
    if k <= NDEG: tanh_l[k] = tanh[k]
comp = [Fraction(0)] * (NB + 1)
cur = [Fraction(0)] * (NB + 1); cur[0] = Fraction(1)
for k in range(1, NB + 1):
    new = [Fraction(0)] * (NB + 1)
    for i2 in range(NB + 1):
        if cur[i2]:
            for j2 in range(1, NB + 1 - i2):
                new[i2 + j2] += cur[i2] * arctanh[j2]
    cur = new
    if tanh_l[k]:
        for i2 in range(NB + 1):
            comp[i2] += tanh_l[k] * cur[i2]
check("A10 tanh(arctanh(x)) = x (deg<=13 exact)", comp[1] == 1 and all(comp[k] == 0 for k in range(2, NB + 1)))

# multiplicative FGL must FAIL all three
Fm = mult_F()
odd_tot, anti, ineg, minus_x = F_props(Fm, "mult")
check("A10 mult FGL x+y+xy fails all three", (not odd_tot) and (not anti) and (not minus_x),
      f"[-1] = {[str(c) for c in ineg[1:5]]} (expect -1,1,-1,1)")
check("A10 mult FGL [-1] = -x+x^2-x^3+... ", ineg[1] == -1 and ineg[2] == 1 and ineg[3] == -1 and ineg[4] == 1)
# tangent FGL passes
Ft = tangent_F()
odd_tot, anti, ineg, minus_x = F_props(Ft, "tangent")
check("A10 tangent FGL (x+y)/(1-xy) passes all three (log=arctan analogue)", odd_tot and anti and minus_x)

# random odd logs pass; random non-odd logs fail (lemma: 3 conditions equivalent)
rng = random.Random(42)
for trial in range(3):
    l = [Fraction(0)] * (NB + 1); l[1] = Fraction(1)
    for k in range(3, 10, 2):
        l[k] = Fraction(rng.randint(-5, 5), rng.randint(1, 7))
    F, _ = biv_from_log(l)
    o, a, ineg2, m = F_props(F, "rndodd")
    check(f"A10 random odd log #{trial} passes all three", o and a and m)
for trial in range(3):
    l = [Fraction(0)] * (NB + 1); l[1] = Fraction(1)
    for k in range(2, 9):
        l[k] = Fraction(rng.randint(-5, 5), rng.randint(1, 7))
    if all(l[k] == 0 for k in range(2, NB + 1, 2)):
        l[2] = Fraction(1, 2)  # force non-odd
    F, _ = biv_from_log(l)
    o, a, ineg2, m = F_props(F, "rndnonodd")
    check(f"A10 random non-odd log #{trial} fails all three", (not o) and (not a) and (not m))

# ============================================================ A11 set identities
print("=" * 72)
print("A11: unit-relabel set identities (H equality checked in census verify script)")
print("=" * 72)
s9 = {(5 * d) % 9 for d in {1, 2, 3, 4}}
check("A11 5*rot_9 = {1,2,5,6} (a tournament set)", s9 == {1, 2, 5, 6} and is_tournament_set(s9, 9), str(sorted(s9)))
s13 = {(7 * d) % 13 for d in {1, 2, 3, 4, 5, 6}}
check("A11 7*rot_13 = {1,2,3,7,8,9} (a tournament set)", s13 == {1, 2, 3, 7, 8, 9} and is_tournament_set(s13, 13), str(sorted(s13)))

# ============================================================ summary
print("=" * 72)
nf = sum(1 for _, ok in CHECKS if not ok)
print(f"TOTAL {len(CHECKS)} checks, {nf} failures")
sys.exit(1 if nf else 0)
