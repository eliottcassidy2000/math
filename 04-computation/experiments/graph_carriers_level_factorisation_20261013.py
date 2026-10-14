"""
graph_carriers_level_factorisation_20261013.py -- S30 (owner: unify the deeper themes; graph theory against every other
area).  Five computations behind the synthesis note graph_theory_as_common_carrier_20261013.md.

Part 1.  Level factorisation of the compressed twisted transfer operator for the maps x/2, (px+1)/2 (p odd), on the
         class means modulo 2^j p^k (THM-4629's operator with p in place of 3):
           (M v)(r, rho) = (u/2) [v(r/2, rho/2) + v(r/2 + 2^(j-1), rho/2)]            (r even)
                         = (v/2) [v(s_r, c(rho)) + v(s_r + 2^(j-1), c(rho))]          (r odd), s_r = (pr+1)/2 mod 2^(j-1)
           (j = 0: both branches with weight 1/2 on the p-adic classes alone), c(rho) = (p rho' + 1)/2 mod p^k with rho'
         the odd representative of rho modulo 2 p^k; u = zeta^alpha e^(-s), v = zeta^(alpha+beta) e^(s(theta_p - 1)).
         Claim (THM-4631): det(I - M) = (1 - u/2 - v/2) prod_(i=1)^k (1 - (u/2)^(L_i))^(phi(p^i)/L_i), L_i = ord_(p^i)(2):
         the graded pieces of the p-adic depth filtration are the Cayley graphs of multiplication by 2^(-1) on (Z/p^i)^x
         (disjoint cycles of length L_i), the odd branch lowers depth (nilpotent), the 2-adic part (de Bruijn) only
         contributes the even density 1/2.  Exact symbolic (sympy, symbols u, v) for sizes <= 18; exact rational
         evaluation at three random rational points (integer Bareiss determinants) for sizes <= 400; numeric complex
         determinants at random (u, v) beyond, including the Wieferich primes 1093 and 3511 at level 1.
Part 2.  ord_(p^i)(2), cycle counts phi(p^i)/L_i at levels 1 and 2 for p < 100; Artin (2 a primitive root <=> one cycle
         per level) and the Wieferich jump at level 2; the Paley split at p = 7 ({1,2,4} | {3,5,6}).
Part 3.  The Farey geodesic of log2 3: partial quotients and convergents (P_n, Q_n) to n = N (exact integers from a
         high-precision continued fraction, prefix validated at two precisions); the residues (Q_n, P_n) mod m for
         m = 2..7 against equidistribution over the primitive pairs; the character of THM-4629 that each of the first 40
         convergents resonates with (m = 2: tau <-> (1,0), omega <-> (0,1), C <-> (1,1)); Gauss-Kuzmin check.
Part 4.  Paley-7 against the Q resonator: eigenvalues of the Paley tournament on 7 vertices, the AR(2) polynomial of
         THM-4630's linearisation, the reciprocal relation, the (n+1)/4 = 2 explanation, Ramanujan-Nagell 2^n - 7 = x^2.
Part 5.  The 2-adic (de Bruijn) part alone: spectrum of E_j (even-class averaging) for j = 1..7: one eigenvalue 1/2, the
         rest 0 (nilpotent: rank of E_j^q), so the 2-adic refinement never adds a mode.
"""
import sys, math, random
import numpy as np
import sympy as sp
from fractions import Fraction
from mpmath import mp, log, floor

random.seed(20261013)

def odd_rep(rho, pk):
    return rho if rho % 2 == 1 else rho + pk            # odd representative of rho modulo 2 p^k

def build_M(p, k, j, u, v, exact):
    """Dense matrix of M on functions of (r mod 2^j, rho mod p^k), index r * p^k + rho (j = 0: rho only).
    exact=True: object dtype (sympy / Fraction entries); else complex."""
    pk = p ** k; two_j = 2 ** j; n = two_j * pk
    inv2 = pow(2, -1, pk) if k > 0 else 0
    M = np.zeros((n, n), dtype=object) if exact else np.zeros((n, n), dtype=complex)
    if exact:
        for i in range(n):
            for jj in range(n):
                M[i, jj] = 0
    for r in range(two_j):
        for rho in range(pk):
            row = r * pk + rho
            rho_h = (rho * inv2) % pk if k > 0 else 0
            rho_o = ((p * odd_rep(rho, pk) + 1) // 2) % pk if k > 0 else 0
            if j == 0:
                M[row, rho_h] += u / 2
                M[row, rho_o] += v / 2
            else:
                half = two_j // 2
                if r % 2 == 0:
                    r2 = (r // 2) % half
                    for lift in (r2, r2 + half):
                        M[row, (lift % two_j) * pk + rho_h] += u / 2
                else:
                    s_r = ((p * r + 1) // 2) % half if half > 1 else 0
                    for lift in (s_r, s_r + half):
                        M[row, (lift % two_j) * pk + rho_o] += v / 2
    return M

def ord_mod(a, n):
    if n == 1: return 1
    x, e = a % n, 1
    while x != 1:
        x = (x * a) % n; e += 1
    return e

def phi_pk(p, i): return (p - 1) * p ** (i - 1)

def predicted(p, k, u, v):
    out = 1 - u / 2 - v / 2
    for i in range(1, k + 1):
        L = ord_mod(2, p ** i)
        out = out * (1 - (u / 2) ** L) ** (phi_pk(p, i) // L)
    return out

def exact_rational_det_I_minus_M(p, k, j, uq, vq):
    """det(I - M) exactly for rational uq, vq: scale by the common denominator and use integer Bareiss."""
    M = build_M(p, k, j, Fraction(uq), Fraction(vq), exact=True)
    n = M.shape[0]
    den = 1
    for i in range(n):
        for jj in range(n):
            den = den * M[i, jj].denominator // math.gcd(den, M[i, jj].denominator)
    A = sp.zeros(n, n)
    for i in range(n):
        for jj in range(n):
            A[i, jj] = int((Fraction(int(i == jj)) - M[i, jj]) * den)
    return sp.Rational(A.det(method="bareiss"), den ** n)

print("# Part 1: level factorisation det(I - M_(j,k)) for the maps x/2, (px+1)/2")
u, v = sp.symbols("u v")
symbolic_cases = [(3,1,0),(3,1,1),(3,1,2),(3,2,0),(3,2,1),(5,1,0),(5,1,1),(7,1,0),(7,1,1),(11,1,0),(13,1,0),(17,1,0)]
for (p, k, j) in symbolic_cases:
    M = sp.Matrix(build_M(p, k, j, u, v, exact=True))
    det = sp.factor(sp.expand((sp.eye(M.shape[0]) - M).det(method="berkowitz")))
    pred = sp.factor(sp.expand(predicted(p, k, u, v)))
    ok = sp.expand(det - pred) == 0
    Ls = [ord_mod(2, p ** i) for i in range(1, k + 1)]
    print(f"   [symbolic] p={p:2d} k={k} j={j} (size {2**j * p**k:3d}): L_i = {Ls}, cycles per level = {[phi_pk(p,i)//ord_mod(2,p**i) for i in range(1,k+1)]}: det == predicted: {ok}   det = {det}", flush=True)
print()
rational_cases = [(3,1,3),(3,2,2),(3,3,0),(3,3,1),(3,4,0),(5,1,2),(5,1,3),(5,2,0),(5,2,1),(5,3,0),(7,1,2),(7,2,0),(7,2,1),(7,2,2),(7,3,0),(11,1,1),(11,1,2),(11,2,0),(11,2,1),(13,1,1),(13,2,0),(17,1,1),(17,2,0),(19,1,2),(23,2,0)]
for (p, k, j) in rational_cases:
    n = 2 ** j * p ** k
    allok = True
    for trial in range(3):
        uq = Fraction(random.randint(-9, 9), random.randint(1, 7)); vq = Fraction(random.randint(-9, 9), random.randint(1, 7))
        det = exact_rational_det_I_minus_M(p, k, j, uq, vq)
        pred = predicted(p, k, Fraction(uq), Fraction(vq))
        allok = allok and (sp.Rational(pred.numerator, pred.denominator) == det)
    Ls = [ord_mod(2, p ** i) for i in range(1, k + 1)]
    print(f"   [exact rational, 3 points] p={p:2d} k={k} j={j} (size {n:3d}): L_i = {Ls}, cycles per level = {[phi_pk(p,i)//ord_mod(2,p**i) for i in range(1,k+1)]}: det == predicted at all three points: {allok}", flush=True)
print()
numeric_cases = [(3,5,0),(5,4,0),(7,3,1),(13,2,1),(31,1,3),(37,2,0),(1093,1,0),(3511,1,0)]
for (p, k, j) in numeric_cases:
    uu = complex(random.uniform(-1, 1), random.uniform(-1, 1)); vv = complex(random.uniform(-1, 1), random.uniform(-1, 1))
    M = build_M(p, k, j, uu, vv, exact=False)
    det = np.linalg.det(np.eye(M.shape[0]) - M)
    pred = predicted(p, k, uu, vv)
    Ls = [ord_mod(2, p ** i) for i in range(1, k + 1)]
    print(f"   [numeric] p={p:4d} k={k} j={j} (size {2**j * p**k:5d}): L_i = {Ls}, cycles per level = {[phi_pk(p,i)//ord_mod(2,p**i) for i in range(1,k+1)]}: |det - predicted|/|predicted| = {abs(det - pred)/abs(pred):.2e}", flush=True)

print("\n# Part 2: ord_(p^i)(2) and cycle counts at levels 1, 2 (p < 100), Artin and Wieferich")
primes = [p for p in range(3, 100) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
rows = []
for p in primes:
    L1, L2 = ord_mod(2, p), ord_mod(2, p * p)
    rows.append((p, L1, (p - 1) // L1, L2, (p * (p - 1)) // L2, L2 == L1))
print("   p   ord_p(2)  cycles@1   ord_(p^2)(2)  cycles@2   ord_(p^2)(2) == ord_p(2)")
for p, L1, c1, L2, c2, w in rows:
    print(f"   {p:3d}   {L1:4d}      {c1:3d}        {L2:5d}        {c2:4d}      {w}")
prim = [p for p, L1, c1, L2, c2, w in rows if c1 == 1]
print(f"   primes < 100 with 2 a primitive root (one cycle per level): {prim}  ({len(prim)} of {len(primes)}; Artin's constant 0.3739558 predicts {0.3739558*len(primes):.1f})")
for p in (1093, 3511):
    L1, L2 = ord_mod(2, p), ord_mod(2, p * p)
    print(f"   Wieferich p={p}: ord_p(2) = {L1}, ord_(p^2)(2) = {L2} (equal), cycles at level 1 = {(p-1)//L1}, at level 2 = {p*(p-1)//L2} (jump by the factor p = {p})")
print("   p = 7: the doubling orbits on (Z/7)^x are {1,2,4} and {3,5,6} (quadratic residues and non-residues), the Paley split; level factor (1 - (u/2)^3)^2.")
print("   orbits of x -> x/2 on (Z/7)^x:", end=" ")
seen = set(); inv2 = pow(2, -1, 7)
for x0 in range(1, 7):
    if x0 in seen: continue
    orb = []; x = x0
    while x not in orb:
        orb.append(x); seen.add(x); x = (x * inv2) % 7
    print(orb, end=" ")
print()

print("\n# Part 3: the Farey geodesic of log2 3 -- convergents, residues mod m, resonant characters")
def cf_log2_3(N, dps):
    mp.dps = dps
    x = log(3) / log(2)
    a = []
    for _ in range(N):
        fl = int(floor(x)); a.append(fl)
        frac = x - fl
        if frac == 0: break
        x = 1 / frac
    return a
N = 2600
a1 = cf_log2_3(N, 4500); a2 = cf_log2_3(N, 5200)
common = 0
while common < min(len(a1), len(a2)) and a1[common] == a2[common]: common += 1
a = a1[:common]
print(f"   partial quotients agreeing at dps 4500 and 5200: {common} (used n <= {common})")
print(f"   first 30: {a[:30]}")
P = [1, a[0]]; Q = [0, 1]
for n in range(1, len(a)):
    P.append(a[n] * P[-1] + P[-2]); Q.append(a[n] * Q[-1] + Q[-2])
conv = list(zip(P[1:], Q[1:]))    # (P_n, Q_n), n = 0, 1, 2, ...
print(f"   convergents 1..14 as (Q, P): {[(q, p) for (p, q) in conv[1:15]]}")
print(f"   digits of Q_n at n = {len(conv)-1}: {len(str(conv[-1][1]))}")
big = sorted(range(len(a)), key=lambda i: -a[i])[:8]
print(f"   largest partial quotients (index, value): {[(i, a[i]) for i in sorted(big)]}")
for m in (2, 3, 4, 5, 6, 7):
    prim_pairs = [(qq, pp) for qq in range(m) for pp in range(m) if math.gcd(math.gcd(qq, pp), m) == 1]
    counts = {pr: 0 for pr in prim_pairs}
    for (pp, qq) in conv[1:]:
        counts[(qq % m, pp % m)] += 1
    tot = sum(counts.values()); exp = tot / len(prim_pairs)
    chi2 = sum((c - exp) ** 2 / exp for c in counts.values())
    worst = max(counts.items(), key=lambda kv: abs(kv[1] - exp))
    print(f"   m={m}: {len(prim_pairs)} primitive pairs (Q,P) mod m, {tot} convergents, expected {exp:.1f} each, chi^2 = {chi2:.1f} (df {len(prim_pairs)-1}), extreme {worst}")
    if m == 2:
        names = {(1, 0): "tau", (0, 1): "omega", (1, 1): "C"}
        print("      m = 2 character of the first 40 convergents (n: character):", " ".join(f"{n}:{names[(qq%2, pp%2)]}" for n, (pp, qq) in enumerate(conv[1:41], start=1)))
        print("      shares:", {names[k]: f"{c/tot:.4f}" for k, c in counts.items()}, "(equidistribution predicts 1/3 each)")
        runs = []
        cur = None; ln = 0
        for (pp, qq) in conv[1:]:
            nm = names[(qq % 2, pp % 2)]
            if nm == cur: ln += 1
            else:
                if cur is not None: runs.append(ln)
                cur, ln = nm, 1
        runs.append(ln)
        print(f"      longest run of one character along the geodesic: {max(runs)}; mean run {np.mean(runs):.3f} (i.i.d. thirds predict 1.5); consecutive convergents never share a character? {max(runs) == 1}")
gk = {k: sum(1 for x in a[1:] if x == k) / (len(a) - 1) for k in (1, 2, 3, 4)}
print(f"   Gauss-Kuzmin: frequencies of a_n = 1,2,3,4: {[round(gk[k],4) for k in (1,2,3,4)]} vs log2(1 + 1/(k(k+2))) = {[round(float(np.log2(1 + 1/(k*(k+2)))),4) for k in (1,2,3,4)]}")

print("\n# Part 4: Paley-7 against the Q resonator (NUMEROLOGY check)")
A = np.zeros((7, 7))
for x in range(7):
    for d in (1, 2, 4):
        A[x, (x + d) % 7] = 1
ev = np.linalg.eigvals(A)
print(f"   Paley tournament P7 (x -> x + {{1,2,4}}): eigenvalues {sorted([complex(round(t.real,6), round(t.imag,6)) for t in ev], key=lambda t: (t.real, t.imag))}")
lam = sp.symbols("lam"); z = sp.symbols("z")
print(f"   minimal polynomial of (-1 + sqrt(-7))/2: {sp.factor(sp.minimal_polynomial((-1 + sp.sqrt(-7))/2, lam))}; the Q linearisation d(n) = -d(n-1)/2 - d(n-2)/2 has characteristic polynomial 2 z^2 + z + 1 with roots {sp.solve(2*z**2 + z + 1, z)}")
print(f"   reciprocal: 2 z^2 + z + 1 at z = 1/lambda for the roots lambda of lambda^2 + lambda + 2: {[sp.simplify(2*(1/r)**2 + (1/r) + 1) for r in sp.solve(lam**2 + lam + 2, lam)]}")
print("   doubly regular tournaments on n = 4t+3 vertices have eigenvalues (-1 +- i sqrt n)/2, roots of lambda^2 + lambda + (n+1)/4: the coincidence with the two-branch, weight-1/2, lags-1-and-2 polynomial needs (n+1)/4 = 2, i.e. n = 7 only.")
print(f"   N((-1+sqrt(-7))/2) = 2: a prime above 2 in Q(sqrt(-7)) (2 splits since -7 = 1 mod 8). Ramanujan-Nagell 2^n - 7 = x^2, n <= 300: {[n for n in range(3, 301) if math.isqrt(2**n - 7) ** 2 == 2**n - 7]}")

print("\n# Part 5: the 2-adic (de Bruijn) part alone: spectrum of E_j v(r) = [r even] (1/2)(v(r/2) + v(r/2 + 2^(j-1)))")
for jj in range(1, 8):
    n = 2 ** jj; E = np.zeros((n, n)); half = n // 2
    for r in range(0, n, 2):
        r2 = (r // 2) % half
        E[r, r2] += 0.5; E[r, (r2 + half) % n] += 0.5
    ev = np.linalg.eigvals(E)
    nz = [complex(round(t.real, 6), round(t.imag, 6)) for t in ev if abs(t) > 1e-9]
    ranks = [np.linalg.matrix_rank(np.linalg.matrix_power(E, q)) for q in range(1, jj + 2)]
    print(f"   j={jj}: nonzero eigenvalues {nz}, zeros {sum(1 for t in ev if abs(t) <= 1e-9)} of {n}; rank(E^q) for q = 1..{jj+1}: {ranks}")
print("done")
