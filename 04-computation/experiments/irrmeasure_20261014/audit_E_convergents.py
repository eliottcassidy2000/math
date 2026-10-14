"""Audit E (2026-10-14) of THM-4633: independent recomputation of the continued fraction of theta = log2 3,
its convergents AND intermediate fractions (semiconvergents) above theta, the Eliahou budgets
m_max(p/q) = q / (3 (p ln2 - q ln3)), the exclusion status at 2^68, 3*2^69 and 2^71, the Legendre threshold,
the error term of the convergent asymptotics claimed in THM-4633 statement 2, and the mu_0 table of statement 6.
Two independent routes: mpmath at 300 digits, and sympy's exact continued_fraction_iterator on log(3)/log(2)."""
from mpmath import mp, mpf, log, floor, sqrt, nstr
mp.dps = 300
theta = log(3) / log(2)
ln2, ln3 = log(2), log(3)

# --- continued fraction, two ways ---
x = theta
a = []
for i in range(80):
    ai = int(floor(x)); a.append(ai); x = 1 / (x - ai)
import sympy
a_sym = []
it = sympy.continued_fraction_iterator(sympy.log(3) / sympy.log(2))
for i, ai in enumerate(it):
    a_sym.append(int(ai))
    if i >= 79: break
print("partial quotients agree (mpmath 300 dps vs sympy exact):", a == a_sym)
print("a_0..a_39 =", a[:40])
print("max partial quotient among a_1..a_40 =", max(a[1:41]), " (a_14 =", a[14], ")")
print("max partial quotient among a_0..a_39 =", max(a[:40]))

# convergents
p = [a[0], a[0] * a[1] + 1]; q = [1, a[1]]
for j in range(2, len(a)):
    p.append(a[j] * p[-1] + p[-2]); q.append(a[j] * q[-1] + q[-2])

def budget(P, Q):
    d = P * ln2 - Q * ln3
    return d, (Q / (3 * d) if d > 0 else None)

M68, M369, M71 = mpf(2) ** 68, 3 * mpf(2) ** 69, mpf(2) ** 71
print("\nM0 = 2^68 = %s, 3*2^69 = %s, 2^71 = %s" % (nstr(M68, 6), nstr(M369, 6), nstr(M71, 6)))
print("\n--- convergents p_j/q_j, sign, budget, status (budget < M0 => excluded) ---")
print("%-3s %-16s %-16s %-6s %-14s %-12s %-6s %-6s %-6s %s" % ("j", "p_j", "q_j", "sign", "delta", "m_max", "2^68", "3*2^69", "2^71", "a_(j+1) q_j^2/(3ln2) and ratio"))
for j in range(0, 26):
    d, b = budget(p[j], q[j])
    if d <= 0:
        print("%-3d %-16d %-16d below" % (j, p[j], q[j])); continue
    approx = a[j + 1] * mpf(q[j]) ** 2 / (3 * ln2)
    st = lambda M: "excl" if b < M else "ALIVE"
    print("%-3d %-16d %-16d above  %-14s %-12s %-6s %-6s %-6s %s  ratio=%s  a_(j+1)=%d" % (
        j, p[j], q[j], nstr(d, 8), nstr(b, 6), st(M68), st(M369), st(M71), nstr(approx, 6), nstr(b / approx, 6), a[j + 1]))

# exact form of the budget at a convergent: |q_j theta - p_j| = 1/(theta_{j+2}^{-1} q_j + q_{j+1})
print("\nCheck of the exact convergent formula m_max = q_j (q_(j+1) + q_j/theta_(j+2)) / (3 ln 2), theta_(j+2) = complete quotient:")
for j in (1, 3, 21, 23):
    # complete quotient theta_{j+2} = [a_{j+2}; a_{j+3}, ...]
    t = theta
    for i in range(j + 2):
        t = 1 / (t - floor(t))
    exact = mpf(q[j]) * (q[j + 1] + q[j] / t) / (3 * ln2)
    d, b = budget(p[j], q[j])
    print("  j=%d: m_max=%s, formula=%s, rel. diff=%s; relative correction (q_(j-1)/q_j + 1/theta_(j+2))/a_(j+1) = %s" % (
        j, nstr(b, 12), nstr(exact, 12), nstr(b / exact - 1, 3), nstr((mpf(q[j - 1]) / q[j] + 1 / t) / a[j + 1], 4)))

# --- intermediate fractions (semiconvergents) above theta ---
print("\n--- intermediate fractions (p_(j-1) + i p_j)/(q_(j-1) + i q_j), 1 <= i < a_(j+1), lying ABOVE theta ---")
print("(these are exactly the other best approximations of the first kind; Eliahou's method must treat them)")
sc = []
for j in range(1, 26):
    for i in range(1, a[j + 1]):
        P, Q = p[j - 1] + i * p[j], q[j - 1] + i * q[j]
        d, b = budget(P, Q)
        if d > 0:
            sc.append((Q, P, j, i, d, b))
sc.sort()
for Q, P, j, i, d, b in sc:
    if Q > 10 ** 13: break
    st = lambda M: "excl" if b < M else "ALIVE"
    print("  K/k = %d/%d  (between p_%d/q_%d and p_%d/q_%d, i=%d)  delta=%s  m_max=%s  2^68:%s  3*2^69:%s  2^71:%s" % (
        P, Q, j - 1, j - 1, j + 1, j + 1, i, nstr(d, 8), nstr(b, 6), st(M68), st(M369), st(M71)))

# the smallest denominator k such that some K/k lies in (theta, theta + 1/(3 M0 ln 2)): brute force over
# best approximations from above (convergents and intermediate fractions), for each M0
def first_survivor(M0):
    eps = 1 / (3 * M0 * ln2)
    cands = [(q[j], p[j]) for j in range(0, 40) if p[j] * ln2 - q[j] * ln3 > 0] + [(Q, P) for Q, P, *_ in sc]
    cands.sort()
    for Q, P in cands:
        if mpf(P) / Q - theta < eps:
            return P, Q
for name, M0 in (("2^68", M68), ("3*2^69", M369), ("2^71", M71), ("2^72", mpf(2)**72), ("2^73", mpf(2)**73)):
    P, Q = first_survivor(M0)
    print("first fraction K/k > theta with K/k - theta < 1/(3 M0 ln2), M0=%s: K=%d, k=%d (K+k=%d)" % (name, P, Q, P + Q))

# Legendre threshold and the self-contained bound
kL = sqrt(mpf(3) / 2 * M71 * ln2)
print("\nLegendre threshold sqrt((3/2) 2^71 ln2) = %s" % nstr(kL, 10))
print("convergents above theta with q_j <= threshold: j =", [j for j in range(0, 30) if q[j] <= kL and p[j] * ln2 - q[j] * ln3 > 0])
print("max budget among them = %s < 2^71 = %s ?" % (nstr(max(budget(p[j], q[j])[1] for j in range(0, 30) if q[j] <= kL and p[j] * ln2 - q[j] * ln3 > 0), 6), nstr(M71, 6)))
print("q_23 =", q[23], "> threshold:", q[23] > kL, "; semiconvergent k = q_21+q_22 =", q[21] + q[22], "> threshold:", q[21] + q[22] > kL)

# mu_0 table of statement 6 (c = 1)
print("\nmu_0 table, k > (3 c M0 ln2)^(1/mu_0) with c = 1, M0 = 2^71:")
for mu in (2, 3, 5.117, 8.616, 4.117 + 1):
    print("  mu_0 = %s: k > %s" % (mu, nstr((3 * M71 * ln2) ** (1 / mpf(mu)), 5)))
print("3*2^71*ln2 =", nstr(3 * M71 * ln2, 8))

# cross-check the table's q_23 = 137528045312 and the Wikipedia-type numbers
print("\np_21+p_22 =", p[21] + p[22], " q_21+q_22 =", q[21] + q[22], " sum =", p[21] + p[22] + q[21] + q[22])
print("p_23 =", p[23], " q_23 =", q[23], " sum =", p[23] + q[23])
print("Eliahou 1993 generators: p_13 =", p[13], " p_15 =", p[15], " p_16 =", p[16])
