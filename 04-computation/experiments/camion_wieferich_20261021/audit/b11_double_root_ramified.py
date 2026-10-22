# Audit b11: counterexamples to the "split or inert" dichotomy of THM-4750 (4) (title and (4e)) and to
# "the triangle class m = 3 is always split" (title), inside the stated range p >= 5.
# Double-root zero classes (chi(Delta) = 0, tag "D" in second_order_3_4000.out): p = 5, triangle class (base m = 3, root
# l = 0) and p = 7, class 2 (base m = 8, root l = 6).  With exact integers, v_p(D(m + l(p-1))) is computed on the whole
# root disc l = r (mod p) for every l with m + l(p-1) <= 800 (all far inside the valid range v < ceil(n/2)).
# By THM-4750 (4e), F_m/p = P U with U a unit on Z_p and P monic quadratic, P = (l - r)^2 mod p, so for valid points
# v_p(D) = 1 + v_p(P(l)).  If the zeros were in Z_p (split), some l in the disc (all residues mod p^2 are sampled)
# would have v_p(D) >= 1 + 2 + 1 = 4; an unramified pair r + p*(conjugates) would give v_p(D) >= 3 on the whole disc.
# A constant valuation 2 means v(P(r + pt)) = 1 for all t: P(r + X) = X^2 + bX + c with v(c) = 1, v(b) >= 1 is
# Eisenstein, so the two zeros are a RAMIFIED conjugate pair (valuation 1/2 from r): neither split nor inert.
# Usage: python3 b11_double_root_ramified.py
NEX = 800
W = [1] + [0] * NEX; row = [1]
for n in range(1, NEX + 1):
    row = [1] + [row[i - 1] + row[i] for i in range(1, n)] + [1]
    W[n] = -sum((row[k] * W[n - k]) << (k * (k - 1) // 2) for k in range(1, n + 1))
s = [0] + [-W[n] for n in range(1, NEX + 1)]
def c(n): return (1 << ((n - 1) // 2)) if n % 2 else 0
def vp(x, p):
    if x == 0: return 99
    v = 0
    while x % p == 0: x //= p; v += 1
    return v
for p, m, r in ((5, 3, 0), (7, 8, 6)):
    eps = 1 if p % 8 in (1, 7) else -1
    A = [eps ** l * (s[m + l * (p - 1)] - c(m + l * (p - 1))) for l in range(3)]
    c0 = (A[0] // p) % p; c1 = ((A[1] - A[0]) // p) % p; c2 = ((A[2] - 2 * A[1] + A[0]) // p) % p
    disc_ = ((c1 - c2 * (p + 1) // 2) ** 2 - 2 * c0 * c2) % p
    pts = [(l, vp(s[m + l * (p - 1)] - c(m + l * (p - 1)), p), -(-(m + l * (p - 1)) // 2)) for l in range(r, (NEX - m) // (p - 1) + 1, p)]
    valid = [(l, v) for l, v, cap in pts if v < cap]
    res2 = sorted({l % (p * p) for l, v in valid})
    print(f"p={p} base m={m}: g = {c0} + {c1} l + {c2} C(l,2) mod {p}, discriminant {disc_} (double root at l={r}); "
          f"q_2 X(m-2) != 0 mod p: {c2 != 0}")
    print(f"   valid points on the disc: {len(valid)}, residues mod p^2 covered: {res2}; valuations: {sorted({v for l, v in valid})}")
    print(f"   => {'RAMIFIED pair (neither split nor inert)' if {v for l, v in valid} == {2} and len(res2) == p else 'inconclusive'}")
