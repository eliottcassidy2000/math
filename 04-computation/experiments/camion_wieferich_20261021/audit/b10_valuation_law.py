# Audit b10: THM-4750 (4e) beyond mod p^2, with exact integers.
# For small p and every zero class m (base in [3, p+1]) the exact valuations v_p(D(m + l(p-1))), l = 0..L, are compared
# with the Weierstrass picture: off the root discs v = 1; on the disc of a simple root (l = l0 mod p) the law
#     v_p(D(m + l(p-1))) = 1 + v_p(l - ell)      (valid while 1 + v_p(l - ell) < ceil(n/2), n = m + l(p-1))
# must hold for ONE p-adic number ell.  ell is reconstructed digit by digit from the data and the law is then checked on
# every l in the disc.  Also records the even-n boundary cases 1 + v = n/2 (where THM-4750's stated range
# "below (n+1)/2" would claim v(D(n)) = n/2 although only v(D(n)) >= n/2 is implied).
# Usage: python3 b10_valuation_law.py
from math import comb
NEX = 800
W = [1] + [0] * NEX; row = [1]
for n in range(1, NEX + 1):
    row = [1] + [row[i - 1] + row[i] for i in range(1, n)] + [1]
    W[n] = -sum((row[k] * W[n - k]) << (k * (k - 1) // 2) for k in range(1, n + 1))
s = [0] + [-W[n] for n in range(1, NEX + 1)]
def c(n): return (1 << ((n - 1) // 2)) if n % 2 else 0
D = [0] + [s[n] - c(n) for n in range(1, NEX + 1)]
def vp(x, p):
    if x == 0: return 10 ** 6
    v = 0
    while x % p == 0: x //= p; v += 1
    return v
for p in (5, 7, 11, 13):
    eps = 1 if p % 8 in (1, 7) else -1; q = p * p
    L = (NEX - p - 1) // (p - 1)
    print(f"p={p}: l = 0..{L}")
    for m in range(3, p + 2):
        if D[m] % p: continue
        A = [eps ** l * D[m + l * (p - 1)] for l in range(3)]
        c0 = (A[0] // p) % p; c1 = ((A[1] - A[0]) // p) % p; c2 = ((A[2] - 2 * A[1] + A[0]) // p) % p
        g = lambda l: (c0 + c1 * l + c2 * (l * (l - 1) // 2)) % p
        roots = [l for l in range(p) if g(l) == 0]
        # derivative of g at a root (simple iff nonzero): g'(l) = c1 + c2 (2l - 1)/2
        inv2 = (p + 1) // 2
        simple = {r: (c1 + c2 * (2 * r - 1) * inv2) % p != 0 for r in roots}
        ls = list(range(L + 1)); vals = {l: vp(D[m + l * (p - 1)], p) for l in ls}
        off = [l for l in ls if g(l % p) != 0 and vals[l] != 1]
        msg = [f"  class base m={m:2d}: g roots mod p {roots} (simple {[simple[r] for r in roots]}); off-disc v!=1 at {off}"]
        for r in roots:
            if not simple[r]: continue
            disc = [l for l in ls if l % p == r]
            # valid l: those with v(D) < ceil(n/2) (then v(D) = v(F)); reconstruct ell greedily from valid data
            valid = [l for l in disc if vals[l] < -(-(m + l * (p - 1)) // 2)]
            ell, k = r, 1
            while True:
                cand = [x for x in (ell + t * p ** k for t in range(p))
                        if all((vals[l] >= 2 + k) == ((l - x) % p ** (k + 1) == 0) for l in valid if vals[l] >= 1 + k)]
                if len(cand) != 1 or k > 6: break
                ell, k = cand[0], k + 1
            bad = [(l, vals[l]) for l in valid if (vp(l - ell, p) < k and vals[l] != 1 + vp(l - ell, p))
                   or (vp(l - ell, p) >= k and vals[l] < 1 + k)]
            # the law must also be consistent with the invalid points: there v(D) >= ceil(n/2) <= 1 + v(l - ell)
            inv = [(l, vals[l]) for l in disc if l not in valid and vp(l - ell, p) < k and 1 + vp(l - ell, p) < -(-(m + l * (p - 1)) // 2)]
            evenb = [l for l in disc if (m + l * (p - 1)) % 2 == 0 and 1 + min(vp(l - ell, p), k) == (m + l * (p - 1)) // 2]
            msg.append(f"      root {r}: ell = {ell} mod {p}^{k}; law violations on valid points {bad}; contradictions at "
                       f"high-valuation points {inv}; even-n boundary points {evenb}")
        print("\n".join(msg))
