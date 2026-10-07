"""A2: exact affine certificates on classes n mod 2^K (Python, Fractions), independent of the C residue search.
For each class: orbit points x_s = (3^a n + c)/2^s; backward words from x_s with b <= a, i < s, explored on exact affine values.
Returns min backward length (depth) of a certificate, and its threshold n* (m(n) < n for n > n*)."""
from fractions import Fraction as Fr
import sys, math

def prefix_affine(rho, K):
    """list of (s, a, c) with x_s = (3^a n + c)/2^s for n = rho mod 2^K"""
    out = [(0, 0, 0)]; a = 0; c = 0; x = rho
    for s in range(K):
        if x % 2 == 0: x //= 2
        else: x = (3 * x + 1) // 2; c = 3 * c + 2 ** s; a += 1
        out.append((s + 1, a, c))
    return out

def certificates(rho, K, maxdepth=None, want_all=False):
    """search branches from each x_s, s<=K; return list of (s, b, i, alpha, gamma, threshold)"""
    found = []
    for (s, a, c) in prefix_affine(rho, K):
        if s == 0: continue
        if 3 ** a < 2 ** s:   # descent
            alpha = Fr(3 ** a, 2 ** s); gamma = Fr(c, 2 ** s)
            found.append((s, 0, 0, alpha, gamma, gamma / (1 - alpha)))
            if not want_all: return found
            continue
        # backward DFS on exact affine value v = A n + B (A, B Fractions); validity: v mod 3 must be decided:
        # v mod 3 is decided iff 3 | (A as 3-adic coefficient)... we track b and require b < a for odd steps.
        stack = [(Fr(3 ** a, 2 ** s), Fr(c, 2 ** s), 0, 0)]
        while stack:
            A, B, b, E = stack.pop()
            i = b + E
            if i > 0:
                if A < 1:
                    found.append((s, b, i, A, B, B / (1 - A) if B > 0 else Fr(-10**9)))
                    if not want_all: return found
            if maxdepth is not None and i >= maxdepth: continue
            if E + 1 < s - a:
                stack.append((2 * A, 2 * B, b, E + 1))
            if b < a:
                # value mod 3: A n + B with A = 3^(a-b) * 2^k, so v = B mod 3 (B has denominator a power of 2)
                Bn, Bd = B.numerator, B.denominator
                vmod3 = (Bn * pow(Bd, -1, 3)) % 3
                if vmod3 == 2:
                    stack.append(((2 * A) / 3, (2 * B - 1) / 3, b + 1, E))
    return found

if __name__ == "__main__":
    K = 10
    for rho in [411, 487, 539, 615, 799, 879]:
        f1 = certificates(rho, K, maxdepth=1)
        fa = certificates(rho, K)
        print(rho, "depth-1 cert:", [(x[0], x[1], x[2], str(x[3]), str(x[4])) for x in f1][:1],
              "| any cert:", [(x[0], x[1], x[2], str(x[3]), str(x[4]), float(x[5])) for x in fa][:1])
    # explicit claim: m = (27n+7)/32 for 539 and (27n+3)/32 for 615; check on many lifts
    def T(x): return x // 2 if x % 2 == 0 else (3 * x + 1) // 2
    def Tn(x, k):
        for _ in range(k): x = T(x)
        return x
    import random
    ok = True
    for rho, cst in [(539, 7), (615, 3)]:
        for j in range(500):
            n = rho + 1024 * random.randint(0, 10**30)
            assert (27 * n + cst) % 32 == 0
            m = (27 * n + cst) // 32
            if not (m < n and Tn(m, 5) == Tn(n, 10)): ok = False
    print("539/615 certificates verified on 500 random lifts each:", ok)
    print("455 orbit:", [Tn(455, k) for k in range(6)], " T^10(539) =", Tn(539, 10))
