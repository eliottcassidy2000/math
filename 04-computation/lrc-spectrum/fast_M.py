"""
Fast EXACT M(S). Integer arithmetic, dedup of candidate t-values.

M(S) = max over candidate t = p/Q of min_v ||v p / Q||, where Q ranges over the
set of binding denominators D = { v_i +- v_j, 2 v_i } (q | one of these), p=1..Q-1.

We represent each candidate t = p/Q. For a fixed Q, ||v p/Q|| = (1/Q)*min(r, Q-r)
with r = (v*p) mod Q. So min_v ||v t|| = (1/Q) * min_v g(v) where g = min(r,Q-r).
We track the best as an exact Fraction by comparing cross-multiplied integers.

To dedup, we reduce each (p,Q) to lowest terms and keep a visited set of (p',Q').
But different Q give different grids; we just compare values exactly.
"""
from fractions import Fraction
from math import gcd


def candidate_Qs(S):
    D = set()
    n = len(S)
    for i in range(n):
        vi = S[i]
        D.add(2 * vi)
        for j in range(i + 1, n):
            vj = S[j]
            D.add(vi + vj)
            d = vi - vj
            if d:
                D.add(abs(d))
    D.discard(0)
    D.discard(1)
    return D


def M_exact_fast(S):
    S = sorted(set(S))
    Qs = candidate_Qs(S)
    best_num, best_den = 0, 1  # best value = best_num/best_den
    best_t = None
    for Q in Qs:
        for p in range(1, Q):
            # min over v of min((v*p)%Q, Q-(v*p)%Q)
            mn = Q  # bigger than any
            for v in S:
                r = (v * p) % Q
                gg = r if r <= Q - r else Q - r
                if gg < mn:
                    mn = gg
                    if mn == 0:
                        break
            # value = mn/Q ; compare to best_num/best_den
            if mn * best_den > best_num * Q:
                best_num, best_den = mn, Q
                best_t = Fraction(p, Q)
    M = Fraction(best_num, best_den)
    return M, best_t


if __name__ == "__main__":
    # validate against AP and skip
    for k in [3, 4, 5, 6, 7]:
        S = list(range(1, k + 1))
        print("AP", k, M_exact_fast(S))
    # cross-check with reference on a few sets
    import sys
    sys.path.insert(0, ".")
    from lrc_maxmin import M_exact
    import random
    random.seed(1)
    bad = 0
    for _ in range(200):
        k = random.randint(3, 7)
        S = random.sample(range(1, 25), k)
        if gcd(*S) if len(S) > 1 else S[0]:
            pass
        a = M_exact_fast(S)[0]
        b = M_exact(S)[0]
        if a != b:
            bad += 1
            print("MISMATCH", S, a, b)
    print("validation mismatches:", bad)
