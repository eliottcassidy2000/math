"""Audit item 1 (independent code): the Ulam / Haar-conditional matrix of S_p(x) = (p x + 1)/2^v on odd residues mod 2^K.

Method (independent of the author's construction):
  * Haar row from the definition: for x = a + 2^K t, t Haar on Z_2, compute the law of S(x) mod 2^K by brute-force
    lift counting over t mod 2^M, separating the 'determined' lifts (v_2(px+1) + K <= K + M, so S(x) mod 2^K is fixed
    by t mod 2^M) from the undetermined mass, whose total is reported (it tends to 0 like 2^-(M-K+1) for a*).
  * Closed-form row (the note's rule): uniform on the 2^v residues = (pa+1)/2^v mod 2^(K-v) if v < K, uniform for a*.
  * Compare exactly (Fractions), check double stochasticity, nnz = K 2^(K-1), the counting matrix C = 2^(K-1) P,
    and the claim 'C(a,b) = #{t mod 2^(K-1): S(a+2^K t) = b mod 2^K}' for non-special rows using integer reps.
"""
from fractions import Fraction
import sys


def v2(n):
    n = abs(n)
    assert n != 0
    return (n & -n).bit_length() - 1


def closed_row(a, K, p):
    """the note's rule, written from scratch"""
    mod = 1 << K
    t = (p * a + 1) % mod
    odd = [b for b in range(1, mod, 2)]
    if t == 0:  # v >= K : the special residue a* = -1/p mod 2^K
        return {b: Fraction(1, len(odd)) for b in odd}, True
    v = v2(t)
    r = (t >> v) % (1 << (K - v))
    sup = [b for b in odd if (b - r) % (1 << (K - v)) == 0]
    assert len(sup) == (1 << v)
    return {b: Fraction(1, len(sup)) for b in sup}, False


def lift_row(a, K, p, M):
    """Haar law of S(x) mod 2^K for x in a + 2^K Z_2, from t mod 2^M.  Returns (determined counts, undetermined count)."""
    mod = 1 << K
    cnt = {}
    undet = 0
    for t in range(1 << M):
        x = a + (t << K)                  # x mod 2^(K+M) is exact; higher digits are Haar
        y = p * x + 1                     # y mod 2^(K+M) exact (p odd unit)
        ymod = y % (1 << (K + M))
        if ymod == 0:
            undet += 1
            continue
        v = v2(ymod)
        # S(x) = y / 2^v ; its residue mod 2^K is determined iff v + K <= K + M  (we know y mod 2^(K+M))
        if v + K <= K + M:
            b = (ymod >> v) % mod
            cnt[b] = cnt.get(b, 0) + 1
        else:
            undet += 1
    return cnt, undet


def main():
    out = []
    for p in (3, 5, 7, 1, -1, 9, 15):
        for K in range(2, 8):
            mod = 1 << K
            odd = list(range(1, mod, 2))
            rows = {}
            special = []
            for a in odd:
                r, sp = closed_row(a, K, p)
                rows[a] = r
                if sp:
                    special.append(a)
            assert len(special) == 1 and (p * special[0] + 1) % mod == 0
            # double stochastic + nnz
            col = {b: Fraction(0) for b in odd}
            nnz = 0
            for a in odd:
                assert sum(rows[a].values()) == 1
                for b, w in rows[a].items():
                    col[b] += w
                    nnz += 1
            assert all(c == 1 for c in col.values()), (p, K)
            assert nnz == K * (1 << (K - 1)), (p, K, nnz)
            # counting matrix integrality
            assert all((w * (1 << (K - 1))).denominator == 1 for a in odd for w in rows[a].values())
            # Haar lift counts with M large: non-special rows exactly, special row: determined part exactly uniform
            M = K + 8
            worst_undet = 0
            for a in odd:
                cnt, undet = lift_row(a, K, p, M)
                tot = 1 << M
                if a == special[0]:
                    # determined part: each odd b should get exactly the same count
                    vals = set(cnt.get(b, 0) for b in odd)
                    assert len(vals) == 1, (p, K, a, vals)
                    worst_undet = undet / tot
                    # undetermined mass is 2^-(M-K+1) exactly (v >= M+1 ... ) : check
                    assert undet * (1 << (M - K + 1)) == tot, (undet, tot)
                else:
                    assert undet == 0
                    got = {b: Fraction(c, tot) for b, c in cnt.items()}
                    assert got == rows[a], (p, K, a)
            # integer-representative counting claim for non-special rows: t mod 2^(K-1)
            for a in odd:
                if a == special[0]:
                    continue
                cnt = {}
                for t in range(1 << (K - 1)):
                    x = a + (t << K)
                    y = p * x + 1
                    s = y >> v2(y) if y != 0 else None
                    b = s % mod
                    cnt[b] = cnt.get(b, 0) + 1
                assert {b: Fraction(c, 1 << (K - 1)) for b, c in cnt.items()} == rows[a], (p, K, a)
            out.append(f'p = {p:3d}, K = {K}: Haar rows (lift counting, M = K + 8) == closed-form rows; doubly stochastic; '
                       f'nnz = {nnz} = K 2^(K-1); a* = {special[0]} (= -1/p mod 2^K); special row: determined lifts exactly '
                       f'uniform, undetermined mass {worst_undet:.2e} = 2^-(M-K+1); integer-rep count over t mod 2^(K-1) '
                       f'exact on the other rows')
    # the special row with integer representatives is NOT uniform (why the note uses the Haar limit)
    K, p = 4, 3
    mod = 1 << K
    a = [x for x in range(1, mod, 2) if (p * x + 1) % mod == 0][0]
    cnt = {}
    for t in range(1 << (K - 1)):
        x = a + (t << K)
        y = p * x + 1
        b = (y >> v2(y)) % mod
        cnt[b] = cnt.get(b, 0) + 1
    out.append(f'K = 4, p = 3, a* = {a}: integer-representative counts over t mod 2^(K-1) = {dict(sorted(cnt.items()))} '
               f'(not uniform: the a* row must be the Haar limit, as the note says)')
    print('\n'.join(out))
    print('A1 ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
