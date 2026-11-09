"""HYP-9590 (C) -> THM-5020 (6). Written after the pre-audit checkpoint 2fe6b181c5, while the audit was running;
the audit found the same argument independently (its m5). No palindrome of length L = n - 2 has P(A) = P(B) = k at
n = 0 or 1 (mod 3).
- n = 3k: P(A) = k forces A to contain the unique maximum packing {0, 3, ..., 3k - 3}, so P(B) <= k - 1 (any word).
- n = 3k + 1 (L = 3k - 1): the k-packings of [0, 3k - 2] are S_j = {0, 3, ..., 3(j - 1)} u {3j + 1, ..., 3k - 2},
  j = 0..k. The only disjoint pair is {S_0, S_k}, and the mirror x -> 3k - 2 - x swaps them. A symmetric colour class
  containing one contains both, so the other class has no k-packing.
The script checks (a) the list of k-packings for k <= 7, and (b) the statement for all palindromes up to n = 31 at
n = 0, 1 (mod 3). It also confirms that (k, k) does occur for palindromes at n = 2 (mod 3).
usage: python palindrome_kk_check.py NMAX
"""
import sys
import itertools


def P(xs):
    c, last = 0, -10
    for x in xs:
        if x - last >= 3:
            c += 1
            last = x
    return c


def main():
    NMAX = int(sys.argv[1])
    for k in range(1, 8):
        L = 3 * k - 1
        packs = [S for S in itertools.combinations(range(L), k) if all(S[i + 1] - S[i] >= 3 for i in range(k - 1))]
        expect = sorted(tuple(list(range(0, 3 * j, 3)) + list(range(3 * j + 1, L, 3))) for j in range(k + 1))
        print(f"k = {k}, L = {L}: k-packings {len(packs)}, equal to the S_j list: {sorted(packs) == expect}")
    for n in range(5, NMAX + 1):
        L, k = n - 2, n // 3
        h = L // 2
        mids = ('',) if L % 2 == 0 else ('0', '1')
        hit = 0
        for code in range(1 << h):
            half = format(code, f'0{h}b') if h else ''
            for mid in mids:
                w = half + mid + half[::-1]
                A = [i for i in range(L) if w[i] == '1']
                B = [i for i in range(L) if w[i] == '0']
                hit += (P(A) == k and P(B) == k)
        print(f"n = {n} (n mod 3 = {n % 3}), k = {k}: palindromes with P(A) = P(B) = k: {hit}")


if __name__ == "__main__":
    main()
