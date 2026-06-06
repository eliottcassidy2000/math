"""SIGNED LRC -- fast confirmation of the two S703 conjectures, pushed further.  monad-explorer-S703c.

(A) AP_n folded sign-orbit = 2^{n-2}  iff  C=2n-1 is prime.   Push to n<=24.
(B) Unfolded cut-spectrum is faithful for EVERY config.        Stress-test n=7 + adversarial sets.

Speedup: a clock value lies in [0, (C-1)//2], a SMALL range, so encode a multiset by its
count-vector (tuple of counts), far faster to build/hash than sorting C(n-1,2) integers.
"""
from itertools import product, combinations
from math import gcd


def is_prime(m):
    if m < 2:
        return False
    d = 2
    while d * d <= m:
        if m % d == 0:
            return False
        d += 1
    return True


def folded_orbit_size(V, C):
    n1 = len(V)
    pairs = [(i, j) for i in range(n1) for j in range(i + 1, n1)]
    half = C // 2
    seen = set()
    for bits in product((0, 1), repeat=n1 - 1):
        col = (0,) + bits
        cnt = [0] * (half + 1)
        for (i, j) in pairs:
            if col[i] != col[j]:
                f = (V[i] + V[j]) % C
                cnt[min(f, C - f)] += 1
            else:
                cnt[abs(V[i] - V[j])] += 1
        seen.add(tuple(cnt))
    return len(seen)


def unfolded_faithful(V):
    """True iff the unfolded cut-spectrum is injective over cuts (up to global swap)."""
    n1 = len(V)
    pairs = [(i, j) for i in range(n1) for j in range(i + 1, n1)]
    maxv = max(V[i] + V[j] for i, j in pairs)
    seen = set()
    total = 0
    for bits in product((0, 1), repeat=n1 - 1):
        col = (0,) + bits
        cnt = [0] * (maxv + 1)
        for (i, j) in pairs:
            cnt[(V[i] + V[j]) if col[i] != col[j] else abs(V[i] - V[j])] += 1
        seen.add(tuple(cnt))
        total += 1
    return len(seen) == total


def main():
    print("=" * 78)
    print("(A) AP_n folded orbit vs 2^{n-2}; conjecture: equal IFF C=2n-1 prime.  n=3..24")
    print("=" * 78)
    print(f"{'n':>3} {'C':>4} {'prime?':>7} {'2^(n-2)':>9} {'folded':>9} {'coll':>6} {'OK?':>4}")
    allok = True
    for n in range(3, 23):
        C = 2 * n - 1
        if n - 2 > 20:   # cap work at 2^20 cuts
            print(f"{n:>3} {C:>4} (skipped: 2^{n-2} too large)")
            continue
        V = list(range(1, n))
        fo = folded_orbit_size(V, C)
        tot = 2 ** (n - 2)
        p = is_prime(C)
        ok = (fo == tot) == p
        allok &= ok
        print(f"{n:>3} {C:>4} {str(p):>7} {tot:>9} {fo:>9} {tot-fo:>6} {'yes' if ok else 'NO!!':>4}")
    print(f"\nConjecture (A) holds for ALL tested n: {allok}")

    print("\n" + "=" * 78)
    print("(B) Unfolded faithfulness stress test: n=7 all configs (speeds<=14) + adversarial sets")
    print("=" * 78)
    # n=7 exhaustive, speeds up to 14
    cnt_faithful = cnt_bad = 0
    bad_examples = []
    for V in combinations(range(1, 15), 6):
        if unfolded_faithful(list(V)):
            cnt_faithful += 1
        else:
            cnt_bad += 1
            if len(bad_examples) < 10:
                bad_examples.append(V)
    print(f" n=7 (speeds in [1,14]): faithful={cnt_faithful}  NON-faithful={cnt_bad}")
    if bad_examples:
        print("  first non-faithful configs:", bad_examples)
    # adversarial: symmetric / AP-rich / Sidon-violating sets at several n
    print("\n adversarial individual configs (symmetric, AP-rich, high additive energy):")
    adv = [
        [1, 2, 4, 5], [1, 2, 3, 5, 6, 7], [1, 2, 4, 8, 16], [1, 3, 5, 7, 9],
        [2, 3, 5, 8, 13, 21], [1, 2, 3, 4, 6, 8, 12], [1, 4, 9, 16, 25],
        [1, 2, 5, 6, 9, 10], [1, 5, 6, 10, 11, 15],  # AP-of-blocks (sum-symmetric)
        [1, 2, 3, 6, 7, 8], [1, 2, 6, 9, 11, 12], list(range(1, 9)),
    ]
    for V in adv:
        f = unfolded_faithful(V)
        print(f"   {str(V):<30} unfolded-faithful={f}")


if __name__ == "__main__":
    main()
