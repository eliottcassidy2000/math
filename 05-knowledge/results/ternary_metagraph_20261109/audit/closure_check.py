"""Auditor's multiplicative-closure check of THM-4950 (3), (4) (own closure code, own piece list for orders <= 10).

Step 1: closure of the strong H_3 values of orders <= 10 (strong_values_le10.txt, from gen_classes.c + ext10.c).
Step 2: add S86's other pieces as the packet defines them (alpha_1 <= 80 inventory; one-apex words of length <= 19;
        Session 16's rays 3 * 2^s and 9 * 2^s - 2), computed here from the packet's formulas, and recompute m*.
Step 3: lane classification of the first missing irreducibles (own factorisation).
usage: python closure_check.py strong_values_le10.txt few_cycle_strong_80.txt M
"""
import sys


def closure(pieces, V):
    reach = bytearray(V + 1)
    reach[1] = 1
    for p in sorted(set(x for x in pieces if 2 <= x <= V)):
        # unbounded multiplicity: x * p for every reachable x, in increasing x
        for x in range(1, V // p + 1):
            if reach[x]:
                reach[x * p] = 1
    return reach


def one_apex(lmax, cap):
    """alpha_1 of the one-apex words 1 u 0 (length <= lmax), as defined in coverage_products.py: sum over i < j with
    w_i = 1, w_j = 0 of c(j - i), c(1) = 1, c(d) = 2^(d-2).  Recomputed by a running sum over the word."""
    vals = set()
    for L in range(2, lmax + 1):
        for u in range(1 << (L - 2)):
            w = [1] + [(u >> (L - 3 - i)) & 1 for i in range(L - 2)] + [0]
            s = 0
            for j in range(L):
                if w[j] == 0:
                    for i in range(j):
                        if w[i] == 1:
                            d = j - i
                            s += 1 if d == 1 else 1 << (d - 2)
                if s > cap:
                    break
            if s <= cap:
                vals.add(s)
    return vals


def factor(v):
    f, d = [], 2
    while d * d <= v:
        while v % d == 0:
            f.append(d)
            v //= d
        d += 1
    if v > 1:
        f.append(v)
    return f


def main():
    strong = [int(x) for x in open(sys.argv[1]).read().split()]
    M = int(sys.argv[3])
    V = 1 + 3 * M
    r1 = closure(strong, V)
    unc1 = [m for m in range(M + 1) if not r1[1 + 3 * m]]
    first1 = [m for m in unc1 if m not in (3, 10)]
    print(f"step 1 (strong orders <= 10 only): uncovered alpha-sums {unc1[:3]}...; least other than 3, 10: "
          f"m = {first1[0]} (value {1 + 3 * first1[0]})")
    # step 2: the other S86 pieces
    inv = []
    for line in open(sys.argv[2]):
        f = line.split()
        a5 = [int(x) for x in f[1:6]]
        H = int(f[6])
        rest = H - 1 - sum(c << (k + 1) for k, c in enumerate(a5))
        a6 = rest // 64
        vec = a5 + [a6]
        inv.append(1 + sum(c * 3 ** (k + 1) for k, c in enumerate(vec)))
    apex = one_apex(19, M)
    rays = {3 * 2 ** s for s in range(1, 20)} | {9 * 2 ** s - 2 for s in range(0, 20)}
    inter = [1 + 3 * a for a in apex | rays if 1 <= a <= M]
    allp = set(strong) | set(inv) | set(inter)
    r2 = closure(allp, V)
    unc2 = [m for m in range(M + 1) if not r2[1 + 3 * m]]
    first2 = [m for m in unc2 if m not in (3, 10)]
    print(f"step 2 (+ alpha_1 <= 80 inventory, {len(apex)} one-apex values, rays): least uncovered other than 3, 10: "
          f"m* = {first2[0]} (value {1 + 3 * first2[0]} = {factor(1 + 3 * first2[0])}); next {first2[1:8]}")
    extra = sorted(v for v in (set(inv) | set(inter)) if v not in set(strong) and v <= V)
    print(f"  piece values from the S86 pieces that are not strong order <= 10 values, below {V}: {len(extra)}; "
          f"smallest {extra[:10]}")
    gain = [m for m in unc1 if r2[1 + 3 * m] and m < first2[0]]
    print(f"  alpha-sums below m* covered only thanks to the S86 pieces: {len(gain)} {gain[:10]}")
    # step 3: is 29,038 or any value below it a piece? lanes of the first missing irreducibles
    print(f"  29038 a piece value: {29038 in allp}; 29035 covered: {bool(r2[29035])}; 9679 in one-apex set: "
          f"{9679 in apex}")


if __name__ == "__main__":
    main()
