#!/usr/bin/env python3
"""
Two triangular towers: additive tower, square tower, and their Pell crossover
families.  monad-explorer-2026-06-12-S712

User seed:
  A_n:  n+1 consecutive integers equal the following n consecutive integers.
  B_n:  n+1 consecutive squares equal the following n consecutive squares.

Write the underlying raw intervals as
  A_L(n) = [n^2, ..., n^2+n],          A_R(n) = [n^2+n+1, ..., n^2+2n]
  B_L(n) = [2n^2+n, ..., 2n^2+2n],    B_R(n) = [2n^2+2n+1, ..., 2n^2+3n]

This script does three things:
  1. verifies the exact closed forms for the tower sums;
  2. records the unique exact shared raw block between the two towers;
  3. classifies the main endpoint-crossover families, which land on Pell curves.

Assumption challenge:
  The useful objects are not the sums alone, but the interval endpoints
  n^2, n(n+1), n(2n+1), ... .  Passing to endpoints preserves overlap / truncation
  structure and reveals Pell families that are invisible at the sum level.

Tournament Analysis note:
  not used here.  The data are symmetric interval equalities / containments rather
  than an oriented pairwise relation, so there is no natural tournament quotient.
"""

from math import isqrt


def T(n):
    return n * (n + 1) // 2


def A_sum(n):
    return n * (n + 1) * (2 * n + 1) // 2


def B_sum(n):
    return n * (n + 1) * (2 * n + 1) * (12 * n * n + 12 * n + 1) // 6


def C_sum(n):
    """A222716: triangular-number analog."""
    return (n - 1) * n * (n + 1) * (3 * n * n - 2) // 6


def A_left(n):
    return n * n, n * n + n


def A_right(n):
    return n * n + n + 1, n * n + 2 * n


def B_left(n):
    return 2 * n * n + n, 2 * n * n + 2 * n


def B_right(n):
    return 2 * n * n + 2 * n + 1, 2 * n * n + 3 * n


def B_left_raw_sum(n):
    a, b = B_left(n)
    return (a + b) * (b - a + 1) // 2


def B_right_raw_sum(n):
    a, b = B_right(n)
    return (a + b) * (b - a + 1) // 2


def interval(lo, hi):
    return tuple(range(lo, hi + 1))


def assert_formulas(limit=20):
    for n in range(1, limit + 1):
        al = interval(*A_left(n))
        ar = interval(*A_right(n))
        bl = interval(*B_left(n))
        br = interval(*B_right(n))

        assert sum(al) == sum(ar) == A_sum(n)
        assert sum(x * x for x in bl) == sum(x * x for x in br) == B_sum(n)
        assert B_left_raw_sum(n) == sum(bl)
        assert B_right_raw_sum(n) == sum(br)
        assert A_sum(n) == 3 * sum(k * k for k in range(1, n + 1))


def exact_shared_block():
    """
    Solve B_left(n) = A_right(m) exactly.

    Lengths force m = n+1, and then the left endpoints give
        2n^2 + n = (n+1)^2 + (n+1) + 1 = n^2 + 3n + 3,
    so n^2 - 2n - 3 = 0, hence n = 3.
    """
    n = 3
    m = 4
    return n, m, B_left(n)


def build_endpoint_maps(m_limit):
    maps = {
        "ALs": {},
        "ALe": {},
        "ARs": {},
        "ARe": {},
    }
    for m in range(1, m_limit + 1):
        maps["ALs"][A_left(m)[0]] = m
        maps["ALe"][A_left(m)[1]] = m
        maps["ARs"][A_right(m)[0]] = m
        maps["ARe"][A_right(m)[1]] = m
    return maps


def endpoint_families(n_limit=5000, m_limit=7000):
    maps = build_endpoint_maps(m_limit)
    fam = {
        "BLs=ALs": [],
        "BLs=ARs": [],
        "BLe=ALe": [],
        "BLe=ARe": [],
        "BRs=ALs": [],
        "BRs=ARs": [],
        "BRe=ALe": [],
        "BRe=ARe": [],
    }
    for n in range(1, n_limit + 1):
        bls, ble = B_left(n)
        brs, bre = B_right(n)
        if bls in maps["ALs"]:
            fam["BLs=ALs"].append((n, maps["ALs"][bls]))
        if bls in maps["ARs"]:
            fam["BLs=ARs"].append((n, maps["ARs"][bls]))
        if ble in maps["ALe"]:
            fam["BLe=ALe"].append((n, maps["ALe"][ble]))
        if ble in maps["ARe"]:
            fam["BLe=ARe"].append((n, maps["ARe"][ble]))
        if brs in maps["ALs"]:
            fam["BRs=ALs"].append((n, maps["ALs"][brs]))
        if brs in maps["ARs"]:
            fam["BRs=ARs"].append((n, maps["ARs"][brs]))
        if bre in maps["ALe"]:
            fam["BRe=ALe"].append((n, maps["ALe"][bre]))
        if bre in maps["ARe"]:
            fam["BRe=ARe"].append((n, maps["ARe"][bre]))
    return fam


def pell_check(label, n, m):
    if label == "BLe=ALe":
        return (2 * m + 1) ** 2 - 2 * (2 * n + 1) ** 2
    if label == "BLs=ALs":
        return (4 * n + 1) ** 2 - 8 * m * m
    if label == "BRs=ALs":
        return 2 * m * m - (2 * n + 1) ** 2
    if label == "BLs=ARs":
        return (4 * n + 1) ** 2 - 2 * (2 * m + 1) ** 2
    if label == "BLe=ARe":
        return (2 * n + 1) ** 2 - 2 * (m + 1) ** 2
    if label == "BRe=ALe":
        return (4 * n + 3) ** 2 - 2 * (2 * m + 1) ** 2
    if label == "BRs=ARs":
        return (2 * m + 1) ** 2 - 2 * (2 * n + 1) ** 2
    if label == "BRe=ARe":
        return None
    raise KeyError(label)


def first_terms(seq_fn, n0, n1):
    return [seq_fn(n) for n in range(n0, n1 + 1)]


def find_sum_crossovers(n_limit=400, m_limit=700):
    hits_left = []
    hits_right = []
    a_map = {A_sum(m): m for m in range(1, m_limit + 1)}
    for n in range(1, n_limit + 1):
        v = B_left_raw_sum(n)
        if v in a_map:
            hits_left.append((n, a_map[v], v))
        v = B_right_raw_sum(n)
        if v in a_map:
            hits_right.append((n, a_map[v], v))
    return hits_left, hits_right


def fmt_pairs(pairs, k=6):
    return ", ".join(f"({n},{m})" for n, m in pairs[:k]) if pairs else "-"


def main():
    assert_formulas()
    fam = endpoint_families()
    raw_left_hits, raw_right_hits = find_sum_crossovers()
    n0, m0, block0 = exact_shared_block()

    print("=" * 78)
    print("TWO TRIANGULAR TOWERS + PELL CROSSOVERS  monad-S712")
    print("=" * 78)
    print("A_n uses consecutive integers; B_n uses consecutive squares on the 4T_n shell.")
    print()

    print("FORMULAS")
    print(f"  A_n   = sum_{'{k=n^2..n^2+n}'} k = sum_{'{k=n^2+n+1..n^2+2n}'} k")
    print("        = n*(n+1)*(2n+1)/2  = 3*sum_{j=1..n} j^2")
    print(f"  B_n   = sum_{'{k=2n^2+n..2n^2+2n}'} k^2 = sum_{'{k=2n^2+2n+1..2n^2+3n}'} k^2")
    print("        = n*(n+1)*(2n+1)*(12n^2+12n+1)/6")
    print("        = A_n * ((12n^2+12n+1)/3) = A_n * ((2n+1)^2 - 2/3)")
    print("  C_n   = triangular-number analog A222716:")
    print("        = sum_{j=n^2-n-1..n^2-1} T_j = sum_{j=n^2..n^2+n-2} T_j")
    print("        = (n-1)*n*(n+1)*(3n^2-2)/6")
    print()

    print("FIRST TERMS")
    print("  A_n:", first_terms(A_sum, 1, 8))
    print("  B_n:", first_terms(B_sum, 1, 6))
    print("  raw B-left sums :", first_terms(B_left_raw_sum, 1, 8))
    print("  raw B-right sums:", first_terms(B_right_raw_sum, 1, 8))
    print("  C_n:", first_terms(C_sum, 1, 8))
    print()

    print("RAW SUM IDENTITIES")
    print("  raw_B_left(n)  = A_n + n*(n+1)^2")
    print("  raw_B_right(n) = A_n + n^2*(n+1)")
    print("  raw_B_left(n) - raw_B_right(n) = n*(n+1)")
    print(f"  bounded raw-sum crossover raw_B_left = A_m: {raw_left_hits[:6]}")
    print(f"  bounded raw-sum crossover raw_B_right = A_m: {raw_right_hits[:6]}")
    print("  In the scanned range, 90 is the only raw left-sum shared by both towers:")
    print("    raw_B_left(3) = 21+22+23+24 = A_4 = 16+17+18+19+20 = 90")
    print()

    print("UNIQUE EXACT SHARED RAW BLOCK")
    print(f"  B_left({n0}) = A_right({m0}) = {block0}")
    print("  Proof: equality of block lengths forces m=n+1; then")
    print("         2n^2+n = m^2+m+1 = (n+1)^2+(n+1)+1, so n^2-2n-3=0 and n=3.")
    print("  So 21+22+23+24 is the unique exact raw interval shared by both towers.")
    print()

    print("ENDPOINT-CROSSOVER FAMILIES")
    print("  These are the infinite truncation / alignment families behind examples like")
    print("  10+11+12 inside 9+10+11+12 and 25+26+27 starting the A_5 left block.")
    for label in [
        "BLe=ALe",
        "BLs=ALs",
        "BRs=ALs",
        "BLs=ARs",
        "BLe=ARe",
        "BRe=ALe",
        "BRs=ARs",
        "BRe=ARe",
    ]:
        pairs = fam[label]
        print(f"  {label:8s}: {fmt_pairs(pairs)}")
        if pairs and label != "BRe=ARe":
            vals = [pell_check(label, n, m) for n, m in pairs[:5]]
            print(f"            Pell constants on first hits: {vals}")
    print()

    print("INTERPRETATION OF THE MAIN FAMILIES")
    print("  BLe=ALe:  B_left(n) is the RIGHT TAIL of A_left(m).")
    print(f"            first hits {fmt_pairs(fam['BLe=ALe'])}")
    print("            equation m(m+1)=2n(n+1), i.e. (2m+1)^2 - 2(2n+1)^2 = -1.")
    print("  BRs=ARs:  B_right(n) is the LEFT TAIL of A_right(m).")
    print(f"            same solution family {fmt_pairs(fam['BRs=ARs'])}")
    print("  BLs=ALs:  B_left(n) starts exactly where A_left(m) starts.")
    print(f"            first hits {fmt_pairs(fam['BLs=ALs'])}")
    print("            equation m^2=n(2n+1), i.e. (4n+1)^2 - 8m^2 = 1.")
    print("  BRs=ALs:  B_right(n) starts exactly where A_left(m) starts.")
    print(f"            first hits {fmt_pairs(fam['BRs=ALs'])}")
    print("            equation m^2=2n^2+2n+1, i.e. 2m^2 - (2n+1)^2 = 1.")
    print("  BLs=ARs:  B_left(n) starts exactly where A_right(m) starts.")
    print(f"            first hits {fmt_pairs(fam['BLs=ARs'])}")
    print("            equation (4n+1)^2 - 2(2m+1)^2 = 7.")
    print("  BRe=ALe:  B_right(n) ends exactly where A_left(m) ends.")
    print(f"            first hits {fmt_pairs(fam['BRe=ALe'])}")
    print("            equation (4n+3)^2 - 2(2m+1)^2 = 7.")
    print()

    print("MAIN BRIDGES")
    print("  1. A_n already hides squares: A_n = 3 * sum_{j<=n} j^2.")
    print("  2. B_n is a multiplicative lift of A_n by the factor (12n^2+12n+1)/3.")
    print("  3. The exact shared block 21..24 is rigid/unique, but one-endpoint bridges")
    print("     persist in infinite Pell families.")
    print("  4. The triangular-number analog C_n shows a third tower exists, but it lives")
    print("     on TRIANGULAR VALUES indexed near n^2, not on a new power shell.")


if __name__ == "__main__":
    main()
