#!/usr/bin/env python3
"""
Two triangular towers: complete taxonomy of endpoint-aligned subset families.
monad-explorer-2026-06-12-S713

Build on S712.  The additive tower partitions all positive integers:

    stage m = [m^2, ..., (m+1)^2-1] = A_L(m) union A_R(m)

where
    A_L(m) = [m^2, ..., m^2+m],   size m+1
    A_R(m) = [m^2+m+1, ..., m^2+2m], size m

The raw square-tower intervals are
    B_L(n) = [2n^2+n, ..., 2n^2+2n],   size n+1
    B_R(n) = [2n^2+2n+1, ..., 2n^2+3n], size n

The user asked for all families like
    10+11+12 inside 9+10+11+12,
    13+14 inside 13+14+15,
    21+22+23+24 exactly,
    25+26+27 inside 25+...+30.

This script shows the complete answer:
  every endpoint-aligned subset pattern is one of 8 Pell-type families:
    B_L/B_R as prefix/suffix of A_L/A_R.

The "when do they show up?" criterion is concise:
  they occur exactly when one endpoint of a B-interval lands on one endpoint of
  an A-half-interval, i.e. when one of
    2n^2+n, 2n^2+2n, 2n^2+2n+1, 2n^2+3n
  equals one of
    m^2, m^2+m, m^2+m+1, m^2+2m.
These reduce to 8 Pell equations.
"""

from math import isqrt


def A_left(m):
    return m * m, m * m + m


def A_right(m):
    return m * m + m + 1, m * m + 2 * m


def B_left(n):
    return 2 * n * n + n, 2 * n * n + 2 * n


def B_right(n):
    return 2 * n * n + 2 * n + 1, 2 * n * n + 3 * n


def size(interval):
    a, b = interval
    return b - a + 1


def locate_in_A(x):
    m = isqrt(x)
    while (m + 1) * (m + 1) <= x:
        m += 1
    while m * m > x:
        m -= 1
    side = "L" if x <= m * m + m else "R"
    return m, side


def raw_shell_gaps(limit=12):
    rows = []
    for n in range(1, limit + 1):
        lo, hi = 2 * n * n + n, 2 * n * n + 3 * n
        next_lo = 2 * (n + 1) * (n + 1) + (n + 1)
        gap = next_lo - hi - 1
        rows.append((n, (lo, hi), 2 * n + 1, gap))
    return rows


def stage_taxonomy(limit=30):
    rows = []
    for n, f, tag in [(n, B_left, "B_L") for n in range(1, limit + 1)] + [
        (n, B_right, "B_R") for n in range(1, limit + 1)
    ]:
        a, b = f(n)
        ma, sa = locate_in_A(a)
        mb, sb = locate_in_A(b)
        rows.append((tag, n, (a, b), (ma, sa), (mb, sb)))
    return rows


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


def endpoint_hits(n_limit=8000, m_limit=12000):
    maps = build_endpoint_maps(m_limit)
    fam = {
        "BL prefix AL": [],
        "BL suffix AL": [],
        "BL prefix AR": [],
        "BL suffix AR": [],
        "BR prefix AL": [],
        "BR suffix AL": [],
        "BR prefix AR": [],
        "BR suffix AR": [],
    }
    for n in range(1, n_limit + 1):
        bl = B_left(n)
        br = B_right(n)
        if bl[0] in maps["ALs"]:
            fam["BL prefix AL"].append((n, maps["ALs"][bl[0]]))
        if bl[1] in maps["ALe"]:
            fam["BL suffix AL"].append((n, maps["ALe"][bl[1]]))
        if bl[0] in maps["ARs"]:
            fam["BL prefix AR"].append((n, maps["ARs"][bl[0]]))
        if bl[1] in maps["ARe"]:
            fam["BL suffix AR"].append((n, maps["ARe"][bl[1]]))
        if br[0] in maps["ALs"]:
            fam["BR prefix AL"].append((n, maps["ALs"][br[0]]))
        if br[1] in maps["ALe"]:
            fam["BR suffix AL"].append((n, maps["ALe"][br[1]]))
        if br[0] in maps["ARs"]:
            fam["BR prefix AR"].append((n, maps["ARs"][br[0]]))
        if br[1] in maps["ARe"]:
            fam["BR suffix AR"].append((n, maps["ARe"][br[1]]))
    return fam


def pell_constant(label, n, m):
    if label == "BL prefix AL":
        return (4 * n + 1) ** 2 - 8 * m * m
    if label == "BL suffix AL":
        return (2 * m + 1) ** 2 - 2 * (2 * n + 1) ** 2
    if label == "BL prefix AR":
        return (4 * n + 1) ** 2 - 2 * (2 * m + 1) ** 2
    if label == "BL suffix AR":
        return (2 * n + 1) ** 2 - 2 * (m + 1) ** 2
    if label == "BR prefix AL":
        return 2 * m * m - (2 * n + 1) ** 2
    if label == "BR suffix AL":
        return (4 * n + 3) ** 2 - 2 * (2 * m + 1) ** 2
    if label == "BR prefix AR":
        return (2 * m + 1) ** 2 - 2 * (2 * n + 1) ** 2
    if label == "BR suffix AR":
        return (4 * n + 3) ** 2 - 2 * (2 * m + 2) ** 2
    raise KeyError(label)


def host_size(label, m):
    if "AL" in label:
        return m + 1
    return m


def block_size(label, n):
    if label.startswith("BL"):
        return n + 1
    return n


def deficit(label, n, m):
    return host_size(label, m) - block_size(label, n)


def describe_family(label, n, m):
    d = deficit(label, n, m)
    if "prefix" in label:
        return f"subset of first {block_size(label, n)} terms; tail deficit {d}"
    return f"subset of last {block_size(label, n)} terms; head deficit {d}"


def print_family(label, pairs):
    first = pairs[:6]
    print(f"  {label}: {first}")
    consts = [pell_constant(label, n, m) for n, m in first[:5]]
    print(f"    Pell constants on first hits: {consts}")
    if first:
        n, m = first[0]
        print(f"    first hit meaning: n={n}, m={m}; {describe_family(label, n, m)}")


def subset_pairs(label, pairs):
    return [(n, m) for n, m in pairs if deficit(label, n, m) >= 0]


def main():
    fam_all = endpoint_hits()
    fam = {label: subset_pairs(label, pairs) for label, pairs in fam_all.items()}

    print("=" * 78)
    print("TWO TRIANGULAR TOWERS: COMPLETE SUBSET-FAMILY TAXONOMY  monad-S713")
    print("=" * 78)
    print("The additive tower covers all positive integers because stage m is")
    print("  [m^2, ..., (m+1)^2-1] = A_L(m) union A_R(m).")
    print("The raw square shell does NOT: its stage-n carrier is")
    print("  [2n^2+n, ..., 2n^2+3n], length 2n+1,")
    print("and the gap before the next carrier has length 2n+2.")
    print()

    print("RAW SHELL COVERAGE")
    for n, interval, length, gap in raw_shell_gaps(8):
        print(f"  n={n}: shell={interval}, length={length}, next gap={gap}")
    print()

    print("THE COMPLETE ENDPOINT-ALIGNED SUBSET FAMILIES")
    print("These are exactly the families where a B-endpoint lands on an A-half-endpoint.")
    print("There are 8 such families, one for each prefix/suffix choice in A_L/A_R.")
    print("  Exception: the Pell family BL prefix AR has a lone startup straddle (n,m)=(1,1),")
    print("  where B_L(1)=[3,4] starts at A_R(1)=[3] but spills one step into A_L(2)=[4].")
    print("  The genuine subset branch of that family starts at (3,4).")
    print()

    order = [
        "BL prefix AL",
        "BL suffix AL",
        "BL prefix AR",
        "BL suffix AR",
        "BR prefix AL",
        "BR suffix AL",
        "BR prefix AR",
        "BR suffix AR",
    ]
    for label in order:
        print_family(label, fam[label])
    print()

    print("THE USER'S EXAMPLES IN THIS TAXONOMY")
    print("  10+11+12 inside 9+10+11+12  -> BL suffix AL, first hit (n,m)=(2,3)")
    print("  13+14 inside 13+14+15      -> BR prefix AR, same Pell family (2,3)")
    print("  21+22+23+24                -> both BL prefix AR and BL suffix AR at (3,4)")
    print("                                so it is the unique exact shared B_L=A_R block")
    print("  25+26+27 inside 25+...+30  -> BR prefix AL, first hit (3,5)")
    print("  36+...+40 inside 36+...+42 -> BL prefix AL, first hit (4,6)")
    print()

    print("CONCISE PREDICTION RULE")
    print("  A family appears exactly when one of")
    print("    2n^2+n, 2n^2+2n, 2n^2+2n+1, 2n^2+3n")
    print("  equals one of")
    print("    m^2, m^2+m, m^2+m+1, m^2+2m.")
    print("  Each choice gives one Pell equation and therefore an infinite sparse family.")
    print("  The host size is m+1 in A_L(m), m in A_R(m);")
    print("  the guest size is n+1 for B_L(n), n for B_R(n).")
    print("  Their difference is the exact prefix/suffix deficit.")
    print()

    print("PELL FORMS")
    print("  BL prefix AL : 2n^2+n = m^2                 -> (4n+1)^2 - 8m^2 = 1")
    print("  BL suffix AL : 2n^2+2n = m^2+m             -> (2m+1)^2 - 2(2n+1)^2 = -1")
    print("  BL prefix AR : 2n^2+n = m^2+m+1            -> (4n+1)^2 - 2(2m+1)^2 = 7")
    print("  BL suffix AR : 2n^2+2n = m^2+2m            -> (2n+1)^2 - 2(m+1)^2 = -1")
    print("  BR prefix AL : 2n^2+2n+1 = m^2             -> 2m^2 - (2n+1)^2 = 1")
    print("  BR suffix AL : 2n^2+3n = m^2+m             -> (4n+3)^2 - 2(2m+1)^2 = 7")
    print("  BR prefix AR : 2n^2+2n+1 = m^2+m+1         -> (2m+1)^2 - 2(2n+1)^2 = -1")
    print("  BR suffix AR : 2n^2+3n = m^2+2m            -> (4n+3)^2 - 2(2m+2)^2 = 1")


if __name__ == "__main__":
    main()
