"""Auditor's check of THM-4950 (2) lane counts, (4) first missing irreducible per lane, the 4,063 count and the
calibration table (m* with strong pieces up to order 8, 9, 10).  Own code: smallest-prime-factor sieve, a direct
irreducibility test in M3 (search for a factorisation x * y with x, y in M3, x, y > 1), and a lane classifier written
from the statement of THM-4950 (2).
usage: python lanes_check.py AUD_CLASS_DIR strong_values_le10.txt few_cycle_strong_80.txt
"""
import sys
import os


def spf_sieve(n):
    spf = list(range(n + 1))
    for i in range(2, int(n ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def factors(v, spf):
    f = []
    while v > 1:
        p = spf[v]
        f.append(p)
        v //= p
    return f


def inM3(v):
    return v % 3 == 1 and v not in (10, 31)


def irreducible(v, divs):
    """v in M3, v > 1: irreducible iff no divisor d with 1 < d < v, d and v/d both in M3"""
    for d in divs:
        if d * d > v:
            break
        if v % d == 0 and inM3(d) and inM3(v // d):
            return False
    return True


def lane(v, f):
    """THM-4950 (2), read literally"""
    if v in (40, 250, 310, 9610, 29791):
        return "sporadic"
    r1 = sorted(p for p in f if p % 3 == 1)
    r2 = sorted(p for p in f if p % 3 == 2)
    if len(f) == 1 and r1 and v != 31:
        return "P"
    if len(f) == 2 and len(r2) == 2 and v != 10:
        return "Q"
    if len(f) == 3 and r2 == [2, 5] and len(r1) == 1 and r1[0] != 31:
        return "10P"
    if len(f) == 2 and 31 in f and len(r1) == 2:
        return "31P"
    if len(f) == 3 and r1 == [31] and len(r2) == 2 and r2 != [2, 5]:
        return "31Q"
    return None


def closure(pieces, V):
    reach = bytearray(V + 1)
    reach[1] = 1
    for p in sorted(set(x for x in pieces if 2 <= x <= V)):
        for x in range(1, V // p + 1):
            if reach[x]:
                reach[x * p] = 1
    return reach


def mstar(pieces, M):
    r = closure(pieces, 1 + 3 * M)
    return next(m for m in range(M + 1) if m not in (3, 10) and not r[1 + 3 * m])


def main():
    X = 200000
    spf = spf_sieve(X)
    m3 = [v for v in range(4, X + 1) if inM3(v)]
    irr = [v for v in m3 if irreducible(v, m3)]
    counts = {}
    bad = []
    for v in irr:
        L = lane(v, factors(v, spf))
        if L is None:
            bad.append(v)
        counts[L] = counts.get(L, 0) + 1
    pred = [v for v in m3 if lane(v, factors(v, spf)) is not None]
    print(f"irreducibles of M3 below {X}: {len(irr)}; by lane {dict(sorted((str(k), c) for k, c in counts.items()))}")
    print(f"irreducibles outside the five lanes and sporadics: {bad[:10]} ({len(bad)}); "
          f"lane elements that are reducible: {sorted(set(pred) - set(irr))[:10]} ({len(set(pred) - set(irr))})")
    # pieces
    aud = sys.argv[1]
    strong = set(int(x) for x in open(sys.argv[2]).read().split())
    by_order = {}
    for n in range(1, 10):
        by_order[n] = set(int(line.split()[7]) for line in open(os.path.join(aud, f"acls_{n}.txt")) if line.split()[1] == "1")
    le8 = set().union(*(by_order[n] for n in range(1, 9)))
    le9 = le8 | by_order[9]
    inv = set()
    for line in open(sys.argv[3]):
        f = line.split()
        a5 = [int(x) for x in f[1:6]]
        H = int(f[6])
        a6 = (H - 1 - sum(c << (k + 1) for k, c in enumerate(a5))) // 64
        inv.add(1 + sum(c * 3 ** (k + 1) for k, c in enumerate(a5 + [a6])))
    # intercyclic pieces as in the packet (one-apex values recomputed by closure_check.py's routine)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from closure_check import one_apex
    M = 40000
    apex = one_apex(19, M)
    rays = {3 * 2 ** s for s in range(1, 20)} | {9 * 2 ** s - 2 for s in range(0, 20)}
    inter = {1 + 3 * a for a in apex | rays if 1 <= a <= M}
    # the order-10 value set (ext10.c, written out as my_t10_values.txt)
    o10 = set(int(x) for x in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "my_t10_values.txt")).read().split())
    s86 = le9 | inv | inter
    print(f"order-10 values not given by S86's pieces (orders <= 9, alpha_1 <= 80, intercyclic): {len(o10 - s86)}")
    print(f"calibration: m* with strong pieces up to order 8 (+ S86 inventory and intercyclic): {mstar(le8 | inv | inter, 2000)}; "
          f"up to order 9: {mstar(le9 | inv | inter, 4000)}; up to order 10: {mstar(le9 | o10 | inv | inter, 12000)}")
    pieces = le9 | o10 | inv | inter
    firsts = {}
    for v in irr:
        if v not in pieces:
            L = lane(v, factors(v, spf))
            firsts.setdefault(L, v)
    print(f"first irreducible not a piece value, by lane: {firsts}")
    print(f"sporadics attained as strong values: {[(s, s in strong) for s in (40, 250, 310, 9610, 29791)]}; "
          f"29791 at order 10: {29791 in o10}")
    print(f"m* = 9679 is {29038 / max(o10):.4f} of the largest order-10 value {max(o10)}")


if __name__ == "__main__":
    main()
