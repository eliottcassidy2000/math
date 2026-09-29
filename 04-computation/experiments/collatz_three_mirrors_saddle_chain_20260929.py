#!/usr/bin/env python3
"""The saddle chain of 1 + 4^m and the Hardy-Ramanujan number 1729 (S23, 2026-09-29).

Proposition 4 of collatz_hedgehog_family_20260927.md: for x = 1 mod 4^j, T^(2j)(x) - 1 = (3/4)^j (x - 1), so the
orbit of 1 + 4^m under T^2 runs down the chain 1 + 3^j 4^(m-j), j = 0..m (the "(3x+1)/4 saddle").  1729 = 1 + 12^3
= 1 + 3^3 4^3 is the balanced point (j = 3) of the chain of m = 6: 4097, 3073, 2305, 1729, 1297, 973, 730.
This script checks the chain identities, the factorisations, where the orbits of 1729 and 27 merge, the position of
1 + 12^m on the chain of 2m, the chain of 17 = 1 + 4^2, the mod-3 residues and odd preimages of the chain points,
and the multiple of 139 = 3^7 - 2^11 (the clock of the -17 cycle) that the chain contains.
Run: python 04-computation/experiments/collatz_three_mirrors_saddle_chain_20260929.py
"""
from __future__ import annotations

from sympy import factorint


def T(n: int) -> int:
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def orbit(n: int, stop=1):
    out = [n]
    while n != stop and n != 1:
        n = T(n)
        out.append(n)
    return out


if __name__ == "__main__":
    print("== the chain of 1 + 4^6 under T^2 ==")
    m = 6
    x = 1 + 4 ** m
    chain = [x]
    for j in range(m):
        x = T(T(x))
        chain.append(x)
    print("   chain:", chain)
    print("   predicted 1 + 3^j 4^(6-j):", [1 + 3 ** j * 4 ** (6 - j) for j in range(7)])
    assert chain == [1 + 3 ** j * 4 ** (6 - j) for j in range(7)]
    for c in chain:
        print(f"   {c:5d} = {dict(factorint(c))}   c mod 3 = {c % 3}, mod 4 = {c % 4}")
    print(f"   1729 = 1 + 12^3: {1 + 12 ** 3}; = 9^3 + 10^3: {9 ** 3 + 10 ** 3}; 730 - 1 = 3^6 = {3 ** 6}; 1729 - 730 = {1729 - 730} = 10^3 - 1")
    print(f"   1297 = 6^4 + 1: {6 ** 4 + 1}; 973 = 7 * 139: {7 * 139}; 139 = 3^7 - 2^11: {3 ** 7 - 2 ** 11}; 7*3^7 - 7*2^11 = 4*3^5 + 1: {7 * 3 ** 7 - 7 * 2 ** 11} = {4 * 3 ** 5 + 1}")
    print("== 1 + 12^m is the midpoint (j = m) of the chain of 1 + 4^(2m) ==")
    for mm in range(1, 7):
        x = 1 + 4 ** (2 * mm)
        for j in range(mm):
            x = T(T(x))
        print(f"   m = {mm}: T^(2m)(1 + 4^(2m)) = {x} = 1 + 12^m = {1 + 12 ** mm}: {x == 1 + 12 ** mm}; factors {dict(factorint(1 + 12 ** mm))}")
    print("== 17 = 1 + 4^2 and its chain 17 -> 13 -> 10 ==")
    print("   orbit of 17:", orbit(17))
    print("   chain of 1 + 4^2:", [1 + 3 ** j * 4 ** (2 - j) for j in range(3)])
    print("== the orbits of 1729 and 27 ==")
    o1729 = orbit(1729)
    o27 = orbit(27)
    print(f"   |orbit 1729| = {len(o1729) - 1} T-steps, max {max(o1729)}; |orbit 27| = {len(o27) - 1}, max {max(o27)}")
    merge = next(v for v in o1729 if v in set(o27))
    print(f"   first common element: {merge} (position {o1729.index(merge)} in the orbit of 1729, {o27.index(merge)} in that of 27); 137 = 1 + 8*17 = {1 + 8 * 17}")
    print("   orbit of 1729 to the merge:", o1729[: o1729.index(merge) + 1])
    print("== odd preimages of the chain points: (2^k a - 1)/3 for even k (a = 1 mod 3) ==")
    for c in chain:
        pre = []
        for k in range(1, 13):
            if (2 ** k * c - 1) % 3 == 0:
                pre.append((k, (2 ** k * c - 1) // 3))
        print(f"   a = {c}: odd preimages (k, (2^k a - 1)/3) = {pre[:5]}")
    print("== the 3x-1 sheet: -17 cycle and its clock ==")
    o = [-17]
    for _ in range(11):
        o.append(T(o[-1]) if o[-1] % 2 == 0 else (3 * o[-1] + 1) // 2)
    print("   3x+1 orbit of -17:", o)
    print("DONE")
