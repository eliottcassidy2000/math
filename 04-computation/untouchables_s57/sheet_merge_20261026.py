#!/usr/bin/env python3
"""Merged sheet statement: centre of the odd multiplication table, reflection centres of the one-step maps, the
unit-shift conjugacy of the root-27 section change, and the classification of the halving-shift alternation
(opus-2026-10-26-S57; completes the exercise left open in kind-pasteur-2026-10-25-S4 Theorem A).

Checks:
 (R) floor(3x/2) commutes with x -> 1 - x, ceil(3x/2) with x -> -1 - x; negation swaps floor and ceil; the unit
     translation conjugates floor(3x/2) into ceil(3z/2) (z = x - 1) but not the full-halving Syracuse maps.
 (C) for odd m = 2k+1: k = -1/2, k+1 = 1/2 (mod m) are the centre labels, k^2 = (k+1)^2 = 1/4, k(k+1) = -1/4
     (2026-10-05 Prop. 1), and 3(1/4) - 1 = -1/4, 3(-1/4) + 1 = 1/4.
 (K) -1 in <2> mod m  <=>  every prime p | m has ord_p(2) even with the same 2-adic valuation  <=>  the two centre
     values +-1/4 lie in one <2>-coset  <=>  some t has E^-_{k+t}(m) = {(-x, y) : (x, y) in E^+_k(m)} for all k.
Usage: python sheet_merge_20261026.py [mmax]
"""
import sys, json, os
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))


def ordmod(a, m):
    x, k = a % m, 1
    while x != 1:
        x = x * a % m
        k += 1
    return k


def v2(n):
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    return c


def prime_factors(n):
    f, p = [], 2
    while p * p <= n:
        if n % p == 0:
            f.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        f.append(n)
    return f


def syr(x, s):
    y = 3 * x + s
    while y % 2 == 0 and y != 0:
        y //= 2
    return y


def tables(m, kmax, sgn):
    """E_k^{sgn}(m) for k = 1..kmax as sets of (x mod m, (3x+sgn) 2^-k mod m) over all residues x."""
    inv2 = pow(2, -1, m)
    out = []
    for k in range(1, kmax + 1):
        ik = pow(inv2, k, m)
        out.append(frozenset((x, (3 * x + sgn) * ik % m) for x in range(m)))
    return out


if __name__ == '__main__':
    mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 2001
    res = {}
    # (R) integer identities
    X = range(-100000, 100001)
    fl = lambda x: (3 * x) // 2
    ce = lambda x: -((-3 * x) // 2)
    assert all(fl(1 - x) == 1 - fl(x) for x in X)
    assert all(ce(-1 - x) == -1 - ce(x) for x in X)
    assert all(-fl(-x) == ce(x) for x in X)
    assert all(ce(x - 1) == fl(x) - 1 for x in X)          # T_+(z) = T_-(z + 1) - 1 with z = x - 1
    bad = [x for x in range(1, 2001, 2) if syr(x - 1, 1) != syr(x, -1) - 1 and x - 1 != 0]
    res['R'] = {'reflection_and_shift_identities_checked_on': '|x| <= 1e5',
                'syracuse_not_conjugated_by_unit_shift_counterexamples_below_2000': len(bad), 'first': bad[:5]}
    # (C) and (K)
    rows = []
    mism = 0
    for m in range(3, mmax + 1, 2):
        k = (m - 1) // 2
        d = pow(4, -1, m)
        assert (2 * k) % m == m - 1 and k * k % m == d and (k + 1) ** 2 % m == d and k * (k + 1) % m == (-d) % m
        assert (3 * d - 1) % m == (-d) % m and (3 * (-d) + 1) % m == d
        H = set()
        x = 1
        while True:
            H.add(x)
            x = 2 * x % m
            if x == 1:
                break
        a = (m - 1) in H                                   # -1 in <2>
        ps = prime_factors(m)
        o = [ordmod(2, p) for p in ps]
        b = all(t % 2 == 0 for t in o) and len(set(v2(t) for t in o)) == 1
        c = ((-d) % m) in {d * h % m for h in H}           # +-1/4 in one <2>-coset
        # table shift (only for moderate m: tables are m-element sets)
        e = None
        if m <= 301:
            o2 = ordmod(2, m)
            Ep = tables(m, 2 * o2, +1)
            Em = tables(m, 3 * o2, -1)
            e = False
            for t in range(o2):
                if all(Em[kk + t] == frozenset(((-x) % m, y) for (x, y) in Ep[kk]) for kk in range(o2)):
                    e = True
                    break
        if not (a == b == c) or (e is not None and e != a):
            mism += 1
        rows.append((m, a, b, c, e))
    alt = [r[0] for r in rows if r[1]]
    res['K'] = {'odd_m_checked': len(rows), 'mismatches': mism, 'alternating_moduli_up_to_101': [m for m in alt if m <= 101],
                'non_alternating_up_to_101': [r[0] for r in rows if not r[1] and r[0] <= 101]}
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(HERE, 'sheet_merge.json'), 'w'), indent=1)
