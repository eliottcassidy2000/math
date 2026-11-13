#!/usr/bin/env python3
"""Independent labelled-orbit census and finite-core dyadic residue checks.

No repository code or class catalogue is imported. All correctness checks are
explicit exceptions, including in optimized Python. H<=21 completeness uses
the cited Busch floor H(strong n)>=25 for n>=7; enumeration alone is n<=6.
"""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
from math import comb, factorial
from pathlib import Path
import json


def need(ok, detail):
    if not ok:
        raise RuntimeError(detail)


def rows_of(code, n):
    rows = [0] * n
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            u, v = (i, j) if code & (1 << bit) else (j, i)
            rows[u] |= 1 << v
            bit += 1
    return rows


def encode(rows, p):
    out = bit = 0
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            out |= ((rows[p[i]] >> p[j]) & 1) << bit
            bit += 1
    return out


def hamilton(rows):
    n = len(rows)
    if not n:
        return 1
    @lru_cache(None)
    def extend(v, remaining):
        if not remaining:
            return 1
        return sum(extend(w, remaining ^ (1 << w))
                   for w in range(n) if remaining & rows[v] & (1 << w))
    allv = (1 << n) - 1
    return sum(extend(v, allv ^ (1 << v)) for v in range(n))


def strong(rows):
    n = len(rows)
    if n == 0:
        return False
    # Independent DFS from each vertex, rather than a matrix closure.
    for root in range(n):
        seen, todo = {root}, [root]
        while todo:
            v = todo.pop()
            for w in range(n):
                if rows[v] & (1 << w) and w not in seen:
                    seen.add(w)
                    todo.append(w)
        if len(seen) != n:
            return False
    return True


def census():
    layers = []
    for n in range(7):
        ps = list(permutations(range(n)))
        universe = 1 << (n * (n - 1) // 2)
        unseen = set(range(universe))
        classes = []
        while unseen:
            code = min(unseen)
            rows = rows_of(code, n)
            orbit = {encode(rows, p) for p in ps}
            need(orbit <= unseen, ("orbits overlap", n, code))
            unseen.difference_update(orbit)
            h = hamilton(rows)
            literal = sum(all(rows[p[i]] & (1 << p[i+1])
                              for i in range(n-1)) for p in ps)
            need(h == literal, ("Hamilton recurrence/permutation", n, code))
            need(h & 1, ("Redei control", n, code))
            classes.append({"code": code, "h": h, "strong": strong(rows),
                            "sc": (universe-1-code) in orbit,
                            "orbit_size": len(orbit),
                            "aut": factorial(n)//len(orbit)})
        need(sum(c["orbit_size"] for c in classes) == universe,
             ("labelled mass", n))
        layers.append(classes)
    need(list(map(len, layers)) == [1,1,1,2,4,12,56], "A000568 control")
    return layers


def poly_mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return c


def parity_pole(numerator, denominator_order):
    p = [a & 1 for a in numerator]
    while p and p[-1] == 0:
        p.pop()
    if not p:
        return 0
    cancelled = 0
    while cancelled < denominator_order and sum(p) % 2 == 0:
        # Exact division by 1+x over F2, constant coefficient first.
        q = [0] * (len(p)-1)
        for i in range(len(q)):
            q[i] = p[i] ^ (q[i-1] if i else 0)
        need(p[-1] == q[-1], "division by 1+x")
        p = q
        cancelled += 1
    return denominator_order-cancelled


def main():
    data = census()
    bound, horizon = 21, 160
    atoms = Counter((r["h"], n) for n, layer in enumerate(data)
                    for r in layer if n > 1 and r["strong"] and r["h"] <= bound)
    centers = Counter((r["h"], n) for n, layer in enumerate(data)
                      for r in layer if n > 1 and r["strong"] and r["sc"]
                      and r["h"] <= bound)
    centers[1,0] = centers[1,1] = 1
    need(dict(atoms) == {(3,3):1,(5,4):1,(9,5):2,(11,5):1,(13,5):1,
                         (15,5):2,(15,6):2,(17,6):2,(19,6):2}, "low atoms")
    # Core words, preserving their order and ordinary class multiplicity.
    cores = {h: Counter() for h in range(1,bound+1,2)}
    def append(product, length, size, weight):
        cores[product][length,size] += weight
        for (q, m), count in atoms.items():
            if product*q <= bound:
                append(product*q, length+1, size+m, weight*count)
    append(1,0,0,1)
    def ordinary(h, n):
        return sum(b*comb(n-m+r,r) for (r,m),b in cores[h].items() if n >= m)
    def fixed(h, n):
        out = 0
        for u in cores:
            if h % (u*u):
                continue
            for (v,c), weight in centers.items():
                if u*u*v == h and n >= c and (n-c)%2 == 0:
                    out += weight*ordinary(u,(n-c)//2)
        return out
    # A separate first-SCC dynamic recurrence does not use core padding.
    recurrence = [Counter({1:1})]
    for n in range(1,horizon+1):
        row = recurrence[n-1].copy()
        for (q,m), weight in atoms.items():
            if m <= n:
                for h,count in recurrence[n-m].items():
                    if h*q <= bound:
                        row[h*q] += weight*count
        recurrence.append(row)
    for n in range(horizon+1):
        for h in cores:
            need(ordinary(h,n) == recurrence[n][h], ("core versus first SCC",n,h))
            need((ordinary(h,n)-fixed(h,n))%2 == 0, ("converse parity",n,h))
    for n, layer in enumerate(data):
        aa = Counter(r["h"] for r in layer)
        cc = Counter(r["h"] for r in layer if r["sc"])
        for h in cores:
            need(ordinary(h,n) == aa[h] and fixed(h,n) == cc[h],
                 ("actual class versus cores",n,h))
    # Every fixed precision and all retained h use the same dyadic clock.
    R = 2  # floor(log_3 21)
    residue_cases = 0
    for k in range(1,5):
        period = 1 << (k-1+R.bit_length())
        for n in range(bound,horizon-period+1):
            for h in cores:
                need((ordinary(h,n+period)-ordinary(h,n))%(1<<k) == 0,
                     ("ordinary modular period",n,h,k))
                need((fixed(h,n+period)-fixed(h,n))%(1<<k) == 0,
                     ("fixed modular period",n,h,k))
                residue_cases += 2
    profiles = {}
    for h, core in cores.items():
        if not core:
            profiles[h] = {"max_core_length": None,"numerator":[],"parity_pole":0,
                           "least_eventual_parity_period":1}
            continue
        rmax = max(r for r,m in core)
        numerator = [0]
        for (r,m),b in core.items():
            factor = [((-1)**j)*comb(rmax-r,j) for j in range(rmax-r+1)]
            term = [0]*m + [b*x for x in factor]
            numerator.extend([0]*max(0,len(term)-len(numerator)))
            for i,x in enumerate(term):
                numerator[i] += x
        pole = parity_pole(numerator,rmax+1)
        period = 1 if pole <= 1 else 1 << (pole-1).bit_length()
        for n in range(bound,horizon-period+1):
            need((ordinary(h,n+period)-ordinary(h,n))%2 == 0,
                 ("reduced pole period",h,n))
        if period > 1:
            need(any((ordinary(h,n+period//2)-ordinary(h,n))%2
                     for n in range(bound,bound+period)), ("least period",h))
        profiles[h] = {"max_core_length":rmax,"numerator":numerator,
                       "parity_pole":pole,"least_eventual_parity_period":period}
    defects = [sum(fixed(h,n)%2 for h in cores) for n in range(horizon+1)]
    for n in range(4,horizon+1):
        need(defects[n] == 2+2*(n%2)+(n//2)%2, ("closed defect",n))
    need(defects[3] == 2 and defects[7] == 5, "onset hostile")
    need(fixed(9,5)==2 and fixed(19,20)==0 and ordinary(19,20)==30,
         "fixed-point/multiplicity hostiles")
    def valuation2(x):
        need(x > 0, "positive valuation input")
        return (x & -x).bit_length()-1
    for n in range(6,horizon+1):
        need(ordinary(19,n)==2*(n-5) and fixed(19,n)==0, ("H19 rank inventory",n))
        need(valuation2(ordinary(19,n))==1+valuation2(n-5), ("H19 exact rank",n))
    for k in range(1,7):
        period = 1 << (k-1)
        for n in range(6,horizon-period+1):
            need(min(k,1+valuation2(n-5)) == min(k,1+valuation2(n+period-5)),
                 ("capped H19 rank period",n,k))
        if period > 1:
            need(any(min(k,1+valuation2(n-5)) != min(k,1+valuation2(n+period//2-5))
                     for n in range(6,6+period)), ("least capped rank period",k))
    # Explicit optimal involution on each actual class set through six.
    pairing_rows = []
    for n,layer in enumerate(data):
        ps = list(permutations(range(n)))
        mask = (1 << (n*(n-1)//2))-1
        codes = {r["code"]:r for r in layer}
        rev = {code:min(encode(rows_of(mask-code,n),p) for p in ps) for code in codes}
        J = {code:rev[code] for code in codes if rev[code]!=code}
        buckets = {}
        for code in codes:
            if rev[code] == code:
                buckets.setdefault(codes[code]["h"],[]).append(code)
        for group in buckets.values():
            group.sort()
            for i in range(0,len(group)-1,2):
                J[group[i]], J[group[i+1]] = group[i+1],group[i]
            if len(group)%2:
                J[group[-1]] = group[-1]
        for code in codes:
            need(J[J[code]] == code and J[rev[code]] == rev[J[code]], "involutions commute")
            need(codes[code]["h"] == codes[J[code]]["h"], "preserve H")
        fixed_codes = [c for c in codes if J[c]==c]
        need(len(fixed_codes)==sum(len(g)%2 for g in buckets.values()), "minimum pairing")
        pairing_rows.append({"n":n,"fixed_H":[codes[c]["h"] for c in fixed_codes]})
    return {"status":"FINITE-EXACT controls; low-H all-order extension uses cited Busch cutoff",
            "universe":"every labelled orientation at n=0..6; full S_n orbit partition; no filters",
            "class_counts":list(map(len,data)),"strong_atoms_le21":[[*q,c] for q,c in sorted(atoms.items())],
            "self_converse_centers_le21":[[*q,c] for q,c in sorted(centers.items())],
            "core_profiles":profiles,"defect_le21_n0_to20":defects[:21],
            "defect_closed_form":"n>=4: 2+2*(n mod 2)+(floor(n/2) mod 2)",
            "H19_maximal_compatible_free_rank":"n>=6: 1+v2(n-5)",
            "H19_rank_n6_to21":[1+valuation2(n-5) for n in range(6,22)],
            "joint_residue_cases":residue_cases,"explicit_pairings":pairing_rows,
            "source_sha256":sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__ == "__main__":
    print(json.dumps(main(),sort_keys=True,indent=2))
