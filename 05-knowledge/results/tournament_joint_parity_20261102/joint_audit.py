#!/usr/bin/env python3
"""Independent labelled census and attachment audit; no repository imports.

Bit k of mask is i->j for lexicographic pairs i<j. Exact universe: every
labelled tournament at orders 0..6. Canonicalization checks all permutations
within equal-score cells; the canonical ordering has nondecreasing outdegree.
Hamilton paths: integer subset DP, independently literal orders through five.
Run with -I -B, and with -I -B -O; no checks rely on Python assert.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial, comb
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(None)
def pairs(n):
    return tuple(combinations(range(n), 2))


def adjacency(n, code):
    out = [0] * n
    for k, (i, j) in enumerate(pairs(n)):
        a, b = (i, j) if (code >> k) & 1 else (j, i)
        out[a] |= 1 << b
    return out


def canon(n, code):
    out = adjacency(n, code)
    cells = defaultdict(list)
    for i, row in enumerate(out):
        cells[row.bit_count()].append(i)
    best, multiplicity = None, 0
    for blocks in product(*(permutations(cells[d]) for d in sorted(cells))):
        p = sum(blocks, ())
        new = sum(((out[p[i]] >> p[j]) & 1) << k
                  for k, (i, j) in enumerate(pairs(n)))
        if best is None or new < best:
            best, multiplicity = new, 1
        elif new == best:
            multiplicity += 1
    return best, multiplicity


def hp(n, code):
    if n == 0:
        return 1
    out = adjacency(n, code)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(full + 1)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, full + 1):
        for i, count in enumerate(dp[mask]):
            if not count:
                continue
            rest = out[i] & (full ^ mask)
            while rest:
                bit = rest & -rest
                j = bit.bit_length() - 1
                dp[mask | bit][j] += count
                rest ^= bit
    return sum(dp[-1])


def hp_literal(n, code):
    out = adjacency(n, code)
    return sum(all((out[a] >> b) & 1 for a, b in zip(p, p[1:]))
               for p in permutations(range(n)))


def extend(n, code, cut):
    out = adjacency(n, code)
    return sum((((out[i] >> j) & 1) if j < n else ((cut >> i) & 1)) << k
               for k, (i, j) in enumerate(pairs(n + 1)))


def histogram(n, code):
    return Counter(hp(n + 1, extend(n, code, cut)) for cut in range(1 << n))


def ordered_covers(n, code):
    out = adjacency(n, code)
    counts = [0] * n
    for p in permutations(range(n)):
        bad = sum(1 - ((out[a] >> b) & 1) for a,b in zip(p,p[1:]))
        for k in range(bad + 1, n + 1):
            counts[k-1] += comb(n-1-bad, k-1-bad)
    return counts


def main():
    expected = [1, 1, 1, 2, 4, 12, 56]
    result = {"universe": "all labelled tournaments, orders 0..6",
              "census": [], "first_H_only_attachment_failure": None}
    for n in range(7):
        fibres = Counter()
        for code in range(1 << len(pairs(n))):
            c, aut = canon(n, code)
            fibres[c] += 1
        check(len(fibres) == expected[n], f"class total {n}")
        F, SC = Counter(), Counter()
        by_h = defaultdict(list)
        for c, size in sorted(fibres.items()):
            _, aut = canon(n, c)
            h = hp(n, c)
            check(aut % 2 == h % 2 == 1, "odd automorphism/path count")
            check(size * aut == factorial(n), "orbit mass")
            if n <= 5:
                check(h == hp_literal(n, c), "DP vs literal orders")
            F[h] += 1
            if canon(n, c ^ ((1 << len(pairs(n))) - 1))[0] == c:
                SC[h] += 1
            by_h[h].append((c, aut))
        check(F[1] == 1, "unique H1")
        check(all((F[h] - SC[h]) % 2 == 0 for h in F), "converse parity by H")
        row = {"n": n, "ordinary": len(fibres), "self_converse": sum(SC.values()),
               "F": dict(sorted(F.items())), "self_converse_F": dict(sorted(SC.items()))}
        if n >= 1:
            check(sum(hp(n,c) // canon(n,c)[1] for c in fibres)
                  == 1 << ((n-1)*(n-2)//2), "Hamilton spine total")
        result["census"].append(row)
        if result["first_H_only_attachment_failure"] is None:
            for h, classes in sorted(by_h.items()):
                base, base_aut = classes[0]
                hist = histogram(n, base)
                for c, aut in classes[1:]:
                    other = histogram(n, c)
                    if hist != other:
                        result["first_H_only_attachment_failure"] = {
                            "n": n, "H": h, "codes": [base, c],
                            "aut_orders": [base_aut, aut],
                            "histograms": [dict(sorted(hist.items())), dict(sorted(other.items()))],
                            "arcs": [[[a,b] if (code>>k)&1 else [b,a]
                                      for k,(a,b) in enumerate(pairs(n))]
                                     for code in (base,c)]}
                        witness = result["first_H_only_attachment_failure"]
                        witness["ordered_path_cover_profiles"] = [ordered_covers(n,code)
                                                                  for code in (base,c)]
                        check(witness["ordered_path_cover_profiles"][0]
                              == witness["ordered_path_cover_profiles"][1], "same covers")
                        witness["first_moment_sums"] = [sum(k*v for k,v in d.items())
                                                        for d in (hist,other)]
                        witness["second_moment_sums"] = [sum(k*k*v for k,v in d.items())
                                                         for d in (hist,other)]
                        for code, expected_hist in ((base,hist),(c,other)):
                            literal = Counter(hp_literal(n+1,extend(n,code,cut))
                                              for cut in range(1<<n))
                            check(literal == expected_hist, "independent literal child orders")
                        break
                if result["first_H_only_attachment_failure"] is not None:
                    break
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
