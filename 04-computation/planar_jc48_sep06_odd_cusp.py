#!/usr/bin/env python3
"""Exact finite controls for the odd-cusp sparse actual-subset theorem.

Standalone Python standard library. Left permutation actions; compose(s,t)
means s after t. No inference of geometric realizability from a passport.
"""

from itertools import permutations
import hashlib
import json


GATES = 0


def check(predicate, label):
    global GATES
    GATES += 1
    if not predicate:
        raise RuntimeError(label)


def compose(s, t):
    return tuple(s[t[i]] for i in range(len(s)))


def power(s, r):
    ans = tuple(range(len(s)))
    for _ in range(r):
        ans = compose(s, ans)
    return ans


def image(s, a):
    return frozenset(s[i] for i in a)


def cycle(d, entries):
    ans = list(range(d))
    for i, x in enumerate(entries):
        ans[x] = entries[(i + 1) % len(entries)]
    return tuple(ans)


def cycles(s):
    todo = set(range(len(s)))
    out = []
    while todo:
        x = min(todo)
        row = []
        while x in todo:
            row.append(x)
            todo.remove(x)
            x = s[x]
        out.append(tuple(row))
    return out


def supported_permutations(d, support):
    for values in permutations(support):
        ans = list(range(d))
        for x, y in zip(support, values):
            ans[x] = y
        yield tuple(ans)


def canonical(q, n=0, extra=0):
    k = (q - 1) // 2
    d = q + n + extra
    aa = frozenset(range(1, k + 1)) | frozenset(range(q, q + n))
    bb = frozenset(range(k + 1, q)) | frozenset(range(q, q + n))
    sigma = cycle(d, (0,) + tuple(range(k + 1, q)))
    tau = cycle(d, (0,) + tuple(range(1, k + 1)))
    return sigma, tau, aa, bb


def describe(q, m, n):
    sigma, tau, aa, bb = canonical(q, n)
    r = (m - 1) // 2
    w = compose(sigma, tau)
    g = power(w, r)
    d = len(sigma)
    check(image(g, aa) == bb, "named reaccess")
    check(compose(g, sigma) == compose(tau, g), "named Artin")
    check(aa <= image(sigma, aa) and bb <= image(tau, bb), "named fixed sets")
    check(all(sigma[x] == x for x in aa), "named A pointwise fixed")
    check(all(tau[x] == x for x in bb), "named B pointwise fixed")
    a = len(aa)
    delta = d - a
    check(d == q + n and delta == (q + 1) // 2, "named degrees")
    check(len(aa & bb) == n and d - len(aa | bb) == 1, "named h,n")
    check(1 - n >= max(0, d - 2 * a), "named N=1 ledger")
    if n == 1:
        check(d == 2 * a and delta == a, "named multi-node ledger")
    return {"q": q, "m": m, "n": n, "d": d, "a": a,
            "delta": delta, "sigma": sigma, "tau": tau,
            "A": sorted(aa), "B": sorted(bb), "product_cycles": cycles(w)}


def main():
    # Exhaust ALL independent permutations with the prescribed fixed A/B.
    # No Artin-relation filter is used before testing actual reaccess.
    census = []
    total_rows = total_survivors = 0
    for k in (1, 2):
        d = 2 * k + 1
        aa = frozenset(range(1, k + 1))
        bb = frozenset(range(k + 1, d))
        sigmas = list(supported_permutations(d, (0,) + tuple(sorted(bb))))
        taus = list(supported_permutations(d, (0,) + tuple(sorted(aa))))
        for r in range(10):
            count = 0
            for sigma in sigmas:
                for tau in taus:
                    total_rows += 1
                    w = compose(sigma, tau)
                    g = power(w, r)
                    retained = image(g, aa) == bb
                    one_sigma = len([c for c in cycles(sigma) if 0 in c][0]) == k + 1
                    one_tau = len([c for c in cycles(tau) if 0 in c][0]) == k + 1
                    predicted = one_sigma and one_tau and (2 * r + 1) % d == 0
                    check(retained == predicted, "complete small permutation classification")
                    if retained:
                        count += 1
                        check(compose(g, sigma) == compose(tau, g), "Artin follows from retained classification")
                        check(len(cycles(w)) == 1, "full odd product cycle")
            expected = (1 if k == 1 else 4) if (2 * r + 1) % d == 0 else 0
            check(count == expected, "labelled survivor count")
            total_survivors += count
            census.append([k, r, count])
    check(total_rows == 400 and total_survivors == 11, "declared small universe")

    # Independent positional transport is checked against tuple permutations.
    canonical_rows = 0
    canonical_survivors = 0
    for q in range(3, 20, 2):
        k = (q - 1) // 2
        for m in range(3, 100, 2):
            r = (m - 1) // 2
            for n in (0, 1):
                canonical_rows += 1
                sigma, tau, aa, bb = canonical(q, n)
                w = compose(sigma, tau)
                g = power(w, r)
                expected_w = tuple((j + 1) % q for j in range(q)) + tuple(range(q, q + n))
                check(w == expected_w, "position-independent product check")
                moved_block = frozenset((j + r) % q for j in range(1, k + 1))
                target_block = frozenset(range(k + 1, q))
                check((moved_block == target_block) == (m % q == 0), "cyclic block period")
                check((image(g, aa) == bb) == (m % q == 0), "divisor reaccess")
                check((compose(g, sigma) == compose(tau, g)) == (m % q == 0), "divisor Artin")
                canonical_survivors += int(m % q == 0)
    check(canonical_rows == 882, "declared canonical universe")

    named = [describe(q, m, n) for q, m, n in
             ((3, 3, 0), (3, 3, 1), (5, 5, 1), (3, 9, 0),
              (3, 9, 1), (9, 9, 0), (9, 9, 1), (5, 15, 0), (5, 15, 1))]

    # Hostile 1: fixed/inertia data and Artin alone lose actual reaccess.
    identity = (0, 1, 2)
    aa, bb = frozenset((1,)), frozenset((2,))
    check(compose(identity, identity) == identity, "identity Artin")
    check(image(identity, aa) != bb, "missing actual reaccess hostile")

    # Hostile 2: an extra unretained fixed sheet preserves all local relations.
    sigma, tau, aa, bb = canonical(3, 1, extra=1)
    g = compose(sigma, tau)
    check(image(g, aa) == bb and compose(g, sigma) == compose(tau, g), "h2 local hostile")
    check(len(sigma) - len(aa | bb) == 2, "h2 outside labels")
    check(1 - len(aa & bb) < max(0, len(sigma) - 2 * len(aa)), "h2 Euler rejection")

    # Hostile 3: two common actual pages violate the Euler budget.
    sigma, tau, aa, bb = canonical(3, 2)
    g = compose(sigma, tau)
    check(image(g, aa) == bb and compose(g, sigma) == compose(tau, g), "n2 local hostile")
    check(len(aa & bb) == 2 and 1 - len(aa & bb) < 0, "n2 Euler rejection")

    # Hostile 4: reversing multiplication without changing transport is wrong.
    sigma, tau, aa, bb = canonical(5, 1)
    check(image(power(compose(sigma, tau), 2), aa) == bb, "correct rightmost-first gauge")
    check(image(power(compose(tau, sigma), 2), aa) != bb, "uncompensated product reversal hostile")

    # Formal q1 is a valid identity passport. Geometry, not permutations,
    # excludes it by delta=1 purity (or the cited degree-two theorem).
    sigma, tau, aa, bb = canonical(1, 1)
    for r in range(10):
        g = power(compose(sigma, tau), r)
        check(image(g, aa) == bb and compose(g, sigma) == compose(tau, g), "q1 survives abstractly")
    check(len(sigma) == 2 and len(aa) == 1, "q1 mapping degree and delta")

    spectra = {}
    for m in (3, 5, 9, 15):
        qq = [q for q in range(3, m + 1, 2) if m % q == 0]
        spectra[m] = {"pre_low_degree_q": qq,
                      "post_cited_low_degree_q": [q for q in qq if q >= 5],
                      "N1_pre_degrees": sorted({q + n for q in qq for n in (0, 1)}),
                      "Nge2_pre_degrees": [q + 1 for q in qq]}
    check(spectra[9]["pre_low_degree_q"] == [3, 9], "m9 two divisor controls")
    check(spectra[9]["post_cited_low_degree_q"] == [9], "m9 cited filter")
    check(spectra[3]["post_cited_low_degree_q"] == [], "ordinary cusp exclusion")

    payload = {"small_census": census, "named": named, "spectra": spectra}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    print("STATUS: FINITE-EXACT controls; analytic proof is in the companion note")
    print("composition: left actions; sigma*tau means sigma after tau")
    print("small_universe: k=1,2; r=0..9; all independent supported permutations")
    print(f"small_rows={total_rows}; reaccess_survivors={total_survivors}; prefilter_Artin=false")
    print("small_counts=" + json.dumps(census, separators=(",", ":")))
    print("canonical_universe: q=3,5,...,19; m=3,5,...,99; n=0,1")
    print(f"canonical_rows={canonical_rows}; divisor_survivors={canonical_survivors}")
    for row in named:
        print("positive_abstract=" + json.dumps(row, sort_keys=True, separators=(",", ":")))
    print("hostiles: inertia-only; h=2; n=2; uncompensated orientation; formal q=1")
    print("spectra=" + json.dumps(spectra, sort_keys=True, separators=(",", ":")))
    print("scope: necessary local+Euler passports only; no Keller realization or general JC closure")
    print(f"semantic_sha256={digest}")
    print(f"PASS gates={GATES}")


if __name__ == "__main__":
    main()
