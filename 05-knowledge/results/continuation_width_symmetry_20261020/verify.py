#!/usr/bin/env python3
"""Exact controls for continuation width and common symmetry support.

Standard library only; every check remains active under -O.
No theorem depends on the finite size of these controls.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb, factorial, gcd
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json

CHECKS = Counter()


def need(value, label):
    CHECKS[label] += 1
    if not value:
        raise RuntimeError(label)


def largest_prime(n):
    p, answer = 2, 1
    while p * p <= n:
        while n % p == 0:
            answer = p
            n //= p
        p += 1
    return max(answer, n)


def primes_to(n):
    return [j for j in range(2, n + 1) if largest_prime(j) == j]


def smooth_part(n, p):
    out = 1
    for q in primes_to(p):
        while n % q == 0:
            out *= q
            n //= q
    return out


def response_dp(n, order, core=False, no_fixed=False):
    lengths = [j for j in range(1, n + 1, 2) if not (no_fixed and j == 1)]
    if order == "prime":
        lengths.sort(key=lambda j: (largest_prime(j), j), reverse=True)
    elif order == "descending":
        lengths.reverse()
    elif order != "ascending":
        raise ValueError(order)
    states = {(0, (0,) * len(lengths)): F(1)}
    rows = []
    for pos, ell in enumerate(lengths):
        future = lengths[pos + 1:]
        new = defaultdict(F)
        transitions = 0
        for (degree, vec), weight in states.items():
            active = [j for j in lengths[pos:] if j <= n - degree]
            responses = dict(zip(active, vec))
            for k in range((n - degree) // ell + 1):
                e = k * (ell - 1) // 2 + ell * comb(k, 2)
                e += k * responses.get(ell, 0)
                factor = weight * F(2**e, ell**k * factorial(k))
                extra = range((n - degree) // ell - k + 1) if core and k else (0,)
                for b in extra:
                    d = degree + ell * (k + b)
                    nxt = tuple(responses[j] + k * gcd(ell, j)
                                for j in future if j <= n - d)
                    signed = (-1)**b * comb(k + b - 1, b) if k else 1
                    new[d, nxt] += factor * signed
                    transitions += 1
        states = {key: val for key, val in new.items() if val}
        rows.append((ell, len(states), transitions))
    values = [F(0)] * (n + 1)
    for (d, vec), value in states.items():
        need(not vec, "DP terminal response empty")
        values[d] += value
    return values, {"peak": max([1] + [x[1] for x in rows]),
                    "transitions": sum(x[2] for x in rows),
                    "rows": rows}


def prime_tests():
    rows = []
    for L in (5, 10, 12, 20, 30, 50):
        ps = [p for p in primes_to(2 * L) if p > L]
        cells = defaultdict(list)
        all_profiles = set()
        for mask in range(1 << len(ps)):
            subset = tuple(p for i, p in enumerate(ps) if mask >> i & 1)
            past = tuple(3 * p for p in subset)
            vec = tuple(sum(gcd(a, p) for a in past) for p in [1] + ps)
            expected = (len(subset),) + tuple(
                len(subset) + (p - 1) * (p in subset) for p in ps)
            need(vec == expected, "prime subset response identity")
            need(3 * sum(subset) + sum(ps) <= 4 * sum(ps),
                 "all prime continuations affordable")
            all_profiles.add(vec)
            cells[len(subset), sum(subset)].append(subset)
        need(len(all_profiles) == 2**len(ps), "prime profiles all distinct")
        cell = max(cells.values(), key=len)
        # Exact LDL pivots of the principal tensor-kernel matrix.
        mat = [[F(2**sum(p - 1 for p in set(a) & set(b)))
                for b in cell] for a in cell]
        pivots = []
        for i in range(len(mat)):
            pivot = mat[i][i]
            need(pivot > 0, "continuation kernel exact positive pivot")
            pivots.append(pivot)
            for j in range(i + 1, len(mat)):
                for k in range(j, len(mat)):
                    mat[j][k] -= mat[j][i] * mat[k][i] / pivot
                    mat[k][j] = mat[j][k]
        rows.append({"L": L, "primes": ps, "horizon": 4 * sum(ps),
                     "distinct_profiles": len(all_profiles),
                     "largest_common_mass_cardinality_cell": len(cell),
                     "cell": cell, "rank": len(pivots)})
    # Verify the old/new smooth-kernel representation at every prime batch.
    for n in range(3, 85, 2):
        for p in primes_to(n):
            if p == 2:
                continue
            future = [j for j in range(1, n + 1, 2) if largest_prime(j) <= p]
            for a in range(3, n + 1, 2):
                if largest_prime(a) > p:
                    c = smooth_part(a, p)
                    need(a > p * c, "old smooth kernel has small mass")
                    need(all(gcd(a, j) == gcd(c, j) for j in future),
                         "smooth kernel preserves every future gcd")
                elif largest_prime(a) == p:
                    need(a == p * (a // p), "current batch exact cofactor")
    return rows


def adjacency(n, mask):
    out = [0] * n
    for i, (a, b) in enumerate(combinations(range(n), 2)):
        if mask >> i & 1:
            a, b = b, a
        out[a] |= 1 << b
    return out


def feedback(adj):
    n = len(adj)
    best, last = [0] * (1 << n), [-1] * (1 << n)
    for mask in range(1, 1 << n):
        opts = [(best[mask ^ (1 << v)] +
                 (adj[v] & (mask ^ (1 << v))).bit_count(), v)
                for v in range(n) if mask >> v & 1]
        best[mask], last[mask] = min(opts)
    order, mask = [], (1 << n) - 1
    while mask:
        order.append(last[mask])
        mask ^= 1 << last[mask]
    return best[-1], tuple(reversed(order))


def automorphisms(adj):
    n = len(adj)
    groups = defaultdict(list)
    for v, row in enumerate(adj):
        groups[row.bit_count()].append(v)
    groups = list(groups.values())
    out = []
    for images in product(*(permutations(g) for g in groups)):
        perm = list(range(n))
        for domain, image in zip(groups, images):
            for a, b in zip(domain, image):
                perm[a] = b
        if all(((adj[a] >> b) & 1) == ((adj[perm[a]] >> perm[b]) & 1)
               for a in range(n) for b in range(a + 1, n)):
            out.append(tuple(perm))
    return out


def check_tournament(adj, brute=False):
    n = len(adj)
    f, order = feedback(adj)
    aut = automorphisms(adj)
    orbits = sorted({tuple(sorted({p[v] for p in aut})) for v in range(n)})
    owner = {v: i for i, orb in enumerate(orbits) for v in orb}
    active = [o for o in orbits if len(o) > 1]
    s = sum(map(len, active))
    c = sum((len(o)**2 - 1) // 8 for o in active)
    delta = 3 * f - s
    need(all(len(o) % 2 for o in orbits), "whole group orbits odd")
    need(f >= c and 3 * c >= s, "common support quadratic budget")
    need(delta not in (1, 2), "symmetry slack one and two absent")
    need(all(len({(adj[v] & sum(1 << x for x in o)).bit_count()
                  for v in o}) == 1 for o in orbits),
         "whole group orbit regularity")
    pos = {v: i for i, v in enumerate(order)}
    backwards = [(a, b) for a in range(n) for b in range(n)
                 if adj[a] >> b & 1 and pos[a] > pos[b]]
    external = [(a, b) for a, b in backwards if owner[a] != owner[b]]
    marked_orbits = {i for i, o in enumerate(orbits) if len(o) >= 5}
    marked_orbits.update(owner[v] for edge in external for v in edge)
    marked = {v for i in marked_orbits for v in orbits[i]}
    need(len(marked) <= 2 * delta, "exceptional skeleton bound")
    quiet = [o for i, o in enumerate(orbits) if i not in marked_orbits and len(o) == 3]
    for o in quiet:
        places = sorted(pos[v] for v in o)
        need(places[-1] - places[0] == 2, "quiet triangle contiguous")
        need(all(len({(adj[v] >> x) & 1 for v in o}) == 1
                 for x in range(n) if x not in o), "quiet triangle module")
    # Reconstruct the entire labelled tournament, including every port identity.
    quiet_id = {v: i for i, o in enumerate(quiet) for v in o}
    rebuilt = [0] * n
    for a, b in combinations(range(n), 2):
        if a in marked and b in marked:
            direction = (adj[a] >> b) & 1
        elif a in quiet_id and quiet_id.get(b) == quiet_id[a]:
            direction = (adj[a] >> b) & 1
        else:
            direction = pos[a] < pos[b]
        if direction:
            rebuilt[a] |= 1 << b
        else:
            rebuilt[b] |= 1 << a
    need(rebuilt == adj, "skeleton rebuild preserves all labelled arcs and ports")
    if delta == 0:
        k = len(quiet)
        need(not marked and f == k and s == 3 * k,
             "zero slack word classification")
        need(len(aut) == 3**k, "zero slack automorphism product")
        fixed_free = sum(all(p[v] != v for v in range(n)) for p in aut)
        need(fixed_free == (2**k if n == 3 * k else 0),
             "zero slack fixed free mass")
    if brute:
        brute_f = min(sum(bool(adj[p[j]] >> p[i] & 1)
                          for i in range(n) for j in range(i + 1, n))
                      for p in permutations(range(n)))
        need(brute_f == f, "feedback independent order enumeration")
    return f, s, delta, len(marked), len(aut)


def tournament_tests():
    counts, slack_counts = [], Counter()
    examples = {}
    for n in range(7):
        total = 1 << comb(n, 2)
        zero = 0
        for mask in range(total):
            result = check_tournament(adjacency(n, mask), brute=n <= 5)
            f, s, delta, marked, aut = result
            slack_counts[n, delta] += 1
            examples.setdefault(delta, {"n": n, "mask": mask, "result": result})
            zero += delta == 0
        expected = sum(comb(n - 2*k, k) * factorial(n) // 3**k
                       for k in range(n // 3 + 1))
        need(zero == expected, "zero slack labelled census matches class series")
        counts.append({"n": n, "labelled_tournaments": total, "zero_slack": zero})
    # Larger positive controls, with actual group computation still cheap because
    # the score sets of ordered components are disjoint.
    for word in [(3,) * k for k in range(1, 5)] + [(1, 3, 1, 3, 1), (5, 3), (3, 5, 1)]:
        n = sum(word)
        adj = adjacency(n, 0)
        start = 0
        for size in word:
            if size > 1:
                for i in range(size):
                    for j in range(i + 1, size):
                        a, b = start+i, start+j
                        if (j-i) % size > size//2:
                            adj[a] ^= 1 << b
                            adj[b] ^= 1 << a
            start += size
        check_tournament(adj)
    return {"universes": counts, "slack_histograms": {
        str(n): {str(d): v for (nn, d), v in slack_counts.items() if nn == n}
        for n in range(7)}, "minimal_examples": examples}


def compose_exponents(word):
    A, C = F(1), F(0)
    for a in word:
        alpha, beta = F(3, 2**a), F(4, 2**a)-1
        A, C = alpha*A, alpha*C+beta
    return A, C


def gap_form(word):
    gaps, bad = [0], []
    for a in word:
        if a == 2:
            gaps[-1] += 1
        else:
            bad.append(a)
            gaps.append(0)
    A = F(3, 4)**sum(gaps)
    for b in bad:
        A *= F(3, 2**b)
    C = Cp = F(0)
    for i, b in enumerate(bad):
        suffix = F(3, 4)**sum(gaps[i+1:])
        for bb in bad[i+1:]:
            suffix *= F(3, 2**bb)
        beta = F(4, 2**b)-1
        C += beta*suffix
        Cp += max(beta, F(0))*suffix
    return A, C, Cp, bad, gaps


def collatz_tests():
    formal = 0
    for ell in range(8):
        for word in product(range(1, 5), repeat=ell):
            A, C = compose_exponents(word)
            AA, CC, Cp, bad, gaps = gap_form(word)
            need((A, C) == (AA, CC), "affine gap coefficients exact")
            need(C <= Cp, "positive carry envelope")
            need(sum(word) == 2*sum(gaps)+sum(bad), "integer clock retained")
            need(F(3,4)**sum(gaps) ==
                 __import__("functools").reduce(
                     lambda x, r: x*F(3,4)**r, gaps, F(1)),
                 "quiet refinement preserves multiplier")
            formal += 1
    block_pass, block_fail = 0, 0
    for s in range(1, 9):
        for r in range(13):
            for a in range(3, 13):
                word = (1,)*s+(2,)*r+(a,)
                A, C = compose_exponents(word)
                u, t = F(3,2)**s, F(3,4)**r
                need(A == 3*t*u/2**a and C == (6*t*(u-1)+4)/2**a-1,
                     "contracting block exact coefficients")
                if 2**a >= 6*t*(u-1)+4:
                    need(A < 1 and C <= 0, "block threshold implies contraction")
                    block_pass += 1
                else:
                    block_fail += 1
    # Actual sources and every first 24 accelerated steps (including terminal 1s).
    authenticated = 0
    for source in range(1, 2002, 2):
        n, word = source, []
        for depth in range(24):
            m = 3*n+1
            a = (m & -m).bit_length()-1
            n = m >> a
            word.append(a)
            A, C, Cp, _, _ = gap_form(word)
            need(A*(source-1)+C == n-1, "actual source terminal preserved")
            if A*(source-1)+Cp < source-1:
                need(n < source, "source-specific margin certifies descent")
            authenticated += 1
    # One defect after a quiet run has an arbitrarily large lawful source.
    # n=2^(2r+1)+1 has exactly r twos, then one one.
    quiet_hostiles = []
    for r in range(1, 31):
        source, n = 2**(2*r+1)+1, 2**(2*r+1)+1
        got = []
        for _ in range(r+1):
            m = 3*n+1
            a = (m & -m).bit_length()-1
            got.append(a)
            n = m >> a
        need(got == [2]*r+[1], "fixed defect count has unbounded quiet clock")
        quiet_hostiles.append({"r":r, "source":source, "terminal":n})
    rotation_hostiles = []
    for L in range(1, 33):
        w = (1,)+(0,)*L
        word = w+w
        shift = L+1
        need(sum(word) == 2 and word[shift:]+word[:shift] == word,
             "two binary defects allow unbounded moving symmetry")
        rotation_hostiles.append((L,len(word)))
    return {"formal_words": formal, "block_threshold_pass":block_pass,
            "block_threshold_fail":block_fail, "authenticated_prefixes":authenticated,
            "quiet_clock_hostiles":quiet_hostiles, "rotation_hostiles":rotation_hostiles}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--bound", type=int, default=100)
    ap.add_argument("--large-bound", type=int, default=0)
    args=ap.parse_args()
    prime=prime_tests()
    graphs=tournament_tests()
    collatz=collatz_tests()
    prior_path=Path(__file__).resolve().parent.parent/"tournament_future_boundary_20261020"/"verify.py"
    spec=importlib.util.spec_from_file_location("prior_future_checker",prior_path)
    prior=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)
    dp={}
    for bound in (30,60,args.bound):
        prime_values, info=response_dp(bound,"prime")
        old_values, old_info=prior.profile_dp(bound)
        need(prime_values==old_values,"prime-order counts equal inherited descending counts")
        if bound<=30:
            need(prime_values==prior.direct_series(bound),"independent direct Burnside sum")
        need(all(v.denominator==1 for v in prime_values),"A000568 output integral")
        dp[str(bound)]={"prime_order":info,"descending_order":old_info,
                        "coefficient_sha256":hashlib.sha256(
                            json.dumps(list(map(int,prime_values))).encode()).hexdigest(),
                        "last_value":str(prime_values[-1])}
    for core,no_fixed in ((True,False),(True,True)):
        val,_=response_dp(24,"prime",core=core,no_fixed=no_fixed)
        ref,_=prior.profile_dp(24,core=core,no_fixed=no_fixed)
        need(val==ref,"prime ordering preserves signed core series")
    if args.large_bound:
        val,info=response_dp(args.large_bound,"prime")
        need(val[:args.bound+1]==prime_values,"large prime run preserves smaller coefficients")
        dp[str(args.large_bound)]={"prime_order":info,"last_value":str(val[-1]),
                                  "last_value_sha256":hashlib.sha256(str(val[-1]).encode()).hexdigest()}
    result={"status":"FINITE-EXACT controls, not all-order proof",
            "checks":dict(CHECKS),"total_checks":sum(CHECKS.values()),
            "prime_separation":prime,"tournaments":graphs,
            "collatz":collatz,"dynamic_programs":dp,
            "reference_checker_sha256":hashlib.sha256(prior_path.read_bytes()).hexdigest()}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
