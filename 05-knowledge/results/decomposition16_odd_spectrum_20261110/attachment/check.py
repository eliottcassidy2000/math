#!/usr/bin/env python3
"""Independent graph-level controls for the H=3 attachment compiler."""
import hashlib
import itertools
import json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from pathlib import Path


def require(ok, message):
    if not ok:
        raise ValueError(message)


@lru_cache(None)
def g(n):
    if n < 2:
        return 1
    q, bit = divmod(n, 2)
    return g(q) if bit else g(q) + q + (q & 1)


def chain_word_count(bits):
    """Even-order increasing paths between a 1 and a 0, by endpoint gaps."""
    f = sum(bits[i] * (1 - bits[j]) * 2 ** (j - i - 1)
            for i in range(len(bits)) for j in range(i + 1, len(bits)))
    d = sum(bits[i] * (1 - bits[i + 1]) for i in range(len(bits) - 1))
    return 1 + f + d


def integer(bits):
    value = 0
    for bit in bits:
        value = 2 * value + bit
    return value


def word(n, width=None):
    if width is None:
        width = n.bit_length()
    return [(n >> i) & 1 for i in reversed(range(width))]


def c(a, b):
    return 12*a*b + 6*a + 3*b + 3*g(a) + 3*g(b) + 2*(a & 1) + (b & 1) - 1


def u(a, b):
    return 12*a*b + 9*a + 3*g(a) + 3*g(b) + 3*(a & 1) - 3


def base(a, b):
    n = a + b + 3
    rows = [sum(1 << j for j in range(i + 1, n)) for i in range(n)]
    rows[a] ^= 1 << (a + 2)
    rows[a + 2] ^= 1 << a
    return tuple(rows)


def attach(rows, bits):
    n = len(rows)
    require(len(bits) == n, "cut size")
    return tuple(rows[i] | ((1 - bits[i]) << n) for i in range(n)) + (
        sum(bit << i for i, bit in enumerate(bits)),)


def canonical(a, b, k=1):
    left = word(a)
    right = [1 - bit for bit in reversed(word(b))]
    bits = left + [int(i < k) for i in range(3)] + right
    return attach(base(len(left), len(right)), bits)


@lru_cache(None)
def hp(rows):
    n = len(rows)
    if n == 0:
        return 1
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            count = dp[mask][v]
            if not count:
                continue
            todo = rows[v] & (full ^ mask)
            while todo:
                bit = todo & -todo
                todo ^= bit
                dp[mask | bit][bit.bit_length()-1] += count
    return sum(dp[-1])


def literal_hp(rows):
    n = len(rows)
    return sum(all((rows[p[i]] >> p[i+1]) & 1 for i in range(n-1))
               for p in itertools.permutations(range(n)))


def tournament(rows):
    n = len(rows)
    return all(not (rows[i] >> i) & 1 for i in range(n)) and all(
        ((rows[i] >> j) & 1) + ((rows[j] >> i) & 1) == 1
        for i in range(n) for j in range(i+1, n))


def strong(rows):
    full = (1 << len(rows))-1
    for i in range(len(rows)):
        reached, frontier = 1 << i, 1 << i
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            new = rows[bit.bit_length()-1] & ~reached
            reached |= new
            frontier |= new
        if reached != full:
            return False
    return True


def predicted(a, b, bits):
    left, tri, right = bits[:a], bits[a:a+3], bits[a+3:]
    aa = integer(left)
    bb = integer([1 - x for x in reversed(right)])
    k = sum(tri)
    if k == 1:
        return c(aa, bb)
    if k == 2:
        return c(bb, aa)
    return u(aa, bb) if k == 0 else u(bb, aa)


def target_pairs(target):
    if target < 5:
        return []
    return [(a, b) for a in range((target-5)//6 + 1)
            for b in range((target-5-6*a)//(12*a+3) + 1)
            if c(a, b) == target]


def histogram(a, b):
    hist = Counter()
    for aa in range(1 << a):
        for bb in range(1 << b):
            left = word(aa, a)
            right = [1 - bit for bit in reversed(word(bb, b))]
            hist[u(aa, bb)] += 1
            hist[u(bb, aa)] += 1
            require(u(aa, bb)-3*g(integer(left+[0]+right)) == 6*aa*(bb+1),
                    "uniform reentry correction")
            require(u(bb, aa)-3*g(integer(left+[1]+right)) == 6*bb*(aa+1),
                    "converse uniform reentry correction")
            hist[c(aa, bb)] += 3
            hist[c(bb, aa)] += 3
    return hist


def main():
    counts = Counter()
    histograms = {}
    for total in range(5):
        for a in range(total+1):
            b = total-a
            old = base(a, b)
            require(tournament(old) and hp(old) == 3, "source")
            actual = Counter()
            for bits in itertools.product((0, 1), repeat=len(old)):
                child = attach(old, bits)
                got = hp(child)
                require(tournament(child), "child matrix")
                require(got == predicted(a, b, list(bits)), "mixed/uniform formula")
                actual[got] += 1
                counts["all_cuts"] += 1
                if len(child) <= 6:
                    require(got == literal_hp(child), "literal paths")
                    counts["literal_children"] += 1
            require(actual == histogram(a, b), "histogram")
            require(sum(actual.values()) == 2**len(old), "cut multiplicity")
            histograms[f"{a},{b}"] = dict(sorted(actual.items()))
            counts["parents"] += 1

    for n in range(1024):
        require(g(n) == chain_word_count(word(n)), "one-apex recurrence")
        counts["independent_g"] += 1

    for a in range(12):
        for b in range(12):
            for k in (1, 2):
                graph = canonical(a, b, k)
                require(tournament(graph) and strong(graph), "strong canonical graph")
                require(hp(graph) == (c(a, b) if k == 1 else c(b, a)), "canonical formula")
                counts["canonical_graphs"] += 1

    # Complete target image <= 501 by the proved inequality (A5), not a
    # vertex-order cutoff. Every candidate is independently counted as a graph.
    cap = 501
    candidates = []
    for a in range((cap-5)//6 + 1):
        for b in range((cap-5-6*a)//(12*a+3) + 1):
            graph = canonical(a, b)
            value = hp(graph)
            require(value == c(a, b), "target graph count")
            require(value >= 12*a*b + 6*a + 3*b + 5, "inverse bound")
            candidates.append((a, b, value))
    counts["inverse_candidates"] = len(candidates)
    by_target = {}
    for target in range(1, cap+1, 2):
        direct = [(a, b) for a, b, value in candidates if value == target]
        require(target_pairs(target) == direct, "inverse completeness")
        if direct:
            by_target[target] = direct
    require(target_pairs(11) == [], "all-order 3 to 11 obstruction")
    require(target_pairs(13) == [(1, 0)], "3 to 13 normal form")
    require(target_pairs(7) == [] and target_pairs(21) == [], "forbidden controls")
    require(hp((30,20,24,18,32,15)) == 27 and hp((6,4,8,3)) == 5,
            "uniform-module collapse hostile")
    hostile_path = (1,2,4,5,0,3)
    hostile_rows = (30,20,24,18,32,15)
    require(all(hostile_rows[hostile_path[i]] >> hostile_path[i+1] & 1
                for i in range(5)), "literal reentry path")

    # All-order one-apex support up to cap: strip trailing one bits using
    # g(2q+1)=g(q); an even representative 2q with count<=cap has q<=cap-1.
    apex = {g(2*q) for q in range(cap) if g(2*q) <= cap}
    uniform = {3*v for v in apex if 3*v <= cap}
    # A=0 gives precisely 3g(B); A>=1 has the stated finite bound.
    for a in range(1, (cap-3)//9 + 1):
        for b in range((cap-3-9*a)//(12*a) + 1):
            value = u(a, b)
            require(hp(canonical(a, b, 0)) == value, "uniform inverse graph")
            counts["uniform_inverse_candidates"] += 1
            if value <= cap:
                uniform.add(value)
    full = set(by_target) | uniform
    require(11 not in full and 13 in full, "family control")
    require([h for h in range(1, 56, 2) if h not in full | apex] == [7, 21, 55],
            "first extra two-parent-lane obstruction")
    witness55 = (54,124,56,113,96,64,5)
    require(tournament(witness55) and hp(witness55) == literal_hp(witness55) == 55,
            "count55 witness")
    deletion55 = []
    for v in range(7):
        keep = [i for i in range(7) if i != v]
        sub = tuple(sum(((witness55[i] >> j) & 1) << k for k, j in enumerate(keep))
                    for i in keep)
        value = hp(sub)
        require(value == literal_hp(sub), "literal55 deletions")
        deletion55.append(value)
    require(deletion55 == [9,25,15,19,29,29,5], "sharp minimum deletion count5")
    diagnostic = {
        "cap": cap,
        "sharp_55_witness": {"rows": witness55, "H": 55, "deletions": deletion55},
        "mixed_support_size": len(by_target),
        "all_H3_source_support_size": len(full),
        "missing_odds_3_to_99": [h for h in range(3, 100, 2) if h not in full],
        "apex_plus_H3_missing_odds_1_to_199": [h for h in range(1, 200, 2)
                                                  if h not in full | apex],
        "selected_pairs": {str(h): target_pairs(h) for h in (5, 7, 9, 11, 13, 17, 19, 21, 23, 25, 29, 31, 37)},
        "strong_core_polynomials_nonmultiple3": {
            str(h): dict(sorted(Counter({degree: 2*sum(
                4+a.bit_length()+b.bit_length() == degree for a, b in target_pairs(h))
                for degree in sorted({4+a.bit_length()+b.bit_length() for a, b in target_pairs(h)})}).items()))
            for h in (5, 11, 13, 17, 19, 23, 25, 29, 31, 37)},
    }
    moments = {}
    for key in ("2,0", "1,1"):
        hist = histograms[key]
        mass = sum(hist.values())
        mean = sum(Fraction(h*w, mass) for h, w in hist.items())
        variance = sum(Fraction(h*h*w, mass) for h, w in hist.items()) - mean*mean
        moments[key] = {"mean": str(mean), "variance": str(variance)}
    require(moments == {"2,0": {"mean": "27/2", "variance": "279/4"},
                        "1,1": {"mean": "27/2", "variance": "351/4"}}, "inherited hostile")
    require(c(1, 0) - 2 == 11 and hp(canonical(1, 0)) == 13,
            "removing parity correction fabricates forbidden source-child11")

    here = Path(__file__).resolve().parent
    root = here.parents[3]
    pins = [
        "01-canon/theorems/THM-002-ocf.md",
        "01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md",
        "05-knowledge/results/session207_degree_three_metagraph_20261106/PROOF.md",
        "05-knowledge/results/tournament_joint_parity_20261102/growth_obstruction.md",
    ]
    output = {"status": "PASS", "checks": dict(sorted(counts.items())),
              "histograms": histograms, "moments": moments, "target_inverse": diagnostic,
              "source_sha256": {p: hashlib.sha256((root/p).read_bytes()).hexdigest() for p in pins},
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "proof_sha256": hashlib.sha256((here/"PROOF.md").read_bytes()).hexdigest()}
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
