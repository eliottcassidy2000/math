"""Finite translated-grid compiler with exact pair-component and gcd sidecars.

Standalone, integer/rational arithmetic only, no producer imports. The full
survivor arrays are intentionally retained in the output, not just maxima.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from pathlib import Path
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
ROOT = Path(__file__).resolve().parents[1]
LIMIT = 97096
CAPS = (90, 30, 9, 4, 2, 1, 1)
PROFILE_PATH = ROOT / "05-knowledge/results/lrc14_joint_shadow_empty_core_next_sep06.json"
PROFILE_SHA = "935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f"
GATES = 0


def need(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise ArithmeticError(message)


def ceiling(n, d):
    return (n + d - 1) // d


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(",", ":")).encode()).hexdigest()


def inert_sum(n):
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            if p % 3 != 2 or e > 2:
                return False
        p += 1
    return n == 1 or n % 3 == 2


def bernoulli(x):
    x -= x.numerator // x.denominator
    return x * x - x + F(1, 6)


def bad(x, u):
    a = x.numerator * u % x.denominator
    return 14 * min(a, x.denominator - a) < x.denominator


def literal_lengths(p, q):
    cuts = {F(0), F(1)}
    for u in (p, q):
        for k in range(u):
            for s in (-1, 1):
                x = F(14 * k + s, 14 * u)
                cuts.add(x - x.numerator // x.denominator)
    cuts = sorted(cuts)
    flags = [bad((a + b) / 2, p) and bad((a + b) / 2, q)
             for a, b in zip(cuts, cuts[1:])]
    pieces = [b - a for a, b, yes in zip(cuts, cuts[1:], flags) if yes]
    if flags[0] and flags[-1]:
        pieces = [pieces[0] + pieces[-1]] + pieces[1:-1]
    return Counter(pieces)


def atlas():
    rows = []
    literal = 0
    for p in range(1, 178):
        for q in range(p + 1, 357 - p):
            if gcd(p, q) != 1 or not inert_sum(p + q):
                continue
            R = (p + q - 1) // 14
            J, D = 2 * R + 1, 14 * p * q
            lengths = Counter({2 * p: 1})
            for r in range(1, R + 1):
                lengths[min(2 * p, p + q - 14 * r)] += 2
            C = tuple(sorted(lengths.items()))
            M = sum(n * mult for n, mult in C)
            mu = F(M, D)
            need(sum(lengths.values()) == J, "all strict open components retained")
            need(mu == F(1, 49) + (bernoulli(F(q - p, 14)) - bernoulli(F(q + p, 14))) / (p * q),
                 "literal lengths agree with independent Bernoulli measure")
            # These three exact inequalities certify the entire lower envelope.
            need(J <= 51, "envelope at x=0")
            need(mu >= F(1, 70), "envelope slope at infinity")
            need(304500 * M >= (4313 + 37 * J) * D, "envelope at x=304500/37")
            if q <= 30 or (p, q) in ((1, 355), (5, 348), (61, 292), (25, 294)):
                claimed = Counter({F(n, D): mult for n, mult in C})
                need(literal_lengths(p, q) == claimed, "independent wall sweep checks full length multiset")
                literal += 1
            rows.append((p, q, J, D, C, M))
    need(len(rows) == 5855, "complete actual inert pair universe")
    need(len({p + q for p, q, *_ in rows}) == 94, "all94 actual inert sums")
    for p, q, J, D, C, M in rows:
        if (p, q) == (1, 10):
            need(F(M, D) == F(1, 70) and J == 1, "low-slope envelope attained")
        if (p, q) == (5, 348):
            need(F(M, D) == F(62, 3045) and J == 51, "high-component envelope attained")
    # Search order changes runtime only; every rejected candidate is paid for.
    rows.sort(key=lambda row: (row[:2] not in ((1, 10), (5, 348)), -row[2], F(row[5], row[3])))
    return rows, literal


def aggregate(t, e):
    return e * min(ceiling(t, 70 * e) - 1, ceiling(62 * t, 3045 * e) - 51)


def pair_aggregate(t, e, row):
    p, q, J, D, C, M = row
    return e * (ceiling((t // e) * M, D) - J)


def component(t, e, row):
    p, q, J, D, C, M = row
    n = t // e
    return e * (sum(mult * ceiling(n * c, D) for c, mult in C) - J)


def bag(t, domain):
    avail = {d: sum(d <= c for c in CAPS) for d in domain if t % d == 0}
    costs = {d: d * ((-(t // d)) % 7) for d in avail}
    for d, cost in costs.items():
        need(cost == 7 * d * ceiling(t, 7 * d) - t, "exact ceiling cost at retained clock divisor")
    entries = sorted(((costs[d], d) for d in avail for _ in range(avail[d])), reverse=True)
    value = sum(w for w, _ in entries[:7])
    need(len(entries) >= 7 and value % 7 == 0, "seven retained costs have integral excess")
    return avail, costs, entries, value // 7


def forced_excess(dp, dq, avail, costs, entries):
    if dp not in avail or dq not in avail or (dp == dq and avail[dp] < 2):
        return None
    used = Counter((dp, dq))
    other = []
    for w, d in entries:
        if used[d]:
            used[d] -= 1
        else:
            other.append(w)
            if len(other) == 5:
                break
    need(len(other) == 5, "five slots remain after actual pair marginal reservation")
    value = costs[dp] + costs[dq] + sum(other)
    need(value % 7 == 0, "forced pair and five remaining costs retain integral excess")
    return value // 7


def compiler(rows, domain, cap, label):
    aggregate_set, component_set, coupled_set = [], [], []
    witnesses = {}
    linear_set = []
    for t in range(1, LIMIT + 1):
        avail, costs, entries, E = bag(t, domain)
        edivs = [e for e in range(cap, 0, -1) if t % e == 0]
        maxe = edivs[0]
        linear_survives = min(F(t, 70) - maxe, F(62 * t, 3045) - 51 * maxe) <= E
        if linear_survives:
            linear_set.append(t)
        es = [e for e in edivs if aggregate(t, e) <= E]
        need(not es or linear_survives, "aggregate ceilings strengthen the linear envelope")
        if not es:
            continue
        aggregate_set.append(t)
        has_component = False
        has_coupled = False
        cache = {}
        for e in es:
            for row in rows:
                A = pair_aggregate(t, e, row)
                if A > E:
                    continue
                C = component(t, e, row)
                need(C >= A, "sum of component ceilings dominates ceiling of total measure")
                if C > E:
                    continue
                has_component = True
                p, q = row[:2]
                dp, dq = e * gcd(t // e, p), e * gcd(t // e, q)
                need(gcd(dp, dq) == e, "forced pair margins retain their exact common gcd")
                key = tuple(sorted((dp, dq)))
                if key not in cache:
                    cache[key] = forced_excess(dp, dq, avail, costs, entries)
                EE = cache[key]
                if EE is None:
                    continue
                need(EE <= E, "pair-coupled bag refines the unconditioned bag")
                if C <= EE:
                    has_coupled = True
                    witnesses[t] = [E, EE, e, p, q, dp, dq, C]
                    break
            if has_coupled:
                break
        if has_component:
            component_set.append(t)
        if has_coupled:
            coupled_set.append(t)
    need(set(coupled_set) <= set(component_set) <= set(aggregate_set) <= set(linear_set),
         "all four retained finite relaxations are nested")
    return {
        label + "_aggregate": aggregate_set,
        label + "_components": component_set,
        label + "_coupled": coupled_set,
    }, {
        "label": label,
        "linear_count": len(linear_set), "linear_max": max(linear_set),
        "coupling_excluded": sorted(set(component_set) - set(coupled_set)),
        "last20_coupled_witnesses": [[t, witnesses[t]] for t in coupled_set[-20:]],
    }


def hostile_controls(rows):
    lookup = {row[:2]: row for row in rows}
    need([aggregate(57995, e) for e in (1, 5, 7)] == [828, 825, 826],
         "rounding is not monotone in the numerical size of a divisor")
    need(component(74550, 30, lookup[(1, 355)]) == 0,
         "open components of exactly one quotient-grid spacing may contain no point")
    need(30 * gcd(74550 // 30, 355) == 10650 > 90,
         "endpoint witness destroyed by its forced pair marginal")
    need(aggregate(27360, 6) == 252, "strict aggregate equality boundary")
    need(component(27360, 6, lookup[(5, 348)]) == 294,
         "per-component ceilings recover credit lost by merging lengths")
    # Literal integer h controls establish the pair-marginal identity beyond h=e.
    for t in range(1, 81):
        for h in range(1, 41):
            e = gcd(t, h)
            need(gcd(t // e, h // e) == 1, "coprime scale quotients")
            for p, q in ((1, 10), (5, 348), (61, 292), (25, 294)):
                need(gcd(t, h * p) == e * gcd(t // e, p), "literal first pair marginal")
                need(gcd(t, h * q) == e * gcd(t // e, q), "literal second pair marginal")
    # Two-line measure envelope cannot replace all shapes after individual ceilings.
    t, e = 23760, 6
    pair = lookup[(25, 294)]
    need(component(t, e, pair) == 246,
         "retained maximum-scale component witness")
    need(component(t, e, pair) < min(component(t, e, lookup[x]) for x in ((1, 10), (5, 348))),
         "full component consumer needs shapes discarded by the two-line measure quotient")


def main():
    raw = PROFILE_PATH.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == PROFILE_SHA, "pinned full profile supplier")
    profile = json.loads(raw)
    allowed = profile["levels"]["6"]["gcds"]
    need(len(allowed) == 42 and max(allowed) == 90 and allowed == sorted(set(allowed)),
         "exact seven-body value projection")
    rows, literal = atlas()
    hostile_controls(rows)
    sets = {}
    summaries = []
    for domain, cap, label in ((list(range(1, 91)), 30, "scalar30"), (allowed, 6, "profile6")):
        out, summary = compiler(rows, domain, cap, label)
        sets.update(out)
        summaries.append(summary)
    expected = {
        "scalar30_aggregate": (34532, 88920),
        "scalar30_components": (32294, 74550),
        "scalar30_coupled": (32272, 74520),
        "profile6_aggregate": (9498, 27360),
        "profile6_components": (8308, 23760),
        "profile6_coupled": (8301, 23760),
    }
    for name, values in sets.items():
        need((len(values), max(values)) == expected[name], "declared complete finite census: " + name)
        need(values == sorted(set(values)) and min(values) == 1, "canonical full survivor array")
    need(summaries[1]["coupling_excluded"] == [12425, 14872, 14910, 15390, 15504, 20520, 21240],
         "all seven additional profile-mode pair-coupling exclusions")
    print("FINITE-EXACT SCOPE: complete clock1..97096 relaxation; survivor does not mean unsafe or realizable.")
    print("PROVED REDUCTIONS: two-line atlas envelope, individual open-component ceilings, exact forced pair marginals.")
    print("PROFILE6 SUPPLIER: third_20260906_decoder_profile_graph.md; some actual edge e<=6, not every edge.")
    print("Pinned profile source SHA256:", PROFILE_SHA)
    print("Atlas:5855 ratios,94 sums; independent literal length sweeps:", literal)
    print("Envelope: min(x/70-1,62*x/3045-51); breakpoint304500/37; attained by(1,10),(5,348).")
    for summary in summaries:
        print("SUMMARY", json.dumps(summary, separators=(",", ":")))
    for name, values in sets.items():
        print("SCALE_SET", json.dumps({"name": name, "count": len(values), "maximum": max(values),
                                       "sha256": digest(values), "survivors": values}, separators=(",", ":")))
    print("Exact gates:", GATES)
    print("PASS")


if __name__ == "__main__":
    main()
