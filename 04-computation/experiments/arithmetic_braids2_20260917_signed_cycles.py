"""Signed 3n+b parameter strata and exact bounded-word cycle census.

Run: python 04-computation/experiments/arithmetic_braids2_20260917_signed_cycles.py
All checks remain active under -O. No bounded census asserts completeness.
"""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import comb, gcd
from pathlib import Path
import argparse
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def valuation(n, p):
    require(n != 0, "valuation of zero")
    n, exponent = abs(n), 0
    while n % p == 0:
        n //= p
        exponent += 1
    return exponent


def step(n, b):
    numerator = 3 * n + b
    require(numerator != 0, "orbit leaves nonzero odd domain")
    exponent = valuation(numerator, 2)
    return numerator // (1 << exponent), exponent


def canonical(values):
    values = list(values)
    start = min(range(len(values)), key=lambda i: (abs(values[i]), values[i]))
    return tuple(values[start:] + values[:start])


def walk_word(n, b, word):
    values = []
    initial = n
    for k in word:
        values.append(n)
        n, actual = step(n, b)
        require(actual == k, "word does not give exact halving exponents")
    require(n == initial, "word does not close")
    # A repeated word may traverse the same primitive cycle several times.
    first_return = values[1:].index(initial) + 1 if initial in values[1:] else len(values)
    primitive = values[:first_return]
    require(len(set(primitive)) == len(primitive), "nonprimitive extracted cycle")
    require(values == (primitive * (len(values) // first_return)), "incorrect period reduction")
    return canonical(primitive)


def word_states(max_length, max_total, word=(), total=0, carry=0):
    # Prefix recursion independently visits every positive composition;
    # the combinatorial total is sum C(max_total,L).
    if word:
        yield word, total, carry
    if len(word) == max_length:
        return
    next_carry = 3 * carry + (1 << total)
    for k in range(1, max_total - total + 1):
        yield from word_states(max_length, max_total, word + (k,), total + k, next_carry)


def affine_fraction_control(word, b):
    # Independent Fraction affine composition, not the integer carry recurrence.
    slope, intercept = Fraction(1), Fraction(0)
    for k in word:
        slope = 3 * slope / (1 << k)
        intercept = (3 * intercept + b) / (1 << k)
    return intercept / (1 - slope)


def cycle_record(cycle, b):
    exponents = [step(n, b)[1] for n in cycle]
    total, carry = 0, 0
    for k in exponents:
        carry = 3 * carry + (1 << total)
        total += k
    delta = (1 << total) - 3 ** len(cycle)
    denominator = abs(delta) // gcd(abs(delta), carry)
    content = gcd(abs(cycle[0]), abs(b))
    require(cycle[0] * delta == b * carry, "cycle equation")
    require(all(gcd(abs(n), abs(b)) == content for n in cycle), "cycle gcd invariant")
    common_content = 0
    for n in cycle:
        common_content = gcd(common_content, abs(n))
    require(common_content == content, "cycle content is not gcd(n,b)")
    require(denominator == abs(b) // content, "rational denominator/content duality")
    if abs(b) == 5 and content == 1:
        require((total - 3 * len(cycle)) % 4 == 0, "primitive mod-five clock")
    require(all(n * cycle[0] > 0 for n in cycle), "mixed-sign cycle")
    return {
        "nodes": cycle,
        "minimal_odd_period": len(cycle),
        "halving_exponents": exponents,
        "total_halving_exponent": total,
        "carry_B": carry,
        "delta_2K_minus_3L": delta,
        "content_d": content,
        "primitive_parameter": b // content,
        "primitive_cycle_nodes": [n // content for n in cycle],
        "rational_cycle_denominator_q": denominator,
    }


def direct_forward_control(b, bound=9999, step_cap=2000, height_cap=10**80):
    # This search uses no word formula and no prescribed period bound.
    resolved = {}
    cycles = set()
    counts = Counter()
    censored = []
    for start in range(-bound, bound + 1, 2):
        path, seen, n = [], {}, start
        for _ in range(step_cap):
            if n in resolved:
                cycle = resolved[n]
                break
            if n in seen:
                cycle = canonical(path[seen[n]:])
                cycles.add(cycle)
                break
            if abs(n) > height_cap:
                cycle = None
                break
            seen[n] = len(path)
            path.append(n)
            n = step(n, b)[0]
        else:
            cycle = None
        if cycle is None:
            censored.append(start)
        else:
            for value in path:
                resolved[value] = cycle
            counts[cycle] += 1
    return {
        "odd_start_absolute_bound": bound,
        "number_of_starts": bound + 1,
        "unresolved_segment_step_cap": step_cap,
        "absolute_height_cap": str(height_cap),
        "cycles": sorted(cycles, key=lambda c: (c[0] < 0, len(c), abs(c[0]))),
        "counts": [{"cycle": c, "starts": counts[c]} for c in sorted(counts)],
        "censored_starts": censored,
    }


def main(output):
    max_length, max_total = 10, 22
    parameters = (-5, -1, 1, 5)
    catalog = {b: set() for b in parameters}
    accepted = Counter()
    word_count = 0
    fraction_checks = 0
    for word, total, carry in word_states(max_length, max_total):
        word_count += 1
        delta = (1 << total) - 3 ** len(word)
        denominator = abs(delta) // gcd(abs(delta), carry)
        if denominator not in (1, 5):
            continue
        accepted[denominator] += 1
        require(delta % 2 and delta % 3, "cycle denominator is not coprime to six")
        for b in parameters:
            if b % denominator:
                continue
            numerator = b * carry
            require(numerator % delta == 0, "denominator divisibility criterion")
            n = numerator // delta
            require(n % 2 != 0, "cycle source not odd")
            require(affine_fraction_control(word, b) == n, "independent affine composition")
            fraction_checks += 1
            catalog[b].add(walk_word(n, b, word))
    require(word_count == sum(comb(max_total, length) for length in range(1, max_length + 1)),
            "ordered-word universe incomplete")
    require(word_count == 1744435, "unexpected word universe")

    # Controls for all odd parameters in a small box, including multiples of3.
    gcd_controls = 0
    for b in range(-25, 26, 2):
        for n in range(-999, 1000, 2):
            if 3 * n + b == 0:
                continue
            successor, _ = step(n, b)
            require(gcd(abs(successor), abs(b)) == gcd(abs(3 * n), abs(b)), "general gcd law")
            for d in (3, 5, 7):
                require(step(d * n, d * b)[0] == d * successor, "odd dilation conjugacy")
            gcd_controls += 1

    braid_controls = 0
    for b in parameters:
        for n in range(-99, 100, 2):
            successor, exponent = step(n, b)
            rn = 4 * n + b
            require(step(rn, b) == (successor, exponent + 2), "inverse-fibre braid")
            for t in range(1, 31):
                rt = 4**t * n + b * (4**t - 1) // 3
                require(valuation(rt - n, 3) == valuation(t, 3), "triadic braid period")
                braid_controls += 1
        require({canonical([-n for n in c]) for c in catalog[b]} == catalog[-b], "sign conjugacy")

    require({canonical([5*n for n in c]) for c in catalog[-1]} ==
            {c for c in catalog[-5] if gcd(abs(c[0]), 5) == 5}, "scaled -1 stratum")
    require({canonical([5*n for n in c]) for c in catalog[1]} ==
            {c for c in catalog[5] if gcd(abs(c[0]), 5) == 5}, "scaled +1 stratum")
    require(step(1, -5)[0] == -1 and step(-1, -5)[0] == -1, "unique sign portal")
    require(all(step(n, -5)[0] > 0 for n in range(3, 1000, 2)), "extra positive sign portal")

    word_catalog = {b: set(catalog[b]) for b in parameters}
    forward = {b: direct_forward_control(b) for b in parameters}
    for b in parameters:
        forward_cycles = set(forward[b]["cycles"])
        require(word_catalog[b] <= forward_cycles, "word cycle absent from independent forward control")
        for cycle in forward_cycles - word_catalog[b]:
            record = cycle_record(cycle, b)
            require(len(cycle) > max_length or record["total_halving_exponent"] > max_total,
                    "forward cycle missing INSIDE the exhaustive word universe")
        require(not forward[b]["censored_starts"], "bounded forward census has unresolved starts")
        catalog[b].update(forward_cycles)
    for b in parameters:
        require({canonical([-n for n in c]) for c in catalog[b]} == catalog[-b], "full observed sign conjugacy")
    for b in (-1, 1):
        require({canonical([5*n for n in c]) for c in catalog[b]} ==
                {c for c in catalog[5*b] if gcd(abs(c[0]), 5) == 5}, "full observed dilation stratum")
    result = {
        "status": "PROVED identities supported by FINITE-EXACT controls; no complete cycle classification",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "word_universe": {
            "parameters_b": parameters,
            "minimal_word_length": 1,
            "maximal_word_length": max_length,
            "positive_word_entries": True,
            "maximal_total_halving_exponent": max_total,
            "all_ordered_words_checked": word_count,
            "node_height_bound": None,
            "sign_filter": None,
            "repeated_words": "included, then primitive cycles deduplicated by cyclic rotation",
        },
        "accepted_ordered_words_by_rational_denominator": dict(sorted(accepted.items())),
        "independent_fraction_compositions": fraction_checks,
        "gcd_and_dilation_controls": gcd_controls,
        "triadic_braid_controls": braid_controls,
        "bounded_word_cycles": {b: [cycle_record(c, b) for c in sorted(word_catalog[b], key=lambda c: (c[0] < 0, len(c), abs(c[0])))]
                                for b in parameters},
        "cycles": {b: [cycle_record(c, b) for c in sorted(catalog[b], key=lambda c: (c[0] < 0, len(c), abs(c[0])))]
                   for b in parameters},
        "independent_direct_forward_census": forward,
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(f"PASS: {word_count:,} ordered words; denominator gates {dict(sorted(accepted.items()))}")
    for b in parameters:
        print(f"b={b}: {len(catalog[b])} cycles; " + str([r['nodes'] for r in result['cycles'][b]]))
    print(f"Independent direct iteration: {sum(r['number_of_starts'] for r in forward.values()):,} signed odd starts; no censored starts")
    print(f"Wrote {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).with_suffix(".json"))
    main(parser.parse_args().output)
