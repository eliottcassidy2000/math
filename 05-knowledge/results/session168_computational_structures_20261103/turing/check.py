#!/usr/bin/env python3
"""Exact finite controls for the S168 profile-invisible machine compiler."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
SEEDS = (8, 10)
PAIRS = list(combinations(range(5), 2))


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def rows_from_mask(n, mask):
    rows = [0] * n
    for bit, (u, v) in enumerate(combinations(range(n), 2)):
        if (mask >> bit) & 1:
            rows[u] |= 1 << v
        else:
            rows[v] |= 1 << u
    return rows


def profile_dp(rows):
    """All-subset endpoint DP, then least-vertex-anchored path-cover DP."""
    n = len(rows)
    total = 1 << n
    ends = [[0] * n for _ in range(total)]
    h = [0] * total
    h[0] = 1
    for v in range(n):
        ends[1 << v][v] = 1
    for mask in range(1, total):
        for v in range(n):
            if not ((mask >> v) & 1):
                continue
            prev = mask ^ (1 << v)
            if prev:
                ends[mask][v] = sum(ends[prev][u] for u in range(n)
                                    if ((prev >> u) & 1) and ((rows[u] >> v) & 1))
        h[mask] = sum(ends[mask])
    covers = [[0] * (n + 1) for _ in range(total)]
    covers[0][0] = 1
    for mask in range(1, total):
        anchor = mask & -mask
        rest = mask ^ anchor
        sub = rest
        while True:
            block = sub | anchor
            complement = mask ^ block
            for j in range(1, n + 1):
                covers[mask][j] += h[block] * covers[complement][j - 1]
            if sub == 0:
                break
            sub = (sub - 1) & rest
    return covers[-1]


def profile_permutations(rows):
    n = len(rows)
    counts = [0] * (n + 1)
    for word in permutations(range(n)):
        for cuts in range(1 << (n - 1)):
            if all(((cuts >> j) & 1) or ((rows[word[j]] >> word[j + 1]) & 1)
                   for j in range(n - 1)):
                counts[1 + cuts.bit_count()] += 1
    for j in range(n + 1):
        require(counts[j] % factorial(j) == 0, "nonintegral unordered cover count")
    return [counts[j] // factorial(j) for j in range(n + 1)]


def falling(r, j):
    value = 1
    for i in range(j):
        value *= r - i
    return value


def slots(profile, r):
    return sum(count * falling(r, j) for j, count in enumerate(profile))


def encode(bits):
    rows = []
    n = 5 * len(bits)
    for j, bit in enumerate(bits):
        later = ((1 << n) - 1) ^ ((1 << (5 * (j + 1))) - 1)
        rows.extend((row << (5 * j)) | later for row in rows_from_mask(5, SEEDS[bit]))
    return rows


def sccs(rows):
    n = len(rows)
    reach = [row | (1 << j) for j, row in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if (reach[i] >> k) & 1:
                reach[i] |= reach[k]
    unseen = set(range(n))
    components = []
    while unseen:
        i = min(unseen)
        part = sorted(j for j in unseen if ((reach[i] >> j) & 1) and ((reach[j] >> i) & 1))
        components.append(part)
        unseen.difference_update(part)
    return components


def arc_changes(left, right):
    require(len(left) == len(right), "order changed during a flip")
    return [(u, v) for u, v in combinations(range(len(left)), 2)
            if ((left[u] >> v) & 1) != ((right[u] >> v) & 1)]


def verify_flip(before, after, coordinate):
    require(sum(x != y for x, y in zip(before, after)) == 1, "not one bit change")
    actual = arc_changes(encode(before), encode(after))
    expected = [(5 * coordinate, 5 * coordinate + 2)]
    require(actual == expected, f"wrong literal flip: {actual}, expected {expected}")
    return list(expected[0])


def digest(bits):
    return hashlib.sha256(json.dumps(encode(bits), separators=(",", ":")).encode()).hexdigest()


def valuation(n, p):
    value = 0
    while n % p == 0:
        n //= p
        value += 1
    return value


def exponent_code(exponents, capacities):
    return [bit for e, cap in zip(exponents, capacities) for bit in ([1] * e + [0] * (cap - e))]


def decode_exponents(bits, capacities):
    pos = 0
    values = []
    for cap in capacities:
        segment = bits[pos:pos + cap]
        e = sum(segment)
        require(segment == [1] * e + [0] * (cap - e), "non-unary register")
        values.append(e)
        pos += cap
    return values


def fractran_step(exponents, capacities, primes, program):
    selected = None
    for i, frac in enumerate(program):
        if all(e >= valuation(frac.denominator, p) for e, p in zip(exponents, primes)):
            selected = i
            break
    if selected is None:
        return {"status": "halt", "source": exponents}
    frac = program[selected]
    delta = [valuation(frac.numerator, p) - valuation(frac.denominator, p) for p in primes]
    target = [e + d for e, d in zip(exponents, delta)]
    if any(e > cap for e, cap in zip(target, capacities)):
        return {"status": "overflow", "source": exponents, "rule": selected, "target": target}
    bits = exponent_code(exponents, capacities)
    micro = []
    cursor = 0
    for register, (e, d, cap) in enumerate(zip(exponents, delta, capacities)):
        positions = range(e, e + d) if d > 0 else range(e - 1, e + d - 1, -1)
        for local in positions:
            next_bits = bits.copy()
            coordinate = cursor + local
            next_bits[coordinate] ^= 1
            arc = verify_flip(bits, next_bits, coordinate)
            decoded = decode_exponents(next_bits, capacities)
            micro.append({"register": register, "coordinate": coordinate, "arc": arc,
                          "exponents": decoded, "graph_sha256": digest(next_bits)})
            bits = next_bits
        cursor += cap
    require(bits == exponent_code(target, capacities), "wrong arithmetic endpoint")
    require(len(micro) == sum(map(abs, delta)), "wrong FRACTRAN flip cost")
    return {"status": "step", "source": exponents, "rule": selected,
            "target": target, "flips": micro}


STATES = ("A", "B", "H")
RULES = {("A", 0): ("B", 1, 1), ("A", 1): ("B", 1, -1),
         ("B", 0): ("A", 1, -1), ("B", 1): ("H", 1, 1)}


def tm_code(tape, head, state):
    bits = []
    for i, symbol in enumerate(tape):
        bits += [symbol] + [int(i == head and state == q) for q in STATES]
    return bits


def tm_step(tape, head, state):
    if (state, tape[head]) not in RULES:
        return {"status": "halt"}
    target_state, write, direction = RULES[state, tape[head]]
    target_head = head + direction
    if not (0 <= target_head < len(tape)):
        return {"status": "overflow"}
    before = tm_code(tape, head, state)
    bits = before.copy()
    operations = []
    if tape[head] != write:
        operations.append(("write", 4 * head))
    operations += [("add_head", 4 * target_head + 1 + STATES.index(target_state)),
                   ("remove_head", 4 * head + 1 + STATES.index(state))]
    micro = []
    for phase, coordinate in operations:
        next_bits = bits.copy()
        next_bits[coordinate] ^= 1
        arc = verify_flip(bits, next_bits, coordinate)
        micro.append({"phase": phase, "coordinate": coordinate, "arc": arc,
                      "graph_sha256": digest(next_bits)})
        bits = next_bits
    target_tape = tape.copy()
    target_tape[head] = write
    require(bits == tm_code(target_tape, target_head, target_state), "wrong TM endpoint")
    require(len(micro) == 2 + int(tape[head] != write), "wrong TM flip count")
    require(len(micro) == sum(x != y for x, y in zip(before, bits)), "nongeodesic TM path")
    return {"status": "step", "tape": target_tape, "head": target_head,
            "state": target_state, "flips": micro}


def main():
    cert = {"scope": "Finite checks, plus structural proof in NOTE.md. Not a universal-program execution."}
    seed_profile = [0, 9, 33, 31, 10, 1]
    seed_checks = []
    invariants = []
    for mask in SEEDS:
        rows = rows_from_mask(5, mask)
        require(profile_dp(rows) == seed_profile, "seed subset-DP profile")
        require(profile_permutations(rows) == seed_profile, "seed permutation profile")
        require(sccs(rows) == [list(range(5))], "seed not strong")
        cycle = [0, 4, 3, 2, 1, 0]
        require(all((rows[u] >> v) & 1 for u, v in zip(cycle, cycle[1:])), "missing seed cycle")
        degrees = [row.bit_count() for row in rows]
        central = degrees.index(2)
        invariant = sorted(degrees[v] for v in range(5) if (rows[central] >> v) & 1)
        invariants.append(invariant)
        seed_checks.append({"mask": mask, "profile": seed_profile,
                            "unique_degree2_outneighbor_degrees": invariant})
    require(invariants[0] != invariants[1], "nonisomorphism witness failed")
    require(arc_changes(rows_from_mask(5, 8), rows_from_mask(5, 10)) == [(0, 2)], "seed flip")
    require(all(slots(seed_profile, r) == r**5 + 6*r**3 + 2*r for r in range(6)), "slot polynomial")
    cert["seed"] = seed_checks
    product_checks = []
    for bits in product((0, 1), repeat=2):
        rows = encode(bits)
        profile = profile_dp(rows)
        require(all(slots(profile, r) == (r**5 + 6*r**3 + 2*r)**2 for r in range(11)),
                "two-slot full-profile product")
        require(sccs(rows) == [list(range(5)), list(range(5, 10))], "SCC product")
        product_checks.append({"bits": list(bits), "profile": profile})
    cert["two_slot_products"] = product_checks
    for left in product((0, 1), repeat=4):
        for right in product((0, 1), repeat=4):
            require(len(arc_changes(encode(left), encode(right))) == sum(x != y for x, y in zip(left, right)),
                    "labelled Hamming isometry")
    cert["four_bit_labelled_distance_pairs"] = 256
    primes, capacities = [2, 3, 5], [3, 3, 3]
    program = [Fraction(3, 2), Fraction(5, 3)]
    statuses = Counter()
    for exponents in product(range(4), repeat=3):
        exponents = list(exponents)
        value = 1
        for p, e in zip(primes, exponents):
            value *= p**e
        choices = [i for i, f in enumerate(program) if value % f.denominator == 0]
        result = fractran_step(exponents, capacities, primes, program)
        statuses[result["status"]] += 1
        require((not choices) == (result["status"] == "halt"), "integer versus vector halt")
        if choices:
            require(result["rule"] == choices[0], "lost priority")
            target = value * program[choices[0]]
            require([valuation(target.numerator, p) for p in primes] == result["target"], "wrong integer result")
    cert["fractran_exhaustive_capacity3"] = dict(statuses)
    trace = []
    exponents = [2, 1, 0]
    while True:
        result = fractran_step(exponents, capacities, primes, program)
        trace.append(result)
        if result["status"] != "step":
            break
        exponents = result["target"]
    require([item["rule"] for item in trace[:-1]] == [0, 0, 1, 1, 1], "demo rule trace")
    cert["fractran_demo"] = {"program": ["3/2", "5/3"], "primes": primes,
                              "capacities": capacities, "integers": [12, 18, 27, 45, 75, 125],
                              "tournament_order": 45, "common_H": 9**9, "trace": trace}
    tm_statuses = Counter()
    for word in product((0, 1), repeat=4):
        for head in range(4):
            for state in STATES:
                result = tm_step(list(word), head, state)
                tm_statuses[result["status"]] += 1
    cert["tm_exhaustive_four_cells"] = dict(tm_statuses)
    tape, head, state = [0] * 4, 2, "A"
    tm_trace = [{"tape": tape.copy(), "head": head, "state": state,
                 "graph_sha256": digest(tm_code(tape, head, state))}]
    all_flips = 0
    for step in range(6):
        result = tm_step(tape, head, state)
        require(result["status"] == "step", "demo prematurely stopped")
        all_flips += len(result["flips"])
        tape, head, state = result["tape"], result["head"], result["state"]
        tm_trace.append(result)
    require(state == "H" and tape == [1, 1, 1, 1] and head == 2, "TM demo endpoint")
    require(tm_step(tape, head, state)["status"] == "halt", "TM demo halt")
    require(all_flips == 16, "TM demo flip total")
    cert["tm_demo"] = {"rules": [{"state": q, "read": a, "new_state": qp, "write": b, "move": d}
                                  for (q, a), (qp, b, d) in RULES.items()],
                       "tape_coordinates": [-2, -1, 0, 1], "tournament_order": 80,
                       "macro_steps": 6, "arc_flips": all_flips, "common_H": 9**16,
                       "trace": tm_trace}
    # The actual 80-vertex tournament has sixteen strong five-vertex SCCs.
    require(sccs(encode(tm_code(tape, head, state))) == [list(range(5*i, 5*i+5)) for i in range(16)],
            "80-vertex decoding SCCs")
    path = ROOT / "certificate.json"
    path.write_text(json.dumps(cert, indent=2) + "\n")
    print(json.dumps({"seed_profile": seed_profile, "two_slot_products": 4,
                      "fractran_states": dict(statuses), "tm_states": dict(tm_statuses),
                      "tm_steps": 6, "tm_flips": all_flips, "certificate": path.name}))


if __name__ == "__main__":
    main()
