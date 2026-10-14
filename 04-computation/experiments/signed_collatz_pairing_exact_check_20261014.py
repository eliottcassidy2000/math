#!/usr/bin/env python3
"""Exact arithmetic audit. Standard library only; checks survive python -O.

No finite orbit census is used as evidence of exhaustion or convergence.
The infinite density statements are proved in the companion note; this script
checks their finite cylinders and exact rational generating sums.
"""
from collections import Counter
from fractions import Fraction
import hashlib
import json

checks = 0


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise ArithmeticError(message)


def valuation(n):
    if n == 0:
        raise ValueError("valuation of zero is outside this audit")
    n = abs(n)
    return (n & -n).bit_length() - 1


def shortcut_odd(n, sign):
    if n % 2 != 1 or sign not in (-1, 1):
        raise ValueError("odd source and sign +/-1 required")
    return (3 * n + sign) // 2


def accelerated(n, sign):
    if n % 2 != 1 or sign not in (-1, 1):
        raise ValueError("odd source and sign +/-1 required")
    return (3 * n + sign) // (1 << valuation(3 * n + sign))


def shortcut(n, sign):
    return shortcut_odd(n, sign) if n % 2 else n // 2


def full(n, sign):
    return 3 * n + sign if n % 2 else n // 2


table = []
for r in range(1, 12, 2):
    row = {
        "source_mod12": r,
        "Splus_mod18": shortcut_odd(r, 1) % 18,
        "Sminus_mod18": shortcut_odd(r, -1) % 18,
        "vplus": ">=2" if r % 4 == 1 else "1",
        "vminus": "1" if r % 4 == 1 else ">=2",
    }
    table.append(row)
    for t in range(-32, 33):
        n = r + 12 * t
        require(shortcut_odd(n, 1) % 18 == row["Splus_mod18"], "mod12 plus")
        require(shortcut_odd(n, -1) % 18 == row["Sminus_mod18"], "mod12 minus")

pair_counts = Counter()
N = 1 << 18
for n in range(1, N, 2):
    a, b = valuation(3 * n + 1), valuation(3 * n - 1)
    u, v = accelerated(n, 1), accelerated(n, -1)
    require(shortcut_odd(n, 1) - shortcut_odd(n, -1) == 1, "adjacent pair")
    require(min(a, b) == 1 and max(a, b) >= 2, "valuation coupling")
    require((1 << a) * u - (1 << b) * v == 2, "lifted pair equation")
    require(u % 3 == (-1) ** a % 3, "plus character")
    require(v % 3 == -((-1) ** b) % 3, "minus character")
    require((u % 3, v % 3) != (1, 2), "forbidden residue pair")
    if n > 1:
        require(min(u, v) < n < max(u, v), "one descent and one ascent")
        selected = u if n % 4 == 1 else v
        require(4 * selected <= 3 * n + 1, "greedy contraction")
    else:
        require(u == v == 1, "root exception")
    if a == 1:
        require(u == (1 << (b - 1)) * v + 1, "plus odd branch")
    else:
        require(v == (1 << (a - 1)) * u - 1, "minus odd branch")
    pair_counts[(u % 3, v % 3)] += 1

for n in range(-2048, 2049):
    for sign in (-1, 1):
        require(full(-n, sign) == -full(n, -sign), "full signed conjugacy")
        require(shortcut(-n, sign) == -shortcut(n, -sign), "shortcut conjugacy")
        if n % 2:
            require(accelerated(-n, sign) == -accelerated(n, -sign), "odd conjugacy")

expected_transitions = {
    1: {0: [0, 2], 1: [2], 2: [1, 2]},
    -1: {0: [0, 1], 1: [1, 2], 2: [1]},
}
for sign in (-1, 1):
    actual = {r: sorted({shortcut(n, sign) % 3 for n in range(6) if n % 3 == r})
              for r in range(3)}
    require(actual == expected_transitions[sign], "shortcut transition relation")

cycle_rows = []
for sign, cycles in ((1, [[1]]), (-1, [[1], [5, 7], [17, 25, 37, 55, 41, 61, 91]])):
    for cycle in cycles:
        vals = [valuation(3 * n + sign) for n in cycle]
        for i, n in enumerate(cycle):
            require(accelerated(n, sign) == cycle[(i + 1) % len(cycle)], "known cycle")
        B, A, m = 0, 0, len(cycle)
        for a in vals:
            B = 3 * B + (1 << A)
            A += a
        require(cycle[0] * ((1 << A) - 3 ** m) == sign * B, "cycle affine equation")
        require(sign * ((1 << A) - 3 ** m) > 0, "cycle sign constraint")
        cycle_rows.append({"sign": sign, "odd_cycle": cycle, "valuations": vals,
                           "A": A, "m": m, "B": B, "denominator": (1 << A) - 3 ** m})

# Exact joint valuation cylinder census: 2^K odd residues modulo 2^(K+1).
K = 14
joint = Counter()
tail = 0
for n in range(1, 1 << (K + 1), 2):
    a, b = valuation(3 * n + 1), valuation(3 * n - 1)
    if max(a, b) <= K:
        joint[(a, b)] += 1
    else:
        tail += 1
for k in range(2, K + 1):
    require(joint[(1, k)] == 1 << (K - k), "minus high exact fibre count")
    require(joint[(k, 1)] == 1 << (K - k), "plus high exact fibre count")
require(tail == 2, "joint omitted tail size")

# Geometric sums prove the three exact limiting pair masses, not just approximations.
even_high = Fraction(1, 4) / (1 - Fraction(1, 4))
odd_high = Fraction(1, 8) / (1 - Fraction(1, 4))
pair_density = {(1, 1): even_high, (2, 2): even_high, (2, 1): 2 * odd_high}
require(set(pair_density.values()) == {Fraction(1, 3)}, "equal-third joint law")

# For each fixed valuation, enumerate its complete CRT fibre through modulus 3^j.
lift_rows = []
for j in range(1, 6):
    q = 3 ** j
    fibres = 0
    for sign in (-1, 1):
        for k in range(1, 9):
            d = 1 << (k + 1)
            r = ((1 << k) - sign) * pow(3, -1, d) % d
            residues = []
            for t in range(3 ** (j - 1)):
                n = r + t * d
                require(valuation(3 * n + sign) == k, "CRT valuation fixed")
                residues.append(accelerated(n, sign) % q)
            target_mod3 = sign * ((-1) ** k) % 3
            require(sorted(residues) == [u for u in range(q) if u % 3 == target_mod3],
                    "uniform lift within mod3 fibre")
            fibres += 1
    lift_rows.append({"modulus": q, "valuation_sign_fibres_checked": fibres,
                      "plus_each_residue_1mod3": str(Fraction(1, q)),
                      "plus_each_residue_2mod3": str(Fraction(2, q))})

# Hostile controls: neither mod3 nor any fixed modest modulus decides U mod3.
modular_witnesses = []
for M in range(1, 65):
    for sign in (-1, 1):
        seen = {}
        witness = None
        for n in range(1, 8 * M + 2, 2):
            r, out = n % M, accelerated(n, sign) % 3
            if r in seen and seen[r][1] != out:
                witness = [seen[r][0], n]
                break
            seen.setdefault(r, (n, out))
        require(witness is not None, "finite-modulus hostile witness")
        if M in (3, 4, 12, 64):
            modular_witnesses.append({"modulus": M, "sign": sign, "sources": witness,
                                       "output_mod3": [accelerated(n, sign) % 3 for n in witness]})

require([n * n % 3 for n in range(3)] == [0, 1, 1], "squaring mod3")
require([n * n % 4 for n in range(4)] == [0, 1, 0, 1], "squaring mod4")

# Optional-sign dynamics: direct checks accompany the elementary descent proof.
max_steps, max_source = 0, 1
for source in range(1, 10000, 2):
    n, steps = source, 0
    while n > 1:
        sign = 1 if n % 4 == 1 else -1
        child = accelerated(n, sign)
        require(1 <= child < n, "greedy strict descent")
        n, steps = child, steps + 1
    if steps > max_steps:
        max_steps, max_source = steps, source

result = {
    "status": "FINITE-EXACT; infinite statements proved in accompanying note",
    "checks": checks,
    "odd_source_universe": {"first": 1, "last": N - 1, "size": N // 2},
    "mod12_table": table,
    "shortcut_mod3_relations": expected_transitions,
    "known_accelerated_cycles_verified_not_exhausted": cycle_rows,
    "joint_pair_counts_finite": {f"{a},{b}": count for (a, b), count in sorted(pair_counts.items())},
    "joint_pair_densities_proved": {f"{a},{b}": str(p) for (a, b), p in sorted(pair_density.items())},
    "valuation_cylinder_audit": {"K": K, "odd_residues": 1 << K, "unresolved_tail_residues": tail},
    "prime_power_lift_checks": lift_rows,
    "modular_hostile_witnesses": modular_witnesses,
    "greedy_sign_check": {"odd_sources_below": 10000, "maximum_steps": max_steps,
                           "first_source_attaining_maximum": max_source},
    "limitations": ["No complete cycle classification.", "No proof for fixed-sign convergence.",
                    "Starting-value density is not orbit frequency.", "Repository prior-art reconciliation is incomplete."],
}
payload = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
result["payload_sha256"] = hashlib.sha256(payload).hexdigest()
print(json.dumps(result, sort_keys=True, indent=2))
