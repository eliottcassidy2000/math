#!/usr/bin/env python3
"""Compile the inherited 65-line support into a positive Mahler equation.
This checks the new arithmetic identities, not tournament membership.
"""
import json

C = tuple(range(1, 130, 2))
threshold = {c: 2 if c == 3 else 3 if c == 5 else 1 for c in C}
seeds = tuple(c << threshold[c] for c in C)
if len(set(seeds)) != 65:
    raise RuntimeError("Source must have exactly 65 distinct seeds")


def member(n):
    if n == 1:
        return True
    if n < 3 or not n & 1 or n in (7, 21):
        return False
    m = n - 1
    c = m // (m & -m)
    return c <= 129


def ray_count(x):
    ans = int(x >= 1) + int(x >= 11)
    for c in C:
        q = (x - 1) // c
        top = q.bit_length() - 1 if q > 0 else -1
        ans += max(0, top - threshold[c] + 1)
    return ans


def compile_prefix(limit):
    b = [0] * (limit + 1)
    source = set(seeds)
    for m in range(1, limit + 1):
        b[m] = int(m in source) + (b[m // 2] if m % 2 == 0 else 0)
    return [int(n == 1) + int(n == 11) + (b[n - 1] if n else 0)
            for n in range(limit + 1)]


def full_source(m):
    # R_*(z) in equation (6).
    return int(m >= 2 and m % 4 == 2) - int(m == 6) - int(m == 10) + int(m == 12) + int(m == 40)


def main():
    limit = 100_000
    coeff = compile_prefix(limit)
    count = 0
    full = [0] * (limit + 1)
    for n in range(limit + 1):
        if coeff[n] not in (0, 1) or bool(coeff[n]) != member(n):
            raise RuntimeError(f"Mahler/direct mismatch at {n}")
        count += coeff[n]
        if ray_count(n) != count:
            raise RuntimeError(f"count mismatch at {n}")
    for m in range(1, limit + 1):
        r = full_source(m)
        if r < 0:
            raise RuntimeError(f"negative full source at {m}")
        full[m] = r + (full[m // 2] if m % 2 == 0 else 0)
        expected = int(m % 2 == 0 and m not in (6, 10, 20))
        if full[m] != expected:
            raise RuntimeError(f"full-target identity mismatch at {m}")
    output = {
        "status": "exact arithmetic controls passed; tournament coverage inherited",
        "checked_through": limit,
        "source_exponents": sorted(seeds),
        "first_values": [n for n in range(1, 70) if coeff[n]],
        "counts_at_powers_of_two": {str(k): ray_count(1 << k) for k in (8, 10, 16, 20, 40, 100)},
        "exceptional_point": 11,
        "forbidden_values": [7, 21],
        "general_chan_pak_status": "OPEN",
    }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
