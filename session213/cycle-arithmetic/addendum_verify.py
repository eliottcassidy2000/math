#!/usr/bin/env python3
"""Exact checks of all four Session212 branches after deleting vertex a."""
import json
from verify import D, apex, count, require, sole_exit


def is_power_of_two(n):
    return n > 0 and n & (n - 1) == 0


def is_three_times_power(n):
    return n % 3 == 0 and is_power_of_two(n // 3)


def main():
    report = {"status": "PASS", "scope": "FINITE-EXACT implementation controls",
              "rows": []}
    total = 0
    successes = {"minus": set(), "plus": set()}
    for m in range(1, 9):
        tested = 0
        for q in range(1 << m):
            bits = [(q >> i) & 1 for i in range(m)]
            d = D(q)
            require((d == 3) == is_three_times_power(q + 1),
                    ("D3 bitword classification", m, q))
            for tau in (0, 1):
                core = apex([0, tau] + bits)
                k = count(core)
                expected = d + (3 - tau) * q + (1 - tau) * (q & 1) + 2 * tau
                require(k == expected, ("core formula", m, q, tau))
                for plus in (False, True):
                    b = sole_exit(core, 1)
                    if plus:
                        a = len(b) - 1
                        z = len(b) - 2
                        b[z] &= ~(1 << a)
                        b[a] |= 1 << z
                    c = q + tau + int(plus)
                    require(count(b, cycles=True) == c, ("cycle count", m, q, tau, plus))
                    if c:
                        if plus:
                            predicted = (tau == 0 and q == 1) or (tau == 1 and d == 3)
                        else:
                            predicted = tau == 1 and is_power_of_two(c)
                        require((k == 2 * c + 1) == predicted,
                                ("complete equality classification", m, q, tau, plus))
                        if predicted:
                            successes["plus" if plus else "minus"].add(c)
                    tested += 1
        total += tested
        report["rows"].append({"m": m, "compiler_cases": tested})
    report["total_compiler_cases"] = total
    report["successful_cycle_counts"] = {k: sorted(v) for k, v in successes.items()}
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
