"""All-four-branch deletion formula and exact target classification oracle."""
import json
from independent_counts import apex, bits, check, count_paths, defect
from audit_comb import comb


def power_two(n):
    return n > 0 and n & (n - 1) == 0


def main():
    cases, d3 = 0, 0
    successes = [set(), set()]
    for m in range(1, 7):
        for q in range(1 << m):
            d = defect(bits(q, m))
            special = (q + 1) % 3 == 0 and power_two((q + 1) // 3)
            check((d == 3) == special, ("D3", m, q))
            d3 += 1
            for tau in range(2):
                actual = count_paths(apex([0, tau] + bits(q, m)))
                expected = d + (3 - tau) * q + (1 - tau) * (q % 2) + 2 * tau
                check(actual == expected, ("deletion", m, q, tau))
                for forward in range(2):
                    c = q + tau + forward
                    if c:
                        good = ((tau == 0 and q == 1) or (tau == 1 and special)) if forward else tau == 1 and power_two(c)
                        check((actual == 2 * c + 1) == good, ("equality", m, q, tau, forward))
                        if good:
                            successes[forward].add(c)
                    cases += 1
    print(json.dumps({"status": "PASS", "branch_cases": cases, "D3_cases": d3,
                      "successful_minus_cycle_counts": sorted(successes[0]),
                      "successful_plus_cycle_counts": sorted(successes[1])}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
