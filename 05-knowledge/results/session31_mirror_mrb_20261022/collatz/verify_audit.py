#!/usr/bin/env python3
"""Exact, dependency-free hostile probes for the supplied Collatz manuscript.

No assertion is used as a validation gate; python -O gives identical output.
All tests concern finite trajectories or the displayed elementary identities.
They do not prove universal Collatz convergence or the native 153N gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from math import gcd
from pathlib import Path


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def v2(n: int) -> int:
    require(type(n) is int and n > 0, "valuation requires a positive integer")
    return (n & -n).bit_length() - 1


def collatz_step(n: int) -> int:
    return n // 2 if n % 2 == 0 else 3 * n + 1


def orbit(n: int, cap: int = 10000) -> list[int]:
    out = [n]
    while n != 1:
        require(len(out) <= cap, f"finite orbit cap reached at start {out[0]}")
        n = collatz_step(n)
        out.append(n)
    return out


def literal_run(x: int) -> dict:
    """Read a maximal increasing odd run through its first >=2-halving event."""
    require(x > 0 and x % 2 == 1, "run starts at a positive odd integer")
    trace = [x]
    valuations = []
    while True:
        peak = 3 * x + 1
        trace.append(peak)
        y = peak
        k = 0
        while y % 2 == 0:
            y //= 2
            k += 1
            trace.append(y)
        valuations.append(k)
        if k >= 2:
            return {"start": trace[0], "end": y, "peak": peak,
                    "run_length": len(valuations), "peak_v2": k,
                    "valuation_word": valuations, "trace": trace}
        x = y


def prescribe_three_parameters(r: int, s: int, t: int) -> dict:
    """r odd steps, terminal peak valuation s+1, next odd run length t."""
    modulus = 3 ** r
    b = ((2 ** s - 1) * pow(2 ** (s + t), -1, modulus)) % modulus
    if b == 0:
        b += modulus
    if b % 2 == 0:
        b += modulus
    a = (2 ** (s + t) * b - 2 ** s + 1) // modulus
    x = a * 2 ** r - 1
    y = 2 ** t * b - 1
    require(a > 0 and a % 2 == 1 and b > 0 and b % 2 == 1,
            "three-parameter construction positivity or parity failed")
    require(a * 3 ** r == 2 ** (s + t) * b - 2 ** s + 1,
            "three-parameter congruence failed")
    literal = literal_run(x)
    require(literal["run_length"] == r, "first run length mismatch")
    require(literal["peak_v2"] == s + 1, "collapse length mismatch")
    require(literal["end"] == y and v2(y + 1) == t,
            "post-collapse trailing-one length mismatch")
    require(literal_run(y)["run_length"] == t, "recharged literal run mismatch")
    return {"r": r, "s": s, "t": t, "a": a, "b": b,
            "x": x, "peak": literal["peak"], "y": y}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("receipt.json"))
    parser.add_argument("--input", type=Path,
                        default=Path(__file__).resolve().parent.parent / "session31-input" /
                        "Collatz_Conjecture_Complete_Proof_English_Version.pdf")
    args = parser.parse_args()
    checks = 0

    # Positive control: unique shifted-valuation grid, including the n=1 branch.
    for n in range(2, 10001, 2):
        e = v2(n + 2)
        a = (n + 2) // 2 ** e
        require(a % 2 == 1 and a * 2 ** e - 2 == n, "grid control")
        checks += 1

    # First invalid transformation in the manuscript.
    require(10 // 2 == 5 and 3 * 2 ** (2 - 1) - 2 == 4, "halving witness")
    require(v2(12 + 2) == 1 and v2(6 + 2) == 3, "base closure witness")
    checks += 2

    # Independent literal first-run verification of the useful peak formula.
    grid_rows = []
    for a in range(1, 102, 2):
        for n in range(2, 33):
            N = a * 2 ** n - 2
            literal = literal_run(N // 2)
            peak = 2 * (a * 3 ** (n - 1) - 1)
            require(literal["peak"] == peak and literal["run_length"] == n - 1,
                    "first peak formula mismatch")
            require(all(collatz_step(u) == v for u, v in
                        zip(literal["trace"], literal["trace"][1:])), "literal edge mismatch")
            grid_rows.append([a, n, peak, literal["end"], literal["peak_v2"]])
            checks += 1

    # Exact envelope inequality and its omitted boundaries.
    envelope_failures = []
    envelope_equalities = []
    for a in range(1, 102, 2):
        for n in range(2, 33):
            p = 2 * (a * 3 ** (n - 1) - 1)
            prev = 2 * (a * 3 ** (n - 2) - 1)
            require(4 * prev - p == 2 * (a * 3 ** (n - 2) - 3), "envelope identity")
            if p > 4 * prev:
                envelope_failures.append([a, n, a * 2 ** n - 2])
            if p == 4 * prev:
                envelope_equalities.append([a, n, a * 2 ** n - 2])
            checks += 1
    require(envelope_failures == [[1, 2, 2]], "envelope failure boundary")
    require(envelope_equalities == [[1, 3, 6], [3, 2, 10]], "envelope equality boundary")

    # Long ascent need not buy a long cascade or any net contraction.
    hostile_family = []
    for n in range(3, 80, 4):
        N = 2 ** n - 2
        literal = literal_run(N // 2)
        require(literal["peak_v2"] == 4, "constant-four-halving family")
        require(literal["end"] == (3 ** (n - 1) - 1) // 8, "family endpoint")
        if n >= 11:
            require(literal["end"] > N, "family expansion")
        hostile_family.append({"n": n, "N": N, "peak": literal["peak"],
                               "end": literal["end"], "peak_v2": 4})
        checks += 1
    n14 = literal_run(7)
    require(n14["end"] == 13 and 4 * 13 > 3 * 14, "three-quarter contraction hostile")
    n2046 = literal_run(1023)
    require(n2046["peak"] == 118096 and n2046["end"] == 7381, "8m cascade hostile")
    checks += 2

    # Arbitrary-length cascade can expose an all-ones integer immediately.
    recharge_rows = []
    for r in range(3, 42, 2):
        for k in range(2, 42, 2):
            y = 2 ** r - 1
            x = (2 ** k * y - 1) // 3
            require(x % 2 == 1 and 3 * x + 1 == 2 ** k * y, "recharge inverse")
            literal = literal_run(x)
            require(literal["run_length"] == 1 and literal["peak_v2"] == k
                    and literal["end"] == y and v2(y + 1) == r, "literal recharge")
            recharge_rows.append([r, k, x, y])
            checks += 1

    # Stronger simultaneous parameter independence, with no density inference.
    parameter_rows = []
    for r in range(1, 13):
        for s in range(1, 13):
            for t in range(1, 13):
                parameter_rows.append(prescribe_three_parameters(r, s, t))
                checks += 1
    parameter_large = [prescribe_three_parameters(20, 20, 20)]
    checks += 1

    # Arbitrarily many whole run/cascade cycles can grow above the baseline.
    rebound_rows = []
    for m in range(1, 41):
        start = 2 ** (3 * m + 1) - 5
        x = start
        endpoints = [x]
        for j in range(m):
            literal = literal_run(x)
            y = literal["end"]
            expected = 9 ** (j + 1) * 2 ** (3 * (m - j - 1) + 1) - 5
            require(literal["run_length"] == 2 and literal["peak_v2"] == 2,
                    "repeated growth-cycle valuation word")
            require(y == expected and y > x, "repeated growth-cycle formula")
            require(min(literal["trace"]) == x, "earlier baseline crossing")
            x = y
            endpoints.append(x)
        rebound_rows.append({"m": m, "start": start, "end": x,
                             "endpoints": endpoints if m <= 7 else None})
        checks += 1

    # Guards here are read from x0 itself, not an unproved native-address map.
    guard_table = [
        {"m_mod4": 0, "x0_mod15": 12, "x0_mod35": 32, "x0_mod105": 102,
         "passes_both_unit_guards": False},
        {"m_mod4": 1, "x0_mod15": 11, "x0_mod35": 11, "x0_mod105": 11,
         "passes_both_unit_guards": True},
        {"m_mod4": 2, "x0_mod15": 3, "x0_mod35": 18, "x0_mod105": 18,
         "passes_both_unit_guards": False},
        {"m_mod4": 3, "x0_mod15": 14, "x0_mod35": 4, "x0_mod105": 74,
         "passes_both_unit_guards": True},
    ]
    require(pow(2, 12, 105) == 1, "guard table periodicity")
    checks += 1
    for m in range(1, 513):
        x0 = 2 ** (3 * m + 1) - 5
        row = guard_table[m % 4]
        require((x0 % 15, x0 % 35, x0 % 105) ==
                (row["x0_mod15"], row["x0_mod35"], row["x0_mod105"]),
                "direct starting-value guard residues")
        require((gcd(x0, 15) == 1 and gcd(x0, 35) == 1) == (m % 2 == 1),
                "direct starting-value guard parity")
        require((gcd(x0, 105) == 1) == (gcd(x0, 15 * 35) == 1) == (m % 2 == 1),
                "lcm/product unit equivalence")
        require(x0 % 15 % 5 == x0 % 35 % 5, "joint guard compatibility")
        checks += 1

    # A non-power-of-two orbit can avoid both proposed gates 3 and 5.
    gate_hostile = orbit(21)
    require(gate_hostile == [21, 64, 32, 16, 8, 4, 2, 1], "gate hostile")
    require(3 not in gate_hostile and 5 not in gate_hostile, "gate omission")
    checks += 1

    # Reproduce the manuscript's finite showcase, not a universal theorem.
    showcase = orbit(14693280766)
    peak = max(showcase)
    peak_index = showcase.index(peak)
    require(peak == 564859072960 and v2(peak) == 6, "showcase peak")
    require(len(showcase) - 1 == 445, "showcase stopping time")
    halving_trace = showcase[peak_index:peak_index + 7]
    require(halving_trace[-1] == 8825923015, "showcase full cascade")
    require(all(z % 8 == 0 for z in halving_trace[:4]), "four overlapping labels")
    checks += 1

    # Exhaustive finite universe 1<=n<=10000, including the below-10000 convention.
    hub1822 = 0
    hit9232 = 0
    maximum9232 = 0
    orbit_digest = hashlib.sha256()
    for start in range(1, 10001):
        seq = orbit(start)
        hub1822 += int(1822 in seq)
        hit9232 += int(9232 in seq)
        maximum9232 += int(max(seq) == 9232)
        orbit_digest.update((str(start) + ":" + ",".join(map(str, seq)) + "\n").encode())
        checks += 1
    require((hub1822, hit9232, maximum9232) == (3997, 4011, 1579), "finite census")
    require(1822 not in orbit(10000) and 9232 not in orbit(10000), "strict bound convention")

    def digest(rows: object) -> str:
        return hashlib.sha256(json.dumps(rows, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

    result = {
        "status": "FINITE-EXACT checks; accompanying proofs are scoped; no global convergence inference",
        "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "checks": checks,
        "grid_formula_cases": len(grid_rows), "grid_formula_sha256": digest(grid_rows),
        "envelope_failures": envelope_failures, "envelope_equalities": envelope_equalities,
        "constant_cascade_family": hostile_family,
        "N14_first_run": n14, "N2046_first_run": n2046,
        "recharge_cases": len(recharge_rows), "recharge_sha256": digest(recharge_rows),
        "three_parameter_cases": len(parameter_rows), "three_parameter_sha256": digest(parameter_rows),
        "three_parameter_large": parameter_large,
        "rebound_cases": rebound_rows, "gate_hostile_21": gate_hostile,
        "guarded_growth_hostile": {
            "guard_coordinate": "ordinary Collatz starting integer x0, not a mapped native source address",
            "checked_m": [1, 512], "periodicity": "2^12 = 1 mod 105",
            "residue_table": guard_table,
            "unbounded_odd_m": "both starting-value unit guards pass; first m whole cycles grow"
        },
        "showcase": {"start": showcase[0], "steps_to_1": len(showcase) - 1,
                     "maximum": peak, "maximum_step": peak_index, "peak_v2": v2(peak),
                     "halving_trace": halving_trace, "orbit_sha256": digest(showcase)},
        "finite_census": {"universe": "integers 1 through 10000 inclusive",
                          "hit_1822": hub1822, "hit_9232": hit9232,
                          "maximum_9232": maximum9232,
                          "same_counts_for_1_through_9999": True,
                          "all_orbits_sha256": orbit_digest.hexdigest()}
    }
    args.out.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"checks": checks, "output": str(args.out),
                      "input_sha256": result["input_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
