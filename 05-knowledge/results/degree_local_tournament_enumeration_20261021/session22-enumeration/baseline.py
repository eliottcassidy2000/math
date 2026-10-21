#!/usr/bin/env python3
"""Packed future-response counting with the minimal common odd denominator.

No core or tournament census is used. The imported previous DP is an
independent implementation control, not a production dependency.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb, factorial, gcd
from pathlib import Path
import argparse
import json
import sys
import time
import tracemalloc

# The declared order-200 exact output exceeds Python's decimal display limit.
sys.set_int_max_str_digits(0)

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "tournament_future_boundary_20261020"
CHECKS = Counter()


def need(ok, label):
    CHECKS[label] += 1
    if not ok:
        raise ValueError(label)


def natural(n):
    if type(n) is not int or n < 0:
        raise ValueError("bound must be a literal nonnegative integer")
    return n


def odd_part(n):
    if type(n) is not int or n < 1:
        raise ValueError("positive literal integer required")
    return n // (n & -n)


def profile_pack(values, bound):
    """Injectively encode a declared finite response list; length is sidecar."""
    natural(bound)
    if type(values) is not tuple or any(type(v) is not int or not 0 <= v <= bound
                                        for v in values):
        raise ValueError("tuple of literal response digits required")
    width = max(1, bound.bit_length())
    return sum(v << (width * i) for i, v in enumerate(values))


def profile_unpack(code, count, bound):
    natural(bound)
    natural(count)
    natural(code)
    width = max(1, bound.bit_length())
    if code >> (width * count):
        raise ValueError("undeclared high response digits")
    mask = (1 << width) - 1
    values = tuple((code >> (width * i)) & mask for i in range(count))
    if any(v > bound for v in values):
        raise ValueError("response digit exceeds declared bound")
    return values


def packed_dp(bound, *, core=False, no_fixed=False, arithmetic="odd"):
    """Return exact coefficients, deterministic work statistics.

    arithmetic='fraction' is a packed rational control; 'odd' is production.
    no_fixed is an observer specialization and can return noninteger values.
    """
    natural(bound)
    if type(core) is not bool or type(no_fixed) is not bool:
        raise ValueError("literal Boolean mode flags required")
    if type(arithmetic) is not str or arithmetic not in ("odd", "fraction"):
        raise ValueError("unknown arithmetic")
    start = 3 if no_fixed else 1
    lengths = list(range(start, bound + 1, 2))[::-1]
    width = max(1, bound.bit_length())
    digit_mask = (1 << width) - 1
    coordinate_count = max(0, (bound - start) // 2 + 1)
    masks = [(1 << (width * c)) - 1 for c in range(coordinate_count + 1)]
    facts = [factorial(k) for k in range(bound + 1)]
    oddfacts = [odd_part(f) for f in facts]
    two_orders = [(f & -f).bit_length() - 1 for f in facts]
    denominator = oddfacts[bound]
    integer = arithmetic == "odd"
    states = {(0, 0): denominator if integer else Fraction(1)}
    rows = []
    division_checks = 0
    for ell in lengths:
        idx = (ell - start) // 2
        row = sum(gcd(ell, j) << (width * i)
                  for i, j in enumerate(range(start, ell, 2)))
        limits = bound // ell
        denoms = [ell ** k * (oddfacts[k] if integer else facts[k])
                  for k in range(limits + 1)]
        bases = [k * (ell - 1) // 2 + ell * k * (k - 1) // 2
                 - (two_orders[k] if integer else 0) for k in range(limits + 1)]
        new = defaultdict(int if integer else Fraction)
        transitions = 0
        response_coordinate_updates = 0
        for (degree, response), weight in states.items():
            charge = (response >> (width * idx)) & digit_mask
            for k in range((bound - degree) // ell + 1):
                exponent = bases[k] + k * charge
                if exponent < 0:
                    raise ValueError("factorial two-cancellation exceeded exponent")
                if integer:
                    factor, remainder = divmod(weight << exponent, denoms[k])
                    if remainder:
                        raise ValueError("nonintegral history transition")
                    division_checks += 1
                else:
                    factor = weight * Fraction(1 << exponent, denoms[k])
                updated = response + k * row
                excesses = range((bound - degree) // ell - k + 1) if core and k else (0,)
                for excess in excesses:
                    newdegree = degree + ell * (k + excess)
                    cap = min(ell - 2, bound - newdegree)
                    count = max(0, (cap - start) // 2 + 1)
                    packed = updated & masks[count]
                    signweight = (-1) ** excess * comb(k + excess - 1, excess) if k else 1
                    new[newdegree, packed] += factor * signweight
                    transitions += 1
                    response_coordinate_updates += count
        states = {key: value for key, value in new.items() if value}
        tuple_cells = sum(max(0, (min(ell - 2, bound - d) - start) // 2 + 1)
                          for d, _ in states)
        rows.append({"cycle_length": ell, "states": len(states),
                     "transitions": transitions,
                     "tuple_coordinate_updates_avoided": response_coordinate_updates,
                     "tuple_response_cells": tuple_cells,
                     "packed_response_bits": sum(r.bit_length() for _, r in states)})
    answer = [Fraction(0)] * (bound + 1)
    for (degree, response), weight in states.items():
        if response:
            raise ValueError("nonempty terminal response")
        answer[degree] += Fraction(weight, denominator) if integer else weight
    return answer, {"bound": bound, "core": core, "no_fixed": no_fixed,
                    "width": width, "odd_denominator_bits": denominator.bit_length(),
                    "precomputed_gcd_calls": coordinate_count * (coordinate_count - 1) // 2,
                    "peak_states": max([1] + [row["states"] for row in rows]),
                    "transitions": sum(row["transitions"] for row in rows),
                    "integer_division_checks": division_checks,
                    "tuple_coordinate_updates_avoided": sum(
                        row["tuple_coordinate_updates_avoided"] for row in rows),
                    "peak_tuple_response_cells": max([0] + [
                        row["tuple_response_cells"] for row in rows]),
                    "peak_packed_response_bits": max([0] + [
                        row["packed_response_bits"] for row in rows]),
                    "rows": rows}


def inherited():
    spec = spec_from_file_location("inherited_future", OLD / "verify.py")
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def encode(sequence):
    return [str(x) for x in sequence]


def rejects(fn, label):
    try:
        fn()
    except (ValueError, TypeError):
        need(True, label)
    else:
        need(False, label)


def controls():
    old = inherited()
    for bound in range(25):
        for core, no_fixed in ((False, False), (True, False), (True, True), (False, True)):
            fast, stats = packed_dp(bound, core=core, no_fixed=no_fixed)
            rational, rational_stats = packed_dp(bound, core=core, no_fixed=no_fixed,
                                                   arithmetic="fraction")
            previous, previous_stats = old.profile_dp(bound, core=core, no_fixed=no_fixed)
            need(fast == previous == rational, "all modes versus inherited rational DP")
            need([(r["states"], r["transitions"]) for r in stats["rows"]] ==
                 [(r["states"], r["transitions"]) for r in previous_stats["rows"]],
                 "identical mathematical states and transitions")
            if bound <= 18:
                need(fast == old.direct_series(bound, core=core, no_fixed=no_fixed),
                     "independent unmerged partition substitution")
            if bound == 24:
                ascending, _ = old.profile_dp(bound, core=core, no_fixed=no_fixed, ascending=True)
                need(fast == ascending, "opposite elimination order")
            need(stats["peak_states"] == rational_stats["peak_states"], "packing arithmetic separation")
    for mass in range(31):
        common = odd_part(factorial(mass))
        for parts in old.odd_partitions(mass):
            weight = old.cycle_weight(parts)
            need(weight.denominator % 2 == 1, "every raw history has odd denominator")
            need(common % weight.denominator == 0, "common denominator covers raw histories")
        fixed = old.cycle_weight((1,) * mass)
        need(fixed.denominator == common, "minimality from fixed-point history")
    for bound in range(1, 41):
        for ell in range(1, bound + 1, 2):
            for k in range(bound // ell + 1):
                two = (factorial(k) & -factorial(k)).bit_length() - 1
                need(ell * k * (k - 1) // 2 >= two, "local factorial two cancellation")
        for count in range(8):
            values = tuple((bound * i + i * i) % (bound + 1) for i in range(count))
            need(profile_unpack(profile_pack(values, bound), count, bound) == values,
                 "exact packed codec and declared length")
    # Smallest bad width: at bound 2, one bit aliases response (2,0) with (0,1).
    need(2 + (0 << 1) == 0 + (1 << 1), "too-small width aliases mathematical states")
    need(profile_pack((2, 0), 2) != profile_pack((0, 1), 2), "correct width repairs alias")
    need(profile_pack((), 5) == profile_pack((0,), 5), "length is retained stage-degree sidecar")
    # The observer is intentionally not rounded into a class count.
    observer, _ = packed_dp(9, core=True, no_fixed=True)
    need(observer[3] == Fraction(2, 3), "fixed-point-free mass is not an integer count")
    # A substitution excess is degree, not an additional cycle response.
    need(profile_pack((1,), 5) != profile_pack((2,), 5), "actual cycles versus substitution excess")
    need(-Fraction(2, 3) * 2 != -Fraction(2, 3) * 4,
         "three-cycle excess then singleton has different interaction if excess counted")
    need(odd_part(factorial(3)) * 2 == 6 and odd_part(factorial(3)) != 1,
         "scaled integer state is not an unlabelled class count")
    core6, _ = packed_dp(6, core=True)
    for n in range(7):
        direct = sum((old.fixed_core_count(parts) for parts in old.odd_partitions(n)), Fraction(0))
        need(core6[n] == direct, "literal invariant orientations and homogeneous-pair test through six")
    for bad in (False, 2.0, -1, "2"):
        rejects(lambda bad=bad: packed_dp(bad), "literal bound guard")
    for name in ("core", "no_fixed"):
        for bad in (0, 1.0, None):
            rejects(lambda name=name, bad=bad: packed_dp(5, **{name: bad}), "literal flag guard")
    for bad in (True, 1.0, "integer"):
        rejects(lambda bad=bad: packed_dp(5, arithmetic=bad), "arithmetic mode guard")
    for bad in ((False,), (1.0,), (-1,), (4,), [1]):
        rejects(lambda bad=bad: profile_pack(bad, 3), "profile source type and digit guard")
    rejects(lambda: profile_unpack(8, 1, 3), "retained length rejects high digits")
    rejects(lambda: profile_unpack(3, 1, 2), "unused binary digit rejects excess")
    return old


def run():
    CHECKS.clear()
    old = controls()
    series = {}
    for name, bound, core, no_fixed in (
        ("tournaments", 200, False, False),
        ("pair_free_cores", 100, True, False),
        ("no_fixed_core_mass", 100, True, True)):
        values, stats = packed_dp(bound, core=core, no_fixed=no_fixed)
        CHECKS["checked exact integer transition divisions"] += stats["integer_division_checks"]
        series[name] = {"coefficients": encode(values), "statistics": stats}
        need(all(v >= 0 for v in values), "final observer coefficients nonnegative")
        if not no_fixed:
            need(all(v.denominator == 1 for v in values), "final class coefficients integral")
        if name == "tournaments":
            for line in (OLD / "a000568-reference.txt").read_text().splitlines():
                if line and not line.startswith("#"):
                    n, value = map(int, line.split())
                    need(values[n] == value, "independent stored A000568 table through 76")
            need(values[:61] == old.direct_series(60), "unmerged odd partitions through 60")
        else:
            previous, previous_stats = old.profile_dp(bound, core=core, no_fixed=no_fixed)
            need(values == previous, "large core and observer inherited control")
            need(stats["peak_states"] == previous_stats["peak_states"] and
                 stats["transitions"] == previous_stats["transitions"], "large inherited work identity")
    previous_receipt = json.loads((OLD / "receipt.json").read_text())
    need(series["tournaments"]["coefficients"] == previous_receipt["large"]["coefficients"],
         "all 201 coefficients agree with frozen inherited 200-order run")
    need(series["pair_free_cores"]["coefficients"] == previous_receipt["core_classes"],
         "all 101 pair-free coefficients agree with frozen receipt")
    need(series["no_fixed_core_mass"]["coefficients"] == previous_receipt["fixed_free_core_mass"],
         "all 101 rational observer coefficients agree with frozen receipt")
    return {"schema": 1, "algorithm": "packed future responses, minimal common odd denominator",
            "checks": sum(CHECKS.values()), "checks_by_kind": dict(sorted(CHECKS.items())),
            "dependencies": {"inherited_script_lf_sha256": digest(OLD / "verify.py"),
                             "inherited_receipt_lf_sha256": digest(OLD / "receipt.json"),
                             "stored_oeis_lf_sha256": digest(OLD / "a000568-reference.txt")},
            "series": series}


def digest(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def benchmark(bound, core=False, memory=False):
    old = inherited()
    functions = (
        ("inherited_tuple_fraction", lambda: old.profile_dp(bound, core=core)),
        ("packed_fraction", lambda: packed_dp(bound, core=core, arithmetic="fraction")),
        ("packed_odd_integer", lambda: packed_dp(bound, core=core)),
    )
    rows = []
    reference = None
    for name, fn in functions:
        if memory:
            tracemalloc.start()
        start = time.perf_counter()
        result, stats = fn()
        elapsed = time.perf_counter() - start
        peak = tracemalloc.get_traced_memory()[1] if memory else None
        if memory:
            tracemalloc.stop()
        if reference is None:
            reference = result
        elif result != reference:
            raise ValueError("benchmark changed the mathematical answer")
        rows.append({"algorithm": name, "seconds": elapsed, "traced_peak_bytes": peak,
                     "peak_states": stats["peak_states"], "transitions": stats["transitions"]})
    return {"bound": bound, "core": core, "memory_tracing": memory,
            "python": sys.version, "runs": rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--benchmark", type=int)
    parser.add_argument("--core", action="store_true")
    parser.add_argument("--memory", action="store_true")
    args = parser.parse_args()
    if args.benchmark is not None:
        print(json.dumps(benchmark(args.benchmark, args.core, args.memory), indent=2))
        return
    result = run()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    path = HERE / "receipt.json"
    if args.verify:
        if path.read_text() != text:
            raise ValueError("saved receipt differs")
    else:
        path.write_text(text, encoding="utf-8", newline="\n")
    print(json.dumps({"checks": result["checks"], "receipt_lf_sha256": digest(path),
                      "verified": args.verify}, sort_keys=True))


if __name__ == "__main__":
    main()
