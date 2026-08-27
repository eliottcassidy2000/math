#!/usr/bin/env python3
"""Finite exact (OS+) boundary audit for the two-reversal sphere.

Left factors: every classified strong fixed-gauge two-reversal presentation
of orders 3..12.  Right factors: every labelled no-sink tournament of orders
3..5.  This is deliberately a presentation census, not a class census.
"""

from __future__ import annotations

from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
DP_PATH = ROOT / "04-computation/tournament_order11_tail_fiber_thm4224.py"
OS_PATH = ROOT / "04-computation/tournament_ordinal_cocycle_parity_thm4184.py"
DP_HASH = "3c9ad5de462dcc71d4312bba12b42086d4e1e43520b5209b4a807f762d059c18"
OS_HASH = "9ab09bf8b70ee5dcf3d86698180a456d67f012655b49a16dfea9903caefbb39c"


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load(name: str, path: Path, expected: str):
    need(sha256(path.read_bytes()).hexdigest() == expected, f"{name} hash")
    spec = spec_from_file_location(name, path)
    need(spec is not None and spec.loader is not None, f"{name} import")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


dp = load("thm4224_dp_os", DP_PATH, DP_HASH)
os = load("thm4184_os", OS_PATH, OS_HASH)


def all_labelled(order: int):
    width = order * (order - 1) // 2
    return tuple(os.parse(f"{word:0{width}b}", order) for word in range(1 << width))


right_by_order = {
    order: tuple(out for out in all_labelled(order) if not os.has_sink(out))
    for order in range(3, 6)
}
need(tuple(map(len, right_by_order.values())) == (2, 32, 704), "right counts")
rights = tuple(
    (order, os.label(out), os.tournament_data(out))
    for order in range(3, 6)
    for out in right_by_order[order]
)
need(len(rights) == 738, "right universe")

digest = sha256()
checks = 0
lefts = 0
print("two_reversal_osplus_audit")
print(f"dependency_sha256_dp={DP_HASH}")
print(f"dependency_sha256_os={OS_HASH}")
print("right_labelled_no_sink=q3:2,q4:32,q5:704,total:738")

for n in range(3, 13):
    minimum = None
    rows = []
    presentations = dp.classified_two_reversal_shell(n)
    need(len(presentations) == (n - 1) ** 2 - 3, f"left count n={n}")
    for reversals in sorted(presentations):
        tournament = dp.from_reversals(n, reversals)
        need(dp.is_strong(tournament), f"left strong n={n}")
        left = os.tournament_data(tournament.out)
        for q, word, right in rights:
            value = os.remainder(left, right)
            need(value > 0, f"OS+ failure n={n},rev={reversals},q={q},word={word}")
            checks += 1
            digest.update(f"{n}:{reversals}:{q}:{word}:{value}\n".encode())
            row = (reversals, left.hamilton, q, word, right.hamilton, value)
            if minimum is None or value < minimum:
                minimum = value
                rows = [row]
            elif value == minimum:
                rows.append(row)
    lefts += len(presentations)
    print(f"n={n} left_presentations={len(presentations)} min_Rplus={minimum} min_mult={len(rows)} min_rows={rows}")

need(lefts == 475 and checks == 475 * 738 == 350550, "OS+ universe")

# Scope control: dropping no-sink is genuinely hostile already at the C3 left.
c3 = os.tournament_data(dp.from_reversals(3, ((0, 1), (1, 2))).out)
sink_controls = []
for q in range(1, 6):
    right = os.tournament_data(os.transitive(q))
    sink_controls.append((q, os.remainder(c3, right)))
need(all(value < 0 for _, value in sink_controls), "sink hostile signs")

print(f"totals_left_presentations={lefts} total_remainders={checks}")
print(f"sink_hostile_C3_right_transitive={sink_controls}")
print(f"audit_digest={digest.hexdigest()}")
print("all_checks=PASS_FINITE_EXACT_ONLY")
