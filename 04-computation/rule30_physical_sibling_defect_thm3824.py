#!/usr/bin/env python3
"""Exact physical sibling-defect companion for THM-3824.

For a fixed scale m, scan actual sibling pairs

    A_(m,t)=(U_m(t),U_m(t+2^m)).

Two pairs have the same exact Q_(2,r) Smith class iff their free differences
agree as integers and their first coordinates agree modulo 2^(r-1).  Since
the physical gap is scale-independent, any output Smith split realizes a
nonzero class of the ambient tariff proved in division_tariff.py.

The companion also checks the proved all-scale mechanism: the packed row has
two leading one bits, whose disjoint magnitude intervals force the exact free
defect to be strictly increasing with nonnegative phase at every fixed scale.
The finite census is only a control for that proof and is not a finite-state
or Rule-30-prize result.
"""

from __future__ import annotations

import hashlib
import json
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


EXPECTED_V = (1, 3, 4, 6, 7, 9, 15, 16, 24, 25, 27, 29, 34)
STOP = 4096
MAX_SCALE = 8
MAX_R = 10


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def phi(value: int) -> int:
    return value ^ ((value << 1) | (value << 2))


def rows() -> tuple[int, ...]:
    out = [1]
    for _ in range(STOP):
        out.append(phi(out[-1]))
    return tuple(out)


def unit(path: tuple[int, ...], m: int, t: int) -> int:
    q = 1 << m
    numerator = path[t + q] - path[t]
    divisor = 1 << EXPECTED_V[m]
    require(numerator % divisor == 0, ("unit divisibility", m, t))
    answer = numerator // divisor
    require(answer & 1, ("unit oddness", m, t))
    return answer


def main() -> None:
    path = rows()
    exact_free_collisions = 0
    smith_collisions = 0
    monotonicity_checks = 0
    tariff_realizations = []
    scale_records = []

    for m in range(MAX_SCALE + 1):
        q = 1 << m
        gap = EXPECTED_V[m + 1] - EXPECTED_V[m]
        phase_stop = STOP - 2 * q
        free_seen: dict[int, list[tuple[int, int, int]]] = {}
        records = []
        for t in range(phase_stop + 1):
            left = unit(path, m, t)
            right = unit(path, m, t + q)
            parent = unit(path, m + 1, t)
            require(left + right == (1 << gap) * parent, ("physical lift", m, t))
            delta = right - left
            records.append((t, left, delta, parent))
            free_seen.setdefault(delta, []).append((t, left, parent))

        repeated_free = sum(len(items) - 1 for items in free_seen.values() if len(items) > 1)
        deltas = tuple(item[2] for item in records)
        require(
            all(deltas[j + 1] > deltas[j] > 0 for j in range(len(deltas) - 1)),
            ("strict free-defect monotonicity", m),
        )
        monotonicity_checks += len(deltas) - 1
        exact_free_collisions += repeated_free
        scale_smith = 0
        scale_tariffs = 0
        for r in range(1, MAX_R + 1):
            modulus = 1 << (r - 1)
            seen: dict[tuple[int, int], tuple[int, int]] = {}
            for t, left, delta, parent in records:
                key = (delta, left % modulus)
                prior = seen.get(key)
                if prior is not None:
                    smith_collisions += 1
                    scale_smith += 1
                    prior_t, prior_parent = prior
                    if prior_parent % modulus != parent % modulus:
                        tariff_realizations.append(
                            {
                                "m": m,
                                "r": r,
                                "gap": gap,
                                "t": (prior_t, t),
                                "delta": delta,
                                "parents_mod": (prior_parent % modulus, parent % modulus),
                            }
                        )
                        scale_tariffs += 1
                else:
                    seen[key] = (t, parent)
        scale_records.append(
            {
                "m": m,
                "gap": gap,
                "phases": len(records),
                "distinct_free": len(free_seen),
                "repeated_free": repeated_free,
                "smith_collisions": scale_smith,
                "tariff_realizations": scale_tariffs,
            }
        )

    semantic = hashlib.sha256(
        json.dumps(
            {
                "strict_monotonicity": True,
                "monotonicity_checks": monotonicity_checks,
                "scales": scale_records,
                "tariffs": tariff_realizations[:20],
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode("ascii")
    ).hexdigest()

    print("RULE30_PHYSICAL_SIBLING_DEFECT_THM3824")
    print("status=PROVED_MONOTONICITY+FINITE-EXACT_CONTROL;no_finite_state;all_rule30_prizes_open")
    print(f"universe=0<=m<={MAX_SCALE};0<=t<=4096-2^(m+1);1<=r<={MAX_R}")
    print(f"exact_free_collisions={exact_free_collisions}")
    print(f"smith_collisions={smith_collisions}")
    print(f"tariff_realizations={len(tariff_realizations)}")
    print(f"strict_monotonicity_checks={monotonicity_checks}")
    print("scale_records=" + repr(tuple(scale_records)).replace(" ", ""))
    print("first_tariff_realizations=" + repr(tuple(tariff_realizations[:20])).replace(" ", ""))
    print("mechanism=top_bits_2t_and_2t-1_force_disjoint_delta_magnitude_intervals")
    print(f"semantic_sha256={semantic}")
    print("PASS")


if __name__ == "__main__":
    main()
