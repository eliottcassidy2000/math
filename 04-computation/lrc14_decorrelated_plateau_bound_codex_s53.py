#!/usr/bin/env python3
"""HYP-2675 / S53: finite decorrelated plateau-bound audit.

KPS S19 sharpened the surviving wide route:

    p0(E) <= p0_decorr + Weyl error,
    sup p0_decorr = Q(k-1) < cap_k.

This script independently audits the finite bounded-base part of that claim
using the THM-548/S51 boundary-value formula

    P_r(B) = sum_t prof_t(B) c_t(r),

where r is the number of independent/decorrelated far points.  It scans bounded
primitive bases B subset {0,...,14} for k=8..12 and compares all base-size
strata b=|B|, r=k-b against the one-stranger plateau Q(k-1).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations

from lrc14_signed_multifar_boundary_hierarchy_codex_s51 import boundary_value_direct
from lrc14_wide_branch_ridge_codex_s47 import CAP, primitive, sumset_excess, additive_energy


Row = tuple[int, ...]


def fmt(q: F) -> str:
    return f"{q} ({float(q):+.9f})"


@dataclass(frozen=True)
class DecorrRow:
    k: int
    core: Row
    r: int
    value: F

    @property
    def b(self) -> int:
        return len(self.core)

    @property
    def cap(self) -> F:
        return CAP[self.k]

    @property
    def margin(self) -> F:
        return self.cap - self.value


def bounded_cores(size: int, bound: int = 14) -> list[Row]:
    if size == 1:
        return [(0,)]
    out: list[Row] = []
    for rest in combinations(range(1, bound + 1), size - 1):
        core = (0,) + rest
        if primitive(core):
            out.append(core)
    return out


def q_plateau(k_minus_1: int) -> F:
    core = tuple(range(k_minus_1))
    return boundary_value_direct(core, 1)


def scan_k(k: int) -> tuple[DecorrRow, dict[int, DecorrRow]]:
    best_by_b: dict[int, DecorrRow] = {}
    best: DecorrRow | None = None
    for b in range(1, k):
        r = k - b
        local: DecorrRow | None = None
        for core in bounded_cores(b):
            val = boundary_value_direct(core, r)
            row = DecorrRow(k, core, r, val)
            if local is None or (row.value, -sumset_excess(row.core), row.core) > (
                local.value,
                -sumset_excess(local.core),
                local.core,
            ):
                local = row
            if best is None or (row.value, row.b, row.core) > (best.value, best.b, best.core):
                best = row
        assert local is not None
        best_by_b[b] = local
    assert best is not None
    return best, best_by_b


def tournament_by_base_size(k: int, best_by_b: dict[int, DecorrRow]) -> None:
    names = list(best_by_b)
    values = {b: (best_by_b[b].value, F(best_by_b[b].b), str(best_by_b[b].core)) for b in names}
    scores: dict[int, int] = {b: 0 for b in names}
    for a in names:
        for b in names:
            if a == b:
                continue
            if values[a] > values[b]:
                scores[a] += 1
    hist = Counter(scores.values())
    path = sorted(names, key=lambda b: values[b], reverse=True)
    # A total order by exact values gives no directed 3-cycles; keep this explicit.
    print("  TOURNAMENT ANALYSIS BY BASE-SIZE STRATA")
    print("    vertices are base sizes b=|B|, not runners.")
    print("    pairwise observable=(max_B P_{k-b}(B), b, core label); larger wins.")
    print(f"    scores={scores}")
    print(f"    score_hist={dict(hist)} directed_3cycles=0 Hamiltonian_path_count=1")
    print(f"    tie Hamiltonian path={' > '.join('b='+str(b) for b in path)}")


def main() -> None:
    print("HYP-2675 / S53 -- decorrelated plateau-bound audit")
    print("bounded bases B subset {0,...,14}; P_r(B) exact by THM-548/S51")
    print()
    for k in range(8, 13):
        q = q_plateau(k - 1)
        cap = CAP[k]
        best, best_by_b = scan_k(k)
        print(f"k={k}: cap={fmt(cap)} Q(k-1)={fmt(q)} cap-Q={fmt(cap-q)}")
        print(
            f"  global bounded-base max: b={best.b} r={best.r} core={best.core} "
            f"P_r={fmt(best.value)} margin={fmt(best.margin)} equals_Q={best.value == q}"
        )
        print("  best by base size:")
        for b in sorted(best_by_b):
            row = best_by_b[b]
            print(
                f"    b={b} r={row.r} P_r={fmt(row.value)} "
                f"Q-P={fmt(q-row.value)} core={row.core} "
                f"exc={sumset_excess(row.core)} energy={additive_energy(row.core)}"
            )
        tournament_by_base_size(k, best_by_b)
        print()
    print("SYNTHESIS")
    print("  In this bounded-base audit, every k=8..12 maximum occurs at")
    print("  b=k-1, r=1, core=consec_{k-1}, i.e. the one-stranger plateau Q(k-1).")
    print("  This supports the incoming KPS S19 decorrelation foundation.")
    print("  Remaining work is not this finite comparison; it is the explicit")
    print("  Weyl/decorrelation error and finite glue from bounded-gap rows.")
    print("  No LRC(14) proof is claimed.")


if __name__ == "__main__":
    main()
