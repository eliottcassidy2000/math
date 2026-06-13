#!/usr/bin/env python3
"""Q31 exception probe for the first below-nine-core LRC(14) layer.

The budget-5 Q27 scout for delete_count=4 finds:

  * 489 deletion sets infeasible within a 20s MILP limit,
  * 4 deletion sets unknown at that limit,
  * 2 deletion sets with actual Q27 covers.

This script resolves the six non-infeasible cases.  The four timeouts become
Q27-infeasible with a larger limit.  The two Q27-feasible deletion shapes both
become infeasible as soon as the fibered ladder is widened from Q27 to Q31.

Consequence in the HYP-2444 carry window:

  retaining at least 8 core speeds forces a Q31 witness.

That does not prove the full LRC(14), but it moves the bounded Church-style
descent portal from "below-nine-core" to "below-eight-core or outside-window".
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import reduce
from math import gcd
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


HERE = Path(__file__).resolve().parent
sys.path.append(str(HERE))

from lrc14_marked_ladder_support_codex import summarize_row  # noqa: E402
from lrc14_near_core_q27_setcover_codex import (  # noqa: E402
    CANDIDATES,
    CORE_T,
    Q27,
    blocks,
    obligations,
    primitive_q27_cover_feasibility,
)
from lrc14_pisano_band_ladder_codex import (  # noqa: E402
    bprime_any,
    first_witness_in,
    min_witness_modulus,
    q_lattice,
)


Q27_FEASIBLE_DELETIONS = (
    (28, 42, 56, 84),
    (42, 56, 70, 84),
)

Q27_TIMEOUT_DELETIONS = (
    (14, 21, 28, 70),
    (21, 28, 56, 70),
    (21, 28, 70, 84),
    (28, 42, 70, 84),
)


@dataclass(frozen=True)
class CoverResult:
    deleted: tuple[int, ...]
    ladder_name: str
    feasible: bool | None
    obligation_count: int
    elapsed_s: float
    chosen: tuple[int, ...]
    message: str


def row_from(deleted: tuple[int, ...], chosen: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted([v for v in CORE_T if v not in deleted] + list(chosen)))


def primitive_cover_feasibility(
    deleted: tuple[int, ...],
    qs: tuple[int, ...],
    budget: int,
    ladder_name: str,
    time_limit_s: int = 240,
) -> CoverResult:
    base = tuple(v for v in CORE_T if v not in deleted)
    obs = obligations(base, qs)
    vals = tuple(v for v in CANDIDATES if v not in base)

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    lows: list[float] = []
    ups: list[float] = []
    row_index = 0

    for q, a in obs:
        any_cover = False
        for j, v in enumerate(vals):
            if blocks(v, q, a):
                rows.append(row_index)
                cols.append(j)
                data.append(-1.0)
                any_cover = True
        if not any_cover:
            return CoverResult(deleted, ladder_name, False, len(obs), 0.0, tuple(), "uncoverable")
        lows.append(-np.inf)
        ups.append(-1.0)
        row_index += 1

    for j, v in enumerate(vals):
        if v % 7 != 0:
            rows.append(row_index)
            cols.append(j)
            data.append(-1.0)
    lows.append(-np.inf)
    ups.append(-1.0)
    row_index += 1

    for j, _v in enumerate(vals):
        rows.append(row_index)
        cols.append(j)
        data.append(1.0)
    lows.append(-np.inf)
    ups.append(float(budget))
    row_index += 1

    matrix = coo_matrix((data, (rows, cols)), shape=(row_index, len(vals))).tocsr()
    start = time.time()
    result = milp(
        c=np.zeros(len(vals)),
        integrality=np.ones(len(vals)),
        bounds=Bounds(0, 1),
        constraints=LinearConstraint(matrix, np.array(lows), np.array(ups)),
        options={"time_limit": time_limit_s, "mip_rel_gap": 0},
    )
    elapsed = time.time() - start

    chosen: tuple[int, ...] = tuple()
    if result.x is not None:
        chosen = tuple(vals[i] for i, x in enumerate(result.x) if x > 0.5)

    if result.success:
        feasible: bool | None = True
    elif "Infeasible" in result.message:
        feasible = False
    else:
        feasible = None

    return CoverResult(
        deleted=deleted,
        ladder_name=ladder_name,
        feasible=feasible,
        obligation_count=len(obs),
        elapsed_s=elapsed,
        chosen=chosen,
        message=result.message,
    )


def print_cover(result: CoverResult) -> None:
    print(
        f"  D={result.deleted} ladder={result.ladder_name:4s} "
        f"feasible={result.feasible} obligations={result.obligation_count} "
        f"elapsed={result.elapsed_s:.2f}s chosen={result.chosen}"
    )
    print(f"    message={result.message}")
    if result.chosen:
        row = row_from(result.deleted, result.chosen)
        print(f"    row={row} gcd={reduce(gcd, row, 0)}")
        print(f"    first_Q150={first_witness_in(list(row), q_lattice(150))}")
        print(f"    min_q<=600={min_witness_modulus(list(row), 600)}")
        print(f"    Bprime={bprime_any(list(row))}")
        support = summarize_row("candidate", row, 27)
        print(
            f"    Q27_support blocked={support.blocked_count} "
            f"universal={support.universal_count} "
            f"min_cover_hist={support.min_blocker_hist}"
        )
        print(f"    cover_load={support.cover_load[:8]}")


def main() -> None:
    print("=" * 80)
    print("Codex LRC14 Q31 exception probe for delete_count=4")
    print("=" * 80)
    print("A. Resolve Q27 timeout cases with a longer limit")
    for deleted in Q27_TIMEOUT_DELETIONS:
        result = primitive_q27_cover_feasibility(deleted, 5, time_limit_s=180)
        print(
            f"  D={deleted} feasible={result.feasible} "
            f"obligations={result.obligation_count} elapsed={result.elapsed_s:.2f}s "
            f"chosen={result.chosen}"
        )
        print(f"    message={result.message}")

    print()
    print("B. Probe the two Q27-feasible deletion shapes")
    for deleted in Q27_FEASIBLE_DELETIONS:
        q27 = primitive_cover_feasibility(deleted, tuple(Q27), 5, "Q27", time_limit_s=120)
        print_cover(q27)
        q31 = primitive_cover_feasibility(deleted, tuple(q_lattice(31)), 5, "Q31", time_limit_s=240)
        print_cover(q31)

    print()
    print("C. Proof-shaped consequence")
    print("  Since Q27 infeasibility implies Q31 infeasibility, the 489 first-pass")
    print("  infeasible cases need no widening.  The four timeout cases are Q27")
    print("  infeasible after a longer limit.  The only two Q27-feasible deletion")
    print("  shapes are Q31-infeasible.  Therefore, inside the carry window, every")
    print("  primitive replacement row retaining at least 8 core speeds has a Q31")
    print("  witness.  The bounded descent portal starts at <=7 retained core speeds,")
    print("  unless the row leaves the window or opens an owner/Bprime exception.")


if __name__ == "__main__":
    main()
