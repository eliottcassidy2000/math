#!/usr/bin/env python3
r"""Exact certificate for the first pure `(3,4)` flood-tail branch.

Fix

    H = {8,9,10,11,12,13,14},  E = H union {3,4},  a = 15.

This script proves that every family

    E union {15,b,c,d},              15 < b < c < d,

is lonely.  It uses the exact THM-741 interval kernel, but implements the
branch traversal and fixed-E2 screen directly rather than calling the flood
driver.  The infinite tail is covered at three levels:

* `b >= V2` is the THM-735 three-needle common-threshold leg;
* for each `b < V2`, the monotone P2/fixed-E2 inequality closes all large c;
* at each remaining c, the exact `m3` fixed-E2 inequality closes all large d,
  and every finite d below its exact rational cap is swept.

The exact swept values use nested sparse subtraction.  Deterministic first,
middle, and last terminal rows in every active b-branch are independently
rebuilt by a full thirteen-comb union through `good_norm`.

Carrier audit.  The faithful vertices are b-indexed proof obligations carrying
their nested survivor bank and certificate type.  A runner, isolated comb,
Fano flag, residue, or wall event loses shared-survivor incidence.  The
periodic danger combs are one-dimensional Kakeya needles only in this precise
adaptive-component sense; no planar Kakeya dimension estimate is invoked.
"""

from __future__ import annotations

import hashlib
import importlib.util
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
H = tuple(range(8, 15))
E = tuple(sorted((*H, 3, 4)))
A = 15


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(sha256(CORE_PATH) == CORE_SHA256, "THM-741 kernel hash changed")
    spec = importlib.util.spec_from_file_location("thm741_j4_34_a15_dependency", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 kernel")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def p2_precloses(core, m: F, r: int, c: int) -> tuple[bool, F]:
    """P2 on c, then the fixed-E2 discrepancy bound on every d>c."""
    denominator_lower = F(5, 7) * m - F(8 * r, 49 * c)
    slack = denominator_lower - core.S2 * r / (7 * (c + 1))
    return denominator_lower > 0 and slack > 0, slack


def first_p2_preclosed(core, m: F, r: int, start: int, stop: int) -> int:
    """First c in [start,stop) whose entire d-tail is P2-preclosed."""
    if start >= stop or not p2_precloses(core, m, r, stop - 1)[0]:
        return stop
    low, high = start, stop - 1
    while low < high:
        middle = (low + high) // 2
        if p2_precloses(core, m, r, middle)[0]:
            high = middle
        else:
            low = middle + 1
    require(p2_precloses(core, m, r, low)[0], "bad P2 cutoff")
    require(low == start or not p2_precloses(core, m, r, low - 1)[0],
            "nonminimal P2 cutoff")
    return low


def tournament_fingerprint(rows: list[dict[str, object]]) -> dict[str, object]:
    """Transitive scheduler on active b-indexed proof obligations.

    Pair observable: larger exact minimum swept clearance, then fewer terminal
    sweeps, then smaller b.  The terminal proof does not use this quotient.
    """
    active = [row for row in rows if int(row["sweeps"]) > 0]
    keys = {
        int(row["b"]): (
            F(row["minimum"]),
            -int(row["sweeps"]),
            -int(row["b"]),
        )
        for row in active
    }
    ordered = tuple(sorted(keys, key=lambda b: keys[b], reverse=True))
    scores = {b: len(ordered) - index - 1 for index, b in enumerate(ordered)}
    score_histogram = Counter(scores.values())
    triangles = 0
    for left, middle, right in combinations(ordered, 3):
        # All edges point forward in `ordered`, hence no directed cycle.
        triangles += int(not (keys[left] > keys[middle] > keys[right]))
    path_digest = hashlib.sha256(
        ",".join(map(str, ordered)).encode("ascii")
    ).hexdigest()
    return {
        "vertices": len(ordered),
        "score_histogram": score_histogram,
        "triangles": triangles,
        "sccs": len(ordered),
        "hamiltonian_paths": 1 if ordered else 0,
        "path": ordered,
        "path_digest": path_digest,
    }


def main() -> None:
    core = load_core()
    require(core.S2 == F(99, 70) and core.S2 * core.S2 > 2,
            "invalid rational sqrt(2) majorant")

    good0, r0, m0 = core.good_norm(E)
    require((r0, m0) == (28, F(433607, 2522520)), "root geometry changed")
    r1, m1, good1 = core.subtract(good0, A)
    sparse_m1 = core.subtract_sparse(good0, A)
    require(m1 == sparse_m1 and m1 > 0, "a=15 subtraction mismatch")

    V2 = core.minV(3, *(F(4) * m1 / (core.S2 * r1)).as_integer_ratio())
    terminal_digest = hashlib.sha256()
    certificate_digest = hashlib.sha256()
    terminal_rows: list[tuple[tuple[int, ...], F]] = []
    b_rows: list[dict[str, object]] = []

    total_c = 0
    p2_preclosed = 0
    exact_m3 = 0
    exact_m3_closed = 0
    fallback = 0
    candidate_nodes = 0
    candidate_sweeps = 0
    minimum: F | None = None
    minimum_family: tuple[int, ...] | None = None

    for b in range(A + 1, V2):
        r2, m2, good2 = core.subtract(good1, b)
        require(m2 == core.subtract_sparse(good1, b) and m2 > 0,
                f"b={b} subtraction mismatch")
        V3 = core.minV(2, *(F(5) * m2 / (core.S2 * r2)).as_integer_ratio())
        c_start = b + 1
        p2_start = first_p2_preclosed(core, m2, r2, c_start, V3)
        row_c = max(0, V3 - c_start)
        row_preclosed = max(0, V3 - p2_start)
        total_c += row_c
        p2_preclosed += row_preclosed

        row_exact = row_closed = row_candidates = row_sweeps = 0
        row_minimum: F | None = None
        for c in range(c_start, p2_start):
            row_exact += 1
            exact_m3 += 1
            m3_sparse = core.subtract_sparse(good2, c)
            require(m3_sparse > 0, f"empty c node b={b},c={c}")
            denominator = m3_sparse - m2 / 7
            if denominator > 0:
                d_cap = core.S2 * r2 / (7 * denominator)
            else:
                r3, m3, _ = core.subtract(good2, c)
                require(m3 == m3_sparse, f"fallback c mismatch b={b},c={c}")
                d_cap = core.S2 * r3 / (6 * m3)
                fallback += 1
            d_max = floor_fraction(d_cap)
            if d_max <= c:
                row_closed += 1
                exact_m3_closed += 1
                continue

            row_candidates += 1
            candidate_nodes += 1
            r3, m3, good3 = core.subtract(good2, c)
            require(m3 == m3_sparse and r3 == len(good3),
                    f"literal c carrier mismatch b={b},c={c}")
            for d in range(c + 1, d_max + 1):
                clearance = core.subtract_sparse(good3, d)
                require(clearance > 0, f"nonpositive terminal {(b,c,d)}")
                family = tuple(sorted(E + (A, b, c, d)))
                row_sweeps += 1
                candidate_sweeps += 1
                terminal_rows.append((family, clearance))
                terminal_digest.update(
                    (f"{b},{c},{d}:{clearance.numerator}/{clearance.denominator}\n").encode(
                        "ascii"
                    )
                )
                if row_minimum is None or clearance < row_minimum:
                    row_minimum = clearance
                if minimum is None or clearance < minimum:
                    minimum = clearance
                    minimum_family = family

        row = {
            "b": b,
            "V3": V3,
            "p2_start": p2_start,
            "c_nodes": row_c,
            "preclosed": row_preclosed,
            "exact": row_exact,
            "closed": row_closed,
            "candidate_nodes": row_candidates,
            "sweeps": row_sweeps,
            "minimum": None if row_minimum is None else str(row_minimum),
        }
        b_rows.append(row)
        certificate_digest.update(
            (
                f"{b},{V3},{p2_start},{row_c},{row_preclosed},{row_exact},"
                f"{row_closed},{row_candidates},{row_sweeps},{row['minimum']}\n"
            ).encode("ascii")
        )

    require(V2 == 189 and len(b_rows) == 173, "b frontier changed")
    require(total_c == 11177 and p2_preclosed == 2849,
            "c/P2 ledger changed")
    require(exact_m3 == 8328 and exact_m3_closed == 7166 and fallback == 0,
            "exact c ledger changed")
    require(candidate_nodes == 1162 and candidate_sweeps == 17198,
            "terminal ledger changed")
    require(minimum == F(32953849, 624660036), "swept minimum changed")
    require(
        minimum_family == (3, 4, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 23),
        "minimum family changed",
    )
    require(
        terminal_digest.hexdigest()
        == "3594fb4a07d9ee79780f7c99cf4cf2427b0a921282f8c2c19249c46c339602b2",
        "terminal ledger digest changed",
    )
    require(
        certificate_digest.hexdigest()
        == "2ce412ad92743b74400e3e4cba57dfb1c1c75a3dc4b919a998d762df050f180a",
        "certificate ledger digest changed",
    )

    # A second path rebuilds first/middle/last terminal rows in every active
    # b branch as a full thirteen-comb union from scratch.
    by_b: dict[int, list[tuple[tuple[int, ...], F]]] = {}
    for family, clearance in terminal_rows:
        b = family[-3] if family[-4] == A else next(
            speed for speed in family if speed > A
        )
        by_b.setdefault(b, []).append((family, clearance))
    sample_manifest = hashlib.sha256()
    sample_count = 0
    for b in sorted(by_b):
        rows = by_b[b]
        for index in sorted({0, len(rows) // 2, len(rows) - 1}):
            family, clearance = rows[index]
            _, _, full_measure = core.good_norm(family)
            require(full_measure == clearance,
                    f"full-union cross-check failed b={b},index={index}")
            sample_count += 1
            sample_manifest.update(
                (
                    f"{','.join(map(str, family))}:{clearance.numerator}/"
                    f"{clearance.denominator}\n"
                ).encode("ascii")
            )

    tournament = tournament_fingerprint(b_rows)
    require(sample_count == 177, "full-union sample census changed")
    require(
        sample_manifest.hexdigest()
        == "623f110f3f6a0f26b275cfe2097e80d01cb11705008228e70ec2734d9c02f4cc",
        "full-union sample manifest changed",
    )
    require(tournament["triangles"] == 0, "scheduler tournament cycle")
    require(
        tournament["score_histogram"]
        == Counter({score: 1 for score in range(int(tournament["vertices"]))}),
        "scheduler tournament score histogram changed",
    )
    require(
        tournament["vertices"] == 61
        and tournament["path_digest"]
        == "ed5c0ab87fa176b82830ac3ed893c5c5529cbe35921b49f19ef078eb6597e621",
        "scheduler tournament path changed",
    )

    print("THM-741 PURE (3,4) FLOOD TAIL: EXACT a=15 BRANCH CERTIFICATE")
    print("=" * 92)
    print(f"dependency_sha256={CORE_SHA256}")
    print(f"body E={E}; fixed a={A}; claim all integers {A}<b<c<d")
    print(f"root r={r0} m={m0}; after-a r={r1} m={m1}; common b cutoff V2={V2}")
    print(
        f"b obligations={len(b_rows)}; c nodes={total_c}; "
        f"P2-preclosed c={p2_preclosed}; exact-m3={exact_m3}"
    )
    print(
        f"exact-m3 closed without d={exact_m3_closed}; fallback={fallback}; "
        f"candidate c nodes/sweeps={candidate_nodes}/{candidate_sweeps}"
    )
    print(f"all terminal sweeps positive; minimum swept={minimum} at {minimum_family}")
    print(f"terminal-ledger SHA256={terminal_digest.hexdigest()}")
    print(f"certificate-ledger SHA256={certificate_digest.hexdigest()}")
    print(
        f"independent full-union samples={sample_count}; "
        f"sample-manifest SHA256={sample_manifest.hexdigest()}"
    )
    print(
        "proof partition: b>=V2 by THM-735 J3; large c by monotone P2/fixed-E2; "
        "large d by exact fixed-E2 cap; every remaining d swept"
    )
    print(
        "Tournament Analysis vertices=active b proof obligations; observable=(minimum swept "
        "clearance,-sweeps,-b); lex switch; no ties"
    )
    print(
        f"tournament vertices={tournament['vertices']}; score histogram 0.."
        f"{int(tournament['vertices'])-1}:1; directed triangles=0; "
        f"SCCs={tournament['sccs']}; Hamiltonian paths={tournament['hamiltonian_paths']}"
    )
    print(f"tournament path SHA256={tournament['path_digest']}")
    print(
        "preserved carrier: b obligations plus nested survivor components, certificate type, "
        "rational caps, and exact terminal margins"
    )
    print(
        "destroyed quotient: the tournament loses interval geometry; runners, isolated combs, "
        "residues, wall events, and Fano flags lose shared-survivor incidence"
    )
    print(
        "Fano/chi7 audit: (3,4) is one edge address only; no nonidentity Fano transport is used"
    )
    print(
        "Kakeya audit: D_w is a translated periodic one-dimensional needle comb; adaptive "
        "component incidence, not dimension, is the proof-bearing statistic"
    )
    print(
        "VERDICT: E_(3,4) union {15,b,c,d} is lonely for every 15<b<c<d; "
        "the other a branches and global THM-741 remain open"
    )
    print(f"source_sha256={sha256(Path(__file__).resolve())}")
    print("ALL EXACT BRANCH CHECKS PASSED")


if __name__ == "__main__":
    main()
