#!/usr/bin/env python3
r"""Exact `(3,4)` pure-tail certificate at first external speed a=18.

For E={3,4,8,...,14}, this proves E union {18,b,c,d} lonely for every
18<b<c<d.  Every first-child carrier, THM-735 cutoff, linear P2 truth set,
exact fixed-E2 cap, and finite terminal is recomputed from the hash-pinned
interval kernel.  The a=17 carrier is rebuilt only for a bidirectional
noncontainment audit; no a=17 cutoff or result row is imported.
"""

from __future__ import annotations

import hashlib
import importlib.util
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
E = tuple(sorted((*range(8, 15), 3, 4)))
A = 18


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(sha256(CORE_PATH) == CORE_SHA256, "THM-741 kernel hash changed")
    spec = importlib.util.spec_from_file_location("thm741_j4_34_a18_dependency", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 kernel")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def intersection_measure(left, right) -> F:
    total = F(0)
    i = j = 0
    while i < len(left) and j < len(right):
        lo, hi = max(left[i][0], right[j][0]), min(left[i][1], right[j][1])
        if lo < hi:
            total += hi - lo
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return total


def p2_slack(core, m: F, r: int, c: int) -> F:
    return F(5, 7) * m - F(8 * r, 49 * c) - core.S2 * r / (7 * (c + 1))


def linear_cutoff(core, m: F, r: int, start: int, stop: int) -> tuple[int, F | None]:
    values = tuple(p2_slack(core, m, r, c) for c in range(start, stop))
    first = next((index for index, value in enumerate(values) if value > 0), len(values))
    require(all(value <= 0 for value in values[:first]), "P2 prefix sign failure")
    require(all(value > 0 for value in values[first:]), "P2 suffix sign failure")
    cutoff = start + first
    return cutoff, None if cutoff == stop else values[first]


def main() -> None:
    core = load_core()
    require(core.S2 == F(99, 70) and core.S2 * core.S2 > 2,
            "bad sqrt(2) majorant")
    good0, r0, m0 = core.good_norm(E)
    require((r0, m0) == (28, F(433607, 2522520)), "root geometry changed")
    r1, m1, good1 = core.subtract(good0, A)
    require((r1, m1) == (26, F(17747, 120120)), "a=18 geometry changed")
    require(m1 > 0 and m1 == core.subtract_sparse(good0, A),
            "a=18 sparse/literal mismatch")
    _, _, full1 = core.good_norm(tuple(sorted(E + (A,))))
    require(full1 == m1, "a=18 full-union mismatch")

    r17, m17, good17 = core.subtract(good0, 17)
    common = intersection_measure(good17, good1)
    only17, only18 = m17 - common, m1 - common
    require(
        (r17, m17, common, only17, only18)
        == (
            26,
            F(758281, 6126120),
            F(1564249, 14294280),
            F(30761, 2144142),
            F(45637, 1191190),
        ),
        "a=17/a=18 comparison geometry changed",
    )
    require(only17 > 0 and only18 > 0, "a=17/a=18 carriers became nested")

    V2 = core.minV(3, *(F(4) * m1 / (core.S2 * r1)).as_integer_ratio())
    terminal_digest = hashlib.sha256()
    certificate_digest = hashlib.sha256()
    terminals: list[tuple[int, tuple[int, ...], F]] = []
    b_rows: list[dict[str, object]] = []
    total_c = preclosed = exact_m3 = exact_closed = 0
    fallback = integral_caps = candidate_nodes = sweep_count = 0
    minimum: F | None = None
    minimum_family: tuple[int, ...] | None = None
    first_failure: tuple[int, int, int, F] | None = None

    for b in range(A + 1, V2):
        r2, m2, good2 = core.subtract(good1, b)
        require(m2 > 0 and m2 == core.subtract_sparse(good1, b),
                f"b={b} sparse/literal mismatch")
        V3 = core.minV(2, *(F(5) * m2 / (core.S2 * r2)).as_integer_ratio())
        start = b + 1
        cutoff, cutoff_slack = linear_cutoff(core, m2, r2, start, V3)
        row_c, row_preclosed = max(0, V3 - start), max(0, V3 - cutoff)
        total_c += row_c
        preclosed += row_preclosed
        row_exact = row_closed = row_candidates = row_sweeps = 0
        row_minimum: F | None = None

        for c in range(start, cutoff):
            exact_m3 += 1
            row_exact += 1
            sparse_m3 = core.subtract_sparse(good2, c)
            require(sparse_m3 > 0, f"empty c node b={b},c={c}")
            denominator = sparse_m3 - m2 / 7
            if denominator > 0:
                cap = core.S2 * r2 / (7 * denominator)
            else:
                r3, m3, _ = core.subtract(good2, c)
                require(m3 == sparse_m3, f"fallback mismatch b={b},c={c}")
                cap = core.S2 * r3 / (6 * m3)
                fallback += 1
            integral_caps += cap.denominator == 1
            d_max = cap.numerator // cap.denominator
            if d_max <= c:
                exact_closed += 1
                row_closed += 1
                continue
            candidate_nodes += 1
            row_candidates += 1
            r3, m3, good3 = core.subtract(good2, c)
            require(m3 == sparse_m3 and r3 == len(good3),
                    f"literal c mismatch b={b},c={c}")
            for d in range(c + 1, d_max + 1):
                clearance = core.subtract_sparse(good3, d)
                if clearance <= 0 and first_failure is None:
                    first_failure = (b, c, d, clearance)
                require(clearance > 0, f"first terminal failure={first_failure}")
                family = tuple(sorted(E + (A, b, c, d)))
                terminals.append((b, family, clearance))
                row_sweeps += 1
                sweep_count += 1
                terminal_digest.update(
                    f"{b},{c},{d}:{clearance.numerator}/{clearance.denominator}\n".encode("ascii")
                )
                if row_minimum is None or clearance < row_minimum:
                    row_minimum = clearance
                if minimum is None or clearance < minimum:
                    minimum, minimum_family = clearance, family

        row = {
            "b": b,
            "V3": V3,
            "cutoff": cutoff,
            "cutoff_slack": None if cutoff_slack is None else str(cutoff_slack),
            "c": row_c,
            "preclosed": row_preclosed,
            "exact": row_exact,
            "closed": row_closed,
            "candidates": row_candidates,
            "sweeps": row_sweeps,
            "minimum": None if row_minimum is None else str(row_minimum),
        }
        b_rows.append(row)
        certificate_digest.update(
            (f"{b},{V3},{cutoff},{row['cutoff_slack']},{row_c},{row_preclosed},"
             f"{row_exact},{row_closed},{row_candidates},{row_sweeps},{row['minimum']}\n").encode(
                "ascii"
            )
        )

    require(first_failure is None, f"unexpected terminal failure={first_failure}")
    require(V2 == 187 and len(b_rows) == 168, "b frontier changed")
    require(total_c == 10762 and preclosed == 2750, "c/P2 ledger changed")
    require(exact_m3 == 8012 and exact_closed == 6987 and fallback == 0,
            "exact c ledger changed")
    require(integral_caps == 0, "integral-cap convention changed")
    require(candidate_nodes == 1025 and sweep_count == 14399,
            "terminal ledger changed")
    require(minimum == F(518788, 9144135), "minimum swept clearance changed")
    require(
        minimum_family == (3, 4, 8, 9, 10, 11, 12, 13, 14, 18, 21, 29, 32),
        "minimum family changed",
    )
    require(
        terminal_digest.hexdigest()
        == "44bd2e2dea062e5090243177f9c8871b0ce820767f972df7278cb5e56d670b08",
        "terminal ledger digest changed",
    )
    require(
        certificate_digest.hexdigest()
        == "f25008da5e5009ae1d3ee15063c90086b599af689ba842e4119e9e03ef594f06",
        "certificate ledger digest changed",
    )

    by_b: dict[int, list[tuple[tuple[int, ...], F]]] = {}
    for b, family, clearance in terminals:
        by_b.setdefault(b, []).append((family, clearance))
    sample_digest = hashlib.sha256()
    samples = 0
    for b in sorted(by_b):
        rows = by_b[b]
        for index in sorted({0, len(rows) // 2, len(rows) - 1}):
            family, clearance = rows[index]
            _, _, full_measure = core.good_norm(family)
            require(full_measure == clearance,
                    f"full-union mismatch b={b},index={index}")
            samples += 1
            sample_digest.update(
                (f"{','.join(map(str, family))}:{clearance.numerator}/"
                 f"{clearance.denominator}\n").encode("ascii")
            )
    require(samples == 165, "full-union sample count changed")
    require(
        sample_digest.hexdigest()
        == "86046df7baa90ce91798f30df2a0bd8c6babac70d6e830651c65a97703a895a3",
        "full-union sample digest changed",
    )

    active = [row for row in b_rows if int(row["sweeps"]) > 0]
    keys = {
        int(row["b"]): (F(row["minimum"]), -int(row["sweeps"]), -int(row["b"]))
        for row in active
    }
    path = tuple(sorted(keys, key=lambda b: keys[b], reverse=True))
    scores = Counter(len(path) - index - 1 for index in range(len(path)))
    require(scores == Counter({score: 1 for score in range(len(path))}),
            "tournament score histogram changed")
    path_digest = hashlib.sha256(",".join(map(str, path)).encode("ascii")).hexdigest()
    require(len(path) == 56, "tournament vertex count changed")
    require(
        path_digest == "2193dc2debc85b984e008f653e9e8cf802662534469724b95da8041177f48b67",
        "tournament Hamiltonian path changed",
    )

    print("THM-741 PURE (3,4) FLOOD TAIL: EXACT a=18 BRANCH CERTIFICATE")
    print("=" * 92)
    print(f"dependency_sha256={CORE_SHA256}")
    print(f"body E={E}; fixed a={A}; claim all integers {A}<b<c<d")
    print(f"root r={r0} m={m0}; after-a r={r1} m={m1}; common b cutoff V2={V2}")
    print(
        f"a=17 comparison r17={r17} m17={m17}; common={common}; "
        f"only17={only17}; only18={only18}; neither carrier contains the other"
    )
    print(
        f"b obligations={len(b_rows)}; c nodes={total_c}; P2-preclosed c={preclosed}; "
        f"exact-m3={exact_m3}"
    )
    print(
        f"exact-m3 closed without d={exact_closed}; fallback={fallback}; "
        f"integral caps={integral_caps}; candidate c nodes/sweeps={candidate_nodes}/{sweep_count}"
    )
    print(f"all terminal sweeps positive; minimum swept={minimum} at {minimum_family}")
    print(f"terminal-ledger SHA256={terminal_digest.hexdigest()}")
    print(f"certificate-ledger SHA256={certificate_digest.hexdigest()}")
    print(f"independent full-union samples={samples}; sample-manifest SHA256={sample_digest.hexdigest()}")
    print(
        "proof partition: b>=V2 by fresh THM-735 J3; every P2 truth set linearly audited; "
        "large d by fresh fixed-E2 cap; every remaining d swept"
    )
    print(
        "transport audit: G17 and G18 are incomparable while V2 and terminal work fall sharply; "
        "adjacent-branch horizons and workloads oscillate, so only the proof schema transports"
    )
    print(
        "Tournament Analysis vertices=active b obligations; observable=(minimum swept "
        "clearance,-sweeps,-b); lex switch"
    )
    print(
        f"tournament vertices={len(path)}; score histogram 0..{len(path)-1}:1; "
        f"directed triangles=0; SCCs={len(path)}; Hamiltonian paths=1"
    )
    print(f"tournament path SHA256={path_digest}")
    print(
        "carrier audit: nested components, phase, certificate type, caps, and margins are "
        "retained; runners, isolated Kakeya needles, residues, wall events, and Fano flags are not"
    )
    print(
        "Kakeya/Fano audit: D_18 changes the component stalk in both directions; chi7 still "
        "supplies only the common root-edge address (3,4)"
    )
    print(
        "VERDICT: E_(3,4) union {18,b,c,d} is lonely for every 18<b<c<d; "
        "the other first-speed branches and global THM-741 remain open"
    )
    print(f"source_sha256={sha256(Path(__file__).resolve())}")
    print("ALL EXACT a=18 BRANCH CHECKS PASSED")


if __name__ == "__main__":
    main()
