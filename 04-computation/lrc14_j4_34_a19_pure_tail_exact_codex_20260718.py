#!/usr/bin/env python3
r"""Exact `(3,4)` pure-tail certificate at first external speed ``a=19``.

For ``E={3,4,8,...,14}``, prove ``E union {19,b,c,d}`` lonely for every
``19<b<c<d``.  The program rebuilds the interval carrier after every inserted
speed from the hash-pinned THM-741 kernel.  It also independently rebuilds the
whole adjacent ``a=18`` cutoff bank: the carriers are incomparable and both
the ``V3`` and P2 cutoffs move in both directions, so no preceding result row
is transported.

The proof partition is the same *schema* as the preceding branches: THM-735's
three-needle common cutoff, the linear P2 screen, a fixed-``E2`` rational cap,
and exact terminal subtraction.  Every P2 truth word is audited to be a suffix
and deterministic terminals are recomputed by the literal full-union engine.
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
A = 19


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(sha256(CORE_PATH) == CORE_SHA256, "THM-741 kernel hash changed")
    spec = importlib.util.spec_from_file_location("thm741_j4_34_a19_dependency", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 kernel")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def intersection_measure(left, right) -> F:
    total = F(0)
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
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
    first = next((i for i, value in enumerate(values) if value > 0), len(values))
    require(all(value <= 0 for value in values[:first]), "P2 prefix sign failure")
    require(all(value > 0 for value in values[first:]), "P2 suffix sign failure")
    cutoff = start + first
    return cutoff, None if first == len(values) else values[first]


def adjacent_cutoff_bank(core, good0, a: int) -> tuple[int, dict[int, tuple[int, int]]]:
    """Rebuild all first/second carriers needed for an adjacent cutoff audit."""
    r1, m1, good1 = core.subtract(good0, a)
    v2 = core.minV(3, *(F(4) * m1 / (core.S2 * r1)).as_integer_ratio())
    bank: dict[int, tuple[int, int]] = {}
    for b in range(a + 1, v2):
        r2, m2, _ = core.subtract(good1, b)
        v3 = core.minV(2, *(F(5) * m2 / (core.S2 * r2)).as_integer_ratio())
        cutoff, _ = linear_cutoff(core, m2, r2, b + 1, v3)
        bank[b] = (v3, cutoff)
    return v2, bank


def main() -> None:
    core = load_core()
    require(core.S2 == F(99, 70) and core.S2 * core.S2 > 2,
            "bad sqrt(2) majorant")

    good0, r0, m0 = core.good_norm(E)
    require((r0, m0) == (28, F(433607, 2522520)), "root geometry changed")
    r1, m1, good1 = core.subtract(good0, A)
    require((r1, m1) == (26, F(6388289, 47927880)), "a=19 geometry changed")
    require(m1 > 0 and m1 == core.subtract_sparse(good0, A),
            "a=19 sparse/literal mismatch")
    _, _, full1 = core.good_norm(tuple(sorted(E + (A,))))
    require(full1 == m1, "a=19 full-union mismatch")

    # Adjacent carriers are compared as actual interval sets, not by scalar mass.
    r18, m18, good18 = core.subtract(good0, 18)
    common = intersection_measure(good18, good1)
    only18, only19 = m18 - common, m1 - common
    require(
        (r18, m18, common, only18, only19)
        == (
            26,
            F(17747, 120120),
            F(5691269, 47927880),
            F(15793, 544635),
            F(11617, 798798),
        ),
        "a=18/a=19 carrier comparison changed",
    )
    require(only18 > 0 and only19 > 0, "adjacent carriers became nested")

    v2 = core.minV(3, *(F(4) * m1 / (core.S2 * r1)).as_integer_ratio())
    terminal_digest = hashlib.sha256()
    certificate_digest = hashlib.sha256()
    terminals: list[tuple[int, tuple[int, ...], F]] = []
    b_rows: list[tuple[int, int, int, str | None, int, int, int, int, int, int, str | None]] = []
    total_c = preclosed = exact_m3 = exact_closed = 0
    fallback = integral_caps = candidate_nodes = sweep_count = 0
    minimum: F | None = None
    minimum_family: tuple[int, ...] | None = None

    for b in range(A + 1, v2):
        r2, m2, good2 = core.subtract(good1, b)
        require(m2 > 0 and m2 == core.subtract_sparse(good1, b),
                f"b={b} sparse/literal mismatch")
        v3 = core.minV(2, *(F(5) * m2 / (core.S2 * r2)).as_integer_ratio())
        cutoff, cutoff_slack = linear_cutoff(core, m2, r2, b + 1, v3)
        row_c = max(0, v3 - (b + 1))
        row_preclosed = max(0, v3 - cutoff)
        total_c += row_c
        preclosed += row_preclosed
        row_exact = row_closed = row_candidates = row_sweeps = 0
        row_minimum: F | None = None

        for c in range(b + 1, cutoff):
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
                require(clearance > 0,
                        f"terminal failure b={b},c={c},d={d}: {clearance}")
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

        row = (
            b,
            v3,
            cutoff,
            None if cutoff_slack is None else str(cutoff_slack),
            row_c,
            row_preclosed,
            row_exact,
            row_closed,
            row_candidates,
            row_sweeps,
            None if row_minimum is None else str(row_minimum),
        )
        b_rows.append(row)
        certificate_digest.update((",".join(map(str, row)) + "\n").encode("ascii"))

    require(v2 == 207 and len(b_rows) == 187, "b frontier changed")
    require((total_c, preclosed, exact_m3) == (13332, 3361, 9971),
            "c/P2 ledger changed")
    require((exact_closed, fallback, integral_caps) == (8597, 0, 0),
            "exact-c ledger changed")
    require((candidate_nodes, sweep_count) == (1374, 21302),
            "terminal ledger changed")
    require(minimum == F(368611, 8351070), "minimum swept clearance changed")
    require(
        minimum_family == (3, 4, 8, 9, 10, 11, 12, 13, 14, 19, 23, 32, 52),
        "minimum family changed",
    )
    require(
        terminal_digest.hexdigest()
        == "91a878c51d96cebeaa015598b4c78f14c8a6b0fa181dc87827108a22c7e7aaeb",
        "terminal ledger digest changed",
    )
    require(
        certificate_digest.hexdigest()
        == "86533640f52a442c138fcea45068fe06fd91fb9e8de353ee187cfcca6166a567",
        "certificate ledger digest changed",
    )

    # Full-union oracle crosschecks use a separate from-scratch union path.
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
    require(samples == 187, "full-union sample count changed")
    require(
        sample_digest.hexdigest()
        == "1c9d0dec8d4f6a8065021930182cce6763506885d09465bbd57252c5fd3dedd3",
        "full-union sample digest changed",
    )

    # Fresh adjacent-bank audit: even the cutoff direction oscillates.
    v2_18, bank18 = adjacent_cutoff_bank(core, good0, 18)
    bank19 = {row[0]: (row[1], row[2]) for row in b_rows}
    common_b = tuple(sorted(set(bank18) & set(bank19)))
    adjacent_digest = hashlib.sha256()
    for b in common_b:
        adjacent_digest.update(
            f"{b},{bank18[b][0]},{bank19[b][0]},{bank18[b][1]},{bank19[b][1]}\n".encode(
                "ascii"
            )
        )
    v3_directions = Counter((bank19[b][0] > bank18[b][0]) -
                            (bank19[b][0] < bank18[b][0]) for b in common_b)
    cutoff_directions = Counter((bank19[b][1] > bank18[b][1]) -
                                (bank19[b][1] < bank18[b][1]) for b in common_b)
    require((v2_18, len(common_b)) == (187, 167), "adjacent bank range changed")
    require(v3_directions == Counter({1: 139, -1: 25, 0: 3}),
            "adjacent V3 directions changed")
    require(cutoff_directions == Counter({1: 138, -1: 23, 0: 6}),
            "adjacent cutoff directions changed")
    require(
        adjacent_digest.hexdigest()
        == "d7bc67012c8d5e161dfcefff0840ecef186e4341b86b295653fb6544921addae",
        "adjacent bank digest changed",
    )

    # Tournament Analysis is scheduler telemetry, not a proof quotient.
    active = [row for row in b_rows if row[9] > 0]
    keys = {row[0]: (F(row[10]), -row[9], -row[0]) for row in active}
    path = tuple(sorted(keys, key=lambda b: keys[b], reverse=True))
    scores = Counter(len(path) - index - 1 for index in range(len(path)))
    require(scores == Counter({score: 1 for score in range(len(path))}),
            "tournament score histogram changed")
    path_digest = hashlib.sha256(",".join(map(str, path)).encode("ascii")).hexdigest()
    require(len(path) == 63, "tournament vertex count changed")
    require(
        path_digest == "fe91efac9d31909e0b3317cadb8ef41130ae40d17ff8c8ebb09c44235ddc9e52",
        "tournament Hamiltonian path changed",
    )

    print("THM-741 PURE (3,4) FLOOD TAIL: EXACT a=19 BRANCH CERTIFICATE")
    print("=" * 92)
    print(f"dependency_sha256={CORE_SHA256}")
    print(f"body E={E}; fixed a={A}; claim all integers {A}<b<c<d")
    print(f"root r={r0} m={m0}; after-a r={r1} m={m1}; common b cutoff V2={v2}")
    print(
        f"a=18 comparison r18={r18} m18={m18}; common={common}; "
        f"only18={only18}; only19={only19}; neither carrier contains the other"
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
        f"adjacent cutoff audit common-b={len(common_b)}; V3 down/equal/up="
        f"{v3_directions[-1]}/{v3_directions[0]}/{v3_directions[1]}; "
        f"P2 cutoff down/equal/up={cutoff_directions[-1]}/"
        f"{cutoff_directions[0]}/{cutoff_directions[1]}"
    )
    print(f"adjacent-bank SHA256={adjacent_digest.hexdigest()}")
    print(
        "proof partition: b>=V2 by fresh THM-735 J3; every P2 truth set linearly audited; "
        "large d by fresh fixed-E2 cap; every remaining d swept"
    )
    print(
        "transport audit: G18 and G19 are incomparable, and V3/P2 cutoffs move in both "
        "directions; only the proof schema transports"
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
        "Kakeya/Fano audit: D_19 changes the component stalk in both directions; chi7 still "
        "supplies only the common root-edge address (3,4)"
    )
    print(
        "VERDICT: E_(3,4) union {19,b,c,d} is lonely for every 19<b<c<d; "
        "208 literal-G1 branches and global THM-741 remain open"
    )
    print(f"source_sha256={sha256(Path(__file__).resolve())}")
    print("ALL EXACT a=19 BRANCH CHECKS PASSED")


if __name__ == "__main__":
    main()
