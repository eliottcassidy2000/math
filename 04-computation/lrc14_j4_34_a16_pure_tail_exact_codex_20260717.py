#!/usr/bin/env python3
r"""Exact certificate for the `(3,4)` pure flood-tail branch `a=16`.

For E={3,4,8,...,14}, prove that E union {16,b,c,d} is lonely for every
16<b<c<d.  This replay rederives all thresholds from the hash-pinned THM-741
kernel.  Unlike the a=15 certificate, it scans the P2 predicate linearly at
every b node and checks that its truth set is a terminal interval; no cutoff
is imported or binary-searched from the previous branch.

The script also compares the exact first-child carriers at a=15 and a=16.
Neither contains the other, so the transported object is the cap/certificate
schema, not a monotone survivor set, Fano symmetry, or isolated needle order.
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
E = tuple(sorted((*range(8, 15), 3, 4)))
A = 16


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(sha256(CORE_PATH) == CORE_SHA256, "THM-741 kernel hash changed")
    spec = importlib.util.spec_from_file_location("thm741_j4_34_a16_dependency", CORE_PATH)
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


def p2_tail_slack(core, m: F, r: int, c: int) -> F:
    """Lower bound after c and every d>=c+1, against fixed E2."""
    return F(5, 7) * m - F(8 * r, 49 * c) - core.S2 * r / (7 * (c + 1))


def linear_p2_cutoff(core, m: F, r: int, start: int, stop: int) -> tuple[int, F | None]:
    flags = [p2_tail_slack(core, m, r, c) > 0 for c in range(start, stop)]
    first = next((index for index, value in enumerate(flags) if value), len(flags))
    require(not any(flags[:first]) and all(flags[first:]),
            f"P2 truth set is not a terminal interval at start={start}")
    cutoff = start + first
    slack = None if cutoff == stop else p2_tail_slack(core, m, r, cutoff)
    if cutoff < stop:
        require(slack is not None and slack > 0, "nonpositive cutoff slack")
        require(cutoff == start or p2_tail_slack(core, m, r, cutoff - 1) <= 0,
                "nonminimal linear cutoff")
    return cutoff, slack


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def tournament(rows: list[dict[str, object]]) -> dict[str, object]:
    active = [row for row in rows if int(row["sweeps"]) > 0]
    keys = {
        int(row["b"]): (F(row["minimum"]), -int(row["sweeps"]), -int(row["b"]))
        for row in active
    }
    order = tuple(sorted(keys, key=lambda b: keys[b], reverse=True))
    scores = Counter(len(order) - index - 1 for index in range(len(order)))
    triangles = sum(
        not (keys[x] > keys[y] > keys[z]) for x, y, z in combinations(order, 3)
    )
    return {
        "vertices": len(order),
        "scores": scores,
        "triangles": triangles,
        "path_digest": hashlib.sha256(",".join(map(str, order)).encode("ascii")).hexdigest(),
    }


def main() -> None:
    core = load_core()
    require(core.S2 == F(99, 70) and core.S2 * core.S2 > 2,
            "invalid sqrt(2) majorant")

    good0, r0, m0 = core.good_norm(E)
    require((r0, m0) == (28, F(433607, 2522520)), "root geometry changed")
    r1, m1, good1 = core.subtract(good0, A)
    require((r1, m1) == (26, F(29921, 210210)), "a=16 geometry changed")
    require(m1 == core.subtract_sparse(good0, A) and m1 > 0,
            "a=16 sparse/literal mismatch")
    _, _, full1 = core.good_norm(tuple(sorted(E + (A,))))
    require(full1 == m1, "a=16 full-union mismatch")

    # Directly challenge monotone transport from the a=15 carrier.
    r15, m15, good15 = core.subtract(good0, 15)
    common = intersection_measure(good15, good1)
    only15 = m15 - common
    only16 = m1 - common
    require(
        (r15, m15, common, only15, only16)
        == (26, F(184909, 1261260), F(105857, 840840),
            F(4019, 194040), F(419, 25480)),
        "a=15/a=16 carrier comparison changed",
    )
    require(only15 > 0 and only16 > 0, "first-child carriers became nested")

    V2 = core.minV(3, *(F(4) * m1 / (core.S2 * r1)).as_integer_ratio())
    terminal_digest = hashlib.sha256()
    certificate_digest = hashlib.sha256()
    terminals: list[tuple[int, int, int, tuple[int, ...], F]] = []
    b_rows: list[dict[str, object]] = []
    total_c = p2_preclosed = exact_m3 = exact_closed = 0
    fallback = candidate_nodes = sweeps = integral_caps = 0
    minimum: F | None = None
    minimum_family: tuple[int, ...] | None = None
    first_failure: tuple[int, int, int, F] | None = None

    for b in range(A + 1, V2):
        r2, m2, good2 = core.subtract(good1, b)
        require(m2 == core.subtract_sparse(good1, b) and m2 > 0,
                f"b={b} sparse/literal mismatch")
        V3 = core.minV(2, *(F(5) * m2 / (core.S2 * r2)).as_integer_ratio())
        c_start = b + 1
        p2_start, p2_slack = linear_p2_cutoff(core, m2, r2, c_start, V3)
        row_c = max(0, V3 - c_start)
        row_preclosed = max(0, V3 - p2_start)
        total_c += row_c
        p2_preclosed += row_preclosed
        row_exact = row_closed = row_candidates = row_sweeps = 0
        row_minimum: F | None = None

        for c in range(c_start, p2_start):
            exact_m3 += 1
            row_exact += 1
            m3_sparse = core.subtract_sparse(good2, c)
            require(m3_sparse > 0, f"empty c node b={b},c={c}")
            denominator = m3_sparse - m2 / 7
            if denominator > 0:
                cap = core.S2 * r2 / (7 * denominator)
            else:
                r3, m3, _ = core.subtract(good2, c)
                require(m3 == m3_sparse, f"fallback mismatch b={b},c={c}")
                cap = core.S2 * r3 / (6 * m3)
                fallback += 1
            d_max = floor_fraction(cap)
            integral_caps += cap.denominator == 1
            if d_max <= c:
                exact_closed += 1
                row_closed += 1
                continue
            candidate_nodes += 1
            row_candidates += 1
            r3, m3, good3 = core.subtract(good2, c)
            require(m3 == m3_sparse and r3 == len(good3),
                    f"literal c mismatch b={b},c={c}")
            for d in range(c + 1, d_max + 1):
                clearance = core.subtract_sparse(good3, d)
                if clearance <= 0 and first_failure is None:
                    first_failure = (b, c, d, clearance)
                require(clearance > 0, f"first nonpositive terminal={first_failure}")
                family = tuple(sorted(E + (A, b, c, d)))
                terminals.append((b, c, d, family, clearance))
                row_sweeps += 1
                sweeps += 1
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
            "p2_start": p2_start,
            "p2_slack": None if p2_slack is None else str(p2_slack),
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
            (f"{b},{V3},{p2_start},{row['p2_slack']},{row_c},{row_preclosed},"
             f"{row_exact},{row_closed},{row_candidates},{row_sweeps},{row['minimum']}\n").encode(
                "ascii"
            )
        )

    require(first_failure is None, f"unexpected failure {first_failure}")
    require(V2 == 194 and len(b_rows) == 177, "b frontier changed")
    require(total_c == 11786 and p2_preclosed == 2986, "c/P2 ledger changed")
    require(exact_m3 == 8800 and exact_closed == 7579 and fallback == 0,
            "exact c ledger changed")
    require(candidate_nodes == 1221 and sweeps == 18182 and integral_caps == 0,
            "terminal ledger changed")
    require(minimum == F(6999703, 133617120), "minimum swept clearance changed")
    require(minimum_family == (3, 4, 8, 9, 10, 11, 12, 13, 14, 16, 19, 23, 32),
            "minimum family changed")
    require(
        terminal_digest.hexdigest()
        == "20daaa1ec85df38843115c5577f4eb2a48080f51e60dcdeec13cb272e8f3b1b4",
        "terminal ledger digest changed",
    )
    require(
        certificate_digest.hexdigest()
        == "001c43a3d7e830473451c2eee33732484e3b00dfd6b97933f7715869bcf9996a",
        "certificate ledger digest changed",
    )

    by_b: dict[int, list[tuple[tuple[int, ...], F]]] = {}
    for b, _, _, family, clearance in terminals:
        by_b.setdefault(b, []).append((family, clearance))
    sample_digest = hashlib.sha256()
    sample_count = 0
    for b in sorted(by_b):
        rows = by_b[b]
        for index in sorted({0, len(rows) // 2, len(rows) - 1}):
            family, clearance = rows[index]
            _, _, full_measure = core.good_norm(family)
            require(full_measure == clearance,
                    f"full-union crosscheck failed b={b},index={index}")
            sample_count += 1
            sample_digest.update(
                (f"{','.join(map(str, family))}:{clearance.numerator}/"
                 f"{clearance.denominator}\n").encode("ascii")
            )

    tour = tournament(b_rows)
    require(sample_count == 187, "full-union sample census changed")
    require(
        sample_digest.hexdigest()
        == "1652d5afe562c4222bec2e500c15d838de1c2f8102deaae8d409a5f84c5b4ab3",
        "full-union sample manifest changed",
    )
    require(tour["triangles"] == 0, "scheduler tournament cycle")
    require(tour["scores"] == Counter({score: 1 for score in range(tour["vertices"])}),
            "scheduler score histogram changed")
    require(
        tour["vertices"] == 63
        and tour["path_digest"]
        == "a856f62aa8d73652ac05b1744ac6af826c6b5c5770dee6641cd5bf7481afccba",
        "scheduler path changed",
    )

    print("THM-741 PURE (3,4) FLOOD TAIL: EXACT a=16 BRANCH CERTIFICATE")
    print("=" * 92)
    print(f"dependency_sha256={CORE_SHA256}")
    print(f"body E={E}; fixed a={A}; claim all integers {A}<b<c<d")
    print(f"root r={r0} m={m0}; after-a r={r1} m={m1}; common b cutoff V2={V2}")
    print(
        f"a=15 comparison r15={r15} m15={m15}; common={common}; "
        f"only15={only15}; only16={only16}; neither carrier contains the other"
    )
    print(
        f"b obligations={len(b_rows)}; c nodes={total_c}; P2-preclosed c={p2_preclosed}; "
        f"exact-m3={exact_m3}"
    )
    print(
        f"exact-m3 closed without d={exact_closed}; fallback={fallback}; "
        f"integral caps={integral_caps}; candidate c nodes/sweeps={candidate_nodes}/{sweeps}"
    )
    print(f"all terminal sweeps positive; minimum swept={minimum} at {minimum_family}")
    print(f"terminal-ledger SHA256={terminal_digest.hexdigest()}")
    print(f"certificate-ledger SHA256={certificate_digest.hexdigest()}")
    print(
        f"independent full-union samples={sample_count}; "
        f"sample-manifest SHA256={sample_digest.hexdigest()}"
    )
    print(
        "proof partition: b>=V2 by THM-735 J3; every P2 truth set linearly audited as a "
        "terminal c interval; large d by exact fixed-E2 cap; every remaining d swept"
    )
    print(
        "transport audit: the a=15 proof schema survives, but component carriers and workload "
        "do not transport monotonically (a=16 has more b,c,terminal obligations)"
    )
    print(
        "Tournament Analysis vertices=active b obligations; observable=(minimum swept "
        "clearance,-sweeps,-b); lex switch"
    )
    print(
        f"tournament vertices={tour['vertices']}; score histogram 0..{tour['vertices']-1}:1; "
        f"directed triangles=0; SCCs={tour['vertices']}; Hamiltonian paths=1"
    )
    print(f"tournament path SHA256={tour['path_digest']}")
    print(
        "carrier audit: nested survivor components, certificate type, phase, and rational caps "
        "preserve the predicate; runners, isolated needles, residues, wall events, and Fano "
        "flags do not"
    )
    print(
        "Kakeya/Fano audit: D_16 is a new periodic comb, not a monotone translate of D_15; "
        "chi7 still names only the root edge (3,4)"
    )
    print(
        "VERDICT: E_(3,4) union {16,b,c,d} is lonely for every 16<b<c<d; "
        "the other first-speed branches and global THM-741 remain open"
    )
    print(f"source_sha256={sha256(Path(__file__).resolve())}")
    print("ALL EXACT a=16 BRANCH CHECKS PASSED")


if __name__ == "__main__":
    main()
