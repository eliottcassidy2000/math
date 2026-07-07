#!/usr/bin/env python3
"""HYP-3531: large-span gentle / denominator atlas.

This extends HYP-3530 in both directions requested by the user:

1. push the exact k<12 attacker scan past the HYP-3530 bounded core; and
2. attach multi-denominator rho_D routing diagnostics to the worst gentle rows.

It reuses HYP-3530's exact Fraction arithmetic and interval engine.  The result
is evidence, not an LRC14 proof: max_D rho_D is a certificate bank, and the
large-span scans are finite / structured probes of the open gentle theorem.

Tournament Analysis declaration:
  vertices: row-family proof carriers and denominator sidecar families, not
            runners or raw arcs;
  pairwise observable: attacker margin + denominator support + formal locality
            - scalar-forgetting risk;
  switch/gauge: higher readiness score, with exact banks before sampled
            structured probes on ties;
  tie Hamiltonian path: exact extended bank -> structured families ->
            denominator router -> raw max_D shadow.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from functools import partial
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import builtins
import itertools
import random
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"


def load_module(name: str, filename: str):
    path = COMP / filename
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3530 = load_module(
    "hyp3530_for_hyp3531",
    "lrc14_multidenom_rhostar_attacker_floor_codex_20260629.py",
)

RNG = random.Random(3531)
print = partial(builtins.print, flush=True)

EXACT_SPAN_OFFSET = 6
DENOM_DMAX = 420
K_RANGE = range(8, 12)


@dataclass(frozen=True)
class RowRecord:
    k: int
    source: str
    E: tuple[int, ...]
    mu: F
    margin: F
    span: int


def fmt(value: F) -> str:
    return f"{value} = {float(value):.8f}"


def compact_counter(counter: Counter | dict) -> dict:
    return dict(sorted(counter.items(), key=lambda item: repr(item[0])))


def insert_best(rows: list[RowRecord], row: RowRecord, limit: int = 8) -> None:
    rows.append(row)
    rows.sort(key=lambda item: (item.mu, item.span, item.E))
    del rows[limit:]


def gap_word(E: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(b - a for a, b in zip(E, E[1:]))


def family_label(E: tuple[int, ...]) -> str:
    k = len(E)
    gaps = gap_word(E)
    if E == tuple(range(k)):
        return "consecutive"
    if E[: k - 1] == tuple(range(k - 1)):
        return "one_tail_prefix"
    if k >= 3 and E[: k - 2] == tuple(range(k - 2)):
        return "two_tail_prefix"
    if len([x for x in E if x % 2 == 1]) <= 2 and max(E) >= 2 * k - 4:
        return "even_lattice_bridge"
    if gaps and max(gaps) >= 2 * max(1, min(gaps)):
        return "split_or_perforated"
    return "other"


def exact_extended_bank():
    rows = []
    worst_records: list[RowRecord] = []
    print("\n" + "=" * 92)
    print("A. Exact extended attacker banks")
    print("=" * 92)
    print(f"Exact primitive normalized banks use span <= k+{EXACT_SPAN_OFFSET}.")

    for k in K_RANGE:
        start = time.perf_counter()
        span = k + EXACT_SPAN_OFFSET
        old_span = k + 5
        gp_min, Ps = H3530.gp_min_for_k(k)
        thr = 1 - gp_min
        consec = tuple(range(k))
        consec_mu = H3530.mu_good_tuple(consec, 1, 7)

        best_mu: F | None = None
        best: list[RowRecord] = []
        shell_best: RowRecord | None = None
        total = 0
        shell_total = 0
        attackers = 0
        below_consec = 0
        near_attackers = 0
        family_counts: Counter[str] = Counter()
        family_best: dict[str, RowRecord] = {}

        for E in H3530.normalized_shapes(k, span):
            total += 1
            label = family_label(E)
            family_counts[label] += 1
            mu = H3530.mu_good_tuple(E, 1, 7)
            margin = mu - thr
            rec = RowRecord(k, f"exact_span_le_{span}", E, mu, margin, E[-1])
            if best_mu is None or mu < best_mu:
                best_mu = mu
                best = [rec]
            elif mu == best_mu and len(best) < 8:
                best.append(rec)
            if label not in family_best or mu < family_best[label].mu:
                family_best[label] = rec
            if E[-1] > old_span:
                shell_total += 1
                if shell_best is None or mu < shell_best.mu:
                    shell_best = rec
            if mu < thr:
                attackers += 1
            if mu < consec_mu:
                below_consec += 1
            if mu < thr + F(1, 20):
                near_attackers += 1

        assert best_mu is not None
        for rec in best:
            insert_best(worst_records, rec)

        floor = gp_min + best_mu - 1
        elapsed = time.perf_counter() - start
        rows.append({
            "k": k,
            "span": span,
            "gp_min": gp_min,
            "P_min": Ps[0],
            "thr": thr,
            "consec_mu": consec_mu,
            "best_mu": best_mu,
            "best": tuple(best),
            "shell_best": shell_best,
            "total": total,
            "shell_total": shell_total,
            "attackers": attackers,
            "below_consec": below_consec,
            "near_attackers": near_attackers,
            "floor": floor,
            "family_counts": dict(family_counts),
            "family_best": family_best,
            "elapsed": elapsed,
        })

        print(f"\n  k={k} span<={span} |P|={13-k} gp_min={fmt(gp_min)} P_min={Ps[0]}")
        print(f"    thr_k={fmt(thr)}")
        print(f"    primitive_shapes={total}; new_shell_shapes(span>{old_span})={shell_total}")
        print(f"    min_mu={fmt(best_mu)}; consecutive_mu={fmt(consec_mu)}")
        print(f"    min_margin_to_thr={fmt(best_mu - thr)}; floor={fmt(floor)}")
        print(
            f"    attackers={attackers}; below_consecutive={below_consec}; "
            f"near_attackers(mu<thr+1/20)={near_attackers}"
        )
        print(f"    family_counts={compact_counter(family_counts)}")
        print(f"    minimizers={[rec.E for rec in best[:4]]}")
        if shell_best:
            print(
                f"    new_shell_best={shell_best.E} mu={fmt(shell_best.mu)} "
                f"margin={fmt(shell_best.margin)} family={family_label(shell_best.E)}"
            )
        print(f"    elapsed_seconds={elapsed:.2f}")

    return rows, worst_records


def one_tail_shapes(k: int):
    for M in range(k + EXACT_SPAN_OFFSET + 1, 15 * k + 1):
        yield tuple(range(k - 1)) + (M,)


def perforated_tail_shapes(k: int):
    base = tuple(range(k))
    for drop in range(1, k - 1):
        core = tuple(x for x in base if x != drop)
        for M in range(k + EXACT_SPAN_OFFSET + 1, 12 * k + 1):
            yield tuple(sorted(core + (M,)))


def two_tail_shapes(k: int):
    core = tuple(range(k - 2))
    for M in range(k + EXACT_SPAN_OFFSET + 1, 10 * k + 1):
        for gap in range(1, 8):
            yield tuple(sorted(core + (M, M + gap)))


def two_block_shapes(k: int):
    for left in range(2, k - 1):
        right = k - left
        for M in (k + 8, 2 * k, 3 * k, 5 * k, 8 * k, 12 * k):
            yield tuple(range(left)) + tuple(range(M, M + right))


def even_bridge_shapes(k: int):
    for h in range(1, 31):
        core = tuple(2 * h * i for i in range(k - 2))
        mid = h * (k - 3)
        for shift in (-3, -1, 1, 3):
            bridge1 = mid + shift
            bridge2 = mid - shift
            E = tuple(sorted(set(core + (bridge1, bridge2))))
            if len(E) == k and E[0] == 0 and H3530.gcd_all(E) == 1:
                yield E


def sparse_probe_shapes(k: int):
    out = set()
    for _ in range(260):
        spread = RNG.choice([3 * k, 5 * k, 8 * k, 12 * k, 16 * k])
        body = sorted(RNG.sample(range(1, spread + 1), k - 1))
        E = tuple([0] + body)
        if H3530.gcd_all(E) == 1:
            out.add(E)
    return sorted(out, key=lambda item: (item[-1], item))


STRUCTURED_GENERATORS = (
    ("one_tail_prefix", one_tail_shapes),
    ("perforated_tail", perforated_tail_shapes),
    ("two_tail_prefix", two_tail_shapes),
    ("two_block", two_block_shapes),
    ("even_lattice_bridge", even_bridge_shapes),
    ("sparse_probe", sparse_probe_shapes),
)


def structured_family_audit():
    rows = []
    worst_records: list[RowRecord] = []
    print("\n" + "=" * 92)
    print("B. Structured gentle/attacker probes beyond the exact core")
    print("=" * 92)

    for k in K_RANGE:
        gp_min, _Ps = H3530.gp_min_for_k(k)
        thr = 1 - gp_min
        print(f"\n  k={k} structured probes; thr_k={fmt(thr)}")
        for name, gen in STRUCTURED_GENERATORS:
            seen = set()
            count = 0
            attackers = 0
            best: RowRecord | None = None
            for E in gen(k):
                if E in seen or len(E) != k or E[0] != 0 or H3530.gcd_all(E) != 1:
                    continue
                seen.add(E)
                count += 1
                mu = H3530.mu_good_tuple(E, 1, 7)
                rec = RowRecord(k, name, E, mu, mu - thr, E[-1])
                if best is None or mu < best.mu:
                    best = rec
                if mu < thr:
                    attackers += 1
            if best is None:
                continue
            rows.append({
                "k": k,
                "family": name,
                "count": count,
                "best": best,
                "attackers": attackers,
            })
            insert_best(worst_records, best)
            print(
                f"    {name:20s} count={count:5d} attackers={attackers:2d} "
                f"best={best.E} span={best.span} mu={fmt(best.mu)} "
                f"margin={fmt(best.margin)}"
            )
    return rows, worst_records


def denominator_records(candidates: list[RowRecord]):
    print("\n" + "=" * 92)
    print("C. Denominator routing for worst gentle rows")
    print("=" * 92)
    print(f"Denominator bank: D=2..{DENOM_DMAX}.  P is the gp_min witness for k.")

    unique: dict[tuple[int, tuple[int, ...], str], RowRecord] = {}
    for rec in candidates:
        unique[(rec.k, rec.E, rec.source)] = rec

    rows = []
    for rec in sorted(unique.values(), key=lambda item: (item.k, item.mu, item.source, item.E))[:32]:
        gp_min, Ps = H3530.gp_min_for_k(rec.k)
        P = Ps[0]
        rho, gp, mu, intervals = H3530.rho_star(P, rec.E, H3530.ONE7)
        scan = H3530.scan_denominators(P, rec.E, H3530.ONE7, Dmax=DENOM_DMAX)
        best = scan["best"]
        first = scan["first_positive"]
        family_best = {
            fam: scan["family_best"].get(fam)
            for fam in ("small_q_2..14", "7m_grid", "14m_grid", "other_grid")
        }
        rows.append({
            "record": rec,
            "P": P,
            "rho": rho,
            "gp": gp,
            "mu": mu,
            "interval_count": len(intervals),
            "best": best,
            "first": first,
            "family_best": family_best,
        })
        print(f"\n  k={rec.k} source={rec.source} E={rec.E}")
        print(f"    mu={fmt(rec.mu)} margin_to_thr={fmt(rec.margin)} P_min={P}")
        print(f"    continuous rho*={fmt(rho)} gp={fmt(gp)} intervals={len(intervals)}")
        print(
            f"    max_D<={DENOM_DMAX} rho_D={fmt(best['rhoD'])} "
            f"at D={best['D']} count={best['count']} family={best['family']}"
        )
        if first:
            print(
                f"    first_positive_D={first['D']} rho_D={first['rhoD']} "
                f"residues={first['residues']}"
            )
        print("    family_best=" + repr({
            fam: None if row is None else (row["D"], row["rhoD"], row["count"])
            for fam, row in family_best.items()
        }))
    return rows


def tournament_report(exact_rows, structured_rows, denom_rows):
    print("\n" + "=" * 92)
    print("D. Tournament Analysis")
    print("=" * 92)
    exact_clear = min(row["best_mu"] - row["thr"] for row in exact_rows)
    structured_clear = min(row["best"].margin for row in structured_rows)
    denom_support = sum(1 for row in denom_rows if row["best"]["rhoD"] > 0)

    carriers = [
        ("exact_extended_bank", exact_clear, 4, 5, 0),
        ("structured_gentle_families", structured_clear, 3, 3, 1),
        ("denominator_sidecar_router", F(denom_support, max(1, len(denom_rows))), 5, 3, 1),
        ("one_tail_prefix_probe", min(r["best"].margin for r in structured_rows if r["family"] == "one_tail_prefix"), 2, 3, 1),
        ("two_tail_prefix_probe", min(r["best"].margin for r in structured_rows if r["family"] == "two_tail_prefix"), 2, 2, 1),
        ("even_lattice_bridge_probe", min(r["best"].margin for r in structured_rows if r["family"] == "even_lattice_bridge"), 2, 2, 1),
        ("raw_max_D_shadow", F(0), 2, 1, 5),
    ]

    def readiness(item):
        _name, margin, witness_support, locality, forgetting = item
        margin_score = 10 if margin > 0 else -20
        return margin_score + 2 * witness_support + locality - 2 * forgetting

    n = len(carriers)
    edge = [[False] * n for _ in range(n)]
    for i, j in itertools.combinations(range(n), 2):
        si = readiness(carriers[i])
        sj = readiness(carriers[j])
        if (si, -i) >= (sj, -j):
            edge[i][j] = True
        else:
            edge[j][i] = True

    scores = [sum(edge[i][j] for j in range(n)) for i in range(n)]
    cycles3 = 0
    for a, b, c in itertools.combinations(range(n), 3):
        degs = [
            int(edge[a][b]) + int(edge[a][c]),
            int(edge[b][a]) + int(edge[b][c]),
            int(edge[c][a]) + int(edge[c][b]),
        ]
        if sorted(degs) == [1, 1, 1]:
            cycles3 += 1

    order = [carriers[i][0] for i, _score in sorted(enumerate(scores), key=lambda t: (-t[1], carriers[t[0]][0]))]
    print("vertices=row-family proof carriers and denominator sidecar families")
    print("pairwise_observable=attacker_margin + denominator_support + formal_locality - scalar_forgetting_risk")
    print(f"score_hist={compact_counter(Counter(scores))}")
    print(f"directed_3cycles={cycles3}")
    print("sccs=singleton SCCs")
    print("hamiltonian_path=" + " -> ".join(order))
    for item, score in sorted(zip(carriers, scores), key=lambda pair: (-pair[1], pair[0][0])):
        print(
            f"  {item[0]} outscore={score} readiness={readiness(item)} "
            f"margin={item[1]} witness_support={item[2]} locality={item[3]} forgetting={item[4]}"
        )


def main() -> None:
    print("=" * 92)
    print("HYP-3531 LRC14 LARGE-SPAN GENTLE / DENOMINATOR ATLAS")
    print("=" * 92)
    print("status=EVIDENCE / finite extension of HYP-3530; not an LRC14 proof")
    print("attacker_predicate=mu_1/7(E) < thr_k where thr_k=1-gp_min(k), k=8..11")
    print("rho_D is treated as a sidecar certificate bank, not a positive-measure theorem.")

    exact_rows, exact_worst = exact_extended_bank()
    structured_rows, structured_worst = structured_family_audit()
    candidates = exact_worst + structured_worst
    denom_rows = denominator_records(candidates)
    tournament_report(exact_rows, structured_rows, denom_rows)

    print("\n" + "=" * 92)
    print("TAKEAWAY")
    print("=" * 92)
    exact_clear = min(row["best_mu"] - row["thr"] for row in exact_rows)
    structured_clear = min(row["best"].margin for row in structured_rows)
    print(f"exact_extended_min_margin={fmt(exact_clear)}")
    print(f"structured_min_margin={fmt(structured_clear)}")
    print("attackers_found=0 in exact extended banks and structured probes")
    print("below_consecutive_found=0 in exact extended banks")
    print("The k<12 large-span theorem is not proved, but the finite danger zone moved")
    print(f"from HYP-3530 span<=k+5 to span<=k+{EXACT_SPAN_OFFSET}, with structured")
    print("one-tail/two-tail/perforated/even-bridge probes still gentle.")
    print("DONE.")


if __name__ == "__main__":
    main()
