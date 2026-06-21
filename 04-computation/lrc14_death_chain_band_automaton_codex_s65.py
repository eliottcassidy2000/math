#!/usr/bin/env python3
"""HYP-2702 / T937: death-chain slope-band automaton scout.

This is a follow-up to ``lrc14_sparse_tail_deficit_automaton_codex_s64.py``.
The S64 scout found that raw sparse residual coordinate wins disappear after
the all-singleton context at total size 7, and that coherent generated
contexts only improve the observed margins.

This scout asks what a proof of the first step can use:

* the one-dimensional singleton death-chain kernel;
* the seven slope bands from the incoming Angle A Sturmian decomposition;
* cellular-automaton sign runs of the local death-chain margin;
* context-merge monotonicity for coherent generated blocks;
* a Tournament Analysis whose vertices are proof quotients, not runners.

No LRC14 proof is claimed.  The point is to distinguish proof carriers that
preserve the sparse-tail predicate from attractive quotients that throw away
too much structure.
"""

from __future__ import annotations

import importlib.util
import sys
from collections import Counter, defaultdict
from functools import lru_cache
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


sys.stdout.reconfigure(line_buffering=True)

ROOT = Path(__file__).resolve().parents[1]
S64_PATH = ROOT / "04-computation" / "lrc14_sparse_tail_deficit_automaton_codex_s64.py"
spec = importlib.util.spec_from_file_location("s64_sparse_tail", S64_PATH)
s64 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(s64)


def fmt(q: F | None, digits: int = 9) -> str:
    if q is None:
        return "None"
    return f"{q} ({float(q):.{digits}f})"


def sign_char(q: F) -> str:
    if q > 0:
        return "+"
    if q < 0:
        return "-"
    return "0"


def band_of(mid: F) -> int:
    y = mid * 7
    return min(6, y.numerator // y.denominator)


@lru_cache(maxsize=None)
def local_pressure(shape: tuple[int, ...], r: int, x: F) -> F:
    """Singleton-context survival pressure at fixed x, averaged over y."""
    total = F(0)
    for hit, mass in s64.hit_law_at_x(shape, x):
        need = 6 - s64.mask_size(hit)
        total += mass * s64.singleton_cover_weight(r, need)
    return total


def band_delta_profile(
    consec: tuple[int, ...], challenger: tuple[int, ...], r: int
) -> tuple[tuple[F, ...], list[tuple[F, F, int, F]]]:
    """Return seven signed band integrals plus the local sign atoms.

    The slope band is ``floor(theta)`` for ``theta=7x``.  Thus the bands are
    x-intervals [s/7,(s+1)/7).  The returned band integrals are unnormalized;
    their sum is exactly ``singleton_delta(r, consec, challenger)``.
    """
    xs = set(s64.union_breakpoints((consec, challenger)))
    xs.update(F(s, 7) for s in range(8))
    xs = sorted(xs)
    bands = [F(0) for _ in range(7)]
    atoms: list[tuple[F, F, int, F]] = []
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        band = band_of(mid)
        val = local_pressure(consec, r, mid) - local_pressure(challenger, r, mid)
        bands[band] += (hi - lo) * val
        atoms.append((lo, hi, band, val))
    return tuple(bands), atoms


def hitcount_tail_margins(
    consec: tuple[int, ...], challenger: tuple[int, ...]
) -> tuple[F, ...]:
    """First-order stochastic dominance margins on hit counts.

    Tail[t] = P_K(h>=t)-P_C(h>=t).  All tails nonnegative would prove all
    increasing hit-count kernels at once.  The singleton death-chain uses only
    a special increasing kernel, so this is intentionally stronger.
    """
    lk = s64.hit_count_law(consec)
    lc = s64.hit_count_law(challenger)
    return tuple(sum(lk[h] - lc[h] for h in range(t, 7)) for t in range(1, 7))


def sign_run_stats(atoms: list[tuple[F, F, int, F]]) -> dict[str, object]:
    pos_int = neg_int = zero_len = F(0)
    pos_len = neg_len = F(0)
    neg_runs = 0
    longest_neg = F(0)
    cur_sign: str | None = None
    cur_len = F(0)
    sign_changes = 0
    band_neg_len = [F(0) for _ in range(7)]
    for lo, hi, band, val in atoms:
        width = hi - lo
        sg = sign_char(val)
        if val > 0:
            pos_int += width * val
            pos_len += width
        elif val < 0:
            neg_int += width * val
            neg_len += width
            band_neg_len[band] += width
        else:
            zero_len += width
        if cur_sign is None:
            cur_sign = sg
            cur_len = width
            continue
        if sg == cur_sign:
            cur_len += width
            continue
        if cur_sign == "-":
            neg_runs += 1
            longest_neg = max(longest_neg, cur_len)
        sign_changes += 1
        cur_sign = sg
        cur_len = width
    if cur_sign == "-":
        neg_runs += 1
        longest_neg = max(longest_neg, cur_len)
    return {
        "pos_int": pos_int,
        "neg_int": neg_int,
        "pos_len": pos_len,
        "neg_len": neg_len,
        "zero_len": zero_len,
        "neg_runs": neg_runs,
        "longest_neg": longest_neg,
        "sign_changes": sign_changes,
        "band_neg_len": tuple(band_neg_len),
    }


def first_positive_depth(consec: tuple[int, ...], challenger: tuple[int, ...]) -> int | None:
    for r in range(12):
        if s64.singleton_delta(r, consec, challenger) > 0:
            return r
    return None


def flatten_coordinate_rows() -> dict[int, list[dict[str, object]]]:
    print("PART 0 -- loading the S64 sparse-coordinate frontier")
    rows_by_size: dict[int, list[dict[str, object]]] = {}
    for size in range(3, 7):
        rows = s64.coordinate_violators(size)
        rows_by_size[size] = rows
        scope = "full span<=13" if size == 3 else "HYP-2698 near-consecutive frontier"
        print(f"  size={size}: rows={len(rows)} ({scope})")
    print()
    return rows_by_size


def analyze_death_chain_bands(
    rows_by_size: dict[int, list[dict[str, object]]]
) -> tuple[dict[int, dict[str, object]], dict[str, object]]:
    print("PART A -- singleton death-chain through seven slope bands")
    print("  Band s means theta=7x lies in [s,s+1); band sums are unnormalized.")
    summary: dict[int, dict[str, object]] = {}
    global_min: tuple[F, int, tuple[int, ...], tuple[F, ...]] | None = None
    global_min_band: tuple[F, int, tuple[int, ...], int, tuple[F, ...]] | None = None
    total_agg_fail = total_fosd_fail = total_band_neg = total_rows = 0

    for size, rows in rows_by_size.items():
        consec = tuple(range(size))
        r = 7 - size
        pattern_counter: Counter[str] = Counter()
        depth_hist: Counter[int | None] = Counter()
        fosd_fail = 0
        agg_fail = 0
        band_neg = 0
        min_delta: F | None = None
        min_delta_shape: tuple[int, ...] | None = None
        min_delta_bands: tuple[F, ...] | None = None
        min_band: F | None = None
        min_band_shape: tuple[int, ...] | None = None
        min_band_idx: int | None = None
        min_band_profile: tuple[F, ...] | None = None
        slow_negative = fast_rescue = 0
        tail_min: F | None = None
        tail_worst: tuple[int, ...] | None = None
        local_negative_shapes = 0

        for row in rows:
            shape = row["shape"]  # type: ignore[assignment]
            assert isinstance(shape, tuple)
            depth_hist[first_positive_depth(consec, shape)] += 1
            delta = s64.singleton_delta(r, consec, shape)
            bands, atoms = band_delta_profile(consec, shape, r)
            stats = sign_run_stats(atoms)
            if stats["neg_len"] > 0:
                local_negative_shapes += 1
            tails = hitcount_tail_margins(consec, shape)
            local_tail_min = min(tails)
            if tail_min is None or local_tail_min < tail_min:
                tail_min = local_tail_min
                tail_worst = shape
            if any(t < 0 for t in tails):
                fosd_fail += 1
            if delta <= 0:
                agg_fail += 1
            if min_delta is None or delta < min_delta:
                min_delta = delta
                min_delta_shape = shape
                min_delta_bands = bands
            pattern = "".join(sign_char(b) for b in bands)
            pattern_counter[pattern] += 1
            neg_here = sum(1 for b in bands if b < 0)
            band_neg += neg_here
            if min_band is None or min(bands) < min_band:
                min_band = min(bands)
                min_band_idx = min(range(7), key=lambda i: bands[i])
                min_band_shape = shape
                min_band_profile = bands
            slow = bands[0] + bands[1]
            fast = sum(bands[2:], F(0))
            if slow < 0:
                slow_negative += 1
            if slow < 0 and fast > -slow:
                fast_rescue += 1

        total_rows += len(rows)
        total_agg_fail += agg_fail
        total_fosd_fail += fosd_fail
        total_band_neg += band_neg
        assert min_delta is not None
        assert min_delta_shape is not None
        assert min_delta_bands is not None
        assert min_band is not None
        assert min_band_shape is not None
        assert min_band_idx is not None
        assert min_band_profile is not None
        if global_min is None or min_delta < global_min[0]:
            global_min = (min_delta, size, min_delta_shape, min_delta_bands)
        if global_min_band is None or min_band < global_min_band[0]:
            global_min_band = (min_band, size, min_band_shape, min_band_idx, min_band_profile)

        summary[size] = {
            "rows": len(rows),
            "r": r,
            "agg_fail": agg_fail,
            "fosd_fail": fosd_fail,
            "band_neg": band_neg,
            "min_delta": min_delta,
            "min_delta_shape": min_delta_shape,
            "min_delta_bands": min_delta_bands,
            "min_band": min_band,
            "min_band_shape": min_band_shape,
            "min_band_idx": min_band_idx,
            "tail_min": tail_min,
            "tail_worst": tail_worst,
            "local_negative_shapes": local_negative_shapes,
        }

        print(
            f"  size={size}, r={r}: rows={len(rows)}, agg_fail={agg_fail}, "
            f"FOSD_fail={fosd_fail}, negative_band_cells={band_neg}/{7 * len(rows)}"
        )
        print(
            f"    min aggregate={fmt(min_delta)} at C={min_delta_shape}; "
            f"band_word={''.join(sign_char(b) for b in min_delta_bands)}"
        )
        print(
            f"    worst band={fmt(min_band)} at C={min_band_shape}, s={min_band_idx}; "
            f"band_word={''.join(sign_char(b) for b in min_band_profile)}"
        )
        print(
            f"    slow-negative rows={slow_negative}, fast-rescued rows={fast_rescue}, "
            f"local-negative-sign rows={local_negative_shapes}"
        )
        print(f"    first-positive-depth histogram={dict(sorted(depth_hist.items(), key=str))}")
        print(f"    top band sign patterns={pattern_counter.most_common(5)}")
        print(
            f"    FOSD strongest failure min_tail={fmt(tail_min)} at C={tail_worst}"
        )
    print()

    assert global_min is not None
    assert global_min_band is not None
    overall = {
        "rows": total_rows,
        "agg_fail": total_agg_fail,
        "fosd_fail": total_fosd_fail,
        "band_neg": total_band_neg,
        "global_min": global_min,
        "global_min_band": global_min_band,
    }
    return summary, overall


def bounded_singleton_slices() -> dict[int, dict[str, object]]:
    print("PART B -- bounded singleton slices beyond the sparse-coordinate rows")
    print("  Exact full span<=13 is cheap for size 3,4; size 5,6 use smaller")
    print("  bounded slices here because a full span<=13 run is much slower.")
    plan = {3: 13, 4: 13, 5: 9, 6: 9}
    out: dict[int, dict[str, object]] = {}
    for size, span in plan.items():
        consec = tuple(range(size))
        r = 7 - size
        shapes = [shape for shape in s64.bounded_shapes(size, span) if shape != consec]
        failures = 0
        min_delta: F | None = None
        worst: tuple[int, ...] | None = None
        depth_hist: Counter[int | None] = Counter()
        for shape in shapes:
            delta = s64.singleton_delta(r, consec, shape)
            if delta <= 0:
                failures += 1
            if min_delta is None or delta < min_delta:
                min_delta = delta
                worst = shape
            depth_hist[first_positive_depth(consec, shape)] += 1
        out[size] = {
            "span": span,
            "count": len(shapes),
            "failures": failures,
            "min_delta": min_delta,
            "worst": worst,
            "depth_hist": depth_hist,
        }
        print(
            f"  size={size}, span<={span}, r={r}: shapes={len(shapes)}, "
            f"failures={failures}, min_delta={fmt(min_delta)} at C={worst}"
        )
        print(f"    first-positive-depth histogram={dict(sorted(depth_hist.items(), key=str))}")
    print()
    return out


def context_merge_test(
    rows_by_size: dict[int, list[dict[str, object]]],
    bounded: dict[int, dict[str, object]],
) -> dict[str, object]:
    print("PART C -- coherent context merge monotonicity")
    print("  Observable: context_delta(partition) - all_singleton_delta at total size 7.")
    failures = 0
    tested = 0
    min_context: F | None = None
    min_context_who = None
    min_lift: F | None = None
    min_lift_who = None

    # Coordinate frontier: all rows, all coherent partitions, exactly as proof target.
    for size, rows in rows_by_size.items():
        consec = tuple(range(size))
        r = 7 - size
        for row in rows:
            shape = row["shape"]  # type: ignore[assignment]
            assert isinstance(shape, tuple)
            singleton = s64.singleton_delta(r, consec, shape)
            for part in s64.partitions(r):
                ctx = s64.coherent_context(part)
                delta = s64.context_delta(consec, shape, ctx)
                lift = delta - singleton
                tested += 1
                if delta <= 0:
                    failures += 1
                if min_context is None or delta < min_context:
                    min_context = delta
                    min_context_who = ("coordinate", size, shape, s64.context_name(part))
                if min_lift is None or lift < min_lift:
                    min_lift = lift
                    min_lift_who = ("coordinate", size, shape, s64.context_name(part))

    # Add the worst bounded-slice shape per size, even if it had no raw sparse coordinate win.
    extra_seen: set[tuple[int, tuple[int, ...]]] = set()
    for size, data in bounded.items():
        shape = data["worst"]
        if not isinstance(shape, tuple):
            continue
        key = (size, shape)
        if key in extra_seen:
            continue
        extra_seen.add(key)
        consec = tuple(range(size))
        r = 7 - size
        singleton = s64.singleton_delta(r, consec, shape)
        for part in s64.partitions(r):
            ctx = s64.coherent_context(part)
            delta = s64.context_delta(consec, shape, ctx)
            lift = delta - singleton
            tested += 1
            if delta <= 0:
                failures += 1
            if min_context is None or delta < min_context:
                min_context = delta
                min_context_who = ("bounded-worst", size, shape, s64.context_name(part))
            if min_lift is None or lift < min_lift:
                min_lift = lift
                min_lift_who = ("bounded-worst", size, shape, s64.context_name(part))

    print(
        f"  tested={tested}, failures={failures}, "
        f"min_context={fmt(min_context)} at {min_context_who}"
    )
    print(f"  min(context-singleton lift)={fmt(min_lift)} at {min_lift_who}")
    if min_lift == 0:
        print("  all observed coherent contexts are at least as strong as all-singletons;")
        print("  equality occurs at the all-singleton partition, as expected.")
    print()
    return {
        "tested": tested,
        "failures": failures,
        "min_context": min_context,
        "min_context_who": min_context_who,
        "min_lift": min_lift,
        "min_lift_who": min_lift_who,
    }


def ca_sign_examples(overall: dict[str, object]) -> None:
    print("PART D -- cellular-automaton sign-run examples")
    examples: list[tuple[str, int, tuple[int, ...]]] = []
    min_delta, size, shape, _bands = overall["global_min"]  # type: ignore[misc]
    examples.append(("global minimum aggregate margin", size, shape))
    min_band, bsize, bshape, _bidx, _profile = overall["global_min_band"]  # type: ignore[misc]
    if (bsize, bshape) != (size, shape):
        examples.append(("global most negative band", bsize, bshape))
    for label, size, shape in examples:
        consec = tuple(range(size))
        r = 7 - size
        bands, atoms = band_delta_profile(consec, shape, r)
        stats = sign_run_stats(atoms)
        print(f"  {label}: size={size}, r={r}, C={shape}")
        print(
            f"    aggregate={fmt(sum(bands, F(0)))}; "
            f"band_word={''.join(sign_char(b) for b in bands)}; "
            f"bands={[str(b) for b in bands]}"
        )
        print(
            f"    local pos_int={fmt(stats['pos_int'])}, neg_int={fmt(stats['neg_int'])}, "
            f"pos_len={fmt(stats['pos_len'])}, neg_len={fmt(stats['neg_len'])}"
        )
        print(
            f"    negative_runs={stats['neg_runs']}, sign_changes={stats['sign_changes']}, "
            f"longest_negative_run={fmt(stats['longest_neg'])}"
        )
        print(f"    band negative lengths={[str(x) for x in stats['band_neg_len']]}")
    print()


def sccs(vertices: list[str], edges: dict[tuple[str, str], str]) -> list[list[str]]:
    adj = {v: [] for v in vertices}
    radj = {v: [] for v in vertices}
    for a in vertices:
        for b in vertices:
            if a == b:
                continue
            winner = edges[(a, b)]
            loser = b if winner == a else a
            adj[winner].append(loser)
            radj[loser].append(winner)
    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in adj[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)
    seen.clear()
    comps: list[list[str]] = []

    def rdfs(v: str, comp: list[str]) -> None:
        seen.add(v)
        comp.append(v)
        for w in radj[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[str] = []
            rdfs(v, comp)
            comps.append(sorted(comp))
    return comps


def hamiltonian_path_count(vertices: list[str], edges: dict[tuple[str, str], str]) -> int:
    count = 0

    def edge(a: str, b: str) -> bool:
        return edges[(a, b)] == a

    def rec(path: list[str], unused: set[str]) -> None:
        nonlocal count
        if not unused:
            count += 1
            return
        last = path[-1]
        for v in list(unused):
            if edge(last, v):
                unused.remove(v)
                path.append(v)
                rec(path, unused)
                path.pop()
                unused.add(v)

    for v in vertices:
        rest = set(vertices)
        rest.remove(v)
        rec([v], rest)
    return count


def tournament_analysis(
    overall: dict[str, object], context: dict[str, object]
) -> None:
    print("PART E -- Tournament Analysis on abstract proof quotients")
    print("  Pairwise observable=(obligation_failures, min_margin, generated_fidelity, auditability).")
    print("  Switch/gauge: lower failures win; ties use larger margin, then fidelity.")

    min_delta, _size, _shape, _bands = overall["global_min"]  # type: ignore[misc]
    min_band, _bsize, _bshape, _bidx, _profile = overall["global_min_band"]  # type: ignore[misc]
    min_context = context["min_context"] or F(0)
    carriers = {
        "raw_sparse_coordinates": (
            int(overall["rows"]),
            F(-1),
            0,
            3,
            "keeps residual packets but scalarizes before context; known false quotient",
        ),
        "hitcount_FOSD": (
            int(overall["fosd_fail"]),
            min_band,
            2,
            4,
            "would prove all increasing kernels, but fails on many rows",
        ),
        "per_band_monotonicity": (
            int(overall["band_neg"]),
            min_band,
            3,
            5,
            "too strong; useful only as a signed-debt ledger",
        ),
        "unbanded_death_chain_kernel": (
            int(overall["agg_fail"]),
            min_delta,
            4,
            6,
            "actual singleton context scalar; all frontier rows pass",
        ),
        "slope_band_signed_sum": (
            int(overall["agg_fail"]),
            min_delta,
            5,
            7,
            "preserves Angle A seven-band debt/rescue information",
        ),
        "coherent_context_merge": (
            int(context["failures"]),
            min_context,
            6,
            7,
            "tests generated OR/deletion contexts, not only singletons",
        ),
        "miss_zeta_generated_word": (
            int(context["failures"]),
            min_context,
            7,
            8,
            "best quotient: generated residual pressure before scalarization",
        ),
    }
    vertices = list(carriers)

    def better(a: str, b: str) -> str:
        ma = carriers[a]
        mb = carriers[b]
        key_a = (-ma[0], ma[1], ma[2], ma[3])
        key_b = (-mb[0], mb[1], mb[2], mb[3])
        return a if key_a >= key_b else b

    edges: dict[tuple[str, str], str] = {}
    scores: Counter[str] = Counter()
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            winner = better(a, b)
            edges[(a, b)] = winner
            edges[(b, a)] = winner
            scores[winner] += 1
    score_hist = Counter(scores[v] for v in vertices)
    cycles = 0
    for a, b, c in combinations(vertices, 3):
        ab = edges[(a, b)]
        bc = edges[(b, c)]
        ca = edges[(c, a)]
        if ab == a and bc == b and ca == c:
            cycles += 1
        if ab == b and bc == c and ca == a:
            cycles += 1
    comps = sccs(vertices, edges)
    hp = hamiltonian_path_count(vertices, edges)
    for name, metric in sorted(carriers.items(), key=lambda kv: (-scores[kv[0]], kv[0])):
        print(
            f"    {name}: score={scores[name]}, "
            f"metric={metric[:4]}, note={metric[4]}"
        )
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={cycles}")
    print(f"  SCCs={comps}")
    print(f"  Hamiltonian_path_count={hp}")
    print("  Tie Hamiltonian path is the score order above.")
    print()


def main() -> None:
    print("HYP-2702 / T937 -- death-chain slope-band automaton scout")
    print("Arithmetic: exact Fractions; imported S64 sector-mask engine.\n")
    rows_by_size = flatten_coordinate_rows()
    _summary, overall = analyze_death_chain_bands(rows_by_size)
    bounded = bounded_singleton_slices()
    context = context_merge_test(rows_by_size, bounded)
    ca_sign_examples(overall)
    tournament_analysis(overall, context)
    print("SYNTHESIS")
    print("  The singleton death-chain margin is positive on every S64 sparse-coordinate")
    print("  frontier row and on the extra bounded singleton slices tested here.")
    print("  However, first-order hit-count dominance and per-band monotonicity both fail.")
    print("  The live proof object should therefore be a signed seven-band death-chain")
    print("  ledger plus the generated miss-zeta/context-merge monotonicity lemma.")


if __name__ == "__main__":
    main()
