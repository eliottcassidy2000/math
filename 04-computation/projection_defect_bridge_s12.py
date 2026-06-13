#!/usr/bin/env python3
"""
opus-2026-05-29-S12

Projection-defect audit across three forgetful maps:

1. directed odd cycles -> vertex supports;
2. deletion / old-coordinate projection at a vertex;
3. tournament orientations -> even graph cycle-space projection (odd n only).

This is a bridge script, not a speed benchmark.  It looks for the same kind of
residue appearing in the complete-Omega, real-root failure, and even-graph
threads.
"""

from __future__ import annotations

import itertools
import sys
from collections import Counter, defaultdict

sys.path.insert(0, "04-computation")

from omega_extreme_fingerprints_s11 import (  # noqa: E402
    alpha_vector,
    delete_vertex,
    gentourng_reps,
    hamiltonian_paths,
    odd_cycle_masks,
    rows_from_nauty_bits,
    score_sequence,
    single_core_signature,
    thm025_rows,
)


def circulant_rows(p: int, connection: set[int]) -> tuple[int, ...]:
    rows = [0] * p
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in connection:
                rows[i] |= 1 << j
    return tuple(rows)


def transitive_rows(n: int) -> tuple[int, ...]:
    return tuple(sum(1 << j for j in range(v + 1, n)) for v in range(n))


def support_tuple(mask: int, n: int) -> tuple[int, ...]:
    return tuple(v for v in range(n) if (mask >> v) & 1)


def odd_cycle_support_stats(rows: tuple[int, ...], include_alpha: bool = True) -> dict[str, object]:
    cycles = odd_cycle_masks(rows)
    masks = [mask for _, mask in cycles]
    by_support = Counter(masks)
    core = (1 << len(rows)) - 1
    for mask in masks:
        core &= mask
    if not masks:
        core = 0
    disjoint_pairs = 0
    complete = True
    for a, b in itertools.combinations(masks, 2):
        if not (a & b):
            disjoint_pairs += 1
            complete = False
            if not include_alpha:
                break
    return {
        "cycles": cycles,
        "masks": masks,
        "alpha": alpha_vector(masks) if include_alpha else None,
        "support_count": len(by_support),
        "support_excess": len(masks) - len(by_support),
        "max_support_mult": max(by_support.values(), default=0),
        "multi_supports": sum(1 for x in by_support.values() if x > 1),
        "core": core,
        "complete": complete,
        "disjoint_pairs": disjoint_pairs,
    }


def deletion_projection_profile(rows: tuple[int, ...]) -> list[dict[str, object]]:
    n = len(rows)
    stats = odd_cycle_support_stats(rows)
    masks = stats["masks"]
    out = []
    for v in range(n):
        old_masks = [mask for mask in masks if not ((mask >> v) & 1)]
        sub = delete_vertex(rows, v)
        out.append(
            {
                "v": v,
                "deg": rows[v].bit_count(),
                "lost": len(masks) - len(old_masks),
                "kept": len(old_masks),
                "loss_frac": (len(masks) - len(old_masks)) / max(len(masks), 1),
                "H_delete": hamiltonian_paths(sub),
                "alpha_delete": alpha_vector([mask for _, mask in odd_cycle_masks(sub)]),
                "score_delete": score_sequence(sub),
            }
        )
    return out


def edge_pairs(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def orientation_bits(rows: tuple[int, ...]) -> list[int]:
    bits = []
    for i, j in edge_pairs(len(rows)):
        bits.append(1 if ((rows[i] >> j) & 1) else 0)
    return bits


def cycle_projection_even_graph(rows: tuple[int, ...]) -> dict[str, object] | None:
    n = len(rows)
    if n % 2 == 0:
        return None
    pairs = edge_pairs(n)
    bits = orientation_bits(rows)
    projected = []
    for e_idx, (a, b) in enumerate(pairs):
        val = bits[e_idx]
        for f_idx, (c, d) in enumerate(pairs):
            if e_idx != f_idx and (a == c or a == d or b == c or b == d):
                val ^= bits[f_idx]
        projected.append(val)
    deg = [0] * n
    for bit, (a, b) in zip(projected, pairs):
        if bit:
            deg[a] += 1
            deg[b] += 1
    return {
        "edges": sum(projected),
        "degrees": tuple(sorted(deg)),
        "degree_even": all(d % 2 == 0 for d in deg),
        "bits": tuple(projected),
    }


def summarize_example(label: str, rows: tuple[int, ...]) -> str:
    n = len(rows)
    stats = odd_cycle_support_stats(rows)
    core = stats["core"]
    deletion = deletion_projection_profile(rows)
    worst = sorted(deletion, key=lambda r: (-r["loss_frac"], r["v"]))[:3]
    even = cycle_projection_even_graph(rows)
    lines = [label, "-" * len(label)]
    lines.append(f"n={n}, H={hamiltonian_paths(rows)}, scores={score_sequence(rows)}")
    lines.append(
        "Omega: cycles={cycles}, supports={supports}, support_excess={excess}, "
        "max_mult={max_mult}, alpha={alpha}, complete={complete}, disjoint_pairs={pairs}".format(
            cycles=len(stats["masks"]),
            supports=stats["support_count"],
            excess=stats["support_excess"],
            max_mult=stats["max_support_mult"],
            alpha=stats["alpha"],
            complete=stats["complete"],
            pairs=stats["disjoint_pairs"],
        )
    )
    lines.append(f"cycle-family core={support_tuple(core, n) if core else ()}")
    if core.bit_count() == 1:
        v = (core & -core).bit_length() - 1
        sig = single_core_signature(rows, v)
        if sig is not None:
            lines.append(f"single-core projection: v={v}, signature={sig[0]}, r_core={sig[1]}")
    lines.append("largest deletion/old-projection losses:")
    for row in worst:
        lines.append(
            "  v={v}, deg={deg}, lost={lost}, kept={kept}, loss={loss_frac:.3f}, "
            "H(T-v)={H_delete}, alpha(T-v)={alpha_delete}".format(**row)
        )
    killers = [row["v"] for row in deletion if row["kept"] == 0]
    lines.append(f"old-projection kill vertices={killers}")
    if even is None:
        lines.append("cycle-space projection: ambiguous at even n (cut ∩ cycle nonzero)")
    else:
        lines.append(
            f"cycle-space projection: even_edges={even['edges']}, "
            f"even_degree_seq={even['degrees']}, degree_even={even['degree_even']}"
        )
    return "\n".join(lines)


def complete_omega_core_census(max_n: int = 8) -> str:
    lines = ["Complete-Omega core census", "=" * 28]
    for n in range(3, max_n + 1):
        reps = gentourng_reps(n)
        core_size_counts: Counter[int] = Counter()
        support_by_core: dict[int, Counter[int]] = defaultdict(Counter)
        target_rows = []
        for cid, bits in enumerate(reps):
            rows = rows_from_nauty_bits(n, bits)
            stats = odd_cycle_support_stats(rows, include_alpha=False)
            if not stats["complete"]:
                continue
            r = len(stats["masks"])
            csize = int(stats["core"].bit_count())
            core_size_counts[csize] += 1
            support_by_core[csize][r] += 1
            if r in {3, 10, 31}:
                target_rows.append((r, cid, csize, hamiltonian_paths(rows), score_sequence(rows)))
        lines.append(f"n={n}: complete={sum(core_size_counts.values())}/{len(reps)}, core sizes={dict(sorted(core_size_counts.items()))}")
        for csize in sorted(support_by_core):
            support = sorted(support_by_core[csize])
            lines.append(f"  core_size={csize}: r_support={support}")
        if target_rows:
            lines.append("  target r rows:")
            for row in target_rows:
                lines.append(f"    r={row[0]}, cid={row[1]}, core_size={row[2]}, H={row[3]}, scores={row[4]}")
        else:
            lines.append("  target r rows: none for r in {3,10,31}")
    return "\n".join(lines)


def target_signature_search(max_m: int = 40) -> str:
    targets = (3, 10, 31, 42, 63)
    max_target = max(targets)
    lines = ["Single-core target persistence", "=" * 30]

    # State = (current r, contribution if next bit is 0, last bit).
    # Since r never decreases, pruning at max_target preserves exact counts
    # for the target values while avoiding a 2^m brute-force pass.
    states: dict[tuple[int, int, int], int] = {(0, 0, 0): 1}
    state_examples: dict[tuple[int, int, int], list[str]] = {(0, 0, 0): [""]}

    for m in range(1, max_m + 1):
        new_states: dict[tuple[int, int, int], int] = defaultdict(int)
        new_examples: dict[tuple[int, int, int], list[str]] = defaultdict(list)
        for (r, add0, last), count in states.items():
            old_next_add = 2 * add0 - last

            r0 = r + add0
            if r0 <= max_target:
                key0 = (r0, old_next_add, 0)
                new_states[key0] += count
                for s in state_examples[(r, add0, last)]:
                    if len(new_examples[key0]) < 3:
                        new_examples[key0].append(s + "0")

            key1 = (r, old_next_add + 1, 1)
            new_states[key1] += count
            for s in state_examples[(r, add0, last)]:
                if len(new_examples[key1]) < 3:
                    new_examples[key1].append(s + "1")

        states = dict(new_states)
        state_examples = dict(new_examples)
        if m < 2:
            continue

        counts = Counter()
        examples: dict[int, list[str]] = {t: [] for t in targets}
        for state, count in states.items():
            r = state[0]
            if r in targets:
                counts[r] += count
                for s in state_examples[state]:
                    if len(examples[r]) < 3:
                        examples[r].append(s)
        summary = ", ".join(
            f"r={t}:{counts[t] if counts[t] else 'absent'}"
            + (f" ex={examples[t]}" if examples[t] else "")
            for t in targets
        )
        lines.append(f"m={m}: {summary}")
    return "\n".join(lines)


def main() -> None:
    print("Projection-defect bridge audit (opus-2026-05-29-S12)")
    print("=" * 72)
    print()
    # Include n=8 because the complete-Omega core strata are the smallest
    # place where the H=63 projection-kill classes appear.
    print(complete_omega_core_census(max_n=8))
    print()
    print(target_signature_search())
    print()

    n8_reps = gentourng_reps(8)
    examples = [
        ("transitive n=7", transitive_rows(7)),
        ("Paley T7 / QR circulant", circulant_rows(7, {1, 2, 4})),
        ("interval T7 / cyclic order", circulant_rows(7, {1, 2, 3})),
        ("THM-344 H=63 cid=2519", rows_from_nauty_bits(8, n8_reps[2519])),
        ("THM-344 H=63 cid=3285", rows_from_nauty_bits(8, n8_reps[3285])),
        ("THM-025 n=9 real-rootedness counterexample", thm025_rows()),
    ]
    for label, rows in examples:
        print(summarize_example(label, rows))
        print()


if __name__ == "__main__":
    main()
