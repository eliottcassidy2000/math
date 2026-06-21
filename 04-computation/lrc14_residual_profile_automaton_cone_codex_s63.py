#!/usr/bin/env python3
"""HYP-2698 / T934 / S63: generated residual-profile cone scout.

HYP-2697 says arbitrary nonnegative residual weights are too strong: the
consecutive block does not dominate every residual coordinate.  The actual
decorrelated LRC context is more rigid.  At fixed slow coordinate x, each
independent cluster has a law on 64 hit masks, and adding a cluster updates the
context by OR-convolution.  Equivalently, residual masks only delete sectors.

This script keeps that generated object:

    x -> w_x(R),  R subset {1,...,6}.

It enumerates coherent-block contexts from integer partitions, records the
integrated residual profiles they generate, and tests the S62 coordinatewise
counterexamples against the actual x-correlated context words.

The broad claim remains open; this is an exact finite scout, not a proof.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd


FULL_MASK = (1 << 6) - 1


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def frac_range(pair: tuple[F, F]) -> str:
    return f"({pair[0]}, {pair[1]})"


def mask_size(mask: int) -> int:
    return mask.bit_count()


def mask_to_sectors(mask: int) -> tuple[int, ...]:
    return tuple(i + 1 for i in range(6) if mask & (1 << i))


def primitive(shape: tuple[int, ...]) -> bool:
    g = 0
    for x in shape:
        g = gcd(g, x)
    return g == 1


def bounded_shapes(size: int, span: int) -> list[tuple[int, ...]]:
    if size == 1:
        return [(0,)]
    out: list[tuple[int, ...]] = []
    for rest in combinations(range(1, span + 1), size - 1):
        shape = (0,) + rest
        if primitive(shape):
            out.append(shape)
    return out


def partitions(n: int, max_part: int | None = None) -> list[tuple[int, ...]]:
    if max_part is None:
        max_part = n
    if n == 0:
        return [()]
    out: list[tuple[int, ...]] = []
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            out.append((first,) + rest)
    return out


def coherent_context(part: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(range(s)) for s in part)


def context_name(part: tuple[int, ...]) -> str:
    return "[" + "+".join(str(x) for x in part) + "]"


@lru_cache(maxsize=None)
def breakpoints(shape: tuple[int, ...]) -> tuple[F, ...]:
    bps = {F(0), F(1)}
    for d in shape:
        if d == 0:
            continue
        for a in range(7 * d + 1):
            bps.add(F(a, 7 * d))
    return tuple(sorted(bps))


def combined_breakpoints(clusters: tuple[tuple[int, ...], ...]) -> tuple[F, ...]:
    bps = {F(0), F(1)}
    for shape in clusters:
        bps.update(breakpoints(shape))
    return tuple(sorted(bps))


@lru_cache(maxsize=None)
def hit_law_at_x(shape: tuple[int, ...], x: F) -> tuple[tuple[int, F], ...]:
    """Law of covered inner-sector mask for one random-phase carrier at fixed x."""
    base = tuple((d * x) % 1 for d in shape)
    cuts = {F(0), F(1)}
    for b in base:
        for s in range(7):
            cuts.add((F(s, 7) - b) % 1)
    cuts = sorted(cuts)
    dist: dict[int, F] = defaultdict(F)
    for lo, hi in zip(cuts, cuts[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        mask = 0
        for b in base:
            sec = int(((b + mid) % 1) * 7)
            if 1 <= sec <= 6:
                mask |= 1 << (sec - 1)
        dist[mask] += hi - lo
    return tuple(sorted(dist.items()))


def or_convolve(a: tuple[tuple[int, F], ...], b: tuple[tuple[int, F], ...]) -> tuple[tuple[int, F], ...]:
    out: dict[int, F] = defaultdict(F)
    for ma, wa in a:
        for mb, wb in b:
            out[ma | mb] += wa * wb
    return tuple(sorted(out.items()))


@lru_cache(maxsize=None)
def hit_union_law_at_x(context: tuple[tuple[int, ...], ...], x: F) -> tuple[tuple[int, F], ...]:
    cur: tuple[tuple[int, F], ...] = ((0, F(1)),)
    for shape in context:
        cur = or_convolve(cur, hit_law_at_x(shape, x))
    return cur


@lru_cache(maxsize=None)
def residual_law_at_x(context: tuple[tuple[int, ...], ...], x: F) -> tuple[tuple[int, F], ...]:
    out: dict[int, F] = defaultdict(F)
    for hit, mass in hit_union_law_at_x(context, x):
        out[FULL_MASK ^ hit] += mass
    return tuple(sorted(out.items()))


def zeta_capacity(hit_law: tuple[tuple[int, F], ...]) -> tuple[F, ...]:
    cap = [F(0) for _ in range(64)]
    for hit, mass in hit_law:
        sub = hit
        while True:
            cap[sub] += mass
            if sub == 0:
                break
            sub = (sub - 1) & hit
    return tuple(cap)


@lru_cache(maxsize=None)
def capacity_at_x(shape: tuple[int, ...], x: F) -> tuple[F, ...]:
    return zeta_capacity(hit_law_at_x(shape, x))


@lru_cache(maxsize=None)
def global_capacity(shape: tuple[int, ...]) -> tuple[F, ...]:
    total = [F(0) for _ in range(64)]
    xs = breakpoints(shape)
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        cap = capacity_at_x(shape, mid)
        for mask, value in enumerate(cap):
            total[mask] += (hi - lo) * value
    return tuple(total)


@lru_cache(maxsize=None)
def integrated_residual_profile(context: tuple[tuple[int, ...], ...]) -> tuple[F, ...]:
    total = [F(0) for _ in range(64)]
    xs = combined_breakpoints(context)
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        for residual, mass in residual_law_at_x(context, mid):
            total[residual] += (hi - lo) * mass
    assert sum(total, F(0)) == 1
    return tuple(total)


def miss_zeta(profile: tuple[F, ...]) -> tuple[F, ...]:
    """Zeta transform: z[A]=Pr(A is contained in the residual mask)."""
    z = [F(0) for _ in range(64)]
    for residual, mass in enumerate(profile):
        sub = residual
        while True:
            z[sub] += mass
            if sub == 0:
                break
            sub = (sub - 1) & residual
    return tuple(z)


def miss_zeta_at_x(context: tuple[tuple[int, ...], ...], x: F) -> tuple[F, ...]:
    profile = [F(0) for _ in range(64)]
    for residual, mass in residual_law_at_x(context, x):
        profile[residual] = mass
    return miss_zeta(tuple(profile))


def cluster_miss_zeta_at_x(shape: tuple[int, ...], x: F) -> tuple[F, ...]:
    profile = [F(0) for _ in range(64)]
    for hit, mass in hit_law_at_x(shape, x):
        profile[FULL_MASK ^ hit] += mass
    return miss_zeta(tuple(profile))


@lru_cache(maxsize=None)
def weighted_delta(
    context: tuple[tuple[int, ...], ...],
    consec: tuple[int, ...],
    challenger: tuple[int, ...],
) -> F:
    """Exact pcover(context+consec)-pcover(context+challenger)."""
    xs = combined_breakpoints(context + (consec, challenger))
    total = F(0)
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        cap_k = capacity_at_x(consec, mid)
        cap_c = capacity_at_x(challenger, mid)
        step = F(0)
        for residual, mass in residual_law_at_x(context, mid):
            step += mass * (cap_k[residual] - cap_c[residual])
        total += (hi - lo) * step
    return total


@lru_cache(maxsize=None)
def raw_averaged_dot(
    context: tuple[tuple[int, ...], ...],
    consec: tuple[int, ...],
    challenger: tuple[int, ...],
) -> F:
    """The tempting but too-coarse averaged profile dot averaged capacity gap."""
    profile = integrated_residual_profile(context)
    cap_k = global_capacity(consec)
    cap_c = global_capacity(challenger)
    return sum((profile[r] * (cap_k[r] - cap_c[r]) for r in range(64)), F(0))


def layer_masses(profile: tuple[F, ...]) -> tuple[F, ...]:
    layers = [F(0) for _ in range(7)]
    for mask, mass in enumerate(profile):
        layers[mask_size(mask)] += mass
    return tuple(layers)


def profile_summary() -> None:
    print("PART A -- coherent contexts generate a small family of residual profiles")
    print("  context size r means r nonzero runners outside the candidate cluster")
    print("  contexts are coherent-block partitions, as in THM-557")
    print(
        f"{'r':>2} {'#ctx':>5} {'support range':>15} {'dominant residual layers':>40} "
        f"{'one-missed range':>28} {'full-missed range':>28}"
    )
    for r in range(0, 9):
        rows = []
        for part in partitions(r):
            context = coherent_context(part)
            profile = integrated_residual_profile(context)
            support = sum(1 for x in profile if x)
            layers = layer_masses(profile)
            dom = max(range(7), key=lambda i: layers[i])
            rows.append((part, support, layers, profile))
        support_range = (min(x[1] for x in rows), max(x[1] for x in rows))
        one_range = (min(x[2][1] for x in rows), max(x[2][1] for x in rows))
        full_range = (min(x[2][6] for x in rows), max(x[2][6] for x in rows))
        dom_hist = Counter(max(range(7), key=lambda i: x[2][i]) for x in rows)
        print(
            f"{r:>2} {len(rows):>5} {str(support_range):>15} "
            f"{dict(sorted(dom_hist.items()))!s:>40} "
            f"{frac_range(one_range):>28} {frac_range(full_range):>28}"
        )
    print()


def factorization_check() -> None:
    print("PART B -- pointwise miss-zeta coordinates factor multiplicatively")
    context = coherent_context((4, 2, 1))
    x = F(5, 37)
    z_context = miss_zeta_at_x(context, x)
    product = [F(1) for _ in range(64)]
    for shape in context:
        z_shape = cluster_miss_zeta_at_x(shape, x)
        for mask in range(64):
            product[mask] *= z_shape[mask]
    assert tuple(product) == z_context
    print(f"  checked context {tuple(len(s) for s in context)} at x={x}")
    for mask in (FULL_MASK, 0b000001, 0b000011, 0b010101, 0b111000):
        print(
            f"    Pr({mask_to_sectors(mask)} all still missed) = "
            f"{z_context[mask]}"
        )
    print("  This is the exact cone coordinate: z_x(A)=prod_i z_{i,x}(A).")
    print("  Averaging over x happens only after this product word is formed.")
    print()


def named_challenger_tests() -> None:
    print("PART C -- S62 coordinatewise counterexamples against generated contexts")
    pairs = [
        ((0, 1, 2), (0, 1, 3)),
        ((0, 1, 2, 3), (0, 1, 2, 4)),
        ((0, 1, 2, 3, 4), (0, 1, 2, 3, 5)),
    ]
    for consec, challenger in pairs:
        cap_k = global_capacity(consec)
        cap_c = global_capacity(challenger)
        bad_masks = tuple(r for r in range(1, 64) if cap_c[r] > cap_k[r])
        print(f"  candidate size={len(consec)} K={consec} C={challenger}")
        print(
            f"    bad residual coordinates={len(bad_masks)}, "
            f"max_bad_gain={max((cap_c[r] - cap_k[r] for r in bad_masks), default=F(0))}"
        )
        for total in range(max(7, len(consec)), 12):
            r = total - len(consec)
            rows = []
            for part in partitions(r):
                context = coherent_context(part)
                actual = weighted_delta(context, consec, challenger)
                raw = raw_averaged_dot(context, consec, challenger)
                profile = integrated_residual_profile(context)
                bad_mass = sum((profile[m] for m in bad_masks), F(0))
                rows.append((actual, raw, bad_mass, part))
            actual_min = min(rows, key=lambda x: (x[0], x[3]))
            raw_min = min(rows, key=lambda x: (x[1], x[3]))
            bad_max = max(rows, key=lambda x: (x[2], x[3]))
            failures = sum(1 for row in rows if row[0] < 0)
            raw_failures = sum(1 for row in rows if row[1] < 0)
            print(
                f"    total_m={total}, context_r={r}, contexts={len(rows):2d}: "
                f"min actual={fmt(actual_min[0])} at {context_name(actual_min[3])}; "
                f"min averaged-dot={fmt(raw_min[1])} at {context_name(raw_min[3])}; "
                f"bad-mass max={fmt(bad_max[2])} at {context_name(bad_max[3])}; "
                f"actual_fail={failures}, averaged_fail={raw_failures}"
            )
        print()


def frontier_challengers(size: int) -> list[tuple[int, ...]]:
    """Near-consecutive challengers that move only a few cells.

    The full span<=8 exact scan is possible but slow because it refines many
    shared-x breakpoint products.  These are the local frontier shapes most
    relevant after the S62 one-gap counterexamples.
    """
    consec = tuple(range(size))
    out = []
    for shape in bounded_shapes(size, size + 2):
        if shape == consec:
            continue
        l1_move = sum(abs(shape[i] - i) for i in range(size))
        max_move = max(abs(shape[i] - i) for i in range(size))
        if l1_move <= 4 and max_move <= 2:
            out.append(shape)
    return out


def bounded_challenger_scan() -> None:
    print("PART D -- exact frontier challenger scan against coherent contexts")
    print("  Shapes of size 3..6 near the consecutive block, total_m=7..11.")
    print("  A negative delta would be a real obstruction to the coherent-context cone.")
    for size in (3, 4, 5, 6):
        consec = tuple(range(size))
        worst = (F(10**9), None, None, None)
        failures = 0
        tested = 0
        challengers = frontier_challengers(size)
        for challenger in challengers:
            for total in range(max(7, size), 12):
                r = total - size
                for part in partitions(r):
                    context = coherent_context(part)
                    delta = weighted_delta(context, consec, challenger)
                    tested += 1
                    if delta < 0:
                        failures += 1
                    if delta < worst[0]:
                        worst = (delta, challenger, total, part)
        print(
            f"  size={size}: challengers={len(challengers)}, tested={tested}, failures={failures}, "
            f"worst_delta={fmt(worst[0])} challenger={worst[1]} "
            f"total_m={worst[2]} context={context_name(worst[3])}"
        )
    print()


def tournament_stats(names: list[str], score) -> None:
    edges = {name: set() for name in names}
    for a, b in combinations(names, 2):
        sa, sb = score(a), score(b)
        if sa > sb:
            edges[a].add(b)
        elif sb > sa:
            edges[b].add(a)
        else:
            edges[a].add(b)
    scores = {name: len(edges[name]) for name in names}
    hist = Counter(scores.values())
    cycles = 0
    for a, b, c in combinations(names, 3):
        for x, y, z in ((a, b, c), (a, c, b)):
            if y in edges[x] and z in edges[y] and x in edges[z]:
                cycles += 1
                break
    hpaths = 0
    for perm in permutations(names):
        if all(perm[i + 1] in edges[perm[i]] for i in range(len(perm) - 1)):
            hpaths += 1
    print(f"    scores={scores}")
    print(f"    score_hist={dict(sorted(hist.items()))} directed_3cycles={cycles}")
    print(f"    Hamiltonian_path_count={hpaths}")


def tournament_analysis() -> None:
    print("PART E -- Tournament Analysis: proof carriers after the automaton scout")
    lanes = {
        "raw_coordinate_weights": (0, 5, 1, 0),
        "averaged_residual_profile": (2, 4, 2, 1),
        "x_correlated_profile_word": (5, 4, 4, 4),
        "miss_zeta_product_word": (5, 5, 5, 5),
        "coherent_context_exact_scan": (4, 3, 4, 3),
        "transfer_tax_label_lift": (4, 4, 5, 4),
        "literal_rule110_analogy": (1, 1, 1, 5),
    }
    print("  vertices: proof-carrier quotients, not runners")
    print("  observable=(exactness, generated-cone fidelity, proof reach, analogy)")
    tournament_stats(list(lanes), lambda n: lanes[n])
    print(
        "  Read: Rule 110 is useful only as a generated-language warning; the live\n"
        "  carrier is the miss-zeta product word plus transfer-tax labels."
    )
    print()


def main() -> None:
    print("HYP-2698 / T934 / S63 -- residual-profile automaton cone scout")
    print("Arithmetic: exact Fractions over shared-x sector-mask breakpoints.")
    print("Rule-110 lens: local finite update generates a constrained language;")
    print("do not replace the generated LRC context word by arbitrary weights.\n")
    profile_summary()
    factorization_check()
    named_challenger_tests()
    bounded_challenger_scan()
    tournament_analysis()
    print("SYNTHESIS")
    print("  Actual contexts are generated by a 64-state OR/deletion automaton.")
    print("  In miss-zeta coordinates, the pointwise context word factors as a")
    print("  product over clusters before x-averaging.")
    print("  The S62 residual-coordinate counterexamples do not produce failures")
    print("  against the coherent generated contexts scanned here; the x-correlated")
    print("  word is stronger than the averaged profile dot averaged capacity gap.")
    print("  No LRC14 proof is claimed; HYP-2698 sharpens the cone coordinates.")


if __name__ == "__main__":
    main()
