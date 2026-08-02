#!/usr/bin/env python3
"""Exact scout for the upper-median boundary star and chord repairs.

This is a research scout, not by itself a theorem promotion.  It works on the
649 bodies remaining after the robust-edge-11 block and checks every injective
assignment from ``m,...,m+D`` at the deterministic upper-median safe cell.
The boundary-star centre is the unique label ``e_*`` satisfying

    e_* j == L/14 (mod L).

For every word missed by its five incident edges, the program checks all ten
non-star chords and computes a minimum chord cover of the missed words.
"""

from __future__ import annotations

import argparse
import hashlib
import random
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation/lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
LOW = ROOT / "04-computation/lrc14_j7_reflected_low_phase_clique_robust_body_closure_thm2941.py"
EDGES = tuple(combinations(range(6), 2))


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def load(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = load("median_star_base", BASE)
LP = load("median_star_low", LOW)


def intersection_mass(first, second) -> F:
    i = j = 0
    total = F(0)
    while i < len(first) and j < len(second):
        total += max(
            F(0),
            min(first[i][1], second[j][1])
            - max(first[i][0], second[j][0]),
        )
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return total


def minimum_covers(masks: dict[tuple[int, int], int], full: int):
    candidates = tuple(edge for edge, mask in masks.items() if mask)
    for size in range(1, len(candidates) + 1):
        rows = []
        for chosen in combinations(candidates, size):
            mask = 0
            for edge in chosen:
                mask |= masks[edge]
            if mask == full:
                rows.append(chosen)
        if rows:
            return tuple(rows)
    return ()


def two_star_mst_credit(
    weights: dict[tuple[int, int], F], centre: int, second: int
) -> F:
    """Maximum spanning-tree weight in ``K_{2,4}`` plus the centre edge."""
    centre_edge = tuple(sorted((centre, second)))
    leaves = tuple(index for index in range(6) if index not in (centre, second))
    arms = {
        leaf: (
            weights[tuple(sorted((centre, leaf)))],
            weights[tuple(sorted((second, leaf)))],
        )
        for leaf in leaves
    }
    through_centre_edge = weights[centre_edge] + sum(
        (max(arms[leaf]) for leaf in leaves), F(0)
    )
    through_leaf = max(
        arms[bridge][0]
        + arms[bridge][1]
        + sum((max(arms[leaf]) for leaf in leaves if leaf != bridge), F(0))
        for bridge in leaves
    )
    return max(through_centre_edge, through_leaf)


def maximum_spanning_tree_credit(
    weights: dict[tuple[int, int], F], edges: tuple[tuple[int, int], ...]
) -> F | None:
    """Kruskal credit, or ``None`` when the supplied intrinsic graph disconnects."""
    parent = list(range(6))

    def root(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    total = F(0)
    used = 0
    for edge in sorted(edges, key=lambda item: (weights[item], item), reverse=True):
        left, right = map(root, edge)
        if left == right:
            continue
        parent[left] = right
        total += weights[edge]
        used += 1
        if used == 5:
            return total
    return None


def body_universe():
    rows = []
    for body in combinations(range(1, 15), 6):
        ruler, _, robust = LP.robust_edges(body)
        if len(robust) <= 10:
            rows.append((body, ruler))
    require(len(rows) == 649, len(rows))
    return tuple(rows)


def scan(
    m: int,
    spread: int,
    inherited: dict[tuple[int, ...], tuple[int, int]] | None = None,
    progress: bool = False,
):
    require(m >= 1 and spread >= 5, (m, spread))
    levels = tuple(range(m, m + spread + 1))
    words = tuple(permutations(levels, 6))
    star_fail_count = 0
    star_fail_bodies = 0
    star_tree_fail_count = 0
    star_tree_fail_bodies = set()
    weakest_star_tree = None
    mst_fail_count = 0
    mst_fail_bodies = set()
    weakest_mst = None
    high_tree_words = 0
    high_tree_fail_count = 0
    weakest_high_tree = None
    two_star_fail_count = 0
    two_star_fail_bodies = 0
    inherited_fail_count = 0
    cover_number_hist = Counter()
    two_star_cover_number_hist = Counter()
    centre_hist = Counter()
    chord_index_hist = Counter()
    chord_label_hist = Counter()
    residual_digest = hashlib.sha256()
    body_rows = []
    chosen_map = {}

    for body_index, (body, L) in enumerate(body_universe()):
        ruler, ranges = R.safe_cell_ranges(body)
        require(ruler == L, (body, ruler, L))
        cells = tuple(j for left, right in ranges for j in range(left, right))
        require(cells, ("empty safe carrier", body))
        cell = cells[len(cells) // 2]
        centres = tuple(
            index
            for index, label in enumerate(body)
            if (label * cell) % L == L // 14
        )
        require(len(centres) == 1, (body, L, cell, centres))
        centre = centres[0]
        centre_hist[(centre, body[centre])] += 1
        second_centre = next(index for index in range(6) if index != centre)
        star = tuple(edge for edge in EDGES if centre in edge)
        chords = tuple(edge for edge in EDGES if centre not in edge)
        second_star_chords = tuple(edge for edge in chords if second_centre in edge)
        two_star = tuple(
            edge for edge in EDGES if centre in edge or second_centre in edge
        )
        require(len(second_star_chords) == 4, (body, centre, second_centre))

        debts = {
            (i, q): F(body[i], 7 * (q * L - body[i]))
            for i in range(6)
            for q in levels
        }
        arcs = {
            (i, q): R.reflected_level_arcs(L, body[i], q, cell)
            for i in range(6)
            for q in levels
        }
        profiles = {
            (i, j, p, q): intersection_mass(arcs[i, p], arcs[j, q])
            for i, j in EDGES
            for p in levels
            for q in levels
            if p != q
        }

        failures = []
        for word in words:
            debt = sum((debts[i, word[i]] for i in range(6)), F(0))
            star_weights = tuple(
                profiles[i, j, word[i], word[j]] for i, j in star
            )
            two_star_weights = {
                (i, j): profiles[i, j, word[i], word[j]] for i, j in two_star
            }
            mst_credit = two_star_mst_credit(
                two_star_weights, centre, second_centre
            )
            mst_margin = mst_credit - debt
            mst_row = (
                mst_margin,
                body,
                L,
                cell,
                centre,
                second_centre,
                word,
                mst_credit,
                debt,
            )
            if weakest_mst is None or mst_row < weakest_mst:
                weakest_mst = mst_row
            if mst_margin <= 0:
                mst_fail_count += 1
                mst_fail_bodies.add(body)
            high_edges = tuple(
                edge
                for edge in EDGES
                if sum(
                    value // gcd(word[edge[0]], word[edge[1]])
                    for value in (word[edge[0]], word[edge[1]])
                )
                >= 8
            )
            high_weights = {
                edge: profiles[edge[0], edge[1], word[edge[0]], word[edge[1]]]
                for edge in high_edges
            }
            high_credit = maximum_spanning_tree_credit(high_weights, high_edges)
            if high_credit is not None:
                high_tree_words += 1
                high_margin = high_credit - debt
                high_row = (
                    high_margin,
                    body,
                    L,
                    cell,
                    word,
                    high_edges,
                    high_credit,
                    debt,
                )
                if weakest_high_tree is None or high_row < weakest_high_tree:
                    weakest_high_tree = high_row
                if high_margin <= 0:
                    high_tree_fail_count += 1
            if max(star_weights) <= debt:
                # A pair-success already forces the star-tree sum above debt,
                # so only pair failures can possibly fail Hunter's star tree.
                tree_margin = sum(star_weights, F(0)) - debt
                tree_row = (
                    tree_margin, body, L, cell, centre, word, debt, star_weights
                )
                if weakest_star_tree is None or tree_row < weakest_star_tree:
                    weakest_star_tree = tree_row
                if tree_margin <= 0:
                    star_tree_fail_count += 1
                    star_tree_fail_bodies.add(body)
                gains = tuple(
                    (profiles[i, j, word[i], word[j]], (i, j))
                    for i, j in chords
                )
                require(
                    max(gain for gain, _ in gains) > debt,
                    ("complete median graph failed", body, m, spread, word, debt),
                )
                failures.append((word, debt, gains))

        if not failures:
            if inherited is not None and body in inherited:
                chosen_map[body] = inherited[body]
            continue

        star_fail_bodies += 1
        star_fail_count += len(failures)
        full = (1 << len(failures)) - 1
        masks = {}
        for edge in chords:
            mask = 0
            for row_index, (_, debt, gains) in enumerate(failures):
                gain = next(value for value, candidate in gains if candidate == edge)
                if gain > debt:
                    mask |= 1 << row_index
            masks[edge] = mask
        covers = minimum_covers(masks, full)
        require(covers, ("no chord cover", body, m, spread))
        cover_number_hist[len(covers[0])] += 1
        second_masks = {edge: masks[edge] for edge in second_star_chords}
        second_union = 0
        for mask in second_masks.values():
            second_union |= mask
        second_failures = len(failures) - second_union.bit_count()
        if second_failures:
            two_star_fail_bodies += 1
            two_star_fail_count += second_failures
        else:
            second_covers = minimum_covers(second_masks, full)
            require(second_covers, ("two-star cover missing", body, m, spread))
            two_star_cover_number_hist[len(second_covers[0])] += 1
        one_chords = tuple(edge for edge in chords if masks[edge] == full)
        if inherited is None:
            chosen = min(one_chords) if one_chords else covers[0][0]
        else:
            chosen = inherited.get(body)
            if chosen is None:
                chosen = min(one_chords) if one_chords else covers[0][0]
            inherited_fail_count += len(failures) - masks.get(chosen, 0).bit_count()
        chosen_map[body] = chosen
        chord_index_hist[(centre, chosen)] += 1
        chord_label_hist[(body[centre], body[chosen[0]], body[chosen[1]])] += 1

        for word, debt, gains in failures:
            winners = tuple(edge for gain, edge in gains if gain > debt)
            residual_digest.update(
                f"{body}|{L}|{cell}|{centre}|{word}|{debt}|{winners}\n".encode()
            )
        body_rows.append(
            (
                body,
                L,
                cell,
                centre,
                body[centre],
                second_centre,
                body[second_centre],
                len(cells),
                len(failures),
                len(covers[0]),
                second_failures,
                one_chords,
                covers[:5],
                chosen,
            )
        )
        if progress and body_index % 50 == 0:
            print("progress", m, spread, body_index, star_fail_count, flush=True)

    summary = {
        "m": m,
        "spread": spread,
        "body_count": 649,
        "words_per_body": len(words),
        "assignment_count": 649 * len(words),
        "star_failures": star_fail_count,
        "star_failure_bodies": star_fail_bodies,
        "star_tree_failures": star_tree_fail_count,
        "star_tree_failure_bodies": len(star_tree_fail_bodies),
        "weakest_star_tree": weakest_star_tree,
        "two_star_mst_failures": mst_fail_count,
        "two_star_mst_failure_bodies": len(mst_fail_bodies),
        "weakest_two_star_mst": weakest_mst,
        "high_tree_words": high_tree_words,
        "high_tree_failures": high_tree_fail_count,
        "weakest_high_tree": weakest_high_tree,
        "two_star_failures": two_star_fail_count,
        "two_star_failure_bodies": two_star_fail_bodies,
        "complete_pair_failures": 0,
        "inherited_chord_failures": inherited_fail_count,
        "cover_number_hist": tuple(sorted(cover_number_hist.items())),
        "two_star_cover_number_hist": tuple(sorted(two_star_cover_number_hist.items())),
        "centre_hist": tuple(sorted(centre_hist.items())),
        "chord_index_hist": tuple(sorted(chord_index_hist.items())),
        "chord_label_hist": tuple(sorted(chord_label_hist.items())),
        "residual_digest": residual_digest.hexdigest(),
    }
    return summary, tuple(body_rows), chosen_map


def sampled_scan(m: int, spread: int, sample_count: int, seed: int = 2941):
    """Deterministic hostile/random probe of the canonical boundary two-star."""
    levels = tuple(range(m, m + spread + 1))
    require(len(levels) >= 6 and sample_count >= 1, (m, spread, sample_count))
    generator = random.Random((seed, m, spread).__repr__())
    words = set()
    structured_subsets = {
        levels[:6],
        levels[-6:],
        (levels[0], levels[1], levels[2], levels[-3], levels[-2], levels[-1]),
        tuple(levels[(index * spread) // 5] for index in range(6)),
    }
    hostile = (2, 3, 4, 6, 8, 12)
    if all(value in levels for value in hostile):
        structured_subsets.add(hostile)
    for subset in structured_subsets:
        require(len(set(subset)) == 6, (m, spread, subset))
        # All permutations of the classical low-phase hostile are important;
        # on generic endpoint/quantile sets a reproducible sample suffices.
        rows = tuple(permutations(subset))
        if subset == hostile:
            words.update(rows)
        else:
            words.update(generator.sample(rows, min(32, len(rows))))
    while len(words) < sample_count:
        words.add(tuple(generator.sample(levels, 6)))
    words = tuple(sorted(words))

    first_failure = None
    weakest = None
    star_failures = 0
    two_star_failures = 0
    full_pair_failures = 0
    mst_failures = 0
    weakest_mst = None
    high_tree_words = 0
    high_tree_failures = 0
    weakest_high_tree = None
    body_checks = 0
    profile_checks = 0
    for body, L in body_universe():
        _, ranges = R.safe_cell_ranges(body)
        cells = tuple(j for left, right in ranges for j in range(left, right))
        cell = cells[len(cells) // 2]
        centres = tuple(
            index
            for index, label in enumerate(body)
            if (label * cell) % L == L // 14
        )
        require(len(centres) == 1, (body, L, cell, centres))
        centre = centres[0]
        second = next(index for index in range(6) if index != centre)
        star = tuple(edge for edge in EDGES if centre in edge)
        two_star = tuple(edge for edge in EDGES if centre in edge or second in edge)
        debt_cache = {
            (i, q): F(body[i], 7 * (q * L - body[i]))
            for i in range(6)
            for q in levels
        }
        arc_cache = {}
        profile_cache = {}

        def profile(i: int, j: int, p: int, q: int) -> F:
            nonlocal profile_checks
            key = (i, j, p, q)
            if key not in profile_cache:
                for index, level in ((i, p), (j, q)):
                    arc_cache.setdefault(
                        (index, level),
                        R.reflected_level_arcs(L, body[index], level, cell),
                    )
                profile_cache[key] = intersection_mass(
                    arc_cache[i, p], arc_cache[j, q]
                )
                profile_checks += 1
            return profile_cache[key]

        for word in words:
            debt = sum((debt_cache[i, word[i]] for i in range(6)), F(0))
            star_gain = max(profile(i, j, word[i], word[j]) for i, j in star)
            if star_gain <= debt:
                star_failures += 1
            gain, edge = max(
                (profile(i, j, word[i], word[j]), (i, j)) for i, j in two_star
            )
            weights = {
                edge0: profile(edge0[0], edge0[1], word[edge0[0]], word[edge0[1]])
                for edge0 in two_star
            }
            mst_credit = two_star_mst_credit(weights, centre, second)
            mst_margin = mst_credit - debt
            mst_row = (
                mst_margin,
                body,
                L,
                cell,
                centre,
                second,
                word,
                mst_credit,
                debt,
            )
            if weakest_mst is None or mst_row < weakest_mst:
                weakest_mst = mst_row
            if mst_margin <= 0:
                mst_failures += 1
            high_edges = tuple(
                edge0
                for edge0 in EDGES
                if sum(
                    value // gcd(word[edge0[0]], word[edge0[1]])
                    for value in (word[edge0[0]], word[edge0[1]])
                )
                >= 8
            )
            high_weights = {
                edge0: profile(
                    edge0[0], edge0[1], word[edge0[0]], word[edge0[1]]
                )
                for edge0 in high_edges
            }
            high_credit = maximum_spanning_tree_credit(high_weights, high_edges)
            if high_credit is not None:
                high_tree_words += 1
                high_margin = high_credit - debt
                high_row = (
                    high_margin,
                    body,
                    L,
                    cell,
                    word,
                    high_edges,
                    high_credit,
                    debt,
                )
                if weakest_high_tree is None or high_row < weakest_high_tree:
                    weakest_high_tree = high_row
                if high_margin <= 0:
                    high_tree_failures += 1
            margin = gain - debt
            row = (margin, body, L, cell, centre, second, word, edge, gain, debt)
            if weakest is None or row < weakest:
                weakest = row
            if margin <= 0:
                two_star_failures += 1
                full_gain, full_edge = max(
                    (profile(i, j, word[i], word[j]), (i, j)) for i, j in EDGES
                )
                if full_gain <= debt:
                    full_pair_failures += 1
                failure = row + (full_edge, full_gain)
                if first_failure is None or failure < first_failure:
                    first_failure = failure
        body_checks += 1

    return {
        "m": m,
        "spread": spread,
        "word_count": len(words),
        "assignment_count": len(words) * 649,
        "body_checks": body_checks,
        "profile_checks": profile_checks,
        "star_failures": star_failures,
        "two_star_failures": two_star_failures,
        "full_pair_failures": full_pair_failures,
        "two_star_mst_failures": mst_failures,
        "weakest_two_star_mst": weakest_mst,
        "high_tree_words": high_tree_words,
        "high_tree_failures": high_tree_failures,
        "weakest_high_tree": weakest_high_tree,
        "weakest_two_star": weakest,
        "first_two_star_failure": first_failure,
    }


def high_channel_audit(max_level: int):
    """Test exact median overlap against universal debt on primitive high channels."""
    channels = tuple(
        (p, q)
        for q in range(2, max_level + 1)
        for p in range(1, q)
        if gcd(p, q) == 1 and p + q >= 8
    )
    require(channels, max_level)
    weakest = None
    failure = None
    checks = 0
    for body, L in body_universe():
        debt = LP.universal_debt(body, L)
        _, ranges = R.safe_cell_ranges(body)
        cells = tuple(j for left, right in ranges for j in range(left, right))
        cell = cells[len(cells) // 2]
        centre = next(
            index
            for index, label in enumerate(body)
            if (label * cell) % L == L // 14
        )
        second = next(index for index in range(6) if index != centre)
        edges = tuple(edge for edge in EDGES if centre in edge or second in edge)
        arc_cache = {}
        for edge in edges:
            i, j = edge
            for p, q in channels:
                for left_level, right_level in ((p, q), (q, p)):
                    for index, level in ((i, left_level), (j, right_level)):
                        arc_cache.setdefault(
                            (index, level),
                            R.reflected_level_arcs(L, body[index], level, cell),
                        )
                    overlap = intersection_mass(
                        arc_cache[i, left_level], arc_cache[j, right_level]
                    )
                    margin = overlap - debt
                    row = (
                        margin,
                        body,
                        L,
                        cell,
                        centre,
                        second,
                        edge,
                        (left_level, right_level),
                        overlap,
                        debt,
                    )
                    if weakest is None or row < weakest:
                        weakest = row
                    if margin <= 0 and (failure is None or row < failure):
                        failure = row
                    checks += 1
    return {
        "max_level": max_level,
        "channel_count": len(channels),
        "body_edge_orientation_checks": checks,
        "weakest": weakest,
        "first_failure": failure,
    }


def projective_low_two_star_rays():
    """All ordered primitive rays on which every two-star edge is low-phase."""
    ratios = frozenset(LP.LOW_RATIOS) | frozenset(
        F(1, 1) / value for value in LP.LOW_RATIOS
    )
    rays = set()
    common_hist = Counter()
    for ratio in ratios:
        common = tuple(
            sorted(
                value
                for value in ratios
                if value not in (F(1), ratio)
                and ratio / value in ratios
            )
        )
        common_hist[(ratio, len(common))] += 1
        for leaves in permutations(common, 4):
            rational = (F(1), ratio, *leaves)
            denominator = lcm(*(value.denominator for value in rational))
            integer = tuple(int(denominator * value) for value in rational)
            divisor = gcd(*integer)
            rays.add(tuple(value // divisor for value in integer))
    ordered = tuple(sorted(rays))
    require(len(ordered) == 11856, len(ordered))
    require(len({value for row in ordered for value in row}) == 27, "ray level alphabet")
    require(max(value for row in ordered for value in row) == 180, "ray height")
    return ordered, tuple(sorted(common_hist))


def projective_ray_audit(body_limit: int = 649, scale: int = 1):
    """Exact all-ray check at one integer dilation scale."""
    require(1 <= body_limit <= 649 and scale >= 1, (body_limit, scale))
    rays, common_hist = projective_low_two_star_rays()
    weakest = None
    first_failure = None
    profile_checks = 0
    checked = 0
    for body, L in body_universe()[:body_limit]:
        _, ranges = R.safe_cell_ranges(body)
        cells = tuple(j for left, right in ranges for j in range(left, right))
        cell = cells[len(cells) // 2]
        centre = next(
            index
            for index, label in enumerate(body)
            if (label * cell) % L == L // 14
        )
        second = next(index for index in range(6) if index != centre)
        leaves = tuple(index for index in range(6) if index not in (centre, second))
        positions = (centre, second, *leaves)
        two_star = tuple(edge for edge in EDGES if centre in edge or second in edge)
        arc_cache = {}
        profile_cache = {}

        def profile(i: int, j: int, p: int, q: int) -> F:
            nonlocal profile_checks
            key = (i, j, p, q)
            if key not in profile_cache:
                for index, level in ((i, p), (j, q)):
                    arc_cache.setdefault(
                        (index, level),
                        R.reflected_level_arcs(L, body[index], level, cell),
                    )
                profile_cache[key] = intersection_mass(
                    arc_cache[i, p], arc_cache[j, q]
                )
                profile_checks += 1
            return profile_cache[key]

        for ray in rays:
            word_list = [0] * 6
            for position, value in zip(positions, ray):
                word_list[position] = scale * value
            word = tuple(word_list)
            debt = sum(
                (F(body[i], 7 * (word[i] * L - body[i])) for i in range(6)),
                F(0),
            )
            gain, edge = max(
                (profile(i, j, word[i], word[j]), (i, j)) for i, j in two_star
            )
            margin = gain - debt
            row = (
                margin,
                body,
                L,
                cell,
                centre,
                second,
                ray,
                word,
                edge,
                gain,
                debt,
            )
            if weakest is None or row < weakest:
                weakest = row
            if margin <= 0 and (first_failure is None or row < first_failure):
                first_failure = row
            checked += 1
        if first_failure is not None:
            break
    return {
        "scale": scale,
        "body_limit": body_limit,
        "ray_count": len(rays),
        "common_neighbour_hist": common_hist,
        "checked": checked,
        "profile_checks": profile_checks,
        "weakest": weakest,
        "first_failure": first_failure,
    }


def hostile_failure_relation(max_level: int):
    """Exact bounded slice of the nine addressed stalk failure relations."""
    body = (1, 2, 3, 4, 6, 12)
    L, ranges = R.safe_cell_ranges(body)
    cells = tuple(j for left, right in ranges for j in range(left, right))
    cell = cells[len(cells) // 2]
    centre = next(
        index for index, label in enumerate(body)
        if (label * cell) % L == L // 14
    )
    second = next(index for index in range(6) if index != centre)
    edges = tuple(edge for edge in EDGES if centre in edge or second in edge)
    debt = LP.universal_debt(body, L)
    arcs = {
        (index, q): R.reflected_level_arcs(L, body[index], q, cell)
        for index in range(6)
        for q in range(1, max_level + 1)
    }
    failures = Counter()
    margins = []
    digest = hashlib.sha256()
    for edge in edges:
        i, j = edge
        for p in range(1, max_level + 1):
            for q in range(1, max_level + 1):
                if p == q:
                    continue
                overlap = intersection_mass(arcs[i, p], arcs[j, q])
                divisor = gcd(p, q)
                key = (edge, p // divisor, q // divisor)
                margin = overlap - debt
                margins.append((margin, edge, (p, q), overlap))
                if margin <= 0:
                    failures[key] += 1
                    digest.update(f"{edge}|{p}|{q}|{overlap}|{debt}\n".encode())
    by_edge = Counter()
    high = 0
    low = 0
    for (edge, p, q), count in failures.items():
        by_edge[edge] += count
        if p + q >= 8:
            high += count
        else:
            low += count
    # A bounded box cannot certify persistence of a primitive channel: near
    # the top of the box only one dilation is visible.  Keep a compact
    # diagnostic sample and avoid attaching an asymptotic meaning to it.
    box_saturated = tuple(
        sorted(
            (count, edge, (p, q))
            for (edge, p, q), count in failures.items()
            if count >= max_level // max(p, q)
        )
    )
    return {
        "body": body,
        "ruler": L,
        "cell": cell,
        "pivots": (centre, second),
        "pivot_labels": (body[centre], body[second]),
        "max_level": max_level,
        "universal_debt": debt,
        "ordered_checks": len(edges) * max_level * (max_level - 1),
        "failure_count": sum(failures.values()),
        "low_failure_count": low,
        "high_failure_count": high,
        "by_edge": tuple(sorted(by_edge.items())),
        "weakest": min(margins),
        "box_saturated_count": len(box_saturated),
        "box_saturated_sample": box_saturated[:12],
        "failure_digest": digest.hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--m", type=int, default=1)
    parser.add_argument("--spread", type=int, default=6)
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--sample", type=int, default=0)
    parser.add_argument("--high-q", type=int, default=0)
    parser.add_argument("--ray-audit", action="store_true")
    parser.add_argument("--body-limit", type=int, default=649)
    parser.add_argument("--scale", type=int, default=1)
    parser.add_argument("--failure-q", type=int, default=0)
    args = parser.parse_args()
    if args.failure_q:
        print(hostile_failure_relation(args.failure_q))
        return
    if args.ray_audit:
        print(projective_ray_audit(args.body_limit, args.scale))
        return
    if args.high_q:
        print(high_channel_audit(args.high_q))
        return
    if args.sample:
        print(sampled_scan(args.m, args.spread, args.sample))
        return
    summary, rows, _ = scan(args.m, args.spread, progress=args.progress)
    print("summary", summary)
    print("body_rows")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
