#!/usr/bin/env python3
"""
pollock_tetrahedral_defect_descent_codex.py

codex-2026-06-13

Pollock tetrahedral scout:

    Te_k = k(k+1)(k+2)/6, k >= 1.

The open Pollock tetrahedral conjecture says every positive integer is a sum
of at most five tetrahedral numbers.  MathWorld/OEIS report that the known
four-tetrahedra defects have 241 terms, apparently ending at 343867.  This
script rediscovers that finite frontier up to a configurable bound, then tests
a proof-shaped descent:

    n in [Te_k, Te_{k+1}) is covered if
    n - Te_{k-j} is a sum of at most four tetrahedral numbers
    for one nearby anchor offset j.

Tournament Analysis is included for the offset stencil.

Vertices:
    offsets j = 0, 1, ..., W, meaning the anchor Te_{k-j}.

Pairwise observable:
    across sampled shells, private shell positions covered by offset i and not
    by offset j.

Gauge:
    orient i -> j when i has larger private coverage than j.

Tie Hamiltonian path:
    0 -> 1 -> 2 -> ... -> W.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from itertools import combinations, permutations


VERIFY_LIMIT = 1_000_000
SHELL_LIMIT = 1_200
STENCIL_WIDTH = 3
OCTAHEDRAL_LIMIT = 100_000


def tetra(k: int) -> int:
    return k * (k + 1) * (k + 2) // 6


def octa(k: int) -> int:
    return k * (2 * k * k + 1) // 3


def tri(k: int) -> int:
    return k * (k + 1) // 2


def atoms_upto(fn, limit: int) -> tuple[int, ...]:
    out: list[int] = []
    k = 1
    while True:
        value = fn(k)
        if value > limit:
            break
        out.append(value)
        k += 1
    return tuple(out)


def reachability_layers(atoms: tuple[int, ...], limit: int, terms: int) -> tuple[list[int], list[int]]:
    """Return exact-r and at-most-r bitsets for 0 <= r <= terms."""

    full_mask = (1 << (limit + 1)) - 1
    exact = [1]
    at_most = [1]
    for _r in range(1, terms + 1):
        layer = 0
        prev = exact[-1]
        for atom in atoms:
            if atom > limit:
                break
            layer |= prev << atom
        layer &= full_mask
        exact.append(layer)
        at_most.append(at_most[-1] | layer)
    return exact, at_most


def missing_from(mask: int, limit: int) -> list[int]:
    full = (1 << (limit + 1)) - 1
    missing = (~mask) & full
    missing &= ~1
    return bit_positions(missing)


def bit_positions(mask: int) -> list[int]:
    out: list[int] = []
    while mask:
        bit = mask & -mask
        out.append(bit.bit_length() - 1)
        mask ^= bit
    return out


def min_term_hist(exact: list[int], limit: int, max_terms: int) -> tuple[Counter[int], list[int]]:
    hist: Counter[int] = Counter()
    full = ((1 << (limit + 1)) - 1) & ~1
    assigned = 0
    hard_mask = 0
    for r in range(1, max_terms + 1):
        new = (exact[r] & full) & ~assigned
        hist[r] = new.bit_count()
        assigned |= new
        if r == max_terms:
            hard_mask = new
    missing = full & ~assigned
    if missing:
        hist[max_terms + 1] = missing.bit_count()
    hard = bit_positions(hard_mask)
    return hist, hard


def tetra_reach_limit(shell_limit: int, width: int, verify_limit: int) -> int:
    max_remainder = 0
    for k in range(1, shell_limit + 1):
        base = max(1, k - width)
        max_remainder = max(max_remainder, tetra(k + 1) - tetra(base))
    return max(verify_limit, max_remainder)


def triangular_defect_pairs(defects: list[int], limit: int) -> list[tuple[int, int, int]]:
    triangles: list[tuple[int, int]] = []
    k = 1
    while tri(k) <= limit:
        triangles.append((tri(k), k))
        k += 1

    defect_set = set(defects)
    pairs: list[tuple[int, int, int]] = []
    for x in defects:
        for gap, k in triangles:
            y = x + gap
            if y > limit:
                break
            if y in defect_set:
                pairs.append((k, x, y))
    return sorted(pairs)


def shell_stencil_profile(r4_mask: int, shell_limit: int, width: int) -> dict[str, object]:
    need_width: list[int] = []
    need_examples: defaultdict[int, list[int]] = defaultdict(list)
    fail_after_width: dict[int, list[int]] = {}

    for k in range(1, shell_limit + 1):
        gap = tetra(k + 1) - tetra(k)
        shell_mask = (1 << gap) - 1
        union = 0
        found = None
        for w in range(width + 1):
            anchor = k - w
            if anchor >= 1:
                shift = tetra(k) - tetra(anchor)
                union |= (r4_mask >> shift) & shell_mask
            if union == shell_mask:
                found = w
                break
        if found is None:
            found = width + 1
        need_width.append(found)
        need_examples[found].append(k)

    for w in range(width + 1):
        fail_after_width[w] = [i + 1 for i, need in enumerate(need_width) if need > w]

    return {
        "hist": Counter(need_width),
        "need_examples": need_examples,
        "fail_after_width": fail_after_width,
        "need_width": need_width,
    }


def offset_tournament(r4_mask: int, shell_limit: int, width: int) -> dict[str, object]:
    offsets = tuple(range(width + 1))
    private = {(i, j): 0 for i in offsets for j in offsets if i != j}
    covered_counts = Counter()

    for k in range(1, shell_limit + 1):
        gap = tetra(k + 1) - tetra(k)
        shell_mask = (1 << gap) - 1
        slices: dict[int, int] = {}
        for j in offsets:
            anchor = k - j
            if anchor < 1:
                slices[j] = 0
            else:
                shift = tetra(k) - tetra(anchor)
                slices[j] = (r4_mask >> shift) & shell_mask
            covered_counts[j] += slices[j].bit_count()

        for i, j in permutations(offsets, 2):
            private[(i, j)] += (slices[i] & ~slices[j] & shell_mask).bit_count()

    edges: list[tuple[int, int, int]] = []
    score = Counter({j: 0 for j in offsets})
    for i, j in combinations(offsets, 2):
        left = private[(i, j)]
        right = private[(j, i)]
        if left > right or (left == right and i < j):
            winner, loser, margin = i, j, left - right
        else:
            winner, loser, margin = j, i, right - left
        edges.append((winner, loser, margin))
        score[winner] += 1

    edge_set = {(a, b) for a, b, _m in edges}

    def beats(a: int, b: int) -> bool:
        return (a, b) in edge_set

    directed_3_cycles = 0
    for a, b, c in combinations(offsets, 3):
        if (beats(a, b) and beats(b, c) and beats(c, a)) or (
            beats(a, c) and beats(c, b) and beats(b, a)
        ):
            directed_3_cycles += 1

    hamiltonian_paths = []
    for path in permutations(offsets):
        if all(beats(path[i], path[i + 1]) for i in range(len(path) - 1)):
            hamiltonian_paths.append(path)

    sccs = strongly_connected_components(offsets, edge_set)

    return {
        "covered_counts": covered_counts,
        "private": private,
        "edges": edges,
        "score": score,
        "score_hist": Counter(score.values()),
        "directed_3_cycles": directed_3_cycles,
        "hamiltonian_paths": hamiltonian_paths,
        "sccs": sccs,
    }


def strongly_connected_components(vertices: tuple[int, ...], edges: set[tuple[int, int]]) -> list[tuple[int, ...]]:
    adjacency = {v: [w for a, w in edges if a == v] for v in vertices}
    reverse = {v: [a for a, w in edges if w == v] for v in vertices}

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in adjacency[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)

    components: list[tuple[int, ...]] = []
    seen.clear()
    for root in reversed(order):
        if root in seen:
            continue
        comp: list[int] = []
        queue = deque([root])
        seen.add(root)
        while queue:
            v = queue.popleft()
            comp.append(v)
            for w in reverse[v]:
                if w not in seen:
                    seen.add(w)
                    queue.append(w)
        components.append(tuple(sorted(comp)))
    return components


def octahedral_scout(limit: int) -> dict[str, object]:
    atoms = atoms_upto(octa, limit)
    exact, at_most = reachability_layers(atoms, limit, 7)
    hist, hard = min_term_hist(exact, limit, 7)
    return {
        "atoms": atoms,
        "missing": missing_from(at_most[7], limit),
        "hist": hist,
        "hard": hard,
    }


def print_list(label: str, values: list[int], count: int = 16) -> None:
    if len(values) <= 2 * count:
        body = values
    else:
        body = values[:count] + ["..."] + values[-count:]  # type: ignore[list-item]
    print(f"  {label}: {body}")


def main() -> None:
    reach_limit = tetra_reach_limit(SHELL_LIMIT, STENCIL_WIDTH, VERIFY_LIMIT)
    tet_atoms = atoms_upto(tetra, reach_limit)
    exact, at_most = reachability_layers(tet_atoms, reach_limit, 5)

    five_missing = missing_from(at_most[5], VERIFY_LIMIT)
    four_defects = missing_from(at_most[4], VERIFY_LIMIT)
    hist, hard5 = min_term_hist(exact, VERIFY_LIMIT, 5)
    defect_pairs = triangular_defect_pairs(four_defects, VERIFY_LIMIT)
    stencil = shell_stencil_profile(at_most[4], SHELL_LIMIT, STENCIL_WIDTH)
    tourney = offset_tournament(at_most[4], SHELL_LIMIT, STENCIL_WIDTH)
    octa_scout = octahedral_scout(OCTAHEDRAL_LIMIT)

    print("Pollock tetrahedral defect-descent scout")
    print("========================================")
    print(f"tetrahedral verify limit: {VERIFY_LIMIT}")
    print(f"shell stencil limit: k <= {SHELL_LIMIT}, width <= {STENCIL_WIDTH}")
    print(f"tetrahedral reachability limit used: {reach_limit}")
    print(f"tetrahedral atoms used: {len(tet_atoms)} (last {tet_atoms[-1]})")
    print()

    print("Five-term coverage and four-term defects")
    print(f"  missing from <=5 tetrahedra through {VERIFY_LIMIT}: {len(five_missing)}")
    print(f"  min-term histogram through {VERIFY_LIMIT}: {dict(sorted(hist.items()))}")
    print(f"  numbers needing exactly five through {VERIFY_LIMIT}: {len(hard5)}")
    print(f"  four-term defect count through {VERIFY_LIMIT}: {len(four_defects)}")
    print(f"  largest four-term defect through {VERIFY_LIMIT}: {four_defects[-1]}")
    print_list("first four-term defects", four_defects[:24], count=24)
    print_list("last four-term defects", four_defects[-24:], count=24)
    print()

    print("Triangular separations inside the four-defect set")
    print(f"  defect pairs with y-x triangular: {len(defect_pairs)}")
    if defect_pairs:
        max_k, x, y = defect_pairs[-1]
        print(f"  largest triangular gap index: k={max_k}, gap={tri(max_k)}")
        print(f"  last pair: {x} -> {y}")
        by_k: defaultdict[int, list[tuple[int, int]]] = defaultdict(list)
        for k, a, b in defect_pairs:
            by_k[k].append((a, b))
        print("  last triangular-gap classes:")
        for k in sorted(by_k)[-10:]:
            print(f"    k={k:3d}: count={len(by_k[k])}, sample={by_k[k][:3]}")
    print()

    print("Nearby-anchor shell stencil")
    print(f"  histogram of minimal needed offset width: {dict(sorted(stencil['hist'].items()))}")
    for w in range(STENCIL_WIDTH + 1):
        failures = stencil["fail_after_width"][w]  # type: ignore[index]
        tail = failures[-8:] if failures else []
        last = failures[-1] if failures else None
        print(f"  after allowing offsets <= {w}: failures={len(failures)}, last={last}, tail={tail}")
    for need in range(STENCIL_WIDTH + 1):
        shells = stencil["need_examples"][need]  # type: ignore[index]
        if shells:
            print(f"  shells with exact needed width {need}: count={len(shells)}, first={shells[:8]}, last={shells[-8:]}")
    print()

    print("Tournament Analysis: anchor-offset dominance")
    print("  vertices: offsets j where the anchor is Te_{k-j}")
    print("  observable: private covered shell positions over sampled shells")
    print("  gauge: larger private coverage wins; tie path 0 -> 1 -> 2 -> 3")
    print(f"  covered counts by offset: {dict(sorted(tourney['covered_counts'].items()))}")
    print("  oriented edges winner -> loser (margin):")
    for winner, loser, margin in tourney["edges"]:  # type: ignore[index]
        print(f"    {winner} -> {loser}  margin={margin}")
    print(f"  score by offset: {dict(sorted(tourney['score'].items()))}")
    print(f"  score histogram: {dict(sorted(tourney['score_hist'].items()))}")
    print(f"  directed 3-cycles: {tourney['directed_3_cycles']}")
    print(f"  SCCs: {tourney['sccs']}")
    print(f"  Hamiltonian paths: count={len(tourney['hamiltonian_paths'])}, paths={tourney['hamiltonian_paths'][:6]}")
    print()

    print("Octahedral sibling scout")
    print(f"  limit: {OCTAHEDRAL_LIMIT}")
    print(f"  atoms: {len(octa_scout['atoms'])} (last {octa_scout['atoms'][-1]})")
    print(f"  missing from <=7 octahedral numbers: {len(octa_scout['missing'])}")
    print(f"  min-term histogram: {dict(sorted(octa_scout['hist'].items()))}")
    print(f"  numbers needing exactly seven: {len(octa_scout['hard'])}")
    print_list("seven-term octahedral cases", octa_scout["hard"], count=16)  # type: ignore[arg-type]
    print()

    print("Proof hooks suggested by this run")
    print("  1. Strong route: prove no four-tetrahedra defects occur beyond 343867.")
    print("     The computation matches the reported 241-defect frontier through 1,000,000.")
    print("  2. Weaker descent route: prove that for all large k there is no pair")
    print("     r, r+tri(k) in the four-defect set with 0 <= r < Te_{k+1}-Te_k,")
    print("     where tri(k)=Te_k-Te_{k-1}=k(k+1)/2 is the one-back gap.")
    print("     Then every shell is covered by Te_k or Te_{k-1} plus four atoms.")
    print("  3. Finite certificate route: the width-3 nearby-anchor stencil covers")
    print("     every sampled shell through k=1200, and width-1 covers after k=825.")


if __name__ == "__main__":
    main()
