"""Tournament perspective ladder and the first A000568 shift failure.

This script compares three related counts:

* U(n): unlabeled tournaments on n vertices, A000568.
* R(n): rooted/node-perspective tournaments on n vertices, equivalently the
  sum of vertex orbits over all unlabeled n-tournaments.
* finite-depth color perspectives inside each tournament, using an in/out
  Weisfeiler-Lehman style refinement.

The goal is not a fast tournament enumerator.  It is a small diagnostic for the
controlled-forgetting question: when does a one-node perspective stop carrying
the same amount of structure as the next unrooted tournament size?
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations, permutations
from math import factorial
from typing import Dict, Iterable, List, Sequence, Tuple


@dataclass(frozen=True)
class CountRow:
    n: int
    unlabeled: int
    rooted: int
    next_unlabeled: int
    equals_next: bool


def integer_partitions(n: int, max_part: int | None = None) -> Iterable[List[int]]:
    if max_part is None or max_part > n:
        max_part = n
    if n == 0:
        yield []
        return
    for part in range(max_part, 0, -1):
        if part <= n:
            for rest in integer_partitions(n - part, part):
                yield [part] + rest


def permutation_from_partition(parts: Sequence[int]) -> List[int]:
    perm: List[int] = []
    start = 0
    for length in parts:
        cycle = list(range(start, start + length))
        for index in range(length):
            perm.append(cycle[(index + 1) % length])
        start += length
    return perm


def conjugacy_class_size(parts: Sequence[int]) -> int:
    counts = Counter(parts)
    centralizer = 1
    for length, multiplicity in counts.items():
        centralizer *= (length**multiplicity) * factorial(multiplicity)
    return factorial(sum(parts)) // centralizer


class ParityDSU:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.parity = [0] * size

    def find(self, item: int) -> Tuple[int, int]:
        if self.parent[item] != item:
            root, parity = self.find(self.parent[item])
            self.parity[item] ^= parity
            self.parent[item] = root
        return self.parent[item], self.parity[item]

    def union(self, left: int, right: int, parity: int) -> bool:
        left_root, left_parity = self.find(left)
        right_root, right_parity = self.find(right)
        if left_root == right_root:
            return (left_parity ^ right_parity) == parity
        self.parent[left_root] = right_root
        self.parity[left_root] = left_parity ^ right_parity ^ parity
        return True


def fixed_tournament_count(parts: Sequence[int]) -> int:
    """Number of labelled tournaments fixed by one permutation of this type."""
    n = sum(parts)
    perm = permutation_from_partition(parts)
    pairs: List[Tuple[int, int]] = []
    pair_index: Dict[Tuple[int, int], int] = {}
    for left in range(n):
        for right in range(left + 1, n):
            pair_index[(left, right)] = len(pairs)
            pairs.append((left, right))

    dsu = ParityDSU(len(pairs))
    for index, (left, right) in enumerate(pairs):
        image_left, image_right = perm[left], perm[right]
        if image_left < image_right:
            image_index = pair_index[(image_left, image_right)]
            flip = 0
        else:
            image_index = pair_index[(image_right, image_left)]
            flip = 1
        if not dsu.union(index, image_index, flip):
            return 0
    components = {dsu.find(index)[0] for index in range(len(pairs))}
    return 2 ** len(components)


def unlabeled_tournament_count(n: int) -> int:
    total = 0
    for parts in integer_partitions(n):
        total += conjugacy_class_size(parts) * fixed_tournament_count(parts)
    return total // factorial(n)


def rooted_tournament_count(n: int) -> int:
    """Burnside count for unlabeled rooted tournaments on n vertices."""
    total = 0
    for parts in integer_partitions(n):
        fixed_roots = parts.count(1)
        total += conjugacy_class_size(parts) * fixed_roots * fixed_tournament_count(parts)
    return total // factorial(n)


def pair_indices(n: int) -> Dict[Tuple[int, int], int]:
    out: Dict[Tuple[int, int], int] = {}
    index = 0
    for left in range(n):
        for right in range(left + 1, n):
            out[(left, right)] = index
            index += 1
    return out


def permuted_bits(bits: int, perm: Sequence[int], index: Dict[Tuple[int, int], int]) -> int:
    n = len(perm)
    out = 0
    new_index = 0
    for left in range(n):
        for right in range(left + 1, n):
            old_left, old_right = perm[left], perm[right]
            if old_left < old_right:
                old_bit = (bits >> index[(old_left, old_right)]) & 1
                new_bit = old_bit
            else:
                old_bit = (bits >> index[(old_right, old_left)]) & 1
                new_bit = 1 - old_bit
            if new_bit:
                out |= 1 << new_index
            new_index += 1
    return out


def canonical_bits(bits: int, n: int, perms: Sequence[Tuple[int, ...]], index: Dict[Tuple[int, int], int]) -> int:
    return min(permuted_bits(bits, perm, index) for perm in perms)


def adjacency_from_bits(bits: int, n: int) -> List[List[bool]]:
    matrix = [[False] * n for _ in range(n)]
    index = 0
    for left in range(n):
        for right in range(left + 1, n):
            if (bits >> index) & 1:
                matrix[left][right] = True
            else:
                matrix[right][left] = True
            index += 1
    return matrix


def wl_perspective_classes(matrix: Sequence[Sequence[bool]], depth: int) -> int:
    n = len(matrix)
    colors = [0] * n
    for _ in range(depth):
        signatures = []
        for vertex in range(n):
            out_colors = sorted(colors[other] for other in range(n) if matrix[vertex][other])
            in_colors = sorted(colors[other] for other in range(n) if matrix[other][vertex])
            signatures.append((tuple(out_colors), tuple(in_colors)))
        palette = {sig: index for index, sig in enumerate(sorted(set(signatures)))}
        colors = [palette[sig] for sig in signatures]
    return len(set(colors))


def rooted_canonical_bits(
    bits: int,
    root: int,
    n: int,
    index: Dict[Tuple[int, int], int],
) -> int:
    others = [vertex for vertex in range(n) if vertex != root]
    return min(permuted_bits(bits, (root,) + perm, index) for perm in permutations(others))


def edge_canonical_bits(
    bits: int,
    edge: Tuple[int, int],
    n: int,
    index: Dict[Tuple[int, int], int],
) -> int:
    tail, tip = edge
    others = [vertex for vertex in range(n) if vertex not in edge]
    return min(permuted_bits(bits, (tail, tip) + perm, index) for perm in permutations(others))


def small_perspective_profile(max_n: int = 6, max_depth: int = 5) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for n in range(1, max_n + 1):
        index = pair_indices(n)
        all_perms = list(permutations(range(n)))
        representatives: Dict[int, List[List[bool]]] = {}
        for bits in range(1 << (n * (n - 1) // 2)):
            canon = canonical_bits(bits, n, all_perms, index)
            if canon not in representatives:
                representatives[canon] = adjacency_from_bits(canon, n)

        rooted_total = 0
        edge_total = 0
        rooted_profile: Counter[int] = Counter()
        wl_totals = [0] * (max_depth + 1)
        for canon, matrix in representatives.items():
            rooted_forms = {rooted_canonical_bits(canon, root, n, index) for root in range(n)}
            rooted_total += len(rooted_forms)
            rooted_profile[len(rooted_forms)] += 1

            edge_forms = set()
            for tail, tip in combinations(range(n), 2):
                if matrix[tail][tip]:
                    edge_forms.add(edge_canonical_bits(canon, (tail, tip), n, index))
                else:
                    edge_forms.add(edge_canonical_bits(canon, (tip, tail), n, index))
            edge_total += len(edge_forms)

            for depth in range(max_depth + 1):
                wl_totals[depth] += wl_perspective_classes(matrix, depth)

        rows.append(
            {
                "n": n,
                "unlabeled": len(representatives),
                "rooted": rooted_total,
                "rooted_profile": dict(sorted(rooted_profile.items())),
                "edge_perspectives": edge_total,
                "wl_depth_totals": wl_totals,
            }
        )
    return rows


def burnside_rows(max_n: int = 10) -> List[CountRow]:
    rows: List[CountRow] = []
    counts = {n: unlabeled_tournament_count(n) for n in range(1, max_n + 2)}
    for n in range(1, max_n + 1):
        rooted = rooted_tournament_count(n)
        rows.append(CountRow(n, counts[n], rooted, counts[n + 1], rooted == counts[n + 1]))
    return rows


def first_shift_failure(rows: Sequence[CountRow]) -> CountRow:
    for row in rows:
        if not row.equals_next:
            return row
    raise ValueError("no failure in supplied range")


def main() -> None:
    rows = burnside_rows(10)
    failure = first_shift_failure(rows)

    print("Tournament perspective ladder (codex-2026-06-26-S211)")
    print("U(n)=A000568 unlabeled tournaments; R(n)=unlabeled rooted/node perspectives")
    print()
    print("[0] Shift-equality table")
    print(" n  U(n)      R(n)      U(n+1)    R(n)=U(n+1)")
    for row in rows:
        print(
            f"{row.n:2d}  {row.unlabeled:<9d} {row.rooted:<9d} "
            f"{row.next_unlabeled:<9d} {row.equals_next}"
        )
    print()
    print(
        "first_failure: compare node perspectives on n-1 with A000568(n); "
        f"n={failure.n + 1}, because R({failure.n})={failure.rooted} "
        f"but U({failure.n + 1})={failure.next_unlabeled}."
    )
    print("defect_at_first_failure=8")
    print()

    print("[1] Burnside obstruction at the first failure")
    print("For U(6), the nonzero symmetry types include fixed-point-free [3,3].")
    for parts in integer_partitions(6):
        fixed = fixed_tournament_count(parts)
        if fixed:
            print(
                f"  cycle_type={parts} class_size={conjugacy_class_size(parts)} "
                f"fixed_tournaments={fixed} fixed_vertices={parts.count(1)}"
            )
    print(
        "readout: rooted perspectives only see symmetries with a named vertex; "
        "the [3,3] term is a rootless cyclic sidecar, so depth on nodes alone "
        "cannot supply the eight missing U(6) classes."
    )
    print()

    print("[2] Finite-depth node-perspective ladder")
    print("WL depth k uses in-neighbor/out-neighbor multisets of depth-(k-1) colors.")
    for row in small_perspective_profile():
        print(
            f"n={row['n']} U={row['unlabeled']} rooted={row['rooted']} "
            f"orbit_profile={row['rooted_profile']} "
            f"edge={row['edge_perspectives']} wl={row['wl_depth_totals']}"
        )
    print()

    print("[3] Controlled-forgetting interpretation")
    print(
        "The equality R(n-1)=U(n) holds through n=5 because one visible node "
        "perspective still behaves like a complete one-step forgetting sidecar. "
        "At n=6 the forgotten coordinate can be a rootless 3-cycle/3-cycle "
        "coupling, not a node.  The ladder therefore needs typed anchors: "
        "node -> edge -> cycle -> clique/extension -> conflict perspectives. "
        "In the matrix-atlas language, this says the missing coordinate must be "
        "kept as an observability sidecar or Schur-complement correction."
    )
    print()

    print("[4] Candidate exact perspective definitions")
    print(
        "node_k(v): recursively color v by the ordered pair of multisets "
        "(out-neighbor node_{k-1}, in-neighbor node_{k-1}); full node perspective "
        "is rooted-tournament isomorphism."
    )
    print(
        "edge_k(a->b): fix tail and tip, then color every external vertex by "
        "its two-bit relation to (a,b): common predator, common prey, transitive "
        "interior a->x->b, or cyclic return b->x->a, recursively decorated by "
        "node/edge colors."
    )
    print(
        "cycle_k(C3): fix a directed 3-cycle up to rotation; color external "
        "vertices by their 3-bit word to the cycle, modulo cyclic rotation, then "
        "iterate on the colored quotient.  This is the first sidecar suggested "
        "by the [3,3] Burnside obstruction."
    )
    print(
        "clique_k(Q): fix an isomorphism class of a rooted subtournament Q "
        "(transitive chain, cyclic triangle, or larger pattern), quotient by "
        "Aut(Q), and color external vertices by their relation word to Q.  "
        "This interpolates between node/edge roots and arbitrary conflict roots."
    )
    print(
        "conflict_k: use any failed quotient fiber as the root object; two "
        "fibers are k-same when their repair carriers have isomorphic colored "
        "incidence to node, edge, cycle, clique, and extension sidecars."
    )


if __name__ == "__main__":
    main()
