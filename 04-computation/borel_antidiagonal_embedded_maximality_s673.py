#!/usr/bin/env python3
"""S673: finite shadows of Borel anti-diagonalization.

Prompt:

    Borel anti-diagonalization, constructive mathematics, tangible
    incompleteness, embedded maximality.

This script builds a deliberately finite model of the selector issue behind
Friedman's Borel diagonalization theme.

Classical diagonalization says: given a named list, construct something outside
it.  Friedman's invariant Borel diagonalization result has the opposite flavor:
if the selector is sufficiently invariant, it cannot uniformly stay outside the
named countable object.  The repo-facing finite shadow is:

    outside selector + invariance => address tax;
    outside selector + outer extension => recursive capture.

The model:
  - U is a finite universe.
  - A state is a named subset A of U.
  - An outside selector must choose y in U \\ A.
  - A symmetry group G acts on U.
  - A selector is possible at A only if some outside y is fixed by every
    automorphism stabilizing A.  Otherwise invariance gives no canonical
    outside point.
  - Naming anchors shrinks G.  The minimum number of anchors needed for every
    non-full state to have an invariant outside selector is the address tax.

This is not Borel mathematics.  It is a finite diagnostic for the same proof
shape: diagonal "missing object" arguments become usable only after the
observer/address coordinate and recursive extension law are named.

Tournament Analysis:
  Vertices are proof/selector lanes, not elements of U.  Pairwise observable is
  (invariance_respect, constructive_witness, needs_address, recursion_rank,
  embedded_maximality_fit, lrc_transfer, overclaim_safety).  The switch is
  majority comparison with listed order as the tie Hamiltonian path.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations, permutations


def banner(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def mask_from_set(xs: tuple[int, ...] | set[int]) -> int:
    out = 0
    for x in xs:
        out |= 1 << x
    return out


def set_from_mask(mask: int, n: int) -> tuple[int, ...]:
    return tuple(i for i in range(n) if (mask >> i) & 1)


def apply_perm_to_mask(mask: int, perm: tuple[int, ...]) -> int:
    out = 0
    for i, j in enumerate(perm):
        if (mask >> i) & 1:
            out |= 1 << j
    return out


def fixed_points(perms: list[tuple[int, ...]], n: int) -> set[int]:
    return {i for i in range(n) if all(p[i] == i for p in perms)}


def symmetric_group(n: int) -> list[tuple[int, ...]]:
    return [tuple(p) for p in permutations(range(n))]


def cyclic_group(n: int) -> list[tuple[int, ...]]:
    return [tuple((i + k) % n for i in range(n)) for k in range(n)]


def dihedral_group(n: int) -> list[tuple[int, ...]]:
    perms = set()
    for k in range(n):
        perms.add(tuple((i + k) % n for i in range(n)))
        perms.add(tuple((k - i) % n for i in range(n)))
    return sorted(perms)


def path_reflection_group(n: int) -> list[tuple[int, ...]]:
    return [tuple(range(n)), tuple(n - 1 - i for i in range(n))]


def trivial_group(n: int) -> list[tuple[int, ...]]:
    return [tuple(range(n))]


def stabilizer_of_mask(group: list[tuple[int, ...]], mask: int) -> list[tuple[int, ...]]:
    return [p for p in group if apply_perm_to_mask(mask, p) == mask]


def group_with_anchors(
    group: list[tuple[int, ...]], anchors: tuple[int, ...]
) -> list[tuple[int, ...]]:
    anchor_set = set(anchors)
    return [p for p in group if all(p[a] == a for a in anchor_set)]


def selectable(mask: int, n: int, group: list[tuple[int, ...]]) -> bool:
    if mask == (1 << n) - 1:
        return False
    stab = stabilizer_of_mask(group, mask)
    fixed = fixed_points(stab, n)
    complement = {i for i in range(n) if not ((mask >> i) & 1)}
    return bool(fixed & complement)


def selectable_point(mask: int, n: int, group: list[tuple[int, ...]]) -> int | None:
    """Return the least invariantly selectable outside point, if one exists."""
    if mask == (1 << n) - 1:
        return None
    stab = stabilizer_of_mask(group, mask)
    fixed = fixed_points(stab, n)
    for i in range(n):
        if i in fixed and not ((mask >> i) & 1):
            return i
    return None


def selector_audit(n: int, group: list[tuple[int, ...]]) -> dict[str, object]:
    selectable_states = 0
    blocked_states = []
    by_size: dict[int, Counter[str]] = {k: Counter() for k in range(n)}
    for mask in range((1 << n) - 1):
        key = "selectable" if selectable(mask, n, group) else "blocked"
        by_size[mask.bit_count()][key] += 1
        if key == "selectable":
            selectable_states += 1
        else:
            blocked_states.append(mask)
    return {
        "n": n,
        "group_size": len(group),
        "states": (1 << n) - 1,
        "selectable": selectable_states,
        "blocked": len(blocked_states),
        "blocked_examples": [set_from_mask(m, n) for m in blocked_states[:6]],
        "by_size": {k: dict(v) for k, v in by_size.items() if v},
    }


def address_tax(n: int, group: list[tuple[int, ...]]) -> tuple[int, tuple[int, ...]]:
    """Minimum pointwise-named anchors making all non-full states selectable."""
    all_masks = range((1 << n) - 1)
    for k in range(n + 1):
        for anchors in combinations(range(n), k):
            anchored_group = group_with_anchors(group, anchors)
            if all(selectable(mask, n, anchored_group) for mask in all_masks):
                return k, anchors
    raise AssertionError("n anchors should trivialize every finite group action")


def closure_sequence(n: int, group: list[tuple[int, ...]]) -> tuple[tuple[int, ...], ...]:
    """Adversarial outer-extension closure A -> A union {selector(A)}."""
    mask = 0
    seq = []
    while mask != (1 << n) - 1:
        y = selectable_point(mask, n, group)
        if y is None:
            break
        seq.append((set_from_mask(mask, n), y))
        mask |= 1 << y
    return tuple(seq)


@dataclass(frozen=True)
class Lane:
    name: str
    vector: tuple[int, int, int, int, int, int, int]
    note: str

    @property
    def total(self) -> int:
        return sum(self.vector)


LANES = [
    Lane(
        "embedded_address_tax",
        (5, 5, 5, 5, 5, 5, 5),
        "Names the least coordinate needed before outside-selection is usable.",
    ),
    Lane(
        "borel_invariant_antiselector",
        (5, 3, 5, 4, 5, 4, 5),
        "Invariant selectors cannot simply choose a new outside point.",
    ),
    Lane(
        "ph_bad_branch_rank",
        (4, 4, 4, 5, 4, 5, 5),
        "Bad colorings are killed by extension rank rather than raw density.",
    ),
    Lane(
        "lrc_owner_carry_rank",
        (4, 4, 5, 5, 5, 5, 4),
        "Owner/carry address should become a rank drop under +27 extension.",
    ),
    Lane(
        "a000568_half_filter_rank",
        (4, 5, 4, 4, 4, 4, 5),
        "Deleted-card L/M/U owner address repairs quotient leakage.",
    ),
    Lane(
        "constructive_named_order",
        (2, 5, 2, 3, 3, 3, 5),
        "A least-outside rule works after an order is named, then gets captured.",
    ),
    Lane(
        "raw_cantor_diagonal",
        (1, 4, 1, 2, 2, 2, 3),
        "Builds an outside object but ignores invariance and extension purity.",
    ),
    Lane(
        "raw_count_shadow",
        (1, 2, 1, 1, 1, 1, 4),
        "Counts states but does not know which selector survives extension.",
    ),
]


def tournament_analysis() -> dict[str, object]:
    names = [lane.name for lane in LANES]
    wins = Counter({name: 0 for name in names})
    adj: dict[tuple[str, str], str] = {}
    for i, a in enumerate(LANES):
        for j, b in enumerate(LANES):
            if i >= j:
                continue
            votes_a = sum(x > y for x, y in zip(a.vector, b.vector))
            votes_b = sum(y > x for x, y in zip(a.vector, b.vector))
            if votes_a > votes_b or (votes_a == votes_b and i < j):
                winner, loser = a.name, b.name
            else:
                winner, loser = b.name, a.name
            wins[winner] += 1
            adj[(winner, loser)] = winner

    directed_3cycles = 0
    for a, b, c in combinations(names, 3):
        ab = (a, b) in adj
        bc = (b, c) in adj
        ca = (c, a) in adj
        ba = (b, a) in adj
        cb = (c, b) in adj
        ac = (a, c) in adj
        if (ab and bc and ca) or (ba and cb and ac):
            directed_3cycles += 1

    return {
        "score_hist": dict(sorted(Counter(wins.values()).items())),
        "scores": dict(sorted(wins.items(), key=lambda kv: (-kv[1], kv[0]))),
        "directed_3cycles": directed_3cycles,
        "hamiltonian_path": [name for name, _ in sorted(wins.items(), key=lambda kv: -kv[1])],
    }


def main() -> None:
    banner("Finite invariant outside-selector audit")
    group_builders = [
        ("ordered/trivial", trivial_group),
        ("path_reflection", path_reflection_group),
        ("cyclic_rotations", cyclic_group),
        ("dihedral_cycle", dihedral_group),
        ("full_symmetric", symmetric_group),
    ]
    for n in range(3, 8):
        print(f"\nn={n}")
        for name, builder in group_builders:
            group = builder(n)
            audit = selector_audit(n, group)
            tax, anchors = address_tax(n, group)
            print(
                f"  {name:17s} | |G|={len(group):5d} | "
                f"selectable={audit['selectable']:3d}/{audit['states']:3d} | "
                f"blocked={audit['blocked']:3d} | address_tax={tax} anchors={anchors}"
            )
            if n == 6 and audit["blocked"]:
                print(f"    blocked examples: {audit['blocked_examples']}")

    banner("Constructive selector captured by recursive outer extension")
    n = 6
    for name, builder in [
        ("ordered/trivial", trivial_group),
        ("path_reflection + anchor 0", lambda m: group_with_anchors(path_reflection_group(m), (0,))),
        ("dihedral + anchors 0,1", lambda m: group_with_anchors(dihedral_group(m), (0, 1))),
        ("symmetric + anchors 0..4", lambda m: group_with_anchors(symmetric_group(m), (0, 1, 2, 3, 4))),
    ]:
        group = builder(n)
        seq = closure_sequence(n, group)
        print(f"\n{name}: closure length {len(seq)}")
        for state, y in seq:
            print(f"  A={state} -> select {y}; outer extension names {y}")

    banner("Tournament Analysis over selector/proof lanes")
    print("coordinate order:")
    print(
        "  invariance_respect, constructive_witness, needs_address, "
        "recursion_rank, embedded_maximality_fit, lrc_transfer, overclaim_safety"
    )
    for lane in LANES:
        print(f"  {lane.name:30s} vector={lane.vector} total={lane.total}")
        print(f"    {lane.note}")
    ta = tournament_analysis()
    print("\nTournament fingerprints:")
    print(f"  score_hist={ta['score_hist']}")
    print(f"  directed_3cycles={ta['directed_3cycles']}")
    print(f"  hamiltonian_path={ta['hamiltonian_path']}")

    banner("Repo-facing claims")
    print(
        "1. Raw diagonalization is a missing-object move; Borel-style "
        "anti-diagonalization is an invariant-selector obstruction."
    )
    print(
        "2. Constructive mathematics enters as named finite procedures: "
        "least outside point, bounded search, and explicit closure sequences."
    )
    print(
        "3. Tangible incompleteness is the uniform version: finite stages are "
        "transparent, but the extension bound/rank can outrun the base theory."
    )
    print(
        "4. Embedded maximality is the repair: state which ambient extensions "
        "are allowed and which address makes the selector/maximal object stable."
    )
    print(
        "5. LRC14 transfer: owner-private deletion is an address tax already paid; "
        "the next theorem is a recursive rank drop on coherent +27 carry lifts."
    )


if __name__ == "__main__":
    main()
