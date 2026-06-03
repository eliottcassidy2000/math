#!/usr/bin/env python3
"""S589: vertex-transitive trienerments, polygon rigidity, and cascade.

The user's proposed lens is:

    vertex-transitive trienerment <=> regular polygon point-set.

This script separates three meanings of "vertex-transitive":

1. Cyclic/rotation-transitive point-set: rigid, forces a regular polygon.
2. Dihedral-transitive point-set: allows imprimitive regular bracelets.
3. Abstract vertex-transitive tournament: Cayley rigidity propagates the local
   star globally, but the group need not contain a cyclic polygon spine.

The point is not to refute the slogan, but to locate the exact rigidity:
cyclic sharp transitivity gives a polygon; dihedral or nonabelian symmetry gives
local isomorphism plus a block/relator cascade.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from itertools import combinations, permutations
from math import gcd


def circular_gaps(points: tuple[int, ...], modulus: int) -> tuple[int, ...]:
    pts = tuple(sorted(points))
    return tuple((pts[(i + 1) % len(pts)] - pts[i]) % modulus for i in range(len(pts)))


def least_period(word: tuple[int, ...]) -> int:
    for p in range(1, len(word) + 1):
        if len(word) % p == 0 and all(word[i] == word[i % p] for i in range(len(word))):
            return p
    return len(word)


def is_regular_polygon(points: tuple[int, ...], modulus: int) -> bool:
    gaps = circular_gaps(points, modulus)
    return len(set(gaps)) == 1


def dihedral_maps(modulus: int):
    for shift in range(modulus):
        yield (1, shift)
        yield (-1, shift)


def apply_map(x: int, eps: int, shift: int, modulus: int) -> int:
    return (shift + eps * x) % modulus


def stabilizer(points: tuple[int, ...], modulus: int) -> tuple[tuple[int, int], ...]:
    P = frozenset(points)
    out = []
    for eps, shift in dihedral_maps(modulus):
        if frozenset(apply_map(x, eps, shift, modulus) for x in P) == P:
            out.append((eps, shift))
    return tuple(out)


def orbits(points: tuple[int, ...], maps: tuple[tuple[int, int], ...], modulus: int) -> tuple[tuple[int, ...], ...]:
    P = set(points)
    unseen = set(points)
    outs = []
    while unseen:
        start = unseen.pop()
        q = deque([start])
        orb = {start}
        while q:
            x = q.popleft()
            for eps, shift in maps:
                y = apply_map(x, eps, shift, modulus)
                if y in P and y not in orb:
                    orb.add(y)
                    unseen.discard(y)
                    q.append(y)
        outs.append(tuple(sorted(orb)))
    return tuple(sorted(outs, key=lambda o: (len(o), o)))


def distance_profile(points: tuple[int, ...], modulus: int, p: int) -> tuple[int, ...]:
    return tuple(sorted(min((q - p) % modulus, (p - q) % modulus) for q in points if q != p))


@dataclass(frozen=True)
class PointSetAudit:
    modulus: int
    points: tuple[int, ...]
    gaps: tuple[int, ...]
    gap_period: int
    stabilizer_size: int
    orbit_count: int
    distance_profile_equal: bool
    regular: bool

    @property
    def dihedral_vt(self) -> bool:
        return self.orbit_count == 1

    @property
    def imprimitive_vt(self) -> bool:
        return self.dihedral_vt and not self.regular


def audit_point_set(points: tuple[int, ...], modulus: int) -> PointSetAudit:
    maps = stabilizer(points, modulus)
    profiles = [distance_profile(points, modulus, p) for p in points]
    return PointSetAudit(
        modulus=modulus,
        points=points,
        gaps=circular_gaps(points, modulus),
        gap_period=least_period(circular_gaps(points, modulus)),
        stabilizer_size=len(maps),
        orbit_count=len(orbits(points, maps, modulus)),
        distance_profile_equal=len(set(profiles)) == 1,
        regular=is_regular_polygon(points, modulus),
    )


def point_set_report(max_modulus: int = 18) -> list[str]:
    audits: list[PointSetAudit] = []
    for N in range(3, max_modulus + 1):
        for r in range(3, N + 1):
            for rest in combinations(range(1, N), r - 1):
                audits.append(audit_point_set((0,) + rest, N))

    dihedral_vt = [a for a in audits if a.dihedral_vt]
    regular = [a for a in dihedral_vt if a.regular]
    imprimitive = [a for a in dihedral_vt if not a.regular]
    local_equal = [a for a in audits if a.distance_profile_equal]
    local_not_vt = [a for a in local_equal if not a.dihedral_vt]

    first_imprimitive = min(imprimitive, key=lambda a: (a.modulus, len(a.points), a.points))
    first_local_not_vt = min(local_not_vt, key=lambda a: (a.modulus, len(a.points), a.points)) if local_not_vt else None

    gap_periods = Counter(a.gap_period for a in imprimitive)
    size_counts = Counter(len(a.points) for a in imprimitive)

    lines = [
        "POINT-SET / TRIENERMENT RIGIDITY AUDIT",
        f"  searched: subsets P of Z/N with 0 in P, 3<=N<={max_modulus}",
        f"  dihedral_vertex_transitive={len(dihedral_vt)}; regular_polygons={len(regular)}; imprimitive_bracelets={len(imprimitive)}",
        f"  local_distance_profile_equal={len(local_equal)}; local_equal_but_not_dihedral_VT={len(local_not_vt)}",
        "",
        "  first imprimitive dihedral-VT point-set:",
        f"    N={first_imprimitive.modulus}, P={first_imprimitive.points}, gaps={first_imprimitive.gaps}, gap_period={first_imprimitive.gap_period}, stabilizer_size={first_imprimitive.stabilizer_size}",
        "    reading: every vertex has the same local trienerment view, but the set is a regular polygon of blocks, not a regular polygon of points.",
    ]
    if first_local_not_vt:
        lines += [
            "",
            "  first local-profile-equal but not globally dihedral-VT point-set:",
            f"    N={first_local_not_vt.modulus}, P={first_local_not_vt.points}, gaps={first_local_not_vt.gaps}, orbit_count={first_local_not_vt.orbit_count}",
            "    reading: local fixed-point data alone can be too weak; a symmetry action is the cascade mechanism.",
        ]
    lines += [
        "",
        f"  imprimitive gap-period histogram: {dict(sorted(gap_periods.items()))}",
        f"  imprimitive size histogram: {dict(sorted(size_counts.items()))}",
    ]
    return lines


class CyclicGroup:
    def __init__(self, n: int):
        self.n = n
        self.elements = tuple(range(n))
        self.identity = 0

    def mul(self, x: int, y: int) -> int:
        return (x + y) % self.n

    def inv(self, x: int) -> int:
        return (-x) % self.n

    def order(self, x: int) -> int:
        if x == 0:
            return 1
        return self.n // gcd(self.n, x)

    def label(self, x: int) -> str:
        return str(x)


class Frobenius21:
    """C7 semidirect C3, with b*a*b^-1 = a^2."""

    def __init__(self):
        self.n = 21
        self.elements = tuple((a, b) for b in range(3) for a in range(7))
        self.identity = (0, 0)

    def mul(self, x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
        a, b = x
        c, d = y
        return ((a + pow(2, b, 7) * c) % 7, (b + d) % 3)

    def inv(self, x: tuple[int, int]) -> tuple[int, int]:
        a, b = x
        neg_b = (-b) % 3
        return ((-pow(2, neg_b, 7) * a) % 7, neg_b)

    def order(self, x: tuple[int, int]) -> int:
        y = self.identity
        for k in range(1, 22):
            y = self.mul(y, x)
            if y == self.identity:
                return k
        raise RuntimeError("order not found")

    def label(self, x: tuple[int, int]) -> str:
        return f"a{x[0]}b{x[1]}"


def inverse_pairs(group) -> list[tuple[object, object]]:
    seen = {group.identity}
    pairs = []
    for x in group.elements:
        if x in seen:
            continue
        y = group.inv(x)
        pairs.append((x, y))
        seen.add(x)
        seen.add(y)
    return pairs


def default_connection(group, mode: str) -> frozenset[object]:
    pairs = inverse_pairs(group)
    chosen = []
    for i, (x, y) in enumerate(pairs):
        if mode == "interval" and isinstance(group, CyclicGroup):
            chosen.append(x if 1 <= x <= group.n // 2 else y)
        elif mode == "alternating":
            chosen.append(x if i % 2 == 0 else y)
        else:
            chosen.append(min((x, y), key=str))
    return frozenset(chosen)


def cayley_arc(group, S: frozenset[object], u: object, v: object) -> bool:
    return group.mul(group.inv(u), v) in S


def tournament_fingerprint(group, S: frozenset[object]) -> dict[str, object]:
    elems = group.elements
    scores = [sum(1 for v in elems if v != u and cayley_arc(group, S, u, v)) for u in elems]
    c3 = 0
    for a, b, c in combinations(elems, 3):
        ab = cayley_arc(group, S, a, b)
        bc = cayley_arc(group, S, b, c)
        ca = cayley_arc(group, S, c, a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            c3 += 1
    out0 = tuple(sorted((group.label(x) for x in S)))
    orders = Counter(group.order(x) for x in elems)
    return {
        "scores": Counter(scores),
        "c3": c3,
        "root_star_size": len(S),
        "root_star_sample": out0[:8],
        "element_order_hist": dict(sorted(orders.items())),
        "max_element_order": max(orders),
        "has_order_n_element": group.n in orders,
    }


def cayley_report() -> list[str]:
    examples = [
        ("C7 interval circulant", CyclicGroup(7), "interval"),
        ("C21 interval circulant", CyclicGroup(21), "interval"),
        ("F21 alternating Cayley", Frobenius21(), "alternating"),
    ]
    lines = [
        "ABSTRACT VT TOURNAMENT CASCADE AUDIT",
        "  Cayley tournaments are vertex-transitive because left translations propagate the root star to every vertex.",
        "  The question is whether that propagation is cyclic-polygonal or group-relator driven.",
    ]
    for name, group, mode in examples:
        S = default_connection(group, mode)
        fp = tournament_fingerprint(group, S)
        lines += [
            f"  {name}:",
            f"    |G|={group.n}; root_star={fp['root_star_size']}; score_hist={dict(fp['scores'])}; c3={fp['c3']}",
            f"    element_order_hist={fp['element_order_hist']}; has_order_|G|_element={fp['has_order_n_element']}",
            f"    root_star_sample={fp['root_star_sample']}",
        ]
    lines += [
        "  reading: cyclic Cayley rigidity is a regular polygon spine; F21 is still vertex-transitive, but no element has order 21, so the local star cascades through a nonabelian relator mesh.",
        "  warning: vertex-transitive tournament => regular polygon is true only after adding a cyclic sharp-transitivity hypothesis.",
    ]
    return lines


@dataclass(frozen=True)
class Lens:
    name: str
    local_rigidity: int
    global_cascade: int
    primitive_polygon: int
    handles_fixed_point: int
    lrc_payload: int

    def key(self) -> tuple[int, int, int, int, int]:
        return (
            self.global_cascade,
            self.local_rigidity,
            self.lrc_payload,
            self.handles_fixed_point,
            self.primitive_polygon,
        )


LENSES = [
    Lens("cyclic_regular_polygon_spine", 5, 5, 5, 3, 4),
    Lens("dihedral_imprimitive_bracelet", 5, 4, 2, 4, 3),
    Lens("nonabelian_Cayley_VT_mesh", 4, 5, 1, 3, 4),
    Lens("source_root_threshold_payload", 4, 4, 1, 5, 5),
    Lens("local_distance_profile_only", 3, 1, 1, 4, 1),
    Lens("regular_score_sequence_only", 1, 1, 0, 1, 1),
]


def tournament_analysis_report() -> list[str]:
    names = [lens.name for lens in LENSES]
    edge = {}
    for a, b in combinations(LENSES, 2):
        winner, loser = (a, b) if a.key() > b.key() else (b, a)
        edge[(winner.name, loser.name)] = True
        edge[(loser.name, winner.name)] = False
    scores = Counter(sum(1 for v in names if u != v and edge[(u, v)]) for u in names)
    c3 = 0
    for triple in combinations(names, 3):
        for a, b, c in permutations(triple):
            if edge[(a, b)] and edge[(b, c)] and edge[(c, a)]:
                c3 += 1
                break
    hpaths = 0
    first = None
    for path in permutations(names):
        if all(edge[(path[i], path[i + 1])] for i in range(len(path) - 1)):
            hpaths += 1
            if first is None:
                first = path
    ranking = [lens.name for lens in sorted(LENSES, key=lambda x: x.key(), reverse=True)]
    return [
        "TOURNAMENT ANALYSIS",
        "  vertices: rigidity lenses, not tournament vertices.",
        "  observable: (global cascade, local rigidity, LRC payload, fixed-point handling, primitive polygonality).",
        f"  ranking: {ranking}",
        f"  score_histogram: {dict(sorted(scores.items()))}",
        f"  directed_3_cycles: {c3}",
        f"  Hamiltonian_path_count: {hpaths}",
        f"  first_Hamiltonian_path: {first}",
    ]


def main() -> None:
    print("S589 vertex-transitive trienerment / polygon rigidity audit")
    print()
    for block in (
        point_set_report(),
        [""],
        cayley_report(),
        [""],
        tournament_analysis_report(),
        [""],
        [
            "SYNTHESIS",
            "  cyclic VT point-set <=> regular polygon, once the transitive group is a single rotation cycle.",
            "  dihedral VT point-set => regular polygon or imprimitive bracelet; local views match, but block structure remains.",
            "  abstract VT tournament => Cayley cascade from one root star; it may be polygonal (cyclic) or relator-mesh (nonabelian).",
            "  LRC rigidity should keep the fixed basepoint/source payload: local isomorphism alone does not say which observer clock survives.",
        ],
    ):
        for line in block:
            print(line)


if __name__ == "__main__":
    main()
