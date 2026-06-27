#!/usr/bin/env python3
"""
HYP-3110 / S263 scout:
De Moivre quintic, Jacobi theta, and crystallographic sidecars for the LRC14
proof frontier.

This is not a proof of LRC14.  It is a proof-carrier audit:

* verify the De Moivre quintic fold as exact Laurent-polynomial cancellation;
* record the finite crystallographic counts requested by the prompt;
* run Tournament Analysis on sidecars, not runners or raw groups.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations
from typing import Dict, Iterable, List, Tuple

Monomial = Tuple[int, int]  # u_power, a_power
Laurent = Dict[Monomial, int]


def clean(poly: Laurent) -> Laurent:
    return {k: v for k, v in poly.items() if v}


def add(*polys: Laurent) -> Laurent:
    out: Counter[Monomial] = Counter()
    for poly in polys:
        out.update(poly)
    return clean(dict(out))


def scale(c: int, poly: Laurent) -> Laurent:
    return clean({m: c * v for m, v in poly.items()})


def mul(p: Laurent, q: Laurent) -> Laurent:
    out: Counter[Monomial] = Counter()
    for (up, ap), cv in p.items():
        for (uq, aq), dw in q.items():
            out[(up + uq, ap + aq)] += cv * dw
    return clean(dict(out))


def pow_poly(p: Laurent, n: int) -> Laurent:
    out: Laurent = {(0, 0): 1}
    for _ in range(n):
        out = mul(out, p)
    return out


def mul_a_power(power: int, p: Laurent) -> Laurent:
    return {(u, a + power): c for (u, a), c in p.items()}


def fmt_poly(poly: Laurent) -> str:
    if not poly:
        return "0"
    parts = []
    for (u, a), c in sorted(poly.items(), key=lambda item: (-item[0][0], -item[0][1])):
        coeff = "" if abs(c) == 1 else str(abs(c))
        sign = "-" if c < 0 else "+"
        u_part = "" if u == 0 else f"u^{u}"
        a_part = "" if a == 0 else f"a^{a}"
        body = coeff + "*".join(x for x in [u_part, a_part] if x)
        if not body:
            body = "1"
        parts.append((sign, body))
    first_sign, first_body = parts[0]
    rendered = ("-" if first_sign == "-" else "") + first_body
    for sign, body in parts[1:]:
        rendered += f" {sign} {body}"
    return rendered


def de_moivre_quintic_identity() -> Tuple[Laurent, Laurent, Laurent, int]:
    """Return lhs, rhs, difference, and pre-cancellation term count.

    For x = u - a/u:

      x^5 + 5 a x^3 + 5 a^2 x = u^5 - a^5/u^5.

    Thus x^5 + 5a x^3 + 5a^2 x + b = 0 reduces to
    y^2 + b y - a^5 = 0 with y = u^5.  HYP-3110 uses this only as a
    normal-form detector for finite-depth degree-5 cancellations.
    """

    x: Laurent = {(1, 0): 1, (-1, 1): -1}
    raw_terms: List[Laurent] = [
        pow_poly(x, 5),
        scale(5, mul_a_power(1, pow_poly(x, 3))),
        scale(5, mul_a_power(2, x)),
    ]
    lhs = add(*raw_terms)
    rhs: Laurent = {(5, 0): 1, (-5, 5): -1}
    diff = add(lhs, scale(-1, rhs))
    pre_cancel = sum(len(term) for term in raw_terms)
    return lhs, rhs, diff, pre_cancel


WALLPAPER_GROUPS = [
    "p1",
    "p2",
    "pm",
    "pg",
    "cm",
    "pmm",
    "pmg",
    "pgg",
    "cmm",
    "p4",
    "p4m",
    "p4g",
    "p3",
    "p3m1",
    "p31m",
    "p6",
    "p6m",
]

CRYSTAL_SYSTEM_SPACE_GROUP_COUNTS = {
    "triclinic": 2,
    "monoclinic": 13,
    "orthorhombic": 59,
    "tetragonal": 68,
    "trigonal": 25,
    "hexagonal": 27,
    "cubic": 36,
}

BRAVAIS_3D_COUNT = 14
JACOBI_THETA_CHANNELS = ["theta1_odd", "theta2_shifted", "theta3_lattice", "theta4_alternating"]


@dataclass(frozen=True)
class Carrier:
    name: str
    exactness: int
    frontier_adjacency: int
    preserved_payload: int
    destroyed_debt: int
    finite_catalog: int
    theorem_readiness: int
    role: str

    def primary_key(self) -> Tuple[int, int, int, int, int, int, str]:
        return (
            self.theorem_readiness,
            self.frontier_adjacency,
            self.preserved_payload,
            self.exactness,
            self.finite_catalog,
            -self.destroyed_debt,
            self.name,
        )

    def novelty_key(self) -> Tuple[int, int, int, int, int, int, str]:
        return (
            self.finite_catalog,
            self.exactness,
            self.preserved_payload,
            self.frontier_adjacency,
            self.theorem_readiness,
            -self.destroyed_debt,
            self.name,
        )


def carriers() -> List[Carrier]:
    return [
        Carrier(
            "finite_address_exit",
            exactness=5,
            frontier_adjacency=5,
            preserved_payload=5,
            destroyed_debt=0,
            finite_catalog=4,
            theorem_readiness=5,
            role="terminal HYP-3107 packet interface",
        ),
        Carrier(
            "observer_gluing_certificate",
            exactness=4,
            frontier_adjacency=5,
            preserved_payload=5,
            destroyed_debt=1,
            finite_catalog=3,
            theorem_readiness=4,
            role="chart overlap certificate producer",
        ),
        Carrier(
            "jacobi_theta_tail",
            exactness=4,
            frontier_adjacency=4,
            preserved_payload=5,
            destroyed_debt=2,
            finite_catalog=3,
            theorem_readiness=4,
            role="signed support-six residue-cusp tail",
        ),
        Carrier(
            "lee_yang_root_curve",
            exactness=4,
            frontier_adjacency=4,
            preserved_payload=4,
            destroyed_debt=2,
            finite_catalog=3,
            theorem_readiness=3,
            role="HYP-3109 zero-locus and ear sidecar",
        ),
        Carrier(
            "de_moivre_quintic_fold",
            exactness=5,
            frontier_adjacency=3,
            preserved_payload=3,
            destroyed_debt=2,
            finite_catalog=2,
            theorem_readiness=3,
            role="degree-5 finite-depth cancellation normal form",
        ),
        Carrier(
            "space_group_230_orbifold",
            exactness=5,
            frontier_adjacency=3,
            preserved_payload=4,
            destroyed_debt=3,
            finite_catalog=5,
            theorem_readiness=2,
            role="3D crystallographic quotient audit",
        ),
        Carrier(
            "wallpaper_17_orbifold",
            exactness=5,
            frontier_adjacency=2,
            preserved_payload=3,
            destroyed_debt=3,
            finite_catalog=5,
            theorem_readiness=2,
            role="2D crystallographic quotient audit",
        ),
        Carrier(
            "raw_scalar_shadow",
            exactness=2,
            frontier_adjacency=1,
            preserved_payload=1,
            destroyed_debt=5,
            finite_catalog=1,
            theorem_readiness=0,
            role="count/root/moment scalar after legal sidecars are dropped",
        ),
    ]


def tournament_edges(nodes: List[Carrier], key_name: str = "primary") -> Dict[Tuple[str, str], str]:
    edges: Dict[Tuple[str, str], str] = {}
    key = Carrier.primary_key if key_name == "primary" else Carrier.novelty_key
    for a, b in combinations(nodes, 2):
        winner, loser = (a, b) if key(a) > key(b) else (b, a)
        edges[(winner.name, loser.name)] = winner.name
    return edges


def out_neighbors(nodes: List[Carrier], edges: Dict[Tuple[str, str], str]) -> Dict[str, List[str]]:
    names = [n.name for n in nodes]
    out = {name: [] for name in names}
    for a, b in combinations(names, 2):
        if (a, b) in edges:
            out[a].append(b)
        elif (b, a) in edges:
            out[b].append(a)
        else:
            raise AssertionError((a, b))
    return out


def directed_three_cycles(nodes: List[Carrier], edges: Dict[Tuple[str, str], str]) -> int:
    names = [n.name for n in nodes]

    def beats(a: str, b: str) -> bool:
        return (a, b) in edges

    count = 0
    for a, b, c in combinations(names, 3):
        if (beats(a, b) and beats(b, c) and beats(c, a)) or (
            beats(a, c) and beats(c, b) and beats(b, a)
        ):
            count += 1
    return count


def scc_sizes(nodes: List[Carrier], edges: Dict[Tuple[str, str], str]) -> List[int]:
    out = out_neighbors(nodes, edges)
    names = [n.name for n in nodes]

    def reach(start: str) -> set[str]:
        seen = {start}
        stack = [start]
        while stack:
            cur = stack.pop()
            for nxt in out[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        return seen

    remaining = set(names)
    sizes = []
    while remaining:
        start = next(iter(remaining))
        comp = {v for v in remaining if start in reach(v) and v in reach(start)}
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes)


def hamiltonian_path_count(nodes: List[Carrier], edges: Dict[Tuple[str, str], str]) -> int:
    names = [n.name for n in nodes]

    def beats(a: str, b: str) -> bool:
        return (a, b) in edges

    count = 0
    # Eight vertices: 40320 permutations, acceptable and exact.
    for perm in permutations(names):
        if all(beats(perm[i], perm[i + 1]) for i in range(len(perm) - 1)):
            count += 1
    return count


def score_histogram(nodes: List[Carrier], edges: Dict[Tuple[str, str], str]) -> Dict[int, int]:
    out = out_neighbors(nodes, edges)
    return dict(sorted(Counter(len(v) for v in out.values()).items()))


def edge_flip_count(nodes: List[Carrier]) -> int:
    primary = tournament_edges(nodes, "primary")
    novelty = tournament_edges(nodes, "novelty")
    names = [n.name for n in nodes]
    flips = 0
    for a, b in combinations(names, 2):
        if ((a, b) in primary) != ((a, b) in novelty):
            flips += 1
    return flips


def main() -> None:
    lhs, rhs, diff, pre_cancel = de_moivre_quintic_identity()
    nodes = carriers()
    edges = tournament_edges(nodes)

    print("HYP-3110 / S263 De Moivre-Jacobi-crystallographic frontier scout")
    print("=" * 78)
    print()
    print("Status: LRC14 is NOT proved here.")
    print("Purpose: identify which new sidecars can be proof-carrier interfaces.")
    print()

    print("[1] De Moivre quintic fold")
    print(f"  x = u - a/u")
    print(f"  lhs = x^5 + 5*a*x^3 + 5*a^2*x = {fmt_poly(lhs)}")
    print(f"  rhs = {fmt_poly(rhs)}")
    print(f"  difference = {fmt_poly(diff)}")
    print(f"  identity_verified = {diff == {}}")
    print(f"  raw Laurent terms before cancellation = {pre_cancel}")
    print("  consequence: x^5+5a*x^3+5a^2*x+b=0 reduces to y^2+b*y-a^5=0 with y=u^5.")
    print("  LRC use: a degree-5 interaction can collapse only with this retained normal-form sidecar.")
    print()

    print("[2] Crystallographic finite catalogs")
    print(f"  wallpaper_group_count = {len(WALLPAPER_GROUPS)}")
    print(f"  wallpaper_groups = {', '.join(WALLPAPER_GROUPS)}")
    print(f"  space_group_3d_count = {sum(CRYSTAL_SYSTEM_SPACE_GROUP_COUNTS.values())}")
    print(f"  space_group_counts_by_crystal_system = {CRYSTAL_SYSTEM_SPACE_GROUP_COUNTS}")
    print(f"  bravais_3d_count = {BRAVAIS_3D_COUNT}")
    print(f"  jacobi_theta_channels = {', '.join(JACOBI_THETA_CHANNELS)}")
    print("  LRC use: these are finite quotient audits, not scalar proof shortcuts.")
    print()

    print("[3] Carrier sidecars")
    for c in nodes:
        print(
            "  {name:30s} exact={e} frontier={f} payload={p} debt={d} finite={fc} ready={r} :: {role}".format(
                name=c.name,
                e=c.exactness,
                f=c.frontier_adjacency,
                p=c.preserved_payload,
                d=c.destroyed_debt,
                fc=c.finite_catalog,
                r=c.theorem_readiness,
                role=c.role,
            )
        )
    print()

    print("[4] Tournament Analysis")
    print("  vertex_set = proof-carrier sidecars, not runners/arcs/group names")
    print("  pairwise_observable = theorem readiness, frontier adjacency, preserved payload, exactness, finite catalog, and destroyed-coordinate debt")
    print("  switch/gauge = primary proof-carrier lexicographic gauge")
    print(f"  score_hist = {score_histogram(nodes, edges)}")
    print(f"  directed_3cycles = {directed_three_cycles(nodes, edges)}")
    print(f"  scc_sizes = {scc_sizes(nodes, edges)}")
    print(f"  hamiltonian_path_count = {hamiltonian_path_count(nodes, edges)}")
    print(f"  edge_flips_vs_novelty_first_gauge = {edge_flip_count(nodes)}")
    tie_path = sorted(nodes, key=lambda c: c.primary_key(), reverse=True)
    print("  tie_hamiltonian_path =")
    print("    " + " > ".join(c.name for c in tie_path))
    print()

    print("[5] Proof-frontier consequences")
    print("  A. De Moivre is a detector for exact finite-depth degree-5 folds, not a general quintic/A5 wall.")
    print("  B. Jacobi theta belongs on the HYP-2614 signed residue-cusp tail after low-height wall deletion.")
    print("  C. Wallpaper/space groups give finite orbifold quotient audits; every use must name stabilizer, translation lattice, glide/screw/torsion, preserved predicate, and destroyed coordinate.")
    print("  D. HYP-3109 root curves stay upstream: crystallographic quotients must not forget the zero-collision sidecar.")
    print("  E. The next proof-facing artifact is a HYP-2963 residual row schema with theta/orbifold columns feeding ObserverGluingCertificate or FiniteAddressBranchPacket.")


if __name__ == "__main__":
    main()
