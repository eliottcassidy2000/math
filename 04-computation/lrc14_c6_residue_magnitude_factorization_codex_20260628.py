#!/usr/bin/env python3
"""Exact CRT/C6 residue-magnitude scout for the LRC14 AP/GW skeleton.

The script does not attempt to prove LRC14.  It records the exact algebraic
ledger behind the current proof route:

* binding residues are precisely the units modulo 14;
* their antipodal binder pairs form one C3 orbit inside the C6 unit group;
* nonunits split by CRT into even cover residues plus the ramified apex 7;
* the 12->24 Goddyn-Wong hinge is a magnitude-side doubling in the even cover
  branch, not a pure residue quotient.

Tournament Analysis uses proof carriers/sidecars as vertices, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations
from math import gcd


MOD = 14
UNITS = tuple(a for a in range(1, MOD) if gcd(a, MOD) == 1)
NONUNITS = tuple(a for a in range(1, MOD) if gcd(a, MOD) != 1)
EVEN_COVER = tuple(a for a in NONUNITS if a % 2 == 0)
APEX7 = 7
QR7 = {1, 2, 4}
NQR7 = {3, 5, 6}


def inv_mod(a: int, m: int = MOD) -> int:
    for x in range(1, m):
        if (a * x) % m == 1:
            return x
    raise ValueError(f"{a} is not invertible modulo {m}")


def complement_pair(a: int, m: int = MOD) -> tuple[int, int]:
    a %= m
    b = (-a) % m
    return tuple(sorted((a, b)))


def crt_label(a: int) -> tuple[int, int]:
    return (a % 2, a % 7)


def chi7(a: int) -> int:
    r = a % 7
    if r == 0:
        return 0
    return 1 if r in QR7 else -1


def v2(n: int) -> int:
    if n == 0:
        return 99
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def unit_power_cycle(generator: int = 3) -> list[int]:
    out = [1]
    x = 1
    while True:
        x = (x * generator) % MOD
        if x == 1:
            return out
        out.append(x)


def c3_pair_orbit(generator: int = 3) -> list[tuple[int, int]]:
    start = complement_pair(1)
    orbit = [start]
    cur = start
    while True:
        cur = complement_pair(cur[0] * generator)
        if cur == start:
            return orbit
        orbit.append(cur)


def unit_contact_rows() -> list[dict[str, object]]:
    rows = []
    for a in UNITS:
        inv = inv_mod(a)
        binder = complement_pair(inv)
        rows.append(
            {
                "unit_time": f"{a}/14",
                "a": a,
                "antipode": (-a) % MOD,
                "inverse": inv,
                "binding_pair": binder,
                "binding_pair_mod7": tuple(x % 7 for x in binder),
                "chi7_pair": tuple(chi7(x) for x in binder),
            }
        )
    return rows


def residue_layer_rows() -> list[dict[str, object]]:
    rows = []
    for a in range(1, MOD):
        if a in UNITS:
            layer = "unit_binding"
        elif a == APEX7:
            layer = "ramified_apex7_cover"
        elif a % 2 == 0:
            layer = "even_cover"
        else:
            layer = "other_nonunit"
        rows.append(
            {
                "residue": a,
                "crt_mod2_mod7": crt_label(a),
                "gcd_with_14": gcd(a, MOD),
                "chi7": chi7(a),
                "layer": layer,
            }
        )
    return rows


def hinge_12_to_24() -> dict[str, object]:
    a, b = 12, 24
    return {
        "from_speed": a,
        "to_speed": b,
        "from_mod14": a % MOD,
        "to_mod14": b % MOD,
        "from_crt": crt_label(a),
        "to_crt": crt_label(b),
        "from_v2": v2(a),
        "to_v2": v2(b),
        "v2_delta": v2(b) - v2(a),
        "from_odd_part": a >> v2(a),
        "to_odd_part": b >> v2(b),
        "from_layer": "even_cover",
        "to_layer": "even_cover",
        "mod14_preserved": (a % MOD) == (b % MOD),
        "unit_contact_killer": (b % MOD) == 0,
        "readout": (
            "doubling stays in the even covering branch and raises v2 by 1; "
            "it changes residue 12 mod 14 to 10 mod 14, so it is a magnitude "
            "sidecar rather than a residue-only quotient"
        ),
    }


def preserves_unit_binders(row: tuple[int, ...]) -> bool:
    row_set = set(x % MOD for x in row)
    return all(pair[0] in row_set and pair[1] in row_set for pair in c3_pair_orbit())


def named_rows() -> list[dict[str, object]]:
    ap = tuple(range(1, 14))
    gw = tuple(list(range(1, 12)) + [13, 24])
    rows = []
    for name, row in (("AP", ap), ("Goddyn-Wong_12_to_24", gw)):
        residues = tuple(sorted({x % MOD for x in row}))
        rows.append(
            {
                "name": name,
                "row": row,
                "residue_set_mod14": residues,
                "preserves_all_unit_binders": preserves_unit_binders(row),
                "contains_speed_0_mod14": any(x % MOD == 0 for x in row),
                "nonunit_residue_set": tuple(x for x in residues if x in NONUNITS),
            }
        )
    return rows


FEATURE_WEIGHTS = {
    "unit_contact_certificate": 12,
    "c6_unit_group": 10,
    "c3_transport": 10,
    "c2_conjugation": 7,
    "qr_nqr_character": 7,
    "seven_adic_residue": 9,
    "two_adic_magnitude": 9,
    "apex7_ramification": 8,
    "jacobsthal_hinge": 9,
    "covering_floor_route": 11,
    "offgrid_bulk_floor": 10,
    "finite_chamber_classifier": 11,
    "observability_glue": 12,
    "controlled_forgetting_guardrail": 10,
    "exact_runner_names": 2,
}


@dataclass(frozen=True)
class Carrier:
    name: str
    preserves: frozenset[str]
    destroys: frozenset[str]
    priority: int

    def score(self) -> int:
        kept = sum(FEATURE_WEIGHTS[x] for x in self.preserves)
        lost = sum(max(1, FEATURE_WEIGHTS[x] // 4) for x in self.destroys)
        return kept - lost


def tournament_carriers() -> list[Carrier]:
    return [
        Carrier(
            "observability_morse_glue",
            frozenset(
                {
                    "unit_contact_certificate",
                    "c6_unit_group",
                    "c3_transport",
                    "c2_conjugation",
                    "seven_adic_residue",
                    "two_adic_magnitude",
                    "apex7_ramification",
                    "jacobsthal_hinge",
                    "covering_floor_route",
                    "offgrid_bulk_floor",
                    "finite_chamber_classifier",
                    "observability_glue",
                    "controlled_forgetting_guardrail",
                }
            ),
            frozenset({"exact_runner_names"}),
            90,
        ),
        Carrier(
            "c6_unit_group",
            frozenset(
                {
                    "unit_contact_certificate",
                    "c6_unit_group",
                    "c3_transport",
                    "c2_conjugation",
                    "qr_nqr_character",
                    "seven_adic_residue",
                    "controlled_forgetting_guardrail",
                }
            ),
            frozenset(
                {
                    "two_adic_magnitude",
                    "apex7_ramification",
                    "jacobsthal_hinge",
                    "offgrid_bulk_floor",
                }
            ),
            80,
        ),
        Carrier(
            "seven_adic_residue_skeleton",
            frozenset(
                {
                    "unit_contact_certificate",
                    "c3_transport",
                    "c2_conjugation",
                    "qr_nqr_character",
                    "seven_adic_residue",
                    "apex7_ramification",
                }
            ),
            frozenset({"two_adic_magnitude", "jacobsthal_hinge", "offgrid_bulk_floor"}),
            75,
        ),
        Carrier(
            "c3_binding_slot_orbit",
            frozenset({"unit_contact_certificate", "c3_transport", "seven_adic_residue"}),
            frozenset(
                {
                    "c2_conjugation",
                    "qr_nqr_character",
                    "two_adic_magnitude",
                    "apex7_ramification",
                    "jacobsthal_hinge",
                    "offgrid_bulk_floor",
                }
            ),
            70,
        ),
        Carrier(
            "two_adic_magnitude_layer",
            frozenset(
                {
                    "two_adic_magnitude",
                    "jacobsthal_hinge",
                    "covering_floor_route",
                    "finite_chamber_classifier",
                    "controlled_forgetting_guardrail",
                }
            ),
            frozenset(
                {
                    "c3_transport",
                    "c2_conjugation",
                    "qr_nqr_character",
                    "seven_adic_residue",
                    "unit_contact_certificate",
                }
            ),
            65,
        ),
        Carrier(
            "jacobsthal_doubling_hinge_12_24",
            frozenset(
                {
                    "two_adic_magnitude",
                    "jacobsthal_hinge",
                    "covering_floor_route",
                    "finite_chamber_classifier",
                }
            ),
            frozenset(
                {
                    "c3_transport",
                    "c2_conjugation",
                    "qr_nqr_character",
                    "seven_adic_residue",
                    "exact_runner_names",
                }
            ),
            60,
        ),
        Carrier(
            "ramified_apex7_cover",
            frozenset(
                {
                    "seven_adic_residue",
                    "apex7_ramification",
                    "covering_floor_route",
                    "offgrid_bulk_floor",
                }
            ),
            frozenset(
                {
                    "unit_contact_certificate",
                    "c3_transport",
                    "c2_conjugation",
                    "two_adic_magnitude",
                    "jacobsthal_hinge",
                }
            ),
            55,
        ),
        Carrier(
            "c2_qr_nqr_conjugation",
            frozenset({"c2_conjugation", "qr_nqr_character", "seven_adic_residue"}),
            frozenset(
                {
                    "c3_transport",
                    "two_adic_magnitude",
                    "apex7_ramification",
                    "jacobsthal_hinge",
                    "offgrid_bulk_floor",
                }
            ),
            50,
        ),
        Carrier(
            "raw_runner_partition",
            frozenset({"exact_runner_names"}),
            frozenset(
                {
                    "unit_contact_certificate",
                    "c6_unit_group",
                    "c3_transport",
                    "c2_conjugation",
                    "qr_nqr_character",
                    "seven_adic_residue",
                    "two_adic_magnitude",
                    "apex7_ramification",
                    "jacobsthal_hinge",
                    "covering_floor_route",
                    "offgrid_bulk_floor",
                    "finite_chamber_classifier",
                    "observability_glue",
                    "controlled_forgetting_guardrail",
                }
            ),
            10,
        ),
    ]


def orient(a: Carrier, b: Carrier) -> tuple[str, str]:
    ka = (a.score(), -len(a.destroys), a.priority, a.name)
    kb = (b.score(), -len(b.destroys), b.priority, b.name)
    return (a.name, b.name) if ka >= kb else (b.name, a.name)


def adjacency(carriers: list[Carrier]) -> dict[str, set[str]]:
    adj = {c.name: set() for c in carriers}
    for a, b in combinations(carriers, 2):
        u, v = orient(a, b)
        adj[u].add(v)
    return adj


def directed_3cycles(adj: dict[str, set[str]]) -> list[tuple[str, str, str]]:
    cycles = []
    names = sorted(adj)
    for a, b, c in combinations(names, 3):
        edges = {
            (a, b): b in adj[a],
            (b, a): a in adj[b],
            (a, c): c in adj[a],
            (c, a): a in adj[c],
            (b, c): c in adj[b],
            (c, b): b in adj[c],
        }
        if edges[(a, b)] and edges[(b, c)] and edges[(c, a)]:
            cycles.append((a, b, c))
        elif edges[(a, c)] and edges[(c, b)] and edges[(b, a)]:
            cycles.append((a, c, b))
    return cycles


def scc_sizes(adj: dict[str, set[str]]) -> list[int]:
    names = list(adj)
    rev = {n: set() for n in names}
    for u, outs in adj.items():
        for v in outs:
            rev[v].add(u)

    def reach(start: str, graph: dict[str, set[str]]) -> set[str]:
        seen = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for v in graph[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        return seen

    remaining = set(names)
    sizes = []
    while remaining:
        n = next(iter(remaining))
        comp = reach(n, adj) & reach(n, rev)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_paths(adj: dict[str, set[str]]) -> list[tuple[str, ...]]:
    names = tuple(adj)
    out = []
    for perm in permutations(names):
        if all(perm[i + 1] in adj[perm[i]] for i in range(len(perm) - 1)):
            out.append(perm)
    return out


def print_table(title: str, rows: list[dict[str, object]], fields: list[str]) -> None:
    print(f"\n## {title}")
    print(" | ".join(fields))
    print(" | ".join("---" for _ in fields))
    for row in rows:
        print(" | ".join(str(row[f]) for f in fields))


def main() -> None:
    print("# LRC14 C6 Residue-Magnitude Factorization Scout")
    print()
    print(f"modulus = {MOD} = 2 * 7")
    print(f"units_mod14 = {UNITS}")
    print(f"nonunits_mod14 = {NONUNITS}")
    print(f"even_cover = {EVEN_COVER}")
    print(f"apex7 = {APEX7}")
    print(f"QR_mod7 = {tuple(sorted(QR7))}, NQR_mod7 = {tuple(sorted(NQR7))}")
    print(f"C6 generator 3 mod 14 powers = {tuple(unit_power_cycle(3))}")
    print(f"C3 antipodal binder orbit under *3 = {tuple(c3_pair_orbit(3))}")
    print()
    print("Cyclotomic readout:")
    print("- Gal(Q(zeta_7)/Q) ~= (Z/7)^* ~= (Z/14)^* is C6.")
    print("- quotient by conjugation +/-1 gives the C3 action on binder slots.")
    print("- the order-3 subgroup fixes the quadratic field Q(sqrt(-7)).")
    print("- therefore the C3 slot action and Q(sqrt(-7)) are complementary")
    print("  projections of the same C6 cyclotomic package, not rival scalars.")

    print_table(
        "CRT residue layer",
        residue_layer_rows(),
        ["residue", "crt_mod2_mod7", "gcd_with_14", "chi7", "layer"],
    )

    print_table(
        "Unit contact/binder ledger",
        unit_contact_rows(),
        [
            "unit_time",
            "antipode",
            "inverse",
            "binding_pair",
            "binding_pair_mod7",
            "chi7_pair",
        ],
    )

    print_table(
        "Named tight rows",
        named_rows(),
        [
            "name",
            "residue_set_mod14",
            "preserves_all_unit_binders",
            "contains_speed_0_mod14",
            "nonunit_residue_set",
        ],
    )

    print("\n## 12->24 hinge")
    for k, v in hinge_12_to_24().items():
        print(f"{k}: {v}")

    print("\n## Proof-route theorem target")
    targets = [
        "prove one binding-pair lemma, then C3-transport it to all three slots",
        "split the covering layer into even cover residues and the ramified apex 7",
        "classify magnitude-side flex inside the covering branch; 12->24 is the live hinge",
        "glue residue and magnitude layers through the HYP-3300 observability/Morse packet",
        "route failures to strict off-grid bulk, covering floor, finite chamber discharge, or named debt",
    ]
    for i, target in enumerate(targets, 1):
        print(f"{i}. {target}")

    carriers = tournament_carriers()
    adj = adjacency(carriers)
    scores = {c.name: c.score() for c in carriers}
    score_hist = dict(sorted(Counter(scores.values()).items()))
    cycles = directed_3cycles(adj)
    scc = scc_sizes(adj)
    paths = hamiltonian_paths(adj)
    priority_path = max(
        paths,
        key=lambda p: tuple((scores[name], -len(next(c for c in carriers if c.name == name).destroys)) for name in p),
    )

    print("\n## Tournament Analysis")
    print("vertices_are = proof carriers / sidecar columns, not runners")
    print("pairwise_observable = retained LRC predicate coordinates minus destroyed sidecars")
    print("switch_gauge = higher weighted retained payload; ties by fewer destroyed sidecars and priority")
    print(f"vertex_count = {len(carriers)}")
    print(f"score_hist = {score_hist}")
    print(f"directed_3cycles = {len(cycles)}")
    print(f"scc_sizes = {scc}")
    print(f"hamiltonian_path_count = {len(paths)}")
    print("priority_hamiltonian_path =")
    for name in priority_path:
        carrier = next(c for c in carriers if c.name == name)
        print(
            f"  {name}: score={carrier.score()}, "
            f"keeps={sorted(carrier.preserves)}, destroys={sorted(carrier.destroys)}"
        )

    print("\n## Assumption challenge")
    print("Candidate vertices considered:")
    for item in (
        "runners",
        "gaps",
        "fixed circle sections",
        "section boundaries",
        "wall-crossing events",
        "residue/CRT fibers",
        "cover arcs",
        "Fourier/cyclotomic modes",
        "matroid contact circuits",
        "proof obligations / sidecar columns",
    ):
        print(f"- {item}")
    print(
        "Chosen vertices are proof obligations because the quotient preserves "
        "the LRC predicate only after unit contacts, covering floor, 2-adic "
        "magnitude, apex ramification, and hinge data are all either retained "
        "or explicitly discharged."
    )


if __name__ == "__main__":
    main()
