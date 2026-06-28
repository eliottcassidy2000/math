#!/usr/bin/env python3
"""HYP-3311: CRT/Galois sidecar audit of the LRC(14) census picture.

This extends the S84 note.  The executable point is deliberately small:

* nonzero classes mod 14 split as U union 2U union {7};
* U carries the cyclotomic C6 = C2 x C3 field data;
* quotienting U by +/- gives the C3 real-cubic binding-pair skeleton;
* the quadratic Q(sqrt(-7)) character cuts across every binding pair;
* 2U is the even covering layer, and 7 is the ramified covering singleton.

The script records the resulting proof-packet and a Tournament Analysis over
carriers.  It is not an LRC proof; it localizes what a proof must still show.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations, permutations
from math import gcd


MOD = 14
P = 7


def legendre7(a: int) -> int:
    r = a % P
    if r == 0:
        return 0
    return 1 if r in {1, 2, 4} else -1


def unit_pair_id(s: int) -> tuple[int, int] | None:
    if gcd(s, MOD) != 1:
        return None
    return tuple(sorted((s % MOD, (-s) % MOD)))


def c3_pair_order() -> list[tuple[int, int]]:
    start = (1, 13)
    out = []
    current = start
    for _ in range(3):
        out.append(current)
        current = tuple(sorted(((3 * current[0]) % MOD, (3 * current[1]) % MOD)))
    return out


PAIR_ORDER = c3_pair_order()
PAIR_INDEX = {pair: idx for idx, pair in enumerate(PAIR_ORDER)}


def unit_pair_slot(s: int) -> int | None:
    pair = unit_pair_id(s)
    if pair is None:
        return None
    return PAIR_INDEX[pair]


def crt_tag(s: int) -> str:
    if gcd(s, MOD) == 1:
        return "unit / binding / 7-residue skeleton"
    if s % P == 0:
        return "apex7 / ramified covering singleton"
    if s % 2 == 0:
        return "even / 2-adic magnitude shadow"
    return "nonunit"


def c6_table() -> list[dict[str, object]]:
    rows = []
    for a in range(1, P):
        pair = tuple(sorted((a, (-a) % P)))
        rows.append(
            {
                "sigma_a": a,
                "quadratic_C2": legendre7(a),
                "real_cubic_C3_pair": pair,
                "fixes_Qsqrt_minus7": a in {1, 2, 4},
                "fixes_real_cubic": a in {1, 6},
            }
        )
    return rows


def layer_table() -> list[dict[str, object]]:
    rows = []
    for s in range(1, MOD):
        rows.append(
            {
                "s": s,
                "mod2": s % 2,
                "mod7": s % 7,
                "legendre7": legendre7(s),
                "unit_pair": unit_pair_id(s),
                "c3_slot": unit_pair_slot(s),
                "double_of_unit": None if s % 2 else (s // 2 if s // 2 in range(1, 7) else (s + MOD) // 2),
                "tag": crt_tag(s),
            }
        )
    return rows


def unit_to_even_shadow() -> dict[int, int]:
    return {u: (2 * u) % MOD for u in range(1, MOD) if gcd(u, MOD) == 1}


def check_factorization() -> dict[str, object]:
    units = [s for s in range(1, MOD) if gcd(s, MOD) == 1]
    evens = [s for s in range(1, MOD) if s % 2 == 0]
    nonunits = [s for s in range(1, MOD) if gcd(s, MOD) != 1]
    shadow = unit_to_even_shadow()
    return {
        "units": units,
        "evens": evens,
        "apex": 7,
        "covering": nonunits,
        "unit_to_even_shadow": shadow,
        "shadow_is_bijection_to_evens": sorted(shadow.values()) == evens,
        "nonzero_partition": sorted(units + evens + [7]) == list(range(1, MOD)),
        "c3_pair_order": PAIR_ORDER,
        "each_pair_has_one_QR_and_one_NQR": all(
            sorted(legendre7(x) for x in pair) == [-1, 1]
            for pair in PAIR_ORDER
        ),
    }


@dataclass(frozen=True)
class Carrier:
    name: str
    preserves_lrc_predicate: int
    detects_binding_layer: int
    detects_covering_layer: int
    field_lattice_payload: int
    rigidity_target: int
    proof_readiness: int

    def score(self) -> int:
        return (
            3 * self.preserves_lrc_predicate
            + 3 * self.detects_binding_layer
            + 3 * self.detects_covering_layer
            + 2 * self.field_lattice_payload
            + 2 * self.rigidity_target
            + 2 * self.proof_readiness
        )


def tournament_fingerprint(carriers: list[Carrier]) -> dict[str, object]:
    def beats(a: Carrier, b: Carrier) -> bool:
        if a.score() != b.score():
            return a.score() > b.score()
        return a.name < b.name

    cycle_count = 0
    for a, b, c in combinations(carriers, 3):
        ab = beats(a, b)
        bc = beats(b, c)
        ca = beats(c, a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            cycle_count += 1

    hpaths = 0
    for order in permutations(carriers):
        if all(beats(order[i], order[i + 1]) for i in range(len(order) - 1)):
            hpaths += 1

    ordered = sorted(carriers, key=lambda item: (-item.score(), item.name))
    return {
        "vertices": len(carriers),
        "score_hist": dict(sorted(Counter(c.score() for c in carriers).items())),
        "directed_3cycles": cycle_count,
        "hamiltonian_path_count": hpaths,
        "selected_path": [c.name for c in ordered],
    }


def proof_carrier_tournament() -> dict[str, object]:
    carriers = [
        Carrier("crt_packet_U_plus_2U_plus_ramified7", 5, 5, 5, 4, 5, 5),
        Carrier("c6_subfield_lattice_Qzeta7", 4, 5, 3, 5, 4, 4),
        Carrier("c3_real_cubic_binding_pair_orbit", 4, 5, 2, 5, 4, 5),
        Carrier("quadratic_Qsqrt_minus7_transverse_character", 3, 4, 3, 5, 3, 4),
        Carrier("two_adic_unit_to_even_shadow_2U", 4, 3, 5, 2, 5, 4),
        Carrier("apex7_ramified_covering_switch", 4, 3, 5, 4, 4, 4),
        Carrier("hinge_12_to_24_height_lift", 3, 2, 5, 1, 5, 4),
        Carrier("raw_mod14_residue_table", 2, 3, 3, 2, 2, 5),
        Carrier("slogan_binding_units_covering_evens", 1, 2, 2, 1, 1, 5),
    ]
    return tournament_fingerprint(carriers)


def main() -> None:
    facts = check_factorization()
    print("HYP-3311 CRT/Galois sidecar audit of LRC14 census")
    print("source=mac-mini-2026-06-28-S84 + codex extension")
    print()

    print("1. CRT FACTORIZATION OF NONZERO CLASSES MOD 14")
    print(f"units U={facts['units']}")
    print(f"evens 2U={facts['evens']}")
    print(f"apex={facts['apex']}")
    print(f"covering nonunits={facts['covering']}")
    print(f"unit_to_even_shadow={facts['unit_to_even_shadow']}")
    print(f"shadow_is_bijection_to_evens={facts['shadow_is_bijection_to_evens']}")
    print(f"nonzero_partition_U_2U_apex={facts['nonzero_partition']}")
    print()

    print("2. CRT TABLE")
    for row in layer_table():
        pair = "-" if row["unit_pair"] is None else row["unit_pair"]
        slot = "-" if row["c3_slot"] is None else row["c3_slot"]
        print(
            f"  s={row['s']:>2}: (mod2={row['mod2']}, mod7={row['mod7']}) "
            f"chi7={row['legendre7']:>2} pair={pair} c3_slot={slot} :: {row['tag']}"
        )
    print()

    print("3. C6 = C2 x C3 SUBFIELD LATTICE INSIDE Q(zeta_7)")
    print("C2 coordinate=quadratic character chi7, fixed field Q(sqrt(-7)).")
    print("C3 coordinate=quotient by complex conjugation +/-1, fixed field Q(cos(2pi/7)).")
    for row in c6_table():
        print(
            f"  sigma_{row['sigma_a']}: chi7={row['quadratic_C2']:>2} "
            f"real_pair={row['real_cubic_C3_pair']} "
            f"fixes_Qsqrt_minus7={row['fixes_Qsqrt_minus7']} "
            f"fixes_real_cubic={row['fixes_real_cubic']}"
        )
    print(f"c3_pair_order_under_times3_mod14={facts['c3_pair_order']}")
    print(f"each_binding_pair_has_one_QR_and_one_NQR={facts['each_pair_has_one_QR_and_one_NQR']}")
    print("readout=the C3 binding-pair quotient and the Q(sqrt(-7)) quadratic character are transverse sidecars.")
    print()

    print("4. TWO-PRIME PROOF PACKET")
    print("binding_layer = U, quotient by +/- gives the C3 orbit seeded by HYP-2909")
    print("covering_layer = 2U plus ramified apex7; 2U is the 2-adic shadow of U")
    print("hinge = 12 -> 24 is a height/doubling lift: it contains the old 12-kill lattice and adds interleaved 24-kills")
    print("guardrail = residue equality alone is false proof currency; same-residue height lifts can be strict positive rows")
    print("proof_target = prove unit-contact rigidity, then prove the covering-flex manifold has only AP/GW integer points")
    print()

    print("5. TOURNAMENT ANALYSIS OVER PROOF CARRIERS")
    fp = proof_carrier_tournament()
    print("pairwise_observable=preserved LRC predicate, binding detection, covering detection, field-lattice payload, rigidity target, proof readiness")
    print("binary_gauge=A->B iff weighted proof-carrier score(A)>score(B); tie Hamiltonian path is lexicographic")
    for key, value in fp.items():
        if key == "selected_path":
            print(f"{key}={' -> '.join(value)}")
        else:
            print(f"{key}={value}")

    print()
    print("6. NEXT PROOF OBLIGATIONS")
    print("A. Formalize unit-contact rigidity: tight rows retaining all six unit contacts force the C3 unit skeleton.")
    print("B. Add the transverse Q(sqrt(-7)) character as a sidecar, not a replacement for the C3 quotient.")
    print("C. Prove the covering layer 2U+{7} has a one-dimensional height flex and that its only integer tight lift is 12->24.")
    print("D. Combine with HYP-3265/HYP-3300: killed unit contacts must route to off-unit covering/Morse chamber witnesses.")


if __name__ == "__main__":
    main()
