#!/usr/bin/env python3
"""Exact companion for the THM-2887 Arf/Q8 semantic lift.

The THM-2884 seam selector is the nonzero indicator on V4.  This script
identifies it as the Arf-one quadratic refinement, constructs its canonical
central extension Q8, compares it with THM-2779's Arf-zero D8 boundary,
and tests the distinguished carry triangle.

The Q8 lift is a structured candidate for the missing origin sign.  No
physical lift of a horn sheet to Q8, no q11/q7 current, and no LRC(14)
conclusion is constructed.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import permutations, product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
THEOREMS = ROOT / "01-canon/theorems"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINNED = {
    THEOREMS
    / "THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate.md":
        "21ac7ec9b19b8ed0abe4763e0b7e13ebc1e5eb776168c8e0088540f29daabccb",
    COMP / "lrc14_bockstein_symplectic_heisenberg_gate_thm2779.py":
        "4c6a58c80ddd4be0fd9bdd297b310df054bbc08996eb223d519d3cce6b8ed13a",
    RESULTS / "lrc14_bockstein_symplectic_heisenberg_gate_thm2779.out":
        "f7c96259777a3ab4a5e46cac8666181ae77a3be2e440cee8785997507706791a",
    COMP / "lrc14_macro_semantic_diagonal_horn_carrier_thm2884.py":
        "b739be20e741d5c061e0febcc8aef9b0f58f4ae8a648aa803610e0dad991929f",
    RESULTS / "lrc14_macro_semantic_diagonal_horn_carrier_thm2884.out":
        "8c3829b1052a641ca08a5e5bda86d9d5e8bd1584f5b2911c57c9fad9da41d4b6",
    COMP / "lrc14_stepped_origin_v4_provenance_transport_thm2886.py":
        "1f2cbffb8151c0c74bb22beff58e0ace5715eba3d4a9afc59481e7ab5e0d6dc9",
    RESULTS / "lrc14_stepped_origin_v4_provenance_transport_thm2886.out":
        "60ba517f1b8fa92fdaeba48b68d10da97f23f00a9fb1526de9e39f1311e49929",
}
for path, expected in PINNED.items():
    require(lf_hash(path) == expected, f"pinned dependency changed: {path.name}")


F2 = (0, 1)
V = tuple(product(F2, repeat=2))
ZERO = (0, 0)
QA = (1, 0)
QB = (0, 1)
QAB = (1, 1)


def add(left, right):
    return left[0] ^ right[0], left[1] ^ right[1]


def determinant(left, right):
    return (left[0] & right[1]) ^ (left[1] & right[0])


def q_d8(vector):
    """THM-2779's Arf-zero refinement."""
    return vector[0] & vector[1]


def q_q8(vector):
    """THM-2884's nonzero indicator / Arf-one refinement."""
    left, right = vector
    return left ^ right ^ (left & right)


def polar(quadratic, left, right):
    return quadratic(left) ^ quadratic(right) ^ quadratic(add(left, right))


def cocycle_d8(left, right):
    return left[1] & right[0]


def cocycle_q8_forward(left, right):
    """Q8 section gauge whose QA-then-QB product pays the central sign."""
    return (
        (left[0] & right[0])
        ^ (left[1] & right[1])
        ^ (left[0] & right[1])
    )


def cocycle_q8_reverse(left, right):
    """Gauge-equivalent Q8 section with the opposite ordered triangle."""
    return (
        (left[0] & right[0])
        ^ (left[1] & right[1])
        ^ (left[1] & right[0])
    )


def multiply(left, right, cocycle):
    vector_left = left[:2]
    vector_right = right[:2]
    vector = add(vector_left, vector_right)
    central = left[2] ^ right[2] ^ cocycle(vector_left, vector_right)
    return vector[0], vector[1], central


def power(element, exponent, cocycle):
    result = (0, 0, 0)
    for _index in range(exponent):
        result = multiply(result, element, cocycle)
    return result


def order(element, cocycle):
    for candidate in range(1, 9):
        if power(element, candidate, cocycle) == (0, 0, 0):
            return candidate
    raise RuntimeError("order exceeded eight")


def inverse(element, cocycle):
    return power(element, order(element, cocycle) - 1, cocycle)


def compose_permutations(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def permutation_order(permutation):
    identity = tuple(range(len(permutation)))
    current = identity
    for candidate in range(1, 25):
        current = compose_permutations(permutation, current)
        if current == identity:
            return candidate
    raise RuntimeError("automorphism order exceeded 24")


def automorphisms(elements, cocycle):
    index = {element: position for position, element in enumerate(elements)}
    table = tuple(
        tuple(index[multiply(left, right, cocycle)] for right in elements)
        for left in elements
    )
    identity_index = index[(0, 0, 0)]
    rows = []
    for tail in permutations(
        tuple(position for position in range(len(elements))
              if position != identity_index)
    ):
        permutation = [None] * len(elements)
        permutation[identity_index] = identity_index
        for position, image in zip(
            (
                position for position in range(len(elements))
                if position != identity_index
            ),
            tail,
        ):
            permutation[position] = image
        permutation = tuple(permutation)
        if all(
            permutation[table[left][right]]
            == table[permutation[left]][permutation[right]]
            for left in range(len(elements))
            for right in range(len(elements))
        ):
            rows.append(permutation)
    return tuple(rows), index, table


def main() -> None:
    # The two quadratic refinements have the same determinant polarization
    # and opposite Arf invariant.
    d8_values = tuple(q_d8(vector) for vector in V)
    q8_values = tuple(q_q8(vector) for vector in V)
    require(
        d8_values == (0, 0, 0, 1)
        and q8_values == (0, 1, 1, 1)
        and all(
            polar(q_d8, left, right) == determinant(left, right)
            and polar(q_q8, left, right) == determinant(left, right)
            for left in V for right in V
        )
        and (q_d8(QA) & q_d8(QB)) == 0
        and (q_q8(QA) & q_q8(QB)) == 1,
        "Arf-zero/Arf-one quadratic anatomy changed",
    )

    # Each multiplication law is associative.  Squares in the central
    # extension recover the corresponding quadratic refinement.
    elements = tuple(product(F2, repeat=3))
    associativity_checks = {}
    for name, quadratic, cocycle in (
        ("D8", q_d8, cocycle_d8),
        ("Q8", q_q8, cocycle_q8_forward),
    ):
        checks = 0
        for left in elements:
            for middle in elements:
                for right in elements:
                    require(
                        multiply(
                            multiply(left, middle, cocycle), right, cocycle
                        )
                        == multiply(
                            left, multiply(middle, right, cocycle), cocycle
                        ),
                        f"{name} multiplication lost associativity",
                    )
                    checks += 1
        require(
            all(
                power((vector[0], vector[1], central), 2, cocycle)
                == (0, 0, quadratic(vector))
                for vector in V for central in F2
            ),
            f"{name} square refinement changed",
        )
        associativity_checks[name] = checks

    d8_order_census = Counter(order(element, cocycle_d8) for element in elements)
    q8_order_census = Counter(
        order(element, cocycle_q8_forward) for element in elements
    )
    require(
        d8_order_census == {1: 1, 2: 5, 4: 2}
        and q8_order_census == {1: 1, 2: 1, 4: 6},
        "D8/Q8 order census changed",
    )

    # The two Q8 cocycles differ by the coboundary of q_Q8.  The forward
    # gauge makes the distinguished ordered semantic triangle pay -1.
    gauge_checks = 0
    for left in V:
        for right in V:
            delta_q = polar(q_q8, left, right)
            require(
                cocycle_q8_reverse(left, right)
                == cocycle_q8_forward(left, right) ^ delta_q,
                "Q8 section gauges ceased to differ by delta(q)",
            )
            gauge_checks += 1
    lifted_qa = (*QA, 0)
    lifted_qb = (*QB, 0)
    lifted_qab_negative = (*QAB, 1)
    require(
        multiply(lifted_qa, lifted_qb, cocycle_q8_forward)
        == lifted_qab_negative
        and multiply(lifted_qb, lifted_qa, cocycle_q8_forward)
        == (*QAB, 0),
        "oriented quaternion triangle changed",
    )

    # Exhaust all normalized F2-valued 2-cocycles on V4.  Among those with
    # determinant commutator there are four quadratic refinements and two
    # normalized cocycle-table representatives per refinement: three
    # Arf-zero D8 forms and one Arf-one Q8 form.  The eight normalized
    # 1-cochains induce only two distinct coboundaries, because the four
    # linear characters form the kernel.  Hence the all-nonzero square
    # law singles out the Q8 cohomology class uniquely.
    nonzero_vectors = V[1:]
    normalized_pairs = tuple(
        (left, right)
        for left in nonzero_vectors for right in nonzero_vectors
    )
    normalized_cocycles = []
    for mask in range(1 << len(normalized_pairs)):
        row = {
            pair: (mask >> index) & 1
            for index, pair in enumerate(normalized_pairs)
        }

        def candidate_cocycle(left, right):
            if left == ZERO or right == ZERO:
                return 0
            return row[left, right]

        if all(
            candidate_cocycle(middle, right)
            ^ candidate_cocycle(add(left, middle), right)
            ^ candidate_cocycle(left, add(middle, right))
            ^ candidate_cocycle(left, middle)
            == 0
            for left in V for middle in V for right in V
        ):
            quadratic_values = tuple(
                candidate_cocycle(vector, vector) for vector in V
            )
            determinant_commutator = all(
                candidate_cocycle(left, right)
                ^ candidate_cocycle(right, left)
                == determinant(left, right)
                for left in V for right in V
            )
            normalized_cocycles.append(
                (
                    tuple(row[pair] for pair in normalized_pairs),
                    quadratic_values,
                    determinant_commutator,
                )
            )
    determinant_quadratic_census = Counter(
        quadratic_values
        for _row, quadratic_values, commutator in normalized_cocycles
        if commutator
    )
    normalized_coboundaries = set()
    for mask in range(1 << len(nonzero_vectors)):
        one_cochain = {
            ZERO: 0,
            **{
                vector: (mask >> index) & 1
                for index, vector in enumerate(nonzero_vectors)
            },
        }
        normalized_coboundaries.add(tuple(
            one_cochain[left]
            ^ one_cochain[right]
            ^ one_cochain[add(left, right)]
            for left, right in normalized_pairs
        ))
    require(
        len(normalized_cocycles) == 16
        and sum(
            int(commutator)
            for _row, _quadratic, commutator in normalized_cocycles
        ) == 8
        and determinant_quadratic_census
        == {
            (0, 1, 0, 0): 2,
            (0, 0, 1, 0): 2,
            (0, 0, 0, 1): 2,
            (0, 1, 1, 1): 2,
        }
        and len(normalized_coboundaries) == 2,
        "normalized V4 central-extension classification changed",
    )

    # The local carry event equals the gauge-invariant polarization of the
    # two semantic edge directions.
    semantic_triangle_polarization = polar(q_q8, QA, QB)
    carry_triangle = (8 + 9) // 13
    require(
        add(QA, QB) == QAB
        and semantic_triangle_polarization == carry_triangle == 1,
        "semantic Arf polarization stopped matching the carry triangle",
    )

    # The source mark resolves the apparent S3-orientation obstruction.
    # Conjugation by the distinguished lift QB=j fixes the Q0/QB axes and
    # pays the central sign on QA/QAB.  Its sign character is exactly
    # h |-> det(QB,h), whose values on the transported q3/q11/q7 seam
    # (Q0,QA,QAB) are THM-2886's (0,1,1).
    source_conjugation_sign = {}
    lifted_qb_inverse = inverse(lifted_qb, cocycle_q8_forward)
    for vector in V:
        lifted_vector = (*vector, 0)
        conjugated = multiply(
            multiply(lifted_qb, lifted_vector, cocycle_q8_forward),
            lifted_qb_inverse,
            cocycle_q8_forward,
        )
        require(
            conjugated[:2] == vector
            and conjugated[2] == determinant(QB, vector),
            "fixed-source Q8 conjugation stopped realizing det(QB,h)",
        )
        source_conjugation_sign[vector] = conjugated[2]
    seam_source_sign = tuple(
        source_conjugation_sign[vector] for vector in (ZERO, QA, QAB)
    )
    require(
        seam_source_sign == (0, 1, 1)
        and source_conjugation_sign[QB] == 0,
        "fixed-source seam orientation stopped matching THM-2886",
    )

    # This coincidence cannot be globalized by a base-independent map from
    # C13 increments to H.  Modulo two, the carry cocycle is the unique
    # normalized coboundary of r(h)=h mod2.  Since q_Q8 is anisotropic,
    # any map ell with delta(q ell)=carry must send every even h to zero,
    # contradicting ell(8)=QA and ell(4)=QAB on the semantic triangle.
    def carry(h, k):
        return (h + k) // 13

    def delta(values, h, k):
        return values[h] ^ values[k] ^ values[(h + k) % 13]

    parity_values = tuple(h % 2 for h in range(13))
    require(
        all(
            delta(parity_values, h, k) == carry(h, k)
            for h in range(13) for k in range(13)
        ),
        "carry parity stopped being delta(h mod2)",
    )
    carry_primitives = []
    for mask in range(1 << 12):
        values = (0,) + tuple(
            (mask >> (h - 1)) & 1 for h in range(1, 13)
        )
        if all(
            delta(values, h, k) == carry(h, k)
            for h in range(13) for k in range(13)
        ):
            carry_primitives.append(values)
    require(
        carry_primitives == [parity_values]
        and q_q8(QA) != parity_values[8]
        and q_q8(QAB) != parity_values[4],
        "global carry/semantic quadratic no-go changed",
    )

    # Aut(Q8) has the exact S4/V4=S3 anatomy recurrent in the repo.
    q8_automorphisms, element_index, _table = automorphisms(
        elements, cocycle_q8_forward
    )
    q8_automorphism_orders = Counter(
        permutation_order(permutation) for permutation in q8_automorphisms
    )
    quotient_actions = set()
    kernel_actions = []
    for permutation in q8_automorphisms:
        action = tuple(
            elements[permutation[element_index[(*vector, 0)]]][:2]
            for vector in V
        )
        quotient_actions.add(action)
        if action == V:
            kernel_actions.append(permutation)

    inner_actions = set()
    for element in elements:
        inverse_element = inverse(element, cocycle_q8_forward)
        permutation = tuple(
            element_index[
                multiply(
                    multiply(element, target, cocycle_q8_forward),
                    inverse_element,
                    cocycle_q8_forward,
                )
            ]
            for target in elements
        )
        inner_actions.add(permutation)
    require(
        len(q8_automorphisms) == 24
        and q8_automorphism_orders == {1: 1, 2: 9, 3: 8, 4: 6}
        and len(quotient_actions) == 6
        and len(kernel_actions) == 4
        and set(kernel_actions) == inner_actions,
        "Aut(Q8)=S4, Inn(Q8)=V4, Out(Q8)=S3 anatomy changed",
    )

    # Q8's central sign cannot be seen faithfully on fewer than eight
    # points.  Every nonidentity cyclic subgroup contains the center, so
    # every nonfree orbit kills the central sign.
    q8_center = (0, 0, 1)
    require(
        all(
            q8_center in {
                power(element, exponent, cocycle_q8_forward)
                for exponent in range(order(element, cocycle_q8_forward))
            }
            for element in elements
            if element != (0, 0, 0)
        ),
        "a nontrivial Q8 subgroup stopped containing the center",
    )
    q8_minimal_faithful_permutation_degree = 8

    # There is no S3-equivariant choice of endpoint sign on the three
    # semantic directions when S3 acts on signs through its sign character.
    nonzero = (QA, QB, QAB)
    nonzero_index = {vector: position for position, vector in enumerate(nonzero)}
    quotient_permutations = set()
    for action in quotient_actions:
        quotient_permutations.add(
            tuple(nonzero_index[action[V.index(vector)]] for vector in nonzero)
        )
    require(len(quotient_permutations) == 6, "nonzero H action stopped being S3")
    source_character_covariance_checks = 0
    source_stabilizer_size = 0
    for action in quotient_actions:
        acted_qb = action[V.index(QB)]
        if acted_qb == QB:
            source_stabilizer_size += 1
        for source in V:
            for vector in V:
                require(
                    determinant(
                        action[V.index(source)], action[V.index(vector)]
                    )
                    == determinant(source, vector),
                    "pointed determinant character lost S3 covariance",
                )
                source_character_covariance_checks += 1
    require(
        source_stabilizer_size == 2,
        "the stabilizer of the marked source QB stopped having order two",
    )

    def sign_of_three(permutation):
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(3) for right in range(left + 1, 3)
        )
        return inversions % 2

    sign_equivariant_maps = []
    trivial_equivariant_maps = []
    for values in product(F2, repeat=3):
        if all(
            values[permutation[index]]
            == (values[index] ^ sign_of_three(permutation))
            for permutation in quotient_permutations
            for index in range(3)
        ):
            sign_equivariant_maps.append(values)
        if all(
            values[permutation[index]] == values[index]
            for permutation in quotient_permutations
            for index in range(3)
        ):
            trivial_equivariant_maps.append(values)
    require(
        sign_equivariant_maps == []
        and trivial_equivariant_maps == [(0, 0, 0), (1, 1, 1)],
        "S3-equivariant semantic orientation boundary changed",
    )

    # No order-thirteen action can twist Q8: Aut(Q8) has order 24.
    # Thus the semidirect-product direction C169 -> Aut(Q8) is direct.
    # This does not classify arbitrary extensions or all groups of order 1352.
    c13_eligible_q8_automorphisms = tuple(
        permutation
        for permutation in q8_automorphisms
        if 169 % permutation_order(permutation) == 0
    )
    require(
        len(c13_eligible_q8_automorphisms) == 1
        and c13_eligible_q8_automorphisms[0] == tuple(range(8))
        and 169 * 8 == 1352,
        "C169/Q8 direct-product invoice changed",
    )

    thm2779 = (
        THEOREMS
        / "THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate.md"
    ).read_text(encoding="utf-8").replace("\r\n", "\n")
    thm2884_output = (
        RESULTS / "lrc14_macro_semantic_diagonal_horn_carrier_thm2884.out"
    ).read_text(encoding="utf-8").replace("\r\n", "\n")
    thm2886_output = (
        RESULTS
        / "lrc14_stepped_origin_v4_provenance_transport_thm2886.out"
    ).read_text(encoding="utf-8").replace("\r\n", "\n")
    require(
        "quadratic refinement `q(x,y)=xy`" in thm2779
        and "H_nonzero_indicator=(1, 0, 1, 1)" in thm2884_output
        and "edges=(3to11=(1, 1, 0),11to7=(1, 0, 1),"
        "3to7=(0, 1, 1))" in thm2884_output,
        "pinned D8/V4 input statements changed",
    )
    require(
        "H_seeded_transported_selector={3: (0, 0), 11: (0, 1), "
        "7: (1, 0)}" in thm2886_output
        and "parity=(0,1,1)=1_XOR_delta_q3_on_the_seam"
        in thm2886_output,
        "pinned THM-2886 source-marked seam statement changed",
    )

    print("THM-2887 QUATERNIONIC SELECTOR-LIFT AUDIT")
    print(
        f"quadratic_values=(THM2779_D8_Arf0:{d8_values},"
        f"THM2884_Q8_Arf1:{q8_values}); "
        "common_polarization=determinant"
    )
    print(
        f"group_order_census=(D8:{tuple(sorted(d8_order_census.items()))},"
        f"Q8:{tuple(sorted(q8_order_census.items()))}); "
        f"associativity_checks={associativity_checks}; "
        f"Q8_gauge_checks={gauge_checks}"
    )
    print(
        f"normalized_V4_cocycles=(all={len(normalized_cocycles)},"
        "determinant_commutator=8,quadratic_census="
        f"{tuple(sorted(determinant_quadratic_census.items()))},"
        f"coboundaries={len(normalized_coboundaries)}); "
        "all_three_nonzero_squares_select_unique_Arf1_Q8_class"
    )
    print(
        "semantic_Q8=QA*i_then_QB*j_equals_minus_QAB*k_in_forward_gauge; "
        "reverse_order_has_plus_k; gauges_differ_by_delta(q_Arf1)"
    )
    print(
        f"local_triangle=polar(QA,QB)={semantic_triangle_polarization}"
        f"=carry_epsilon(8,9)={carry_triangle}; "
        "local_Arf/carry_match_is_gauge-invariant"
    )
    print(
        "source_marked_character=det(QB,h); "
        f"seam_(Q0,QA,QAB)={seam_source_sign}; "
        "Ad_QB=(Q0->+Q0,QB->+QB,QA->-QA,QAB->-QAB); "
        "fixed_source_breaks_S3_and_recovers_THM2886_selector_pair_XOR_parity"
    )
    print(
        f"global_no_go=carry_mod2_has_unique_normalized_primitive_"
        f"{carry_primitives[0]}; even_h_forces_zero_semantic_direction; "
        "conflicts_with_ell8=QA_and_ell4=QAB"
    )
    print(
        f"AutQ8=(size={len(q8_automorphisms)},"
        f"order_census={tuple(sorted(q8_automorphism_orders.items()))}); "
        f"InnQ8={len(inner_actions)}=V4; quotient_actions="
        f"{len(quotient_actions)}=S3; hence_AutQ8_is_S4; "
        f"mu_perm(Q8)="
        f"{q8_minimal_faithful_permutation_degree}"
    )
    print(
        f"S3_orientation_maps=(sign_action:{sign_equivariant_maps},"
        f"trivial_action:{trivial_equivariant_maps}); "
        "semantic_P1(F2)_does_not_canonically_choose_endpoint_sign; "
        f"pointed_det_character_covariance_checks="
        f"{source_character_covariance_checks}; "
        f"Stab(QB)_size={source_stabilizer_size}"
    )
    print(
        f"typed_semidirect_candidate=C169xQ8_direct_product_size={169 * 8}; "
        "the_C169_to_AutQ8_action_is_trivial; this_does_not_classify_all_"
        "order1352_extensions; Q8_is_the_unique_Arf1_central-sign_"
        "extension_of_the_semantic_V4"
    )
    print(
        "VERDICT=Q8_is_the_natural_structured_central-sign_extension_of_"
        "THM2884_support_parity_and_a_candidate_origin-sign_lift; the_Arf_"
        "polarization_detects_the_distinguished_carry_triangle; the_marked_"
        "QB_source_recovers_exactly_THM2886s_"
        "selector-pair_XOR_parity_by_inner_conjugation_but_not_its_origin_"
        "occupancy_or_signed_orientation; the_match_cannot_globalize_to_"
        "C13_via_such_a_base-independent_ell_and_no_physical_Q8_lift_or_"
        "canonical_Q8-valued_current_is_constructed"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
