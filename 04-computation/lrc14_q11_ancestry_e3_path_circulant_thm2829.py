#!/usr/bin/env python3
"""Exact path-space and circulant boundary companion for THM-2829."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT / "04-computation/"
    "lrc14_q11_semantic_reselection_ancestry_thm2829.py"
)
SOURCE_OUTPUT = (
    ROOT / "05-knowledge/results/"
    "lrc14_q11_semantic_reselection_ancestry_thm2829.out"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


require(
    lf_hash(SOURCE)
    == "6c29f820bea97a00a5fae3a808b0fb741afed5468e12b289ffae2e9f95c40f98"
    and lf_hash(SOURCE_OUTPUT)
    == "8bf43f8570d1fded3fd02257246f465e3adbd3c513ce4302477c13b84aa19712",
    "q11 ancestry source pair changed",
)

P = 13
P4 = P**4
DEPTH = P**5
ADDRESS_MODULUS = P**6
A0 = 66099
A7 = 65652
C7 = 65612
E_SUPPORT = (0, 3, 11)
ARROW_PAIRS = (
    (0, 0), (0, 7), (3, 10), (3, 4), (11, 2), (11, 9),
)
SAME_SECTION_COMMON = (66099, 65612, 66099, 65612, 66099, 65612)
NATURAL_LIFTED_COMMON = (66099, 65612, 65579, 65612, 65579, 65098)


def convolution(left, right):
    return tuple(
        sum(
            left[index] * right[(target - index) % P]
            for index in range(P)
        )
        for target in range(P)
    )


def affine_lift(t, address):
    """THM-2813 relative lift A_t on Z/13^6."""
    return (
        address
        + t * DEPTH * ((address - 7) % P)
    ) % ADDRESS_MODULUS


def top_digit_wall(residue, t, address):
    """Formal depth-six analogue of one THM-2788 carry-wall column."""
    if address % P == residue:
        return (address + t * DEPTH) % ADDRESS_MODULUS
    return address


def lifted_future_arrow(ancestry, residue, phase):
    """Add one future phase in Z/13^6, retaining its sheet carry."""
    total = residue + phase
    return (
        (ancestry + total // P) % DEPTH,
        total % P,
    )


def main():
    ancestry = [0] * P
    ancestry[0] = A0
    ancestry[7] = A7
    target = [0] * P
    for index in E_SUPPORT:
        target[index] = 1

    # Since 7 is a unit modulo 13,
    #
    # (a+b z^7) sum_j (-b)^j a^(12-j) z^(7j)
    #     =a^13+b^13.
    inverse_numerator = [0] * P
    for exponent in range(P):
        inverse_numerator[7 * exponent % P] += (
            (-A7) ** exponent * A0 ** (P - 1 - exponent)
        )
    denominator = A0**P + A7**P
    delta = (denominator,) + (0,) * (P - 1)
    require(
        convolution(ancestry, inverse_numerator) == delta,
        "two-point ancestry inverse identity changed",
    )

    intertwiner_numerator = convolution(target, inverse_numerator)
    require(
        convolution(ancestry, intertwiner_numerator)
        == tuple(denominator * value for value in target),
        "unique circulant ancestry-to-E3 map changed",
    )
    sign_pattern = "".join(
        "+" if value > 0 else "-" if value < 0 else "0"
        for value in intertwiner_numerator
    )
    require(
        sign_pattern == "+++++++------"
        and sum(value > 0 for value in intertwiner_numerator) == 7
        and sum(value < 0 for value in intertwiner_numerator) == 6
        and Fraction(sum(intertwiner_numerator), denominator)
        == Fraction(1, 43917),
        "signed circulant boundary changed",
    )

    translate_pairs = tuple(
        (shift, (shift + 7) % P)
        for shift in range(P)
    )
    contained_pairs = tuple(
        pair for pair in translate_pairs
        if set(pair).issubset(E_SUPPORT)
    )
    require(
        not contained_pairs,
        "a positive ancestry translate entered E3 support",
    )

    # The strict shared-q fibre product is empty away from q=0.  Retaining
    # a path phase h instead gives a positive homotopy/path fibre product.
    path_coupling = tuple(
        tuple(target[q] * ancestry[(q + h) % P] for h in range(P))
        for q in range(P)
    )
    path_support = tuple(
        (q, h, path_coupling[q][h])
        for q in range(P)
        for h in range(P)
        if path_coupling[q][h]
    )
    expected_path_support = (
        (0, 0, A0),
        (0, 7, A7),
        (3, 4, A7),
        (3, 10, A0),
        (11, 2, A0),
        (11, 9, A7),
    )
    q_marginal = tuple(sum(path_coupling[q]) for q in range(P))
    r_marginal = tuple(
        sum(path_coupling[q][(r - q) % P] for q in range(P))
        for r in range(P)
    )
    require(
        path_support == expected_path_support
        and q_marginal
        == tuple((A0 + A7) * target[q] for q in range(P))
        and r_marginal
        == tuple(3 * ancestry[r] for r in range(P)),
        "positive path-space coupling or its marginals changed",
    )
    common_ancestry = list(ancestry)
    common_ancestry[7] = C7
    common_path_coupling = tuple(
        tuple(
            target[q] * common_ancestry[(q + h) % P]
            for h in range(P)
        )
        for q in range(P)
    )
    common_path_support = tuple(
        (q, h, common_path_coupling[q][h])
        for q in range(P)
        for h in range(P)
        if common_path_coupling[q][h]
    )
    require(
        tuple((q, h) for q, h, _value in common_path_support)
        == tuple((q, h) for q, h, _value in expected_path_support)
        and common_path_coupling[11][9] == C7,
        "common-ancestry path support changed",
    )
    strict_nonzero = tuple(
        q for q in range(1, P)
        if target[q] and ancestry[q]
    )
    require(
        not strict_nonzero,
        "strict nonzero shared-q fibre product became positive",
    )

    # THM-2820's conditional normal decoder is beta(q)=2q because the
    # physical successor/transvection commutator has q=7t.  Its graph meets
    # the nontrivial ancestry section h=7-q uniquely at q=11.
    normal_phase = tuple((2 * q) % P for q in range(P))
    normal_section = tuple(
        (q, normal_phase[q], path_coupling[q][normal_phase[q]])
        for q in E_SUPPORT
        if path_coupling[q][normal_phase[q]]
    )
    require(
        normal_section == ((0, 0, A0), (11, 9, A7))
        and (2 * 3) % P not in (4, 10)
        and (2 * 11) % P == (7 - 11) % P == 9,
        "normal-decoder path selection changed",
    )
    common_normal_section = tuple(
        (
            q,
            normal_phase[q],
            common_path_coupling[q][normal_phase[q]],
        )
        for q in E_SUPPORT
        if common_path_coupling[q][normal_phase[q]]
    )
    require(
        common_normal_section == ((0, 0, A0), (11, 9, C7)),
        "common normal-decoder path selection changed",
    )
    carry_defect = tuple(
        lifted - same
        for lifted, same in zip(
            NATURAL_LIFTED_COMMON, SAME_SECTION_COMMON
        )
    )
    require(
        carry_defect == (0, 0, -520, 0, -520, -514)
        and tuple(value % P for value in carry_defect)
        == (0, 0, 0, 0, 0, 6)
        and gcd(abs(carry_defect[-1]), 91) == 1
        and tuple(
            lifted_future_arrow(17, q, h)
            for q, h in ARROW_PAIRS
        )
        == (
            (17, 0), (17, 7), (18, 0),
            (17, 7), (18, 0), (18, 7),
        ),
        "borrow-aware normal-path selector changed",
    )
    for start_residue in range(P):
        state = (23, start_residue)
        for _step in range(P):
            state = lifted_future_arrow(*state, 1)
        require(
            state == (24, start_residue),
            "thirteen future generators ceased to make one deck shift",
        )
    u = 7
    nontrivial_ancestry_residue = 7
    selected_t = (
        pow(u + 1, -1, P) * nontrivial_ancestry_residue
    ) % P
    selected_q = (u * selected_t) % P
    require(
        (selected_t, selected_q) == (9, 11)
        and (u * selected_t + selected_t) % P
        == nontrivial_ancestry_residue,
        "general affine normal-path formula changed",
    )

    # Exact role trace.  A_9 has already produced the physical q=11
    # displacement on the residue-eight source.  It fixes the residue-seven
    # target, so it cannot also be the missing target-relative ancestry
    # phase.  The formal one-column wall has the same role asymmetry.
    source_base = 6716
    target_base = 6715
    source_lifted = 3348353
    target_lifted = 3348352
    require(
        affine_lift(selected_t, source_base) == source_lifted
        and affine_lift(selected_t, target_base) == target_base
        and target_lifted == source_lifted - 1
        and target_lifted != affine_lift(selected_t, target_base)
        and top_digit_wall(8, selected_t, source_base) == source_lifted
        and top_digit_wall(8, selected_t, target_base) == target_base
        and (
            7 * ((source_lifted - source_base) // DEPTH)
        ) % P == 11,
        "A9/carry-wall role trace changed",
    )

    # The existing representative gauge translates the coarse ancestry
    # label a by s*13^4 and leaves q fixed.  A_t itself is the identity
    # modulo 13^5.  Neither operation creates the fine path phase h.
    require(
        DEPTH == P * P4
        and gcd(P4, DEPTH) == P4
        and all(
            len({(a + s * P4) % DEPTH for a in range(P)})
            == P
            for s in range(P)
        )
        and all(
            affine_lift(t, a) % DEPTH == a
            for t in range(P)
            for a in range(DEPTH)
        )
        and 11 * P4 == 314171,
        "coarse ancestry-gauge or depth reduction changed",
    )

    print("ANCESTRY-TO-E3 CIRCULANT BOUNDARY PROBE")
    print("status=THM-2829 VERIFIED-EXACT sidecar; no physical map or LRC14")
    print(
        f"ancestry={tuple(ancestry)}; "
        f"E3_indicator={tuple(target)}"
    )
    print(
        f"inverse_denominator=a^13+b^13={denominator}; "
        f"intertwiner_numerator={intertwiner_numerator}"
    )
    print(
        f"intertwiner_sign_pattern={sign_pattern}; "
        "positive=7; negative=6; "
        f"augmentation={Fraction(sum(intertwiner_numerator), denominator)}"
    )
    print(
        f"translated_ancestry_support_pairs={translate_pairs}; "
        f"pairs_contained_in_E3={contained_pairs}"
    )
    print(
        f"strict_nonzero_shared_q={strict_nonzero}; "
        f"path_support={path_support}"
    )
    print(
        f"common_ancestry_path_support={common_path_support}"
    )
    print(
        f"path_q_marginal={q_marginal}; "
        f"path_r_marginal={r_marginal}"
    )
    print(
        f"normal_beta_graph_on_E3={normal_section}; "
        f"common_normal_beta_graph={common_normal_section}; "
        "general_formula=(u+1)t=r; "
        f"(u,r,t,q)=({u},{nontrivial_ancestry_residue},"
        f"{selected_t},{selected_q})"
    )
    print(
        f"natural_lifted_common={NATURAL_LIFTED_COMMON}; "
        f"same_section_common={SAME_SECTION_COMMON}; "
        f"carry_defect_lift_minus_same={carry_defect}; "
        f"carry_defect_mod13={tuple(value % P for value in carry_defect)}; "
        "beta_arrow_abs_defect=514_gcd91_1; "
        "thirteen_generators=one_positive_deck_shift"
    )
    print(
        "address_role_trace="
        f"(A9({source_base})={source_lifted},"
        f"A9({target_base})={affine_lift(selected_t, target_base)},"
        f"adjacent_target={target_lifted}); "
        "formal_C8^9_has_same_source_only_action"
    )
    print(
        "gauge_boundary="
        "ell+sW_reindexes_a_by_s*13^4_with_q_fixed; "
        "A_t_is_identity_mod13^5; "
        "q11_ancestor_shift_reindexes_a_by_11*13^4=314171"
    )
    print(
        "CONCLUSION: the ancestry profile is a rational C13 "
        "group-algebra unit, so a unique circulant coefficient "
        "intertwiner to E3 exists.  It necessarily has both signs.  "
        "No nonzero nonnegative circulant kernel can land inside E3, "
        "because each positive source coefficient creates a translate "
        "of {0,7}, while E3={0,3,11} contains no such pair.  A physical "
        "repair must therefore be signed, noncirculant, or retain an "
        "extra path fibre.  Retaining h gives a positive six-point path "
        "coupling.  Its residue-marginal same-section weights are not the "
        "natural lifted weights on wrapped arrows.  The conditional normal "
        "label beta(q)=2q selects q=11, where the +1 borrow leaves 65098 "
        "common sheets and creates the unique nonzero mod13 carry response "
        "6 (absolute defect 514, a 91-unit).  But A9 already "
        "produces the physical q11 displacement and fixes the target role; "
        "neither A9, the representative gauge, nor THM-2788's carry-wall "
        "shape supplies the new target-relative U_Q action h=beta."
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
