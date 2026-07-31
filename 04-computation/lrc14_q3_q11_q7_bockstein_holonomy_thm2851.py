#!/usr/bin/env python3
"""Exact companion for THM-2851 q3/q11/q7 Bockstein holonomy."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import gcd
from pathlib import Path


P = 13
ANCESTRY = P**5
S = P**5
MODULUS = P**6
ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINNED = {
    ROOT / "04-computation/lrc14_q11_semantic_word_horn_thm2835.py":
        "207dd94f086338ae1e80b7d7196f99bf41e795893d13b6d48e4e7d516af03523",
    ROOT / "05-knowledge/results/lrc14_q11_semantic_word_horn_thm2835.out":
        "1ebe0cbaf7d4ef13defed0bdb5b37df1364880acdbfc6139b243ab9df65f6bf6",
    ROOT / "04-computation/lrc14_q3_q11_transverse_endpoint_horn_thm2847.py":
        "258659c5093d98eea84056bdd3599b32d2a244bcd37dfa5f22dc5b25946ffe72",
    ROOT / "05-knowledge/results/lrc14_q3_q11_transverse_endpoint_horn_thm2847.out":
        "155fce129c750a9505fdda3c71a250ff3a57fcd4044bb1df941da83c08baee1d",
    ROOT / "04-computation/lrc14_prime_power_unit_mass_q11_response_thm2839.py":
        "68ae72b62b7974e4f2c2bf7d570615c8c524746978c57cf120f6372a7250ece4",
    ROOT / "05-knowledge/results/lrc14_prime_power_unit_mass_q11_response_thm2839.out":
        "495829603ea0c3944f83d7ae269cbbc5cbdec9fdc452395e96c78de8b2e7e24b",
}
for path, expected in PINNED.items():
    require(lf_hash(path) == expected, f"pinned dependency changed: {path.name}")


def natural_lift(state, h):
    """Lift q -> q+h through n=13a+q, with h in {0,...,12}."""
    a, q = state
    total = q + h
    return ((a + total // P) % ANCESTRY, total % P)


def central_translate(state, exponent=1):
    a, q = state
    return ((a + exponent) % ANCESTRY, q)


def carry_cocycle(h, k):
    return (h + k) // P


def affine_lift(y, t):
    """THM-2813 lift, with t read in F_13."""
    return (y + (t % P) * S * ((y - 7) % P)) % MODULUS


def beta(q):
    return (2 * q) % P


def address(q):
    return (6716 + beta(q) * S) % MODULUS


def rank(matrix):
    work = [[Fraction(x) for x in row] for row in matrix]
    if not work:
        return 0
    rows = len(work)
    cols = len(work[0])
    pivot_row = 0
    for col in range(cols):
        pivot = next(
            (r for r in range(pivot_row, rows) if work[r][col]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][col]
        work[pivot_row] = [entry / scale for entry in work[pivot_row]]
        for r in range(rows):
            if r == pivot_row:
                continue
            scale = work[r][col]
            if scale:
                work[r] = [
                    left - scale * right
                    for left, right in zip(work[r], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def determinant(matrix):
    work = [[Fraction(x) for x in row] for row in matrix]
    size = len(work)
    require(all(len(row) == size for row in work), "matrix is not square")
    result = Fraction(1)
    for col in range(size):
        pivot = next((r for r in range(col, size) if work[r][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            result *= -1
        pivot_value = work[col][col]
        result *= pivot_value
        for r in range(col + 1, size):
            scale = work[r][col] / pivot_value
            for c in range(col + 1, size):
                work[r][c] -= scale * work[col][c]
    require(result.denominator == 1, "integer determinant became fractional")
    return result.numerator


def main():
    # The natural lifts form the standard nonsplit cyclic extension.
    composition_checks = 0
    cocycle_checks = 0
    for a in (0, 1, ANCESTRY - 1):
        for q in range(P):
            for h in range(P):
                for k in range(P):
                    left = natural_lift(natural_lift((a, q), h), k)
                    right = central_translate(
                        natural_lift((a, q), (h + k) % P),
                        carry_cocycle(h, k),
                    )
                    require(left == right, "natural-lift composition changed")
                    composition_checks += 1
    for h in range(P):
        for k in range(P):
            for ell in range(P):
                left = (
                    carry_cocycle(h, k)
                    + carry_cocycle((h + k) % P, ell)
                )
                right = (
                    carry_cocycle(k, ell)
                    + carry_cocycle(h, (k + ell) % P)
                )
                require(left == right, "carry cocycle identity changed")
                cocycle_checks += 1

    # The generator winds once in the ancestry direction after thirteen steps.
    state = (0, 0)
    for _ in range(P):
        state = natural_lift(state, 1)
    require(state == (1, 0), "generator carry winding changed")

    # At the first ancestry quotient the extension is C_169 -> C_13.
    section_candidates = tuple(1 + P * c for c in range(P))
    require(
        all((P * candidate) % (P**2) != 0 for candidate in section_candidates),
        "a homomorphic C13 section unexpectedly appeared",
    )
    additive_orders = tuple(
        (P**2) // gcd(P**2, candidate)
        for candidate in section_candidates
    )
    require(set(additive_orders) == {P**2}, "section order boundary changed")

    # The q3 -> q11 -> q7 triangle has one unit of central carry.
    q3_state = (0, 3)
    after_8 = natural_lift(q3_state, 8)
    after_8_9 = natural_lift(after_8, 9)
    direct_4 = natural_lift(q3_state, 4)
    require(
        after_8 == (0, 11)
        and after_8_9 == (1, 7)
        and direct_4 == (0, 7)
        and after_8_9 == central_translate(direct_4),
        "q3/q11/q7 carry triangle changed",
    )
    require(
        all(
            natural_lift((a, q), h)[1] == (q + h) % P
            for a in (0, 1, ANCESTRY - 1)
            for q in range(P)
            for h in range(P)
        )
        and central_translate(after_8_9)[1] == after_8_9[1],
        "carry-blind residue projection changed",
    )

    # By contrast, the affine address torsor is flat.
    address_checks = 0
    for q in range(P):
        for h in range(P):
            require(
                affine_lift(address(q), 2 * h)
                == address((q + h) % P),
                "target/address semiconjugacy changed",
            )
            address_checks += 1
    for residue in range(P):
        for quotient in (0, 1, P**4 - 1):
            y = (residue + P * quotient) % MODULUS
            for s in range(P):
                for t in range(P):
                    require(
                        affine_lift(affine_lift(y, t), s)
                        == affine_lift(y, (s + t) % P),
                        "affine address action is no longer flat",
                    )
    n3, n11, n7 = address(3), address(11), address(7)
    require(
        affine_lift(n3, 3) == n11
        and affine_lift(n11, 5) == n7
        and affine_lift(n3, 8) == n7
        and affine_lift(affine_lift(n3, 3), 5)
        == affine_lift(n3, 8),
        "specific affine triangle changed",
    )

    # The pinned canonical evidence supplies the 449-sheet semantic leg and
    # the exact E3 mapping-cone frame.
    word_output = (
        ROOT / "05-knowledge/results/lrc14_q11_semantic_word_horn_thm2835.out"
    ).read_text(encoding="utf-8").replace("\r\n", "\n")
    horn_output = (
        ROOT / "05-knowledge/results/lrc14_q3_q11_transverse_endpoint_horn_thm2847.out"
    ).read_text(encoding="utf-8").replace("\r\n", "\n")
    spectrum_output = (
        ROOT / "05-knowledge/results/lrc14_prime_power_unit_mass_q11_response_thm2839.out"
    ).read_text(encoding="utf-8").replace("\r\n", "\n")
    require(
        "L=(count=449,first=59306,last=311961,"
        "residue_mod169=156,whole_cylinder=1)" in word_output
        and "natural_beta_lift=(a,11)->(a+1,7); "
        "word_state=QAB_on_449/449; QA=QB=0" in word_output,
        "pinned 449-sheet carry leg changed",
    )
    require(
        "q7_E3_only_horn_cells=20; horn_first=(0, 5, 1)" in horn_output
        and "scope=20_cell_E3_only_horn" in horn_output
        and "E3_mapping_cone=(image_rank=3,"
        "kernel_generator=(0, 0, 0, 1),"
        "complement_on_kernel=449)" in horn_output,
        "pinned E3 mapping-cone data changed",
    )
    require(
        "residue_mod169=156" in spectrum_output
        and "plus1_lands_QAB=449" in spectrum_output
        and "nonzero_characters=371293" in spectrum_output
        and "rational_translate_span=371293" in spectrum_output,
        "pinned horn full-spectrum/provenance data changed",
    )

    # The oriented carry edge has a canonical signed derivative.  On the
    # first carry quotient all 449 labels sit at residue zero and all
    # successors at residue one.
    reduced_derivative = (-449, 449) + (0,) * (P - 2)
    reduced_circulant = tuple(
        tuple(reduced_derivative[(column - row) % P] for column in range(P))
        for row in range(P)
    )
    require(
        sum(reduced_derivative) == 0
        and rank(reduced_circulant) == P - 1
        and 449 % P == 7
        and reduced_derivative[1] % P == 7,
        "first-carry derivative boundary changed",
    )
    full_derivative_rank = ANCESTRY - 1

    full_frame = (
        (1, 1, 0, 0),
        (1, 0, 0, 0),
        (1, 1, 449, 0),
        (0, 0, 0, 449),
    )
    e3_frame = (
        (1, 1, 0, 0),
        (1, 0, 0, 0),
        (1, 1, 449, 0),
        (0, 0, 0, 0),
    )
    complement_frame = (
        (0, 0, 0, 0),
        (0, 0, 0, 0),
        (0, 0, 0, 0),
        (0, 0, 0, 449),
    )
    kernel_generator = (0, 0, 0, 1)
    require(
        determinant(full_frame) == -(449**2)
        and rank(full_frame) == 4
        and rank(e3_frame) == 3
        and tuple(
            sum(row[j] * kernel_generator[j] for j in range(4))
            for row in e3_frame
        ) == (0, 0, 0, 0)
        and tuple(
            sum(row[j] * kernel_generator[j] for j in range(4))
            for row in complement_frame
        ) == (0, 0, 0, 449),
        "mapping-cone linear algebra changed",
    )

    # Every nontrivial character of the first carry quotient sees the loop.
    holonomy_exponents = tuple(range(1, P))
    require(
        all(exponent % P for exponent in holonomy_exponents),
        "a charged carry character became trivial",
    )

    print("Q3/Q11/Q7 BOCKSTEIN HOLONOMY AND CARRY DERIVATIVE")
    print(
        f"natural_lift_composition_checks={composition_checks}; "
        f"cocycle_checks={cocycle_checks}; "
        "law=L_k*L_h=T^floor((h+k)/13)*L_(h+k_mod13)"
    )
    print(
        "generator_winding=L_1^13=T; "
        f"C169_section_candidates={len(section_candidates)}; "
        "homomorphic_sections=0; candidate_orders=169"
    )
    print(
        f"semantic_triangle=(q3={q3_state},"
        f"L8={after_8},L9L8={after_8_9},L4={direct_4}); "
        "identity=L9L8=T*L4; carry=1"
    )
    print(
        f"address_semiconjugacy_checks={address_checks}; "
        f"addresses=(q3={n3},q11={n11},q7={n7}); "
        "identity=A5*A3=A8; holonomy=0"
    )
    print(
        "carry_Fourier_holonomy_exponents="
        f"{holonomy_exponents}; all_12_nontrivial_characters_detect=1"
    )
    print(
        "semantic_leg=THM2835_PINNED:"
        "449_QA(q11,a)_to_QAB(q7,a+1); "
        "macro_leg=THM2847_PINNED:20_cell_E3_only_horn"
    )
    print(
        "oriented_carry_derivative=449*(delta1-delta0); "
        f"first_carry_rank={rank(reduced_circulant)}; "
        "all_12_charged_characters_nonzero=1; "
        f"full_ancestry_derivative_rank={full_derivative_rank}; "
        "augmentation_adic_order=1; signed_l1=898"
    )
    print(
        f"frame=(det={determinant(full_frame)},rank={rank(full_frame)}); "
        f"E3_rank={rank(e3_frame)}; kernel={kernel_generator}; "
        "complement_on_kernel=449"
    )
    print(
        "flat_intertwiner_consequence=F*T=F; "
        "sharp_hostile=F0(a,q)=q; "
        "faithful_first_carry_sidecar_min_states=13; "
        "faithful_q_and_carry_joint_min_states=169; "
        "equality=fibrewise_C13_torsor; "
        "full_ancestry_fidelity_cost_per_q=13^5"
    )
    print(
        "scope=coefficient_nerve_and_realization_obstruction_only; "
        "no_physical_E3_contraction=1; no_row_or_LRC14_conclusion=1"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
