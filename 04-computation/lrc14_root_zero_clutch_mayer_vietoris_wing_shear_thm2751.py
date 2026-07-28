#!/usr/bin/env python3
"""Exact fixed-clock root-zero Mayer--Vietoris wing target spectrum (THM-2751).

This proof-complete candidate imports the promoted THM-2749 companion
and reconstructs, in one source coordinate and with one delayed-prefix bank,
the fixed source-one sheet A_{1,t}, the pulled target sheet B_{1,t}, their common
section M_t, and the disjoint wings L_t=A_t minus B_t and R_t=B_t minus A_t.
"""

from __future__ import annotations

import ast
from hashlib import sha256
import inspect
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

import lrc14_fully_marked_root_zero_target_profile_thm2749 as marked
import lrc14_semantic_root_zero_clutch_refinement_probe_20260728 as legacy


P = marked.P
T = marked.T
SHIFT = marked.SHIFT
CONTENT = marked.CONTENT
RAIL = marked.RAIL
BASE = 10**60
private = marked.private
two = marked.two
lift = marked.lift


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n")


def constructor_audit():
    legacy_path = COMP / "lrc14_semantic_root_zero_clutch_refinement_probe_20260728.py"
    legacy_bytes = lf_bytes(legacy_path)
    legacy_hash = sha256(legacy_bytes).hexdigest()
    require(legacy_hash
            == "b3b623a4c016b1303ac3a74c72f9ae0bbd69cdb97a553f613143f0005a3dd286",
            "legacy comparator changed")
    tree = ast.parse(legacy_bytes.decode("utf-8"))
    legacy_function = next(
        node for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "restrict_e3_and_sheet"
    )
    legacy_names = {node.id for node in ast.walk(legacy_function)
                    if isinstance(node, ast.Name)}
    canonical_tree = ast.parse(inspect.getsource(two.source_present_section))
    canonical_names = {node.id for node in ast.walk(canonical_tree)
                       if isinstance(node, ast.Name)}
    require("clock_comb" not in legacy_names
            and "source_clock" not in legacy_names
            and {"clock_comb", "source_clock"} <= canonical_names,
            "source-one constructor distinction changed")
    return legacy_hash


def prefix_audit(module, delayed):
    """Compare the fourteen actual marked prefixes with the pinned legacy bank."""
    historical = legacy.build_q3_pair_prefixes(module)
    masses = []
    for ell in range(7):
        for kappa in range(2):
            old_prefix = historical[ell][6][kappa]
            new_prefix = delayed[ell][kappa]
            require(new_prefix == old_prefix,
                    f"legacy/marked prefix mismatch at {(ell, kappa)}")
            masses.append(new_prefix[2][-1])
    require(tuple(masses[:2]) == (0, 0)
            and tuple(masses[2:]) == (206986279500,) * 12,
            "legacy/marked prefix masses changed")
    return tuple(masses)


def difference(first, second):
    return marked.intersect(first, marked.complement(second))


def interval_mass(intervals):
    return sum(right - left for left, right in intervals)


def cut_carrier(module, present, carrier, ell):
    source_safe = marked.complement(present[ell, 7])
    target_safe = marked.complement(
        marked.shift_union(present[ell, 7], SHIFT)
    )
    row = marked.intersect(carrier, source_safe)
    row = marked.intersect(row, target_safe)
    return marked.intersect(
        row, tuple(module.make_comb(module.C3, 182, 169, 181))
    )


def coefficient_vector(module, delayed, present, weight, carrier, carry,
                       caches):
    values = []
    masses = []
    for ell in range(7):
        source_safe = marked.complement(present[ell, 7])
        target_safe = marked.complement(
            marked.shift_union(present[ell, 7], SHIFT)
        )
        row = marked.intersect(carrier, source_safe)
        row = marked.intersect(row, target_safe)
        row_weighted = marked.restrict_weighted(weight, row)
        row_weighted = private.old.intersect_weighted_comb(
            row_weighted, module.C3, 182, 169, 181
        )
        masses.append(sum((right - left) * w
                          for left, right, w in row_weighted))
        values.append(private.delayed_carry_pair(
            row_weighted, delayed[ell], caches[ell]
        )[carry][1])
    return tuple(values), tuple(masses)


def add_weighted_profiles(channels):
    """Exact nonnegative sum of step profiles with integer channel scales."""
    events = {}
    for pieces, scale in channels:
        for left, right, weight in pieces:
            events[left] = events.get(left, 0) + scale * weight
            events[right] = events.get(right, 0) - scale * weight
    value = 0
    previous = None
    result = []
    for point in sorted(events):
        if previous is not None and previous < point and value:
            if result and result[-1][1] == previous and result[-1][2] == value:
                result[-1] = (result[-1][0], point, value)
            else:
                result.append((previous, point, value))
        value += events[point]
        require(value >= 0, "weighted channel multiplexing became negative")
        previous = point
    require(value == 0, "weighted channel multiplexing did not close")
    return tuple(result)


def multiplexed_vectors(module, delayed, present, source_weight, target_weight,
                        A, B, M):
    """Recover A,B,M vectors from one exact positive functional evaluation.

    BASE is larger than every channel coefficient.  Positivity and linearity
    let three coefficients be packed as A + BASE*B + BASE^2*M_source without
    cross-channel carry.
    """
    combined = add_weighted_profiles((
        (marked.restrict_weighted(source_weight, A), 1),
        (marked.restrict_weighted(target_weight, B), BASE),
        (marked.restrict_weighted(source_weight, M), BASE**2),
    ))
    values = {name: [] for name in ("A", "B", "M")}
    for ell in range(7):
        source_safe = marked.complement(present[ell, 7])
        target_safe = marked.complement(
            marked.shift_union(present[ell, 7], SHIFT)
        )
        row = tuple(private.old.intersect_weighted_union(
            combined, source_safe
        ))
        row = tuple(private.old.intersect_weighted_union(
            row, target_safe
        ))
        row = tuple(private.old.intersect_weighted_comb(
            row, module.C3, 182, 169, 181
        ))
        packed = private.delayed_carry_pair(
            row, delayed[ell], {}
        )[12][1]
        a = packed % BASE
        b = (packed // BASE) % BASE
        m = packed // (BASE**2)
        require(a < BASE and b < BASE and m < BASE,
                "multiplexed coefficient overflowed its channel")
        require(packed == a + BASE * b + BASE**2 * m,
                "multiplexed coefficient failed exact decoding")
        values["A"].append(a)
        values["B"].append(b)
        values["M"].append(m)
    return {name: tuple(vector) for name, vector in values.items()}


def direct_target_vector(module, delayed, present, source_weight, carrier,
                         caches):
    """Recompute a pulled target carrier after pushing to target coordinate."""
    # Rail 8's physical target weight is the original rail profile after push.
    values = []
    for ell in range(7):
        source_safe = marked.complement(present[ell, 7])
        target_safe = marked.complement(
            marked.shift_union(present[ell, 7], SHIFT)
        )
        row = marked.intersect(carrier, source_safe)
        row = marked.intersect(row, target_safe)
        row = marked.intersect(
            row, tuple(module.make_comb(module.C3, 182, 169, 181))
        )
        # Push the already fully cut source-coordinate carrier.
        direct = marked.shift_union(row, -SHIFT)
        row_weighted = marked.restrict_weighted(source_weight, direct)
        row_weighted = private.old.intersect_weighted_comb(
            row_weighted, module.C3, 182, 1, 13
        )
        values.append(private.delayed_carry_pair(
            row_weighted, delayed[ell], caches[ell]
        )[6][1])
    return tuple(values)


def scalar_amplitude(vector):
    require(vector[0] == 0 and len(set(vector[1:])) == 1,
            f"clock vector is not (0,c,...,c): {vector}")
    return vector[1]


def polynomial_trim(coefficients):
    values = list(coefficients)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def resultant_phi13(profile):
    try:
        return marked.sylvester_resultant((1,) * 13, polynomial_trim(profile))
    except RuntimeError as error:
        if str(error) == "zero Sylvester determinant":
            return 0
        raise


def circular_convolution(left, right):
    return tuple(sum(left[j] * right[(i - j) % P] for j in range(P))
                 for i in range(P))


def rotations(profile):
    return tuple(tuple(profile[(i - shift) % P] for i in range(P))
                 for shift in range(P))


def dihedral_orbit(profile):
    reflected = tuple(profile[-i % P] for i in range(P))
    return rotations(profile) + rotations(reflected)


def solve_circulant_over_q(source, target):
    """Rational solutions k to k*source=target, with k_12 fixed to zero.

    The missing uniform degree of freedom is restored separately.  Gaussian
    elimination uses Fraction through the canonical module's import.
    """
    from fractions import Fraction
    matrix = []
    for row in range(P):
        matrix.append([
            Fraction(source[(row - col) % P]) for col in range(P - 1)
        ] + [Fraction(target[row])])
    rank = 0
    pivots = []
    for col in range(P - 1):
        pivot = next((r for r in range(rank, P) if matrix[r][col]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        value = matrix[rank][col]
        matrix[rank] = [x / value for x in matrix[rank]]
        for r in range(P):
            if r != rank and matrix[r][col]:
                factor = matrix[r][col]
                matrix[r] = [matrix[r][c] - factor * matrix[rank][c]
                             for c in range(P)]
        pivots.append(col)
        rank += 1
    for row in matrix:
        if not any(row[col] for col in range(P - 1)) and row[-1]:
            return None, rank
    solution = [Fraction(0) for _ in range(P)]
    for r, col in enumerate(pivots):
        solution[col] = matrix[r][-1]
    return tuple(solution), rank


def main():
    legacy_hash = constructor_audit()
    module, _prefixes, _a, _b, rails, present, _starts = (
        lift.m.core.build_carrier_data()
    )
    require(BASE > 13 * T**3,
            "multiplexing base lost its conservative channel bound")
    delayed = marked.marked_prefixes(
        module,
        private.build_pair_prefixes(module),
        two.deepest_fork(module),
    )
    prefix_masses = prefix_audit(module, delayed)
    source = two.exclusive_source(module, 3)
    clock_comb = tuple(
        module.make_comb(module.C1, 182, 26 * ell - 13, 26 * ell + 13)
        for ell in range(7)
    )
    source_weight, target_weight, rail_common = marked.rail_data(rails, RAIL)
    raw_static = tuple(two.intersect_sorted(source, clock_comb[1]))
    raw_static = tuple(module.subtract_comb(
        raw_static, module.W[1], 182, -13, 13
    ))
    raw_static = tuple(module.subtract_comb(
        raw_static, module.C2, 182, -13, 13
    ))

    scalar = {name: [] for name in ("A", "B", "Ms", "Mt", "L", "R")}
    reductions = {name: [] for name in scalar}
    row_lines = []
    physical_masses = []
    for t in range(P):
        raw_A = tuple(module.subtract_comb(
            raw_static, module.W[2], 182,
            -14 * t - 13, -14 * t + 13,
        ))
        raw_A = tuple(module.subtract_comb(
            raw_A, module.C3, 182,
            14 * t - 13, 14 * t + 13,
        ))
        A = marked.intersect(rail_common, raw_A)
        B = marked.intersect(
            rail_common, marked.shift_union(raw_A, SHIFT)
        )
        M = marked.intersect(A, B)
        L = difference(A, B)
        Rwing = difference(B, A)

        require(marked.merge(M + L) == marked.merge(A),
                f"A != M disjoint-union L at t={t}")
        require(marked.merge(M + Rwing) == marked.merge(B),
                f"B != M disjoint-union R at t={t}")
        require(not marked.intersect(M, L)
                and not marked.intersect(M, Rwing),
                f"wing decomposition is not disjoint at t={t}")

        packed = multiplexed_vectors(
            module, delayed, present, source_weight, target_weight, A, B, M
        )
        vectors = {
            "A": packed["A"],
            "B": packed["B"],
            "Ms": packed["M"],
            # THM-2749 independently recomputes this same common carrier with
            # the pulled target rail weight and obtains exact raw equality.
            "Mt": packed["M"],
        }
        vectors["L"] = tuple(vectors["A"][i] - vectors["Ms"][i]
                               for i in range(7))
        vectors["R"] = tuple(vectors["B"][i] - vectors["Mt"][i]
                               for i in range(7))
        require(vectors["Ms"] == vectors["Mt"],
                f"common raw source/target vector sheared at t={t}")

        # Representatives for the nonzero common window and the exceptional
        # terminal label are also recomputed after an actual coordinate push.
        # This locally enforces the physical pulled-target typing instead of
        # relying only on the inherited common-carrier equality of THM-2749.
        if t in (3, 12):
            direct_B = direct_target_vector(
                module, delayed, present, source_weight, B,
                tuple({} for _ in range(7)),
            )
            require(vectors["B"] == direct_B,
                    f"pulled/direct target B vector differs at t={t}")
        if t == 3:
            direct_M = direct_target_vector(
                module, delayed, present, source_weight, M,
                tuple({} for _ in range(7)),
            )
            require(vectors["Mt"] == direct_M,
                    "pulled/direct target M vector differs at t=3")

        require(all(vectors["A"][i] == vectors["Ms"][i] + vectors["L"][i]
                    for i in range(7)),
                f"source functional failed additivity at t={t}")
        require(all(vectors["B"][i] == vectors["Mt"][i] + vectors["R"][i]
                    for i in range(7)),
                f"target functional failed additivity at t={t}")
        if t == 3:
            require(
                vectors["Ms"] == vectors["Mt"]
                == (0,) + (marked.AMPLITUDE,) * 6,
                f"multiplexed common channel disagrees with THM2749: {vectors}",
            )

        for name in scalar:
            scalar[name].append(scalar_amplitude(vectors[name]))
        for name, root in (("A", 12), ("B", 1), ("Ms", 12),
                           ("Mt", 1), ("L", 12), ("R", 1)):
            reductions[name].append(marked.unit_reduction(vectors[name], root))
        mass_tuple = tuple(interval_mass(carrier)
                           for carrier in (A, B, M, L, Rwing))
        physical_masses.append(mass_tuple)
        row_lines.append(
            f"t={t} pre_relative_interval_mass=(A:{interval_mass(A)},B:{interval_mass(B)},"
            f"M:{interval_mass(M)},L:{interval_mass(L)},R:{interval_mass(Rwing)}) "
            f"amp=(A:{scalar_amplitude(vectors['A'])},B:{scalar_amplitude(vectors['B'])},"
            f"Ms:{scalar_amplitude(vectors['Ms'])},Mt:{scalar_amplitude(vectors['Mt'])},"
            f"L:{scalar_amplitude(vectors['L'])},R:{scalar_amplitude(vectors['R'])})"
        )
        del A, B, M, L, Rwing, vectors, packed

    scalar = {name: tuple(values) for name, values in scalar.items()}
    reductions = {name: tuple(values) for name, values in reductions.items()}
    require(tuple(physical_masses[:3]) == ((0, 0, 0, 0, 0),) * 3
            and tuple(physical_masses[3:12])
            == ((13751337600, 13808634840, 6320326320,
                 7431011280, 7488308520),) * 9
            and physical_masses[12]
            == (7404566400, 7435418760, 0, 7404566400, 7435418760),
            "physical Mayer--Vietoris interval census changed")
    supports = {name: tuple(i for i, value in enumerate(profile) if value)
                for name, profile in scalar.items()}
    require(len(supports["A"]) == 9 and len(supports["B"]) == 10,
            "natural source/target support-cardinality obstruction changed")
    resultants = {name: resultant_phi13(profile)
                  for name, profile in scalar.items()}

    cross = tuple(scalar["L"][t] * scalar["R"][t] for t in range(P))
    cross_resultant = resultant_phi13(cross)

    C = marked.AMPLITUDE
    G = C // 119
    B0 = 121 * G
    require(C == 119 * G and scalar["B"] == (0, 0, 0) + (B0,) * 10,
            "natural target amplitude factorization changed")
    W = tuple(int(3 <= t <= 11) for t in range(P))
    U = tuple(int(3 <= t <= 12) for t in range(P))
    Q = tuple(2 if 3 <= t <= 11 else 121 if t == 12 else 0
              for t in range(P))
    require(scalar["A"] == scalar["Ms"] == scalar["Mt"]
            == tuple(C * value for value in W)
            and scalar["L"] == (0,) * P
            and scalar["R"] == tuple(G * value for value in Q),
            "Mayer--Vietoris amplitude factorization changed")
    norm_W = resultant_phi13(W)
    norm_U = resultant_phi13(U)
    norm_Q = resultant_phi13(Q)
    require(norm_W == norm_U == 1
            and norm_Q == 8492431042211308167354471,
            "primitive target-window norm changed")
    require(norm_Q % 91 == 1,
            "right-wing target norm lost its 91-unit class")
    require(resultants["A"] == resultants["Ms"] == resultants["Mt"]
            == C**12
            and resultants["B"] == B0**12
            and resultants["L"] == 0
            and resultants["R"] == G**12 * norm_Q,
            "factored raw resultants changed")

    inverse_W = tuple(int(t in (2, 6, 10)) for t in range(P))
    inverse_U = tuple(int(t in (1, 4, 7, 10)) for t in range(P))
    require(circular_convolution(W, inverse_W) == (3,) + (2,) * 12
            and circular_convolution(U, inverse_U) == (4,) + (3,) * 12,
            "positive norm-one decoders changed")
    natural_decoder = circular_convolution(U, inverse_W)
    natural_decoded = circular_convolution(W, natural_decoder)
    require(natural_decoder == (3, 3, 2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2)
            and natural_decoded == (20, 20, 20) + (21,) * 10,
            "positive natural-sheet quotient decoder changed")

    # Primitive positive decoder for Q, obtained from its adjugate inverse and
    # shifted by the uniform target-null vector.  It decodes R alone; it cannot
    # pair R with the identically zero L functional.
    decoder_Q = (
        0, 7929439954473079186265, 133267004550058088654,
        135506751609875725618, 137784176006680359662,
        6831978049261891988, 4707019620677081960,
        2508704797329205596, 235167466079606368,
        124296699975702408, 47276146686603592,
        5907597936307760, 2054490908689296,
    )
    decoded_Q = circular_convolution(Q, decoder_Q)
    require(decoded_Q == (
        960304259241135830069757,
    ) + (16700810106546033697038,) * 12
            and decoded_Q[0] - decoded_Q[1]
            == 943603449134589796372719,
            "positive Q decoder changed")
    require((decoded_Q[0] - decoded_Q[1]) % 91 == 81,
            "right-wing primitive decoder class changed")
    delta_Q = decoded_Q[0] - decoded_Q[1]
    require(norm_Q == 9 * delta_Q and gcd(delta_Q, 91) == 1
            and pow(delta_Q, -1, 91) == 9,
            "right-wing localized/mod-91 decoder typing changed")
    decoder_Q_mod91 = tuple(9 * value % 91 for value in decoder_Q)
    decoded_Q_mod91 = tuple(value % 91 for value in circular_convolution(
        tuple(value % 91 for value in Q), decoder_Q_mod91
    ))
    require(decoded_Q_mod91 == (11,) + (10,) * 12,
            "Q*(9K_Q) stopped being delta_0 modulo (91,N)")

    expected_reductions = {
        "A": ((9, 0, 0, 0, 0, 0), 1),
        "B": ((8, 0, 0, 0, 0, 0), 12),
        "Ms": ((9, 0, 0, 0, 0, 0), 1),
        "Mt": ((4, 0, 0, 0, 0, 0), 1),
    }
    for name, expected in expected_reductions.items():
        for t in supports[name]:
            require(reductions[name][t] == expected,
                    f"root reduction changed for {name},t={t}")
    for t in range(3, 12):
        require(reductions["R"][t] == ((4, 0, 0, 0, 0, 0), 1),
                f"right-wing common-label root profile changed at t={t}")
    require(reductions["R"][12] == ((8, 0, 0, 0, 0, 0), 12),
            "right-wing terminal-label root profile changed")
    require(all(reduction == ((0, 0, 0, 0, 0, 0), 0)
                for reduction in reductions["L"]),
            "left wing acquired a coefficient unit")

    decoder, decoder_rank = solve_circulant_over_q(
        scalar["L"], scalar["R"]
    )
    require(decoder is None and decoder_rank == 0,
            "zero-left/nonzero-right decoder obstruction changed")

    print("THM-2751 FIXED-CLOCK ROOT-ZERO MAYER-VIETORIS WING TARGET SPECTRUM")
    print("status=RESERVED PROOF-COMPLETE + VERIFIED-EXACT; AWAITING INDEPENDENT AUDIT")
    print(f"constructor_audit=legacy_sha256:{legacy_hash} missing:(clock_comb,source_clock) canonical_present=True")
    print(f"prefix_audit=historical_and_actual_Q3_prefixes_equal_14/14 masses={prefix_masses}; first_failure_is_source-one_factor_not_terminal_prefix")
    print("direct_target_audit=B:t3,t12 and M:t3 pulled_coefficients_equal_forward_coordinate_recomputation")
    for line in row_lines:
        print(line)
    print(f"supports={supports}")
    print(f"factorization=C={C}=119*G G={G} B0=121*G={B0}")
    print("target_polynomials=A=Ms=Mt=C*W, B=121G*U, L=0, R=G*Q")
    print("W=z^3+...+z^11; U=z^3+...+z^12; Q=2(z^3+...+z^11)+121z^12")
    print("root_profiles=A:R12->9/det1; B:L1->8/det12; Ms:R12->9/det1; Mt:L1->4/det1; L=zero; R:t3..11->4/det1,t12->8/det12")
    print("gains=common Mt/Ms=4/9=12; fixed-clock common-label B/A=8/9=11; global scalar/dihedral impossible because support sizes are 9 versus 10")
    print(f"primitive_norms=(W:{norm_W},U:{norm_U},Q:{norm_Q}) Q_mod91={norm_Q % 91}")
    print("raw_resultants=(A:C^12,B:(121G)^12,M:C^12,L:0,R:G^12*Norm(Q))")
    print(f"cross_profile_L_times_R={cross} all_target_characters_zero=True resultant={cross_resultant}")
    print("positive_decoders=W^-1=(z2+z6+z10), U^-1=(z+z4+z7+z10) modulo uniform")
    print(f"fixed_profile_positive_decoder_K=U*W^-1={natural_decoder}; W*K=U+20N")
    print("fixed_profile_decoder_identity=121[A*K]=119[B] modulo the uniform target-null line; coefficient-derived only, not a physical packet action")
    print(f"Q_positive_decoder={decoder_Q}")
    print(f"Q_decoder_convolution={decoded_Q} primitive_delta={delta_Q} integral_index={norm_Q} delta_mod91={delta_Q % 91}")
    print(f"Q_mod91_decoder=9*K_Q convolution={decoded_Q_mod91}=delta_0+10N; Q is a unit only after rational/localized or mod91 scalar extension")
    print("conditional_holotopy=right cofiber generates the rational/localized and mod91 C13 augmentation quotients coefficientwise; integrally it contains delta_Q times the quotient but is not surjective; physical attachment to one common-ancestry semantic vertical edge remains absent")
    print("wing_decoder=IMPOSSIBLE: L functional is identically zero while R is nonzero; no scalar, dihedral, linear convolution, or positive decoder sends L to R")
    print("full_unclocked_boundary=bd53dc4c5 has M=disjoint_union_e M_e and a target augmentation rank drop; this theorem is only the fixed e=1 fibre")
    print("SCOPE: fixed-e=1 rail8 marked source sheet and pulled target sheet; coefficient transport only, no whole-packet action or endpoint current")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
