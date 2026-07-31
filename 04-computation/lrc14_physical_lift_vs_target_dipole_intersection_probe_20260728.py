#!/usr/bin/env python3
"""Exact typing probe for the missing two-index LRC14 carrier.

There are three different finite alphabets near the present frontier:

* ``n mod 13`` is the residue of a THM-2707 physical lift address;
* ``s,t`` are endpoint target-twist parameters;
* ``u,v`` are THM-2334 Fourier multiindices, whose dipole residues are
  selected by Fourier transforming in ``s,t``.

This script proves that the first alphabet cannot implement either nonzero
THM-2350 dipole action.  It also checks the enlarged residue-zero
following-atom macro inside the full THM-2707 packet fibre and gives two
finite non-identifiability controls showing why a one-sided shift table, or
a marginal which forgets ``u``, cannot determine ``eta.(u-v)``.

All packet tests use exact integer/Fraction arithmetic inherited from the
audited THM-2694/2707 companions.  The target-action no-go is a symbolic
divisibility proof; it does not depend on the packet census.
"""

from bisect import bisect_right
from collections import Counter
from fractions import Fraction
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_full_physical_lift_fibre_thm2707 as full


old = full.old
m = full.m

P = 13
R = P**6
H = 1
C1 = P
A = P**3
C = 2 * P**5
K_A = (27, 40, 53, 66)
K_B = 14

HALF_SUPPORTS = (
    (
        (27, (0, 3, 4, 5, 8, 9, 10, 11, 12)),
        (40, (0, 1, 2, 3, 4, 5, 8, 9, 10, 11)),
        (53, (0, 1, 2, 3, 4, 5, 8, 11, 12)),
        (66, (0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12)),
    ),
    (
        (27, (1, 2, 3, 4, 7, 8, 9, 10, 11)),
        (40, (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
        (53, (1, 2, 3, 6, 7, 8, 9, 10, 11)),
        (66, (1, 2, 3, 4, 5, 6, 7, 8, 9)),
    ),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def frac(value):
    return value - value.numerator // value.denominator


def strict_interval_index(value, starts, intervals):
    index = bisect_right(starts, value) - 1
    return index >= 0 and intervals[index][0] < value < intervals[index][1]


def cyclotomic_reduction(coefficients):
    """Reduce degree-at-most-12 coefficients modulo Phi_13."""
    require(len(coefficients) == P, "wrong cyclotomic vector length")
    top = coefficients[-1]
    return tuple(value - top for value in coefficients[:-1])


def primitive_character_reductions(values):
    """Exact reductions of sum_s values[s] zeta^(h*s), h=1,...,12."""
    require(len(values) == P, "wrong finite-function length")
    reductions = []
    for frequency in range(1, P):
        coefficients = [0] * P
        for shift, value in enumerate(values):
            coefficients[(frequency * shift) % P] += value
        reductions.append(cyclotomic_reduction(coefficients))
    return tuple(reductions)


def dipole_intersection_certificate(blocker, graft):
    """Prove diagonal physical translations meet a target dipole trivially.

    A physical translation by k/R moves the two factor phases by
    blocker*k/R and graft*k/R.  A target dipole moves them by +s/13 and
    -s/13.  Adding the congruences gives

        (blocker+graft) k = 0 mod R.

    Here blocker+graft is a 13-adic unit, so k=0 and then s=0.
    """
    pair_sum = blocker + graft
    require(blocker % P == 0 and graft % P == 1,
            "canonical dipole residue types changed")
    require(pair_sum % P == 1 and gcd(pair_sum, R) == 1,
            "paired-speed sum stopped being a 13-adic unit")

    # ALL-depth form.  On C_(13^r) x C_(13^r), the phase map
    #
    #     (x,y) -> (blocker*x+y, graft*x-y)
    #
    # has determinant -(blocker+graft), a 13-adic unit.  It is therefore an
    # automorphism at every r>=1.  Restricting y to the embedded order-13
    # subgroup cannot create a kernel.  The finite depth list below is only a
    # regression control; the gcd argument is the all-r proof.
    determinant = -pair_sum
    require(gcd(determinant, P) == 1,
            "phase-map determinant stopped being a 13-adic unit")
    checked_depths = tuple(range(1, 9))
    checked_kernel_sizes = []
    for depth in checked_depths:
        modulus = P**depth
        require(gcd(determinant, modulus) == 1,
                "phase map lost invertibility at a checked depth")
        inverse = pow(determinant % modulus, -1, modulus)
        require((determinant * inverse) % modulus == 1,
                "phase-map determinant inverse failed")
        checked_kernel_sizes.append(1)

    # The pair-sum proof forces k=0.  Check every target colour at that
    # forced value, rather than scanning 13^6 physical translations.
    pair_solutions = []
    rhs_unit = R // P
    forced_k = 0
    for s in range(P):
        if ((blocker * forced_k - s * rhs_unit) % R == 0
                and (graft * forced_k + s * rhs_unit) % R == 0):
            pair_solutions.append((s, forced_k))
    require(pair_solutions == [(0, 0)],
            "a nontrivial physical/target dipole intersection appeared")

    # Independent neutral-anchor derivation.  The guard coordinate alone
    # requires R|k.  Even after deliberately dropping it, target-neutral c1
    # requires 13^5|k; the blocker phase is then integral, so s=0.
    require(gcd(H, R) == 1, "guard stopped being primitive")
    guard_forced_step = R // gcd(H, R)
    source_forced_step = R // gcd(C1, R)
    require(guard_forced_step == R and source_forced_step == P**5,
            "neutral-anchor divisibility changed")
    source_candidates = tuple(h * source_forced_step for h in range(P))
    require(all((C1 * k) % R == 0 for k in source_candidates),
            "source-neutral candidates changed")
    require(all((blocker * k) % R == 0 for k in source_candidates),
            "source anchor no longer kills the blocker target phase")

    # Positive control: if the oppositely moving graft coordinate is erased,
    # the blocker coordinate alone has gcd(blocker,R) physical solutions for
    # every target colour.  Thus the paired coordinate is load-bearing.
    one_coordinate_counts = []
    divisor = gcd(blocker, R)
    for s in range(P):
        rhs = s * rhs_unit
        require(rhs % divisor == 0,
                "single-coordinate congruence unexpectedly became insoluble")
        one_coordinate_counts.append(divisor)

    return (
        blocker,
        graft,
        pair_sum,
        gcd(pair_sum, R),
        determinant % P,
        tuple(zip(checked_depths, checked_kernel_sizes)),
        tuple(pair_solutions),
        guard_forced_step,
        source_forced_step,
        tuple(one_coordinate_counts),
    )


def full_scc_following_macro_control():
    """Rebuild the full address universe and the residue-zero atom macro."""
    x = Fraction(649039434905733, 1304692766858936)
    z = Fraction(46873542509301, 100360982066072)
    radius = Fraction(1, 1304692766858936)
    require(frac(13 * x) == z, "THM-2694 D endpoint changed")

    module, prefixes, _, _, rails, present, starts = (
        m.core.build_carrier_data()
    )
    rows = m.shard((0, 1))[6][0]
    good_residues = []
    for residue in range(P):
        carry = (2 + 7 * residue) % P
        root = (6 + residue) % P
        if root and m.is_unit(rows[0][0][carry][1][6], root, 26):
            good_residues.append(residue)
    good_residues = tuple(good_residues)
    require(good_residues == (0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12),
            "THM-2707 private unit residue bank changed")

    rail_support = old.merge_support(rails[0][3])
    present_support = tuple(present[1, (-6) % P])
    following = old.d.build_atom(
        module, prefixes, present, starts, rails[0], 1, 6, 1, 1
    )
    following_support = old.merge_support(following["pieces"])
    prefix_starts, prefix_lengths = prefixes[
        following["future"]
    ][following["h"]][:2]
    following_delayed_support = tuple(
        (Fraction(left), Fraction(left + length))
        for left, length in zip(prefix_starts, prefix_lengths)
    )

    denominator = (z * m.T).denominator
    base_numerator = (z * m.T).numerator
    modulus = m.T * denominator
    step = (7 * m.T // R) * denominator
    require(m.T % R == 0 and step * R == 7 * m.T * denominator,
            "THM-2707 address grid changed")

    def scale_support(support):
        intervals = tuple(
            (left * denominator, right * denominator)
            for left, right in support
        )
        return intervals, tuple(left for left, _ in intervals)

    scaled_rail, rail_starts = scale_support(rail_support)
    scaled_present, present_starts = scale_support(present_support)
    scaled_following, following_starts = scale_support(following_support)

    good = []
    following_hits = []
    point = base_numerator
    for address in range(R):
        if (address % P in good_residues
                and strict_interval_index(
                    point, rail_starts, scaled_rail
                )
                and strict_interval_index(
                    point, present_starts, scaled_present
                )):
            good.append(address)
            if strict_interval_index(
                    point, following_starts, scaled_following):
                following_hits.append(address)
        point = (point + step) % modulus

    residue_counts = Counter(address % P for address in good)
    require(len(good) == 3346, "THM-2707 packet census changed")
    require(dict(sorted(residue_counts.items())) == {
        0: 304, 1: 305, 2: 304, 3: 305, 4: 304, 5: 305,
        6: 304, 9: 301, 10: 304, 11: 305, 12: 305,
    }, "THM-2707 residue census changed")

    residue_zero = tuple(address for address in good if address % P == 0)
    require(tuple(following_hits) == residue_zero,
            "following atom is no longer exactly the residue-zero part")

    # The q-address is evaluated after x -> {13x}; hence an x-radius `radius`
    # becomes the base-grid radius below.  Every midpoint hit retains the
    # entire inherited open interval, not just its centre.
    base_radius = 13 * m.T * radius
    whole_following = []
    for address in following_hits:
        q = frac(z + Fraction(7 * address, R)) * m.T
        if old.open_arc_is_contained(
                q, base_radius, following_support, m.T):
            whole_following.append(address)
    require(tuple(whole_following) == residue_zero,
            "a residue-zero following atom lost the whole common interval")

    # The base support above is only one half of the THM-2680 atom.  Its
    # delayed prefix is invariant under q_n={13x}+7n/R because multiplication
    # by R turns the address displacement into the integer 7n.  Check the
    # whole pulled-back x-cylinder against that prefix as well.
    delayed_center = frac(R * z) * m.T
    delayed_radius = 13 * R * m.T * radius
    require(old.open_arc_is_contained(
        delayed_center,
        delayed_radius,
        following_delayed_support,
        m.T,
    ), "the frozen following delayed prefix lost the whole common interval")

    other_residue = tuple(address for address in good if address % P != 0)
    require(len(other_residue) == 3042,
            "other-residue packet count changed")
    require(all((7 * address) % P != 0 for address in other_residue),
            "an alleged nonzero quotient lift became residue-trivial")
    rooted_two_step_loops = len(residue_zero) * len(other_residue)
    require(rooted_two_step_loops == 924768,
            "residue-zero two-step macro count changed")

    return (
        len(good),
        tuple(sorted(residue_counts.items())),
        len(residue_zero),
        len(other_residue),
        rooted_two_step_loops,
        True,
        residue_zero[:5],
        residue_zero[-5:],
    )


def difference_pushforward(table):
    """Push a joint residue table (u,v) to q=u-v."""
    require(len(table) == P and all(len(row) == P for row in table),
            "wrong joint residue table shape")
    result = [0] * P
    for u in range(P):
        for v in range(P):
            result[(u - v) % P] += table[u][v]
    return tuple(result)


def nonidentifiability_controls():
    """Check primal-marginal and dual fixed-left non-identifiability."""
    rows = []
    for event_index, event_rows in enumerate(HALF_SUPPORTS):
        for pivot, support_tuple in event_rows:
            support = set(support_tuple)
            require(0 < len(support) < P,
                    "half-cycle support stopped being proper")
            profile = tuple(int(v in support) for v in range(P))
            primitive = primitive_character_reductions(profile)
            require(all(any(reduction) for reduction in primitive),
                    "a half-cycle primitive moving character vanished")

            # Primal hostile: after *abstractly* forgetting a putative left
            # residue u, the same v-marginal admits a diagonal lift u=v and a
            # fixed-left lift u=0.  These are information-theoretic controls,
            # not claims that the shift support itself is a Fourier residue.
            diagonal = [[0] * P for _ in range(P)]
            fixed_left = [[0] * P for _ in range(P)]
            for v in range(P):
                diagonal[v][v] = profile[v]
                fixed_left[0][v] = profile[v]
            diagonal_marginal = tuple(
                sum(diagonal[u][v] for u in range(P)) for v in range(P)
            )
            fixed_marginal = tuple(
                sum(fixed_left[u][v] for u in range(P)) for v in range(P)
            )
            require(diagonal_marginal == fixed_marginal == profile,
                    "the two residue lifts lost their common v-marginal")
            q_diagonal = difference_pushforward(diagonal)
            q_fixed = difference_pushforward(fixed_left)
            require(q_diagonal[0] == len(support)
                    and all(q_diagonal[q] == 0 for q in range(1, P)),
                    "diagonal residue lift acquired a nonzero difference")
            require(any(q_fixed[q] for q in range(1, P)),
                    "fixed-left residue lift lost every nonzero difference")

            # Dual hostile faithful to what the existing calculation actually
            # fixes: the row at unshifted left target parameter s=0.  Two
            # nonnegative 13x13 twist tables have that identical row.  Filling
            # the otherwise unseen diagonal reverses the primitive-character
            # verdict.  This is why the missing object is a joint endpoint
            # orbit, not another right-shift calculation.
            base = [[0] * P for _ in range(P)]
            filled = [[0] * P for _ in range(P)]
            for t in range(P):
                base[0][t] = profile[t]
                filled[0][t] = profile[t]
            for s in range(1, P):
                filled[s][s] = 1
            require(tuple(base[0]) == tuple(filled[0]) == profile,
                    "dual controls lost their common fixed-left row")
            base_diagonal = tuple(base[s][s] for s in range(P))
            filled_diagonal = tuple(filled[s][s] for s in range(P))
            base_characters = primitive_character_reductions(base_diagonal)
            filled_characters = primitive_character_reductions(
                filled_diagonal
            )
            base_all_nonzero = all(any(value) for value in base_characters)
            filled_all_nonzero = all(
                any(value) for value in filled_characters
            )
            require(base_all_nonzero != filled_all_nonzero,
                    "unseen diagonal stopped reversing the target verdict")

            rows.append((
                event_index,
                pivot,
                support_tuple,
                int(0 in support),
                len(support),
                sum(q_diagonal[q] != 0 for q in range(1, P)),
                sum(q_fixed[q] != 0 for q in range(1, P)),
                base_all_nonzero,
                filled_all_nonzero,
            ))
    require(len(rows) == 8, "half-cycle control-row count changed")
    return tuple(rows)


def main():
    require((P, R, H, C1, A, C, K_A, K_B) == (
        13, 4826809, 1, 13, 2197, 742586,
        (27, 40, 53, 66), 14,
    ), "canonical typed constants changed")

    dipoles = tuple(
        dipole_intersection_certificate(A, pivot) for pivot in K_A
    ) + (dipole_intersection_certificate(C, K_B),)
    macro = full_scc_following_macro_control()
    controls = nonidentifiability_controls()

    print("LRC14 PHYSICAL-LIFT / TARGET-DIPOLE INTERSECTION EXACT PROBE")
    print("status=FINITE-EXACT typing no-go plus positive packet-macro control")
    print(f"universe=physical_addresses n in Z/{R}Z; target_twists "
          f"s,t in F_{P}; abstract projected endpoint residues "
          f"ubar=eta.u,vbar=eta.v in F_{P}")
    print("type_dictionary=packet_address:n; packet_residue:n_mod_13; "
          "carry/root:event labels; s,t:dual target-twist parameters; "
          "u,v:global Fourier multiindices; "
          "ubar,vbar:abstract target projections; "
          "q=ubar-vbar=eta.(u-v)")
    print("pointwise_warning=u,v are global Fourier summation indices, not "
          "measurable labels of x, a runner, a packet vertex, a carry, or a "
          "private root")
    print("dipole_equations=blocker*k/R=s/13 and "
          "graft*k/R=-s/13 mod 1")
    print("pair_sum_derivation=R|(blocker+graft)k and "
          "gcd(blocker+graft,R)=1 force k=0, then s=0")
    print("all_depth_phase_map=(x,y)->(blocker*x+y,graft*x-y) on "
          "C_(13^r)^2 has determinant -(blocker+graft), a 13-adic unit; "
          "hence its kernel is trivial for every r>=1")
    print("kakeya_reading=the physical diagonal-translation needle and the "
          "opposite target-dipole needle are transverse in the two-factor "
          "phase torus; deeper 13-adic toothpicks do not create an "
          "intersection")
    print("neutral_anchor_derivation=guard H=1 forces R|k; even without H, "
          "target-neutral c1=13 forces 13^5|k and then each blocker phase "
          "is integral, forcing s=0")
    print("dipole_certificates=(blocker,graft,sum,gcd,det_mod13," 
          "checked_depth_kernel_sizes,solutions,guard_step,source_step," 
          "single_coordinate_solution_counts)="
          f"{dipoles}")
    print("single_coordinate_control=after deleting the oppositely moving "
          "graft, every target colour has gcd(blocker,R) physical solutions")
    print("full_scc_control=(packet_count,residue_counts,residue0_following," 
          "other_residue_packets,rooted_two_step_loops," 
          "delayed_prefix_whole,first5,last5)="
          f"{macro}")
    print("following_macro_law=each of 304 residue-zero endpoints retains the "
          "whole frozen following atom on I and has a physical two-step loop "
          "through each of 3042 other-residue packets")
    print("macro_target_boundary=every nontrivial macro edge is a global "
          "physical translation outside both nonzero target-dipole "
          "subgroups; the closed two-step sum is target-neutral, not a "
          "definition of eta.u or eta.v")
    print("nonidentifiability_rows=(event,pivot,support,zero_present,size," 
          "diag_nonzero_q_count,fixed_left_nonzero_q_count," 
          "base_row_target_chars_nonzero,filled_diag_target_chars_nonzero)="
          f"{controls}")
    print("two_controls=primal: diagonal u=v versus fixed u=0 have the same "
          "forgotten-u v-marginal but opposite nonzero-q support; dual: two "
          "joint twist tables have the same observed left-zero row but "
          "opposite primitive common-diagonal verdict")
    print("minimal_sidecar=a factorwise THM2350-covariant joint endpoint "
          "orbit on one physical ancestry, a selected pivot, and a fixed "
          "marked triangle X,m,Y (or proved normal/jet weight) whose "
          "bi-Fourier coefficient survives the pair-to-target line sum")
    print("SCOPE=no physical endpoint current, scalar-row exclusion, ledger "
          "decrement, or LRC14 conclusion")


if __name__ == "__main__":
    main()
