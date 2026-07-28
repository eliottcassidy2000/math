#!/usr/bin/env python3
"""Exact probe for the filtered affine-handoff germ groupoid.

This is the exact companion for the THM-2697 candidate.  It combines three
independently proved layers without promoting any implication between them:

* the BS(1,13)-modeled affine normal form on the circle's local inverse
  branches and its distinct forward ``C_R semidirect N`` quotient;
* the intrinsic three-tooth rail-clock support from THM-2689; and
* the frozen physical mixed D/slope witness from the direct referee.

The load-bearing new coordinate is the 13-adic depth of a root translation.
Pushing a translation through D does not erase its root digit: it raises the
digit by one filtration level.  Forgetting that level creates the apparent
root loss.  The smallest nonkernel extension is ``C_(2R)``: its parity
quotient changes the delayed base map from ``13y`` to ``13y+1/2``.
"""

from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys


HERE = Path(__file__).resolve()
ROOT = next(
    candidate for candidate in HERE.parents
    if (candidate / "01-canon").is_dir()
    and (candidate / "04-computation").is_dir()
)
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_affine_three_tooth_full_clock_classification as affine
import lrc14_mixed_d_slope7_two_simplex_independent_referee as mixed
import lrc14_predecessor_carry_private_root_atlas_thm2640 as packet
import lrc14_dilation_reversed_clock_fibre_product_probe as d_edge
import lrc14_odometer_alternating_lift_labelled_tail_scout_20260728 as tail


P = 13
N = 6
R = P ** N
INV2 = 7


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compose(outer, inner):
    """Compose A_l D^m after A_k D^n in forward affine normal form."""
    m, ell = outer
    n, k = inner
    return m + n, (ell + P ** m * k) % R


def affine_coeff(normal):
    """Return slope/intercept of A_k D^n on the universal cover."""
    n, k = normal
    return P ** n, F(k, R)


def depth(k):
    """13-adic depth in Z/13^6, with 6 reserved for zero."""
    k %= R
    if k == 0:
        return N
    j = 0
    while j < N and k % P == 0:
        j += 1
        k //= P
    return j


def root_digit(k, j):
    k %= R
    require(j < N and k % (P ** j) == 0,
            "root digit requested outside its filtration level")
    return (2 * (k // (P ** j))) % P


def shift(j, delta):
    """Generator-linear raw translation numerator with graded root delta."""
    return (INV2 * P ** j * delta) % R


def circle_fraction(value):
    """Canonical representative in [0,1) for an exact rational."""
    return value - (value.numerator // value.denominator)


def apply_normal(normal, value):
    """Apply A_k D^n on the circle."""
    n, k = normal
    return circle_fraction(P ** n * value + F(k, R))


def delayed_base(value):
    """The THM-2693 terminal coordinate pi(x)={13^6 x}."""
    return circle_fraction(R * value)


def base_iterate(value, depth_value):
    """Iterate B_0(y)={13y} exactly."""
    return circle_fraction(P ** depth_value * value)


def circle_distance(value):
    value = circle_fraction(value)
    return min(value, 1 - value)


def translate_integer_union(intervals, offset, period):
    """Translate an exact half-open integer-grid circle union."""
    out = []
    for old_left, old_right in intervals:
        length = old_right - old_left
        left = (old_left + offset) % period
        right = left + length
        if right <= period:
            out.append((left, right))
        else:
            out.append((left, period))
            out.append((0, right - period))
    return tail.merge_intervals(out)


def translate_fraction_union(intervals, offset):
    """Translate an exact half-open Fraction circle union."""
    out = []
    for old_left, old_right in intervals:
        length = old_right - old_left
        left = circle_fraction(old_left + offset)
        right = left + length
        if right <= 1:
            out.append((left, right))
        else:
            out.append((left, F(1)))
            out.append((F(0), right - 1))
    return tail.merge_intervals(out)


def algebra_audit():
    # A_k C = C A_(13k), equivalently D A_k = A_(13k) D.
    bs_checks = 0
    for k in range(-24, 25):
        # Compare universal-cover affine coefficients, not decimal samples.
        left = (F(1, P), F(k, R))       # A_k after contraction C
        right = (F(1, P), F(13 * k, P * R))
        require(left == right, "BS(1,13) local-germ relation failed")
        require(compose((1, 0), (0, k)) == (1, (13 * k) % R),
                "forward D/translation relation failed")
        bs_checks += 1

    assoc_checks = 0
    samples = (-371294, -91, -14, -7, -1, 0, 1, 7, 14, 91, 371294)
    for a in range(3):
        for b in range(3):
            for c in range(3):
                for k in samples:
                    for ell in samples:
                        for q in (-14, 0, 14):
                            x = (a, k % R)
                            y = (b, ell % R)
                            z = (c, q % R)
                            require(compose(z, compose(y, x))
                                    == compose(compose(z, y), x),
                                    "affine normal-form composition lost associativity")
                            assoc_checks += 1

    carry_checks = 0
    for j in range(N):
        for delta in range(1, P):
            k = shift(j, delta)
            require(depth(k) == j and root_digit(k, j) == delta,
                    "graded root section changed")
            if j + 1 < N:
                require((P * k) % R == shift(j + 1, delta),
                        "D stopped raising root depth by one")
                require(root_digit(P * k, j + 1) == delta,
                        "graded root digit changed under D")
            else:
                require((P * k) % R == 0,
                        "top root depth did not expire at depth six")
            for epsilon in range(P):
                total = delta + epsilon
                low = total % P
                carry = total // P
                lhs = (shift(j, delta) + shift(j, epsilon)) % R
                rhs = (shift(j, low)
                       + (shift(j + 1, carry) if j + 1 < N else 0)) % R
                require(lhs == rhs, "same-depth odometer carry law failed")
                carry_checks += 1

    interchange = []
    for delta in range(1, P):
        k = shift(0, delta)
        post = compose((0, k), (1, 0))   # A_k after D
        pre = compose((1, 0), (0, k))    # D after A_k
        require(post == (1, k) and pre == (1, (P * k) % R),
                "pre/post-D normal forms changed")
        require(root_digit(post[1], 0) == delta,
                "post-D shift lost its depth-zero root")
        require(depth(pre[1]) == 1 and root_digit(pre[1], 1) == delta,
                "pre-D shift did not become the same depth-one root")
        commutator = (post[1] - pre[1]) % R
        require(root_digit(commutator, 0) == delta,
                "oriented interchange defect lost the quotient root")
        interchange.append((delta, post[1], pre[1], commutator))

    return bs_checks, assoc_checks, carry_checks, tuple(interchange)


def base_functor_audit():
    """Audit pi o (A_k D^n)=B_0^n o pi and the inherited depth-four zero."""
    samples = (
        F(-17, 19), F(-1, R), F(0), F(1, R), F(5, 17),
        F(649039434905733, 1304692766858936), F(41, 43),
    )
    numerators = (-R - 1, -91, -7, -1, 0, 1, 7, 91, R + 1)
    functor_checks = 0
    for n in range(5):
        for k in numerators:
            for value in samples:
                require(
                    delayed_base(apply_normal((n, k % R), value))
                    == base_iterate(delayed_base(value), n),
                    "affine normal form stopped projecting to its base degree",
                )
                functor_checks += 1

    kernel_checks = 0
    for j in range(N):
        for delta in range(1, P):
            for value in samples:
                require(
                    delayed_base(apply_normal((0, shift(j, delta)), value))
                    == delayed_base(value),
                    "a graded slope/root translation moved the delayed base",
                )
                kernel_checks += 1

    # Replay THM-2693's sparse inherited grammar using only the target tooth
    # D_(13^3) and the speed-14 safe factor.  This is the exact base object
    # through which every integer-lift/slope-decorated chronology word factors.
    module, _, _, _, _, _, _ = packet.core.build_carrier_data()
    require(module.C2 == P ** 3 and module.W[module.UNIT_IDX[0]] == P + 1,
            "THM-2693 sparse delayed grammar changed")
    minimal_base = module.make_comb(module.C2, 182, -13, 13)
    minimal_base = module.subtract_comb(
        minimal_base, module.W[module.UNIT_IDX[0]], 182, -13, 13
    )
    minimal_fraction = [
        (F(left, tail.T), F(right, tail.T))
        for left, right in minimal_base
    ]
    expected = (
        (1, 1886, F(6, 49)),
        (2, 1606, F(22187, 2798978)),
        (3, 1452, F(1452, 2599051)),
        (4, 0, F(0)),
    )
    rows = []
    for event_depth in range(1, 5):
        support = tail.fraction_repeated_tail(minimal_fraction, event_depth)
        row = (event_depth, len(support), tail.interval_mass(support))
        require(row == expected[event_depth - 1],
                "THM-2693 sparse delayed nilpotence replay changed")
        rows.append(row)

    return functor_checks, kernel_checks, tuple(rows)


def parity_sidecar_audit():
    """Audit the C_(2R) parity quotient and its recurrent raw delayed word."""
    samples = (
        F(-17, 19), F(-1, 2 * R), F(0), F(1, 2 * R), F(5, 17),
        F(649039434905733, 1304692766858936), F(41, 43),
    )
    parity_checks = 0
    for n in range(-24, 25):
        epsilon = n % 2
        for value in samples:
            # D tau_n=tau_(13n)D on the circle.
            lhs = circle_fraction(P * circle_fraction(value + F(n, 2 * R)))
            rhs = circle_fraction(P * value + F(13 * n, 2 * R))
            require(lhs == rhs, "half-grid D covariance changed")
            # Post-D tau_n has base map B_epsilon(y)=13y+epsilon/2.
            target = circle_fraction(P * value + F(n, 2 * R))
            require(
                delayed_base(target)
                == circle_fraction(P * delayed_base(value) + F(epsilon, 2)),
                "C_(2R) parity stopped controlling the delayed base map",
            )
            # The central half-circle representative commutes with odd D.
            require(
                circle_fraction(P * circle_fraction(value + F(1, 2)))
                == circle_fraction(P * value + F(1, 2)),
                "central parity representative stopped commuting with D",
            )
            parity_checks += 1

    handoff_checks = 0
    for prior_numerator in (-2 * R - 1, -91, -1, 0, 1, 91, 2 * R + 1):
        for new_numerator in range(-24, 25):
            target_numerator = new_numerator + P * prior_numerator
            for value in samples:
                lhs = circle_fraction(
                    P * circle_fraction(
                        value + F(prior_numerator, 2 * R)
                    ) + F(new_numerator, 2 * R)
                )
                rhs = circle_fraction(
                    P * value + F(target_numerator, 2 * R)
                )
                require(lhs == rhs,
                        "tau_new D tau_prior chronology handoff changed")
                require(target_numerator % 2
                        == (new_numerator + prior_numerator) % 2,
                        "handoff parity stopped adding across chronology")
                handoff_checks += 1

    digit_checks = 0
    for c in range(P):
        for h in range(P):
            for kappa in range(2):
                digit = 2 * h + kappa
                upper = digit // P
                for edge in range(2):
                    absolute = (2 * c + upper + (edge == 0)) % P
                    for k in range(-13, 14):
                        n = 2 * k + 1
                        digit2 = (digit + P) % (2 * P)
                        upper2 = digit2 // P
                        carry2 = (c + k + upper) % P
                        absolute2 = (
                            2 * carry2 + upper2 + (edge == 0)
                        ) % P
                        require(upper2 == 1 - upper,
                                "odd half-grid shift stopped toggling the half digit")
                        require(absolute2 == (absolute + n) % P,
                                "THM-2640 labels stopped reconstructing the half root")
                        digit_checks += 1

    module, _, _, _, _, _, _ = packet.core.build_carrier_data()
    sector_words = tail.prior.sector_words(module)
    raw = tail.merge_intervals(sector_words[0] + sector_words[1])
    fixed_rows = []
    expected_distances = (
        F(1, 24),
        (F(5, 12), F(3, 8), F(1, 3), F(7, 24), F(1, 4)),
        F(1, 12),
    )
    for value in (F(11, 24), F(13, 24)):
        require(circle_fraction(P * value + F(1, 2)) == value,
                "parity delayed fixed point changed")
        require(tail.inside_unweighted(value * tail.T, raw),
                "parity fixed point left the raw delayed word")
        sectors = tuple(
            sector for sector in range(2)
            if tail.inside_unweighted(value * tail.T, sector_words[sector])
        )
        distances = (
            circle_distance(module.C2 * value),
            tuple(circle_distance(module.W[index] * value)
                  for index in module.UNIT_IDX),
            circle_distance(module.C3 * value),
        )
        require(distances == expected_distances,
                "parity fixed-point delayed factor word changed")
        require(distances[0] < F(1, 14)
                and min(distances[1] + (distances[2],)) > F(1, 14),
                "parity fixed point lost strict danger/safety typing")
        fixed_rows.append((value, sectors, distances))

    # Exact endpoint composite of the independently replayed half-odometer
    # two-cycle.  Two odd edges cancel parity and descend to integral D^2
    # macro arrows, but the forced B_0 midpoint leaves W.
    x0 = F(55_232_507, 115_843_416)
    x1 = F(58_313_459, 115_843_416)
    k0, k1 = 1_472_973, 4_502_560
    n0, n1 = 2 * k0 + 1, 2 * k1 + 1
    require(circle_fraction(P * x0 + F(n0, 2 * R)) == x1
            and circle_fraction(P * x1 + F(n1, 2 * R)) == x0,
            "half-odometer endpoint two-cycle changed")
    macro0 = ((n1 + P * n0) // 2) % R
    macro1 = ((n0 + P * n1) // 2) % R
    require((n1 + P * n0) % 2 == 0
            and (n0 + P * n1) % 2 == 0,
            "two odd half edges stopped descending to integral macro arrows")
    require(circle_fraction(P ** 2 * x0 + F(macro0, R)) == x0
            and circle_fraction(P ** 2 * x1 + F(macro1, R)) == x1,
            "integral D^2 endpoint macro loops changed")
    require((macro0, macro1, 2 * macro0 % P, 2 * macro1 % P)
            == (4_343_980, 2_084_552, 8, 4),
            "integral macro numerator/root-step bank changed")
    y0, y1 = delayed_base(x0), delayed_base(x1)
    require((y0, y1) == (F(11, 24), F(11, 24)),
            "two-cycle endpoints stopped sharing the fixed delayed phase")
    midpoint0 = base_iterate(y0, 1)
    reflected_midpoint = base_iterate(F(13, 24), 1)
    require((midpoint0, reflected_midpoint) == (F(23, 24), F(1, 24)),
            "forced B_0 midpoint pair changed")
    require(not tail.inside_unweighted(midpoint0 * tail.T, raw)
            and not tail.inside_unweighted(reflected_midpoint * tail.T, raw),
            "an integral B_0 subdivision midpoint entered W")
    macro_rows = (
        (x0, macro0, 2 * macro0 % P, y0, midpoint0),
        (x1, macro1, 2 * macro1 % P, y1, midpoint0),
    )

    parity_word_rows = []
    parity_support = None
    refined_grid = None
    for bits in product(range(2), repeat=3):
        phases = [F(0)]
        cumulative = F(0)
        for bit in bits:
            cumulative = circle_fraction(P * cumulative + F(bit, 2))
            phases.append(cumulative)
        integer_sets = [
            translate_integer_union(raw, int(-phase * tail.T), tail.T)
            for phase in phases
        ]
        support, grid = tail.integer_itinerary(integer_sets)
        mass = F(sum(right - left for left, right in support), grid)
        parity_word_rows.append((bits, len(support), mass))
        if bits == (1, 1, 1):
            parity_support = support
            refined_grid = grid
            require((len(support), mass)
                    == (23_496, F(64_614, 5_710_115_047)),
                    "four-state parity delayed cylinder census changed")
        else:
            require(not support and mass == 0,
                    "a non-all-odd three-edge parity word acquired support")
    require(parity_support is not None and refined_grid is not None,
            "all-odd parity word was not scanned")

    phases = (F(0), F(1, 2), F(0), F(1, 2))
    integer_support = parity_support
    integer_mass = F(
        sum(right - left for left, right in integer_support), refined_grid
    )

    raw_fraction = [(F(left, tail.T), F(right, tail.T))
                    for left, right in raw]
    fraction_sets = [
        translate_fraction_union(raw_fraction, -phase) for phase in phases
    ]
    fraction_support = fraction_sets[0]
    for index in range(1, 4):
        fraction_support = tail.fraction_pullback_intersection(
            fraction_support, fraction_sets[index], P ** index
        )
    require([
        (int(left * refined_grid), int(right * refined_grid))
        for left, right in fraction_support
    ] == integer_support,
            "Fraction/integer parity delayed replays disagree")

    return (parity_checks, handoff_checks, digit_checks, tuple(fixed_rows),
            tuple(parity_word_rows), len(integer_support), integer_mass,
            macro_rows, (midpoint0, reflected_midpoint))


def intrinsic_affine_audit():
    require(not affine.chain_components(F(0), 1),
            "pure unshifted D unexpectedly acquired a three-event corridor")
    rows = []
    for delta in range(1, P):
        beta = F(shift(0, delta), R)
        components = affine.chain_components(beta, 1)
        require(len(components) == 2,
                "generator-section root lift changed its intrinsic component count")
        require(affine.measure(components) == F(14 * delta, P ** 8),
                "generator-section root-lift corridor mass changed")
        signatures = tuple(
            affine.component_signature(component, 1, beta)
            for component in components
        )
        require(tuple(row["arrivals"] for row in signatures)
                == ((0, 0, 6), (6, 6, 12)),
                "generator-section lift changed its two arrival words")
        require(tuple(row["stored_edges"] for row in signatures)
                == (((0, 3), (3, 0), (0, 3)),
                    ((4, 0), (0, 4), (4, 0))),
                "generator-section lift changed its local clock words")
        require((91 * shift(0, delta)) % R != 0,
                "nonzero root lift became globally seven-clock covariant")
        rows.append((delta, len(components), affine.measure(components),
                     signatures[0]["arrivals"], signatures[1]["arrivals"]))
    return tuple(rows)


def mixed_physical_audit():
    x = F(649039434905733, 1304692766858936)
    z = F(46873542509301, 100360982066072)
    require(mixed.frac(P * x) == z, "frozen D coordinate changed")

    module, prefixes, _, _, rails, present, starts = (
        packet.core.build_carrier_data()
    )
    current = d_edge.build_atom(
        module, prefixes, present, starts, rails[2], 3, 2, 0, 0
    )
    following = d_edge.build_atom(
        module, prefixes, present, starts, rails[0], 1, 6, 1, 1
    )
    require(mixed.atom_membership(current, x, prefixes)
            and mixed.atom_membership(following, z, prefixes),
            "frozen strict D edge lost an endpoint atom")

    pair_prefixes = packet.build_pair_prefixes(module)
    shard = packet.shard((0, 1))
    rows = shard[6][0]
    carry, root = 2, 6
    require(mixed.private_membership(
        module, pair_prefixes, rails, present, z, carry, root
    ), "frozen D target lost the private base configuration")

    successes = []
    failures = []
    retained_following = []
    next_rail_teeth = []

    def rail_tooth(value):
        value = mixed.frac(value)
        if F(0) <= value < F(1, 28):
            return 0
        if F(13, 28) <= value < F(15, 28):
            return 6
        if F(27, 28) <= value < F(1):
            return 12
        return None

    for delta in range(1, P):
        shifted = mixed.frac(z + F(shift(0, delta), R))
        next_rail_teeth.append(
            (delta, rail_tooth(shifted), rail_tooth(P * shifted))
        )
        if mixed.atom_membership(following, shifted, prefixes):
            retained_following.append(delta)
        carry2 = (carry + shift(0, delta)) % P
        root2 = (root + delta) % P
        require(carry2 == (carry + 7 * delta) % P,
                "carry/root section slope changed")
        if root2 == 0:
            failures.append((delta, "target_root_zero"))
            continue
        unit = packet.is_unit(rows[0][0][carry2][1][6], root2, 26)
        physical = mixed.private_membership(
            module, pair_prefixes, rails, present,
            shifted, carry2, root2,
        )
        if unit and physical:
            successes.append(delta)
        else:
            failures.append((delta, "unit" if not unit else "physical"))
    require(tuple(successes) == (1, 2, 3, 4, 5, 6, 9, 10, 11, 12),
            "frozen mixed physical success bank changed")
    require(tuple(failures) == ((7, "target_root_zero"), (8, "unit")),
            "frozen mixed physical hostile bank changed")
    require(not retained_following,
            "a nonzero slope shift unexpectedly retained the exact following D atom")
    require(tuple(next_rail_teeth)
            == tuple((delta, 6, None) for delta in range(1, P)),
            "the frozen mixed target acquired a further D-rail continuation")
    return (x, z, tuple(successes), tuple(failures),
            tuple(retained_following), tuple(next_rail_teeth))


def main():
    bs, assoc, carries, interchange = algebra_audit()
    functor_checks, kernel_checks, sparse_tail = base_functor_audit()
    (parity_checks, parity_handoffs, digit_checks, parity_fixed,
     parity_words, parity_components, parity_mass, parity_macros,
     parity_B0_midpoints) = parity_sidecar_audit()
    intrinsic = intrinsic_affine_audit()
    (x, z, successes, failures, retained_following,
     next_rail_teeth) = mixed_physical_audit()

    print("LRC14 FILTERED AFFINE-HANDOFF GROUPOID SCRATCH PROBE")
    print("ambient=BS(1,13)-modeled local germs; forward circle quotient=C_(13^6) semidirect N with A^R=identity")
    print("local_relations=A_k C = C A_(13k), D A_k = A_(13k) D")
    print(f"exact_algebra_checks=BS:{bs}:associativity:{assoc}:graded_carries:{carries}")
    print("forward_normal_form=A_k D^n; compose((m,l),(n,k))=(m+n,l+13^m*k mod13^6)")
    print("graded_root=rho_j(k)=2*(k/13^j) mod13 on 13^j/13^(j+1)")
    print("D_covariance=D S_(j,delta)=S_(j+1,delta) D; S_(6,delta)=identity")
    print(f"oriented_interchange_rows={interchange}")
    print("interchange_meaning=post-D S_(0,delta) retains depth0 root; pre-D S_(0,delta) becomes the same depth1 root")
    print(f"base_functor_checks=normal_forms:{functor_checks}:kernel_translations:{kernel_checks}")
    print("base_functor=pi(x)={13^6*x}; pi(A_k D^n x)=B_0^n(pi(x)); every integral root/slope translation lies in ker(pi)")
    print(f"sparse_inherited_delayed_tail={sparse_tail}")
    print("kernel_decoration_no_go=four old delayed event states remain empty under arbitrary integral slope/root insertions between three degree-one chronology edges")
    print(f"parity_sidecar_checks=covariance_and_base:{parity_checks}:chronology_handoffs:{parity_handoffs}:edge_typed_digit_law:{digit_checks}")
    print("parity_extension=tau_n(x)=x+n/(2R); post-D tau_n D has raw n, pre-D D tau_n=tau_(13n)D has raised raw13n; epsilon=n mod2 survives; pi(tau_n D x)=B_epsilon(pi(x))")
    print(f"parity_fixed_raw_delayed_points={parity_fixed}")
    print(f"parity_three_edge_word_bank={parity_words}")
    print(f"parity_four_state_raw_delayed_cylinders=components:{parity_components}:mass:{parity_mass}")
    print("parity_verdict=B_1(y)={13y+1/2} has strict fixed old-W loops at 11/24 and 13/24; static (c,h,kappa,e) labels reconstruct the half-root step, but no semantic transition follows")
    print(f"two_odd_edge_integral_D2_macro_loops={parity_macros}")
    print(f"forced_integral_B0_subdivision_midpoints_outside_W={parity_B0_midpoints}")
    print(f"generator_lift_intrinsic_corridors={intrinsic}")
    print("generator_lift_formula=each generator multiple delta has two components, total mass 14*delta/13^8, and no global seven-clock permutation")
    print(f"frozen_mixed_point=x:{x}:z=D(x):{z}")
    print(f"frozen_mixed_physical_successes={successes}")
    print(f"frozen_mixed_physical_failures={failures}")
    print(f"same_exact_following_D_atom_retained_after_slope_shift={retained_following}")
    print(f"frozen_mixed_target_and_next_D_rail_teeth={next_rail_teeth}")
    print("typed_verdict=the affine germ algebra composes exactly, but physical realization is a branch-refined partial mask/span: pure DD is empty, D-then-root-shift has ten strict packet charts, and fibre holonomy cannot fill the inherited delayed base's missing four-state simplex")
    print("scope=THM-2698 supplies the physical local parity cycle; no semantic terminal endpoint bibundle, common owner, row exclusion, or LRC14 conclusion")


if __name__ == "__main__":
    main()
