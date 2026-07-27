#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2531.

The script audits the prime C_13 necklace selector, its Mersenne-score
realization, affine/reflection/complement covariance, the dihedral tie
boundary, and its pointwise compatibility with the THM-2527/2528 Boolean
cut/path imbalance.
"""

from collections import Counter, defaultdict


P = 13
ALL = (1 << P) - 1
MERSENNE = ALL
checks = 0


def require(condition, label):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError("FAILED: " + label)


def bit(mask, r):
    return (mask >> (r % P)) & 1


def word(mask, a, tau=1):
    return tuple(bit(mask, a + j * tau) for j in range(P))


def score(mask, a, tau=1):
    return sum(bit(mask, a + j * tau) << (P - 1 - j) for j in range(P))


def selector(mask, tau=1):
    require(mask not in (0, ALL), "selector domain")
    scores = tuple(score(mask, a, tau) for a in range(P))
    top = max(scores)
    anchors = tuple(a for a, value in enumerate(scores) if value == top)
    require(len(anchors) == 1, "unique lexicographic anchor")
    a = anchors[0]
    w = word(mask, a, tau)
    require(w[0] == 1, "maximal rotation starts occupied")
    q = next(j for j in range(1, P) if w[j] == 0)
    occupied = (a + (q - 1) * tau) % P
    empty = (a + q * tau) % P
    return a, q, occupied, empty, w, top


def support_shift(mask, b):
    return sum(bit(mask, r) << ((r + b) % P) for r in range(P))


def support_scale(mask, u):
    return sum(bit(mask, r) << ((u * r) % P) for r in range(P))


def support_affine(mask, u, b):
    return support_shift(support_scale(mask, u), b)


def support_reflect(mask, c=0):
    return support_affine(mask, -1, c)


def cyclic_run(mask, a, tau=1, value=1):
    length = 0
    while length < P and bit(mask, a + length * tau) == value:
        length += 1
    return length


def autocorrelation(mask):
    return tuple(
        sum(bit(mask, r) * bit(mask, r + u) for r in range(P))
        for u in range(P)
    )


def psi_from_correlation(c):
    return 7 * c[0] - 12 * c[1] + 8 * c[2] - 6 * c[3] + 7 * c[4] - 6 * c[5] + 2 * c[6]


def chi7(s):
    return 1 if s in (1, 2, 4) else -1


def a_kernel(v, tau=1):
    for s in range(1, 7):
        if v % P in ((2 * tau * s) % P, (-2 * tau * s) % P):
            return -chi7(s)
    raise RuntimeError("zero A-kernel direction")


def h_kernel(h, tau=1):
    for s in range(1, 7):
        if h % P == (-2 * tau * s) % P:
            return 1
        if h % P == (2 * tau * s) % P:
            return -1
    raise RuntimeError("zero H-kernel direction")


def is_prime(n):
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


require(is_prime(MERSENNE), "8191 is prime")
require(pow(2, P, MERSENNE) == 1, "Mersenne rotation closes")
require(all(pow(2, k, MERSENNE) != 1 for k in range(1, P)), "order of two is thirteen")

# At t=-4 the two uncollapsed +/- path kernels have 72 routes each.
t_fixed = -4 % P
positive_kernel = [0] * P
negative_kernel = [0] * P
for v in range(1, P):
    for h in range(1, P):
        displacement = (t_fixed + v + h) % P
        if a_kernel(v) * h_kernel(h) > 0:
            positive_kernel[displacement] += 1
        else:
            negative_kernel[displacement] += 1
require(sum(positive_kernel) == sum(negative_kernel) == 72, "balanced path signs")

anchor_counts = Counter()
boundary_counts = Counter()
run_counts = Counter()
weight_counts = Counter()
canonical_scores = set()
reflection_canonical_scores = set()
dihedral_tie_scores = set()
reflection_symmetric_masks = 0
sharp_masks = []
sharp_scores = set()
class_counts = Counter()
class_positive = defaultdict(int)
class_negative = defaultdict(int)
boundary_path_difference = Counter()
path_positive_total = 0
path_negative_total = 0
psi_total = 0
psi_square_total = 0
weighted_numerator_checksum = 0
same_reverse_order = 0
converse_reverse_order = 0

for mask in range(1, ALL):
    n = mask.bit_count()
    weight_counts[n] += 1
    selected = selector(mask, 1)
    a, q, occupied, empty, w, top = selected

    scores = tuple(score(mask, r, 1) for r in range(P))
    require(len(set(scores)) == P, "trivial cyclic stabilizer")
    for r in range(P):
        require(
            score(mask, r + 1, 1) == 2 * score(mask, r, 1) - MERSENNE * bit(mask, r),
            "Mersenne doubling recurrence",
        )
    require(bit(mask, occupied) == 1 and bit(mask, empty) == 0, "selected one-to-zero boundary")
    require(empty == (occupied + 1) % P, "selected boundary follows guard")
    maximum_run = max(cyclic_run(mask, r, 1, 1) for r in range(P) if bit(mask, r))
    require(q == maximum_run, "lex anchor begins a longest occupied run")

    anchor_counts[a] += 1
    boundary_counts[(occupied, empty)] += 1
    run_counts[q] += 1
    canonical_scores.add(top)
    class_counts[(a, q)] += 1

    # Translation covariance and multiplicative/slope covariance generate
    # the full affine law (e,tau) -> (uS+b,u tau).
    for b in range(P):
        shifted = selector(support_shift(mask, b), 1)
        require(
            shifted[:4] == ((a + b) % P, q, (occupied + b) % P, (empty + b) % P),
            "translation-covariant selector",
        )
    for u in range(1, P):
        scaled = selector(support_scale(mask, u), u)
        require(
            scaled[:4] == ((u * a) % P, q, (u * occupied) % P, (u * empty) % P),
            "multiplicative slope covariance",
        )

    # Reflection at fixed printed slope is the opposite-slope selector of
    # the original word; simultaneous reflection and slope reversal is
    # affine covariance, not an untyped invariance assertion.
    reverse_selected = selector(mask, -1)
    reflected = selector(support_reflect(mask), 1)
    require(
        reflected[:4]
        == ((-reverse_selected[0]) % P, reverse_selected[1], (-reverse_selected[2]) % P, (-reverse_selected[3]) % P),
        "fixed-slope reflection law",
    )

    order_forward = tuple(sorted(range(P), key=lambda r: score(mask, r, 1), reverse=True))
    order_reverse = tuple(sorted(range(P), key=lambda r: score(mask, r, -1), reverse=True))
    same_reverse_order += order_reverse == order_forward
    converse_reverse_order += order_reverse == tuple(reversed(order_forward))

    # Complement sends every score L to 8191-L and hence takes the word
    # tournament to its exact global converse.
    complement = ALL ^ mask
    complement_scores = tuple(score(complement, r, 1) for r in range(P))
    require(
        all(complement_scores[r] == MERSENNE - scores[r] for r in range(P)),
        "complement score reversal",
    )
    complement_selected = selector(complement, 1)
    minimum_anchor = min(range(P), key=lambda r: scores[r])
    require(complement_selected[0] == minimum_anchor, "complement selects old minimum")

    # A nonconstant prime necklace has at most one reflection axis.
    axes = tuple(c for c in range(P) if support_reflect(mask, c) == mask)
    require(len(axes) <= 1, "at most one reflection axis")
    if axes:
        reflection_symmetric_masks += 1
        reflection_canonical_scores.add(top)
    reverse_top = reverse_selected[5]
    if top == reverse_top:
        dihedral_tie_scores.add(top)

    c = autocorrelation(mask)
    psi = psi_from_correlation(c)
    require(1 <= psi <= 98, "strict Boolean cut score")
    positive_paths = sum(positive_kernel[d] * c[d] for d in range(P))
    negative_paths = sum(negative_kernel[d] * c[d] for d in range(P))
    require(positive_paths - negative_paths == psi, "uncollapsed path imbalance equals Psi")
    require(positive_paths > negative_paths, "pointwise strict four-arm imbalance")

    path_positive_total += positive_paths
    path_negative_total += negative_paths
    psi_total += psi
    psi_square_total += psi * psi
    class_positive[(a, q)] += positive_paths
    class_negative[(a, q)] += negative_paths
    boundary_path_difference[occupied] += psi

    # Exact finite Boolean threshold expansion.  The rational owner control
    # uses denominator 17 and introduces no new sample-space coordinate.
    require(sum(1 for j in range(1, 99) if psi >= j) == psi, "Psi threshold layers")
    owner_numerator = (3 * n + a + q) % 17
    threshold_count = sum(
        1
        for ell in range(1, 18)
        for j in range(1, 99)
        if owner_numerator >= ell and psi >= j
    )
    require(threshold_count == owner_numerator * psi, "rational owner threshold layers")
    weighted_numerator_checksum += threshold_count

    if 42 * psi == n * (P - n):
        sharp_masks.append(mask)
        sharp_scores.add(top)
        require(q == 4 and n in (6, 7), "sharp cut wall has selected four-run")

require(len(canonical_scores) == 630, "630 nonconstant rotation necklaces")
require(all(value == 630 for value in anchor_counts.values()), "uniform exhaustive anchor census")
require(len(boundary_counts) == P and all(value == 630 for value in boundary_counts.values()), "uniform exhaustive boundary census")
require(reflection_symmetric_masks == 1638, "reflection-symmetric mask census")
require(len(reflection_canonical_scores) == 126, "reflection-symmetric necklace census")
require(dihedral_tie_scores == reflection_canonical_scores, "dihedral ties exactly reflection necklaces")
require(len(sharp_masks) == 52 and len(sharp_scores) == 4, "sharp equality necklace census")
require(len(class_counts) == P * (P - 1), "all anchor/run selector classes occur")

class_differences = {
    key: class_positive[key] - class_negative[key]
    for key in class_counts
}
require(all(value > 0 for value in class_differences.values()), "every selector class retains path sign")
require(path_positive_total - path_negative_total == psi_total, "global path/Psi checksum")
require(all(value == psi_total // P for value in boundary_path_difference.values()), "boundary path census")

# THM-2439's homometric pair proves that the selector cannot be recovered
# after the collision profile has been collapsed to autocorrelation.
mask_a = sum(1 << r for r in (0, 1, 3, 9))
mask_b = sum(1 << r for r in (1, 2, 5, 7))
sel_a = selector(mask_a, 1)
sel_b = selector(mask_b, 1)
require(autocorrelation(mask_a) == autocorrelation(mask_b), "homometric hostile")
require(sel_a[2:4] == (1, 2) and sel_b[2:4] == (2, 3), "homometric boundaries differ")
dihedral_a = {
    support_affine(mask_a, u, b)
    for u in (-1, 1)
    for b in range(P)
}
require(mask_b not in dihedral_a, "homometric pair is non-dihedral")

# On the actual deep comb there is one cyclic occupied run, of size one or
# two.  The ambient lex marker and terminal wall therefore collapse to
# quadratic local predicates.  In the target-anchored chart root zero is
# forbidden, leaving the twelve-vertex path and exactly two endpoint cases
# in each natural deep direction.
deep_masks = tuple(1 << r for r in range(P)) + tuple(
    (1 << r) | (1 << ((r + 1) % P))
    for r in range(P)
)
for mask in deep_masks:
    selected = selector(mask, 1)
    a, _, occupied, _, _, _ = selected
    marker_values = tuple(bit(mask, r) * (1 - bit(mask, r - 1)) for r in range(P))
    terminal_values = tuple(bit(mask, r) * (1 - bit(mask, r + 1)) for r in range(P))
    require(sum(marker_values) == 1 and marker_values[a] == 1, "quadratic deep marker")
    require(sum(terminal_values) == 1 and terminal_values[occupied] == 1, "quadratic deep wall")
    for tau in range(1, P):
        a_tau = selector(mask, tau)[0]
        inverse_step = pow(tau, -1, P)
        neighbour = -1 if inverse_step <= 6 else 1
        marker_tau = tuple(
            bit(mask, r) * (1 - bit(mask, r + neighbour))
            for r in range(P)
        )
        require(sum(marker_tau) == 1 and marker_tau[a_tau] == 1, "all-slope quadratic deep marker")

target_masks = tuple(1 << r for r in range(1, P)) + tuple(
    (1 << r) | (1 << (r + 1))
    for r in range(1, P - 1)
)
plus_target_masks = []
minus_target_masks = []
for mask in target_masks:
    plus_boundary = selector(mask, 1)[2:4]
    minus_boundary = selector(mask, -1)[2:4]
    if plus_boundary == (12, 0):
        plus_target_masks.append(mask)
    if minus_boundary == (1, 0):
        minus_target_masks.append(mask)
    require((plus_boundary == (12, 0)) == bool(bit(mask, 12)), "plus endpoint diagonal")
    require((minus_boundary == (1, 0)) == bool(bit(mask, 1)), "minus endpoint diagonal")
require(
    set(plus_target_masks) == {1 << 12, (1 << 11) | (1 << 12)},
    "plus target-ending masks",
)
require(
    set(minus_target_masks) == {1 << 1, (1 << 1) | (1 << 2)},
    "minus target-ending masks",
)
require(selector(1 << 1, -1)[2:4] == (1, 0), "THM-2529 hostile target wall")

# A deterministic positive weighting checks the exact linear recovery of
# start-marker and terminal-wall masses from the anchored Gram matrix.
weights = {mask: index + 1 for index, mask in enumerate(target_masks)}
gram = [
    [sum(weight * bit(mask, u) * bit(mask, v) for mask, weight in weights.items()) for v in range(P)]
    for u in range(P)
]
gamma = [
    sum(weight * bit(mask, a) * (1 - bit(mask, a - 1)) for mask, weight in weights.items())
    for a in range(P)
]
terminal_mass = [
    sum(weight * bit(mask, a) * (1 - bit(mask, a + 1)) for mask, weight in weights.items())
    for a in range(P)
]
require(all(gamma[a] == gram[a][a] - gram[(a - 1) % P][a] for a in range(P)), "Gram start-marker recovery")
require(all(terminal_mass[a] == gram[a][a] - gram[a][(a + 1) % P] for a in range(P)), "Gram terminal-wall recovery")
require(gamma[0] == 0 and sum(gamma) == sum(weights.values()), "positive target-arc partition")
require(any(gamma[a] != gamma[0] for a in range(1, P)), "cyclotomic selector vector is nonconstant")

# There is no complement-dihedral fixed mask at odd length: any such
# identity would force |S|=13-|S|.  We also keep the exact finite count.
complement_dihedral_ties = sum(
    1
    for mask in range(1, ALL)
    if any((ALL ^ mask) == support_affine(mask, u, b) for u in (-1, 1) for b in range(P))
)
require(complement_dihedral_ties == 0, "no complement-dihedral tie")

expected_run_necklaces = {
    1: 40,
    2: 172,
    3: 178,
    4: 114,
    5: 62,
    6: 32,
    7: 16,
    8: 8,
    9: 4,
    10: 2,
    11: 1,
    12: 1,
}
require({q: run_counts[q] // P for q in sorted(run_counts)} == expected_run_necklaces, "run-length necklace census")
require(same_reverse_order == 0 and converse_reverse_order == 0, "slope reversal is neither identity nor converse")

sharp_words = sorted(format(value, "013b") for value in sharp_scores)
run_summary = ",".join(f"{q}:{expected_run_necklaces[q]}" for q in sorted(expected_run_necklaces))
print("THM-2531 exact referee")
print(f"p={P} mersenne_modulus={MERSENNE} prime={is_prime(MERSENNE)} order_of_2={P}")
print(f"mixed_masks={ALL-1} rotation_necklaces={len(canonical_scores)} cyclic_stabilizer_failures=0")
print(f"anchor_counts={','.join(str(anchor_counts[a]) for a in range(P))}")
print(f"boundary_counts={','.join(str(boundary_counts[(a,(a+1)%P)]) for a in range(P))}")
print(f"max_one_run_necklaces={run_summary}")
print(f"reflection_symmetric_masks={reflection_symmetric_masks} reflection_necklaces={len(reflection_canonical_scores)}")
print(f"dihedral_orbits=378 orientation_ties={len(dihedral_tie_scores)} complement_dihedral_ties={complement_dihedral_ties}")
print(f"fixed_mask_slope_reversal_same={same_reverse_order} fixed_mask_slope_reversal_converse={converse_reverse_order}")
print(f"sharp_masks={len(sharp_masks)} sharp_necklaces={len(sharp_scores)} selected_run=4")
print("sharp_canonical_words=" + ",".join(sharp_words))
print(f"psi_range=1..98 psi_sum={psi_total} psi_square_sum={psi_square_total}")
print(f"boolean_layer_sum={psi_total} rational_owner_denominator=17 weighted_numerator_checksum={weighted_numerator_checksum}")
print(f"path_positive_total={path_positive_total} path_negative_total={path_negative_total} difference={path_positive_total-path_negative_total}")
print(f"selector_classes={len(class_counts)} class_imbalance_min={min(class_differences.values())} class_imbalance_max={max(class_differences.values())}")
print(f"boundary_imbalance_each={psi_total//P}")
print(f"homometric_correlation={autocorrelation(mask_a)} boundaries={sel_a[2:4]},{sel_b[2:4]}")
print(f"deep_comb_masks={len(deep_masks)} quadratic_marker_wall=PASS")
print(f"deep_slope_marker_tests={len(deep_masks)*(P-1)}")
print("epsilon_convention=tau_inverse_1..6:+1,tau_inverse_7..12:-1")
print("target_anchored_masks=23 plus_target={12},{11,12} minus_target={1},{1,2}")
print("thm2529_hostile_minus_boundary=(1,0)")
print(f"gram_marker_mass={sum(gamma)} gamma0={gamma[0]} gram_linear_recovery=PASS all_root_colours=NONZERO_BY_PHI13")
print(f"translation_covariance_tests={(ALL-1)*P} slope_covariance_tests={(ALL-1)*(P-1)}")
print(f"checks={checks}")
print("status=PASS")
