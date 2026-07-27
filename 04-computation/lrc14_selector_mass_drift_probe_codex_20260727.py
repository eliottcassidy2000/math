#!/usr/bin/env python3
"""Exact probe for selector-star mass versus THM-2365 target drift.

This is a research companion, not a theorem-ID reservation.  It records the
smallest path-cone controls separating three cell observables:

* the cell mass ``mu``;
* the THM-2531 positive-direction start-marker vector ``gamma``; and
* the THM-2365 relative diagonal ``h_j=H(t+j,s,t)``.

It also resolves the deep-phase endpoint geometry on the denominator-182
mesh and inserts that geometry into THM-2367's lawful full-owner circulant
hostile.  Every calculation is exact integer or ``Fraction`` arithmetic.
"""

from fractions import Fraction


P = 13
MESH = 2 * P * 7  # all endpoints a/13 +/- 1/14 lie on this mesh


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def circular_distance(x):
    x %= 1
    return min(x, 1 - x)


def danger(y, root):
    return circular_distance(y - Fraction(root, P)) < Fraction(1, 14)


def path_coordinates(alpha=None, beta=None):
    """Return (mu,gamma,h) on the anchored path 1--...--12.

    ``alpha[j]`` is the singleton-ray mass at j and ``beta[j]`` is the
    adjacent-pair-ray mass at {j,j+1}.  Missing entries mean zero.
    """
    alpha = dict(alpha or {})
    beta = dict(beta or {})
    require(all(1 <= j <= 12 and value >= 0 for j, value in alpha.items()), "alpha cone")
    require(all(1 <= j <= 11 and value >= 0 for j, value in beta.items()), "beta cone")
    gamma = [Fraction(0) for _ in range(P)]
    h = [Fraction(0) for _ in range(P)]
    for j in range(1, P):
        # A start at j is either singleton {j} or the pair {j,j+1}.
        gamma[j] = alpha.get(j, 0) + beta.get(j, 0)
        # Root j is active in a start at j or in the preceding pair.
        h[j] = gamma[j] + beta.get(j - 1, 0)
        require(
            h[j]
            == alpha.get(j, 0) + beta.get(j - 1, 0) + beta.get(j, 0),
            "selector/run-length reconstruction",
        )
    mu = sum(alpha.values(), Fraction(0)) + sum(beta.values(), Fraction(0))
    require(mu == sum(gamma), "start-marker partition of mass")
    return mu, tuple(gamma), tuple(h)


def vector_energy(vectors):
    """Mean squared distance from the cell average, with no root normalization."""
    count = len(vectors)
    width = len(vectors[0])
    means = [
        sum((vector[j] for vector in vectors), Fraction(0)) / count
        for j in range(width)
    ]
    return sum(
        (vector[j] - means[j]) ** 2
        for vector in vectors
        for j in range(width)
    ) / count


def scalar_energy(values):
    mean = sum(values, Fraction(0)) / len(values)
    return sum((value - mean) ** 2 for value in values) / len(values)


def compare_cells(cells):
    coordinates = [path_coordinates(*cell) for cell in cells]
    masses = [item[0] for item in coordinates]
    gammas = [item[1] for item in coordinates]
    diagonals = [item[2] for item in coordinates]
    return scalar_energy(masses), vector_energy(gammas), vector_energy(diagonals)


# Control A: the selector start and cell mass are identical, but the run
# changes from singleton {2} to adjacent pair {2,3}.  THM-2365's diagonal
# sees root 3 appear, while (mu,gamma) see no change.
selector_blind = compare_cells(
    [
        ({2: Fraction(1)}, {}),
        ({}, {2: Fraction(1)}),
    ]
)
require(selector_blind == (0, 0, Fraction(1, 4)), "selector-blind H drift")

# Control B: one pair versus two singleton rays has the same diagonal but
# different total mass and start-marker vector.  Thus selector or mass drift
# does not imply THM-2365 H-drift, even inside the nonnegative path cone.
selector_false_positive = compare_cells(
    [
        ({}, {2: Fraction(1)}),
        ({2: Fraction(1), 3: Fraction(1)}, {}),
    ]
)
require(
    selector_false_positive == (Fraction(1, 4), Fraction(1, 4), 0),
    "selector false-positive control",
)

# Control C keeps even the total mass fixed: exchange pair {2,3} for its two
# singleton rays while making the opposite exchange at the disjoint edge
# {6,7}.  H and mu are unchanged, but gamma changes.
mass_neutral_false_positive = compare_cells(
    [
        ({6: Fraction(1), 7: Fraction(1)}, {2: Fraction(1)}),
        ({2: Fraction(1), 3: Fraction(1)}, {6: Fraction(1)}),
    ]
)
require(
    mass_neutral_false_positive == (0, Fraction(1, 2), 0),
    "mass-neutral selector false positive",
)


# Resolve one turn of the actual deep comb.  Sampling the 182 elementary
# open cells at their midpoints is exact because every deep endpoint lies on
# the mesh boundary.  The start-marker tooth is Delta_a(1-Delta_{a-1}).
singleton_counts = [0] * P
pair_counts = [0] * P
marker_counts = [0] * P
safe_marker_counts = [0] * P
safe_diagonal_counts = [0] * P
safe_cells = 0

for cell in range(MESH):
    y = Fraction(2 * cell + 1, 2 * MESH)
    bits = [int(danger(y, root)) for root in range(P)]
    support = [root for root, value in enumerate(bits) if value]
    require(len(support) in (1, 2), "one-or-two deep roots")
    if len(support) == 2:
        require(
            (support[1] - support[0]) % P == 1 or support == [0, P - 1],
            "cyclic adjacency",
        )
    starts = [root for root in range(P) if bits[root] and not bits[(root - 1) % P]]
    require(len(starts) == 1, "unique positive-direction start")
    start = starts[0]
    marker_counts[start] += 1
    if len(support) == 1:
        singleton_counts[support[0]] += 1
    else:
        pair_start = start
        pair_counts[pair_start] += 1

    # Target root zero is the excluded deep root in THM-2530/2531.
    if not bits[0]:
        safe_cells += 1
        safe_marker_counts[start] += 1
        for root in support:
            safe_diagonal_counts[root] += 1

require(marker_counts == [14] * P, "each marker tooth has length 1/13")
require(singleton_counts == [2] * P, "each singleton phase has length 1/91")
require(pair_counts == [12] * P, "each adjacent-pair phase has length 6/91")
require(safe_cells == 156, "deep-safe mass is 6/7")
require(
    safe_marker_counts == [0] + [14] * 11 + [2],
    "anchored marker histogram",
)
require(
    safe_diagonal_counts == [0, 14] + [26] * 10 + [14],
    "anchored diagonal histogram",
)

# THM-2367 peels every lower factor from its full nine-factor lawful hostile.
# Its residual deep phase is uniform and is multiplied by
# A=(13/30)(6/7)^3(1/7).  The exact selector histogram below is consequently
# independent of both lawful target cells (s,t), just like its H diagonal.
peeling_constant = Fraction(468, 12005)
hostile_gamma = tuple(peeling_constant * Fraction(count, MESH) for count in safe_marker_counts)
hostile_h = tuple(peeling_constant * Fraction(count, MESH) for count in safe_diagonal_counts)
hostile_mu = peeling_constant * Fraction(safe_cells, MESH)

require(hostile_mu == peeling_constant * Fraction(6, 7), "hostile owner mass")
require(hostile_gamma[0] == 0, "hostile excluded marker")
require(
    hostile_gamma[1:12] == (peeling_constant / 13,) * 11
    and hostile_gamma[12] == peeling_constant / 91,
    "hostile selector normal form",
)
require(
    hostile_h[1] == hostile_h[12] == peeling_constant / 13
    and hostile_h[2:12] == (peeling_constant / 7,) * 10,
    "THM-2367 circulant diagonal",
)


print("LRC14 selector-mass / THM-2365 drift exact probe")
print(f"mesh={MESH} marker_tooth_cells={marker_counts[0]} marker_tooth_mass=1/13")
print(
    "deep_phase_counts="
    f"singleton_each:{singleton_counts[0]},pair_each:{pair_counts[0]},safe:{safe_cells}"
)
print("selector_blind=(mass_energy,gamma_energy,H_energy)=" + ",".join(map(str, selector_blind)))
print(
    "selector_false_positive=(mass_energy,gamma_energy,H_energy)="
    + ",".join(map(str, selector_false_positive))
)
print(
    "mass_neutral_false_positive=(mass_energy,gamma_energy,H_energy)="
    + ",".join(map(str, mass_neutral_false_positive))
)
print(f"thm2367_peeling_constant={peeling_constant}")
print(f"thm2367_selector_mu={hostile_mu}")
print("thm2367_selector_gamma=" + ",".join(map(str, hostile_gamma)))
print("thm2367_relative_H=" + ",".join(map(str, hostile_h)))
print("exact_reconstruction=H_j=gamma_j+beta_(j-1)")
print("minimal_missing_sidecar=run_length_two_mass_beta")
print("profile_only_uniform_test=IMPOSSIBLE-WITHOUT-SCALAR-COVER-SIDECAR")
print("status=PASS")
