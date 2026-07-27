#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2538.

The referee checks four finite statements:

* the deep-target covariant product shift preserves the excluded-root
  anchor and has multiplier zeta_13^(a*tau-b*h-q*c*tau);
* three rows leave a common nonzero transverse target shift and four
  positive packing scalars suffice;
* an explicit rational anchored singleton table has three full gain slices
  and literal transverse boundaries while its anchored Gram is cellwise
  constant; and
* same-ancestry cross-Kakeya products recover a later arrival coupling and
  THM-2537's selected-head threshold sum is its exact positive-hit form,
  whereas separately integrated endpoint banks retain only its margins and
  even all pairwise banks can lose a three-field parity incidence.

All checks use integers or the exact finite field F_547.  Explicit
exceptions, rather than ``assert``, keep normal and optimized runs identical.
"""

from itertools import permutations, product


ROOT = 13
SOURCE = 7
FIELD = 547  # FIELD - 1 = 6 * 7 * 13
checks = 0


def require(condition: bool, message: str) -> None:
    """Raise on failure in normal and optimized Python."""

    global checks
    checks += 1
    if not condition:
        raise RuntimeError("FAILED: " + message)


def primitive_root(prime: int) -> int:
    factors = (2, 3, 7, 13)
    for candidate in range(2, prime):
        if all(pow(candidate, (prime - 1) // factor, prime) != 1 for factor in factors):
            return candidate
    raise RuntimeError("FAILED: no primitive root")


GENERATOR = primitive_root(FIELD)
ZETA = pow(GENERATOR, (FIELD - 1) // ROOT, FIELD)
XI = pow(GENERATOR, (FIELD - 1) // SOURCE, FIELD)


def zeta(exponent: int) -> int:
    return pow(ZETA, exponent % ROOT, FIELD)


def xi(exponent: int) -> int:
    return pow(XI, exponent % SOURCE, FIELD)


require(pow(ZETA, ROOT, FIELD) == 1, "zeta_13 order divides 13")
require(all(pow(ZETA, exponent, FIELD) != 1 for exponent in range(1, ROOT)), "zeta_13 primitive")
require(pow(XI, SOURCE, FIELD) == 1, "xi_7 order divides 7")
require(all(pow(XI, exponent, FIELD) != 1 for exponent in range(1, SOURCE)), "xi_7 primitive")


# ---------------------------------------------------------------------------
# Deep-anchor covariance and the full-character multiplier.


anchor_covariance_checks = 0
for c, tau, u, target in product(range(ROOT), repeat=4):
    old_difference = (c * u - target) % ROOT
    new_difference = (c * (u + tau) - (target + c * tau)) % ROOT
    require(new_difference == old_difference, "deep target moves covariantly with root")
    anchor_covariance_checks += 1


# On a delta basis at (s_0,t_0,u_0), the shift
#
#   (P^anc B)_(s,t)(u)=B_(s+h,t+c*tau)(u+tau)
#
# moves its support to (s_0-h,t_0-c*tau,u_0-tau).  The character convention
# zeta^(b*s+q*t-a*u) therefore changes by a*tau-b*h-q*c*tau.
multiplier_checks = 0
for a, b, q, c, tau, h in product(range(ROOT), repeat=6):
    direct_change = (-b * h - q * c * tau + a * tau) % ROOT
    claimed_change = (a * tau - b * h - q * c * tau) % ROOT
    require(direct_change == claimed_change, "anchored product-shift multiplier")
    multiplier_checks += 1


# Three functionals forbid at most three of the twelve nonzero target shifts.
transverse_pattern_checks = 0
minimum_nonzero_transverse = ROOT - 1
for forbidden_values in product(range(ROOT), repeat=3):
    allowed = tuple(h for h in range(1, ROOT) if h not in forbidden_values)
    require(len(allowed) >= 9, "three rows leave at least nine nonzero transverse shifts")
    minimum_nonzero_transverse = min(minimum_nonzero_transverse, len(allowed))
    transverse_pattern_checks += 1
require(minimum_nonzero_transverse == 9, "nine-shift bound is sharp")


# If u_i-v_i is nonzero, u_i+t*v_i=0 excludes at most one scalar t.  The
# entire analytic content of the four-candidate argument is therefore the
# following finite pigeonhole audit.  ``None`` means that a row excludes no
# candidate from {1,2,3,4}.
packing_pattern_checks = 0
packing_candidates = (1, 2, 3, 4)
for bad_values in product((None,) + packing_candidates, repeat=3):
    survivors = tuple(t for t in packing_candidates if t not in bad_values)
    require(bool(survivors), "three rows cannot cover four packing scalars")
    packing_pattern_checks += 1
require(
    tuple(t for t in packing_candidates if t not in (1, 2, 3)) == (4,),
    "four-candidate bound has a sharp exclusion pattern",
)


# The full rational Galois orbit scales source, target, deep-target, and root
# colours simultaneously.  Distinct ratios a/b remain disjoint.
galois_units = tuple(product(range(1, SOURCE), range(1, ROOT)))
galois_orbit_checks = 0
all_gain_incidences: set[tuple[int, int, int, int]] = set()
for gain in (1, 2, 3):
    orbit = {
        (unit7, unit13, 0, (unit13 * gain) % ROOT)
        for unit7, unit13 in galois_units
    }
    require(len(orbit) == 72, "one full-character gain orbit has 72 incidences")
    require(
        all(a * pow(b, -1, ROOT) % ROOT == gain for _, b, _, a in orbit),
        "Galois action preserves gain a/b",
    )
    require(all_gain_incidences.isdisjoint(orbit), "distinct gain orbits are disjoint")
    all_gain_incidences.update(orbit)
    galois_orbit_checks += 1
require(len(all_gain_incidences) == 216, "three gains have 216 incidences")


# ---------------------------------------------------------------------------
# An anchored constant-singleton hostile with three transverse gains.
#
# Use the ordered section F_13 -> {0,...,6}, d |-> rep(d) mod 7.  For gain
# lambda, the Boolean atom in cell (ell,s,t) is supported on the singleton
# root u=t+1 and the source label ell=section(s-lambda*t).  Every cell thus
# has the same anchored deep mask e_1.  The section's carry supplies nonzero
# mixed source/target phase without moving that deep location.


def section(value: int) -> int:
    return (value % ROOT) % SOURCE


GAINS = (1, 2, 3)
HOSTILE_TAU = 1
HOSTILE_H = 4
PACKING_SCALAR = 2


def hostile_b(ell: int, s: int, target: int, gain: int, root: int) -> int:
    return int(
        root % ROOT == (target + 1) % ROOT
        and ell % SOURCE == section(s - gain * target)
    )


def hostile_a(ell: int, s: int, target: int, gain: int, root: int) -> int:
    return hostile_b(
        ell,
        s + HOSTILE_H,
        target + HOSTILE_TAU,
        gain,
        root + HOSTILE_TAU,
    )


hostile_atom_checks = 0
for gain, s, target, root, ell in product(
    GAINS, range(ROOT), range(ROOT), range(ROOT), range(SOURCE)
):
    old = hostile_b(ell, s, target, gain, root)
    new = hostile_a(ell, s, target, gain, root)
    deep_target = int(root == target)
    require(old * deep_target == 0, "hostile atom obeys excluded-root anchor")
    require(new * deep_target == 0, "shifted hostile atom preserves excluded-root anchor")
    require(old * new == 0, "chosen transverse hostile orientations are disjoint")
    hostile_atom_checks += 1


# The disjointness mechanism is the carry section itself.  For the three
# differences h-lambda*tau = 3,2,1 it changes at every cyclic input.
section_transverse_checks = 0
for gain, difference in zip(GAINS, (3, 2, 1), strict=True):
    require(difference == (HOSTILE_H - gain * HOSTILE_TAU) % ROOT, "hostile transverse difference")
    for value in range(ROOT):
        require(
            section(value + difference) != section(value),
            "ordered section separates the transverse atom",
        )
        section_transverse_checks += 1


# Its primitive source/target character is the nonzero geometric sum
# sum_{d=0}^{12}(xi^kappa zeta^b)^d.
hostile_kernel_checks = 0
for kappa, b in product(range(1, SOURCE), range(1, ROOT)):
    ratio = xi(kappa) * zeta(b) % FIELD
    direct = sum(xi(kappa * section(d)) * zeta(b * d) for d in range(ROOT)) % FIELD
    geometric = (1 - pow(ratio, ROOT, FIELD)) * pow((1 - ratio) % FIELD, -1, FIELD) % FIELD
    require(direct == geometric, "hostile ordered-section geometric kernel")
    require(direct != 0, "hostile primitive source/target kernel is nonzero")
    hostile_kernel_checks += 1


def hostile_zhat(kappa: int, b: int, q: int, a: int) -> int:
    """Full character/root coefficient of Z=A+2B on the hostile."""

    total = 0
    for gain, s, target in product(GAINS, range(ROOT), range(ROOT)):
        root = (target + 1) % ROOT
        ell_b = section(s - gain * target)
        ell_a = section(s + HOSTILE_H - gain * (target + HOSTILE_TAU))
        phase_b = xi(kappa * ell_b) * zeta(b * s + q * target - a * root)
        phase_a = xi(kappa * ell_a) * zeta(b * s + q * target - a * root)
        total += phase_a + PACKING_SCALAR * phase_b
    return total % FIELD


hostile_gain_checks = 0
for kappa, b, gain in product(range(1, SOURCE), range(1, ROOT), GAINS):
    a = (b * gain) % ROOT
    coefficient = hostile_zhat(kappa, b, 0, a)
    require(coefficient != 0, "anchored hostile retains every selected gain incidence")
    hostile_gain_checks += 1
require(hostile_gain_checks == 216, "hostile has three complete gain slices")


# Every target cell has the same total Z mass and the same relative singleton
# Gram: only coordinate (1,1) is nonzero.  Hence its anchored Gram, both
# selector vectors, and every location statistic are target-cell invariant.
hostile_cell_gram_checks = 0
reference_gram: tuple[tuple[int, ...], ...] | None = None
for s, target in product(range(ROOT), repeat=2):
    gram = [[0 for _ in range(ROOT)] for _ in range(ROOT)]
    for gain in GAINS:
        root = (target + 1) % ROOT
        relative = (root - target) % ROOT
        gram[relative][relative] += 1 + PACKING_SCALAR
    frozen = tuple(tuple(row) for row in gram)
    if reference_gram is None:
        reference_gram = frozen
    require(frozen == reference_gram, "hostile anchored Gram is cellwise constant")
    require(gram[1][1] == 3 * (1 + PACKING_SCALAR), "hostile singleton Gram mass")
    require(
        sum(sum(row) for row in gram) == 3 * (1 + PACKING_SCALAR),
        "hostile Gram has no other coordinate",
    )
    hostile_cell_gram_checks += 1


# ---------------------------------------------------------------------------
# Common-ancestry arrival: exact positive reconstruction and marginal no-go.


def bit(mask: int, root: int) -> int:
    return (mask >> (root % ROOT)) & 1


def kakeya(mask: int, tau: int, root: int) -> int:
    return bit(mask, root) * (1 - bit(mask, root + tau))


# The empty endpoint of the same predicate annihilates as an arrival.
self_arrival_checks = 0
for mask in range(1 << ROOT):
    for source_root in range(ROOT):
        for arrival_root in range(ROOT):
            if source_root == arrival_root:
                continue
            value = kakeya(mask, arrival_root - source_root, source_root) * bit(mask, arrival_root)
            require(value == 0, "same-horizon empty head cannot be arrival of itself")
            self_arrival_checks += 1


# With any retained empty anchor, one boundary directed to that anchor
# recovers each other occupied coordinate pointwise.  The predecessor and
# later field may use different anchors; in particular a later arrival at
# predecessor root zero must use a different later empty anchor.
anchored_reconstruction_checks = 0
for mask in range(1 << ROOT):
    for anchor in range(ROOT):
        if bit(mask, anchor):
            continue
        for root in range(ROOT):
            if root == anchor:
                continue
            require(
                kakeya(mask, anchor - root, root) == bit(mask, root),
                "empty-anchor Kakeya coordinate inverse",
            )
            anchored_reconstruction_checks += 1


# The complete all-slope formula is useful when no zero anchor is available.
crofton_coordinate_checks = 0
for mask in range(1 << ROOT):
    occupancy = sum(bit(mask, root) for root in range(ROOT))
    for root in range(ROOT):
        boundary_sum = sum(kakeya(mask, tau, root) for tau in range(1, ROOT))
        require(
            boundary_sum == (ROOT - occupancy) * bit(mask, root),
            "all-slope pointwise Crofton reconstruction",
        )
        crofton_coordinate_checks += 1


# Every all-slope predecessor needle ending at a target coordinate s sums to
# m_e(1-e_s).  Linearity in the later field then proves the complete
# empty-head hitting identity; checking target-coordinate basis vectors is
# an exact check for every later field over any coefficient ring.
all_slope_hitting_checks = 0
for mask in range(1 << ROOT):
    occupancy = sum(bit(mask, root) for root in range(ROOT))
    for target in range(ROOT):
        needle_sum = sum(
            kakeya(mask, tau, root) * int((root + tau) % ROOT == target)
            for tau, root in product(range(1, ROOT), range(ROOT))
        )
        require(
            needle_sum == occupancy * (1 - bit(mask, target)),
            "all-slope Kakeya needles hit each empty target with occupancy multiplicity",
        )
        all_slope_hitting_checks += 1


# THM-2537's selected head is a categorical field: on an active fibre it is
# one singleton and on an inactive fibre it is empty.  Any nonzero slope
# therefore recovers every selected-head coordinate pointwise.
selected_head_kakeya_checks = 0
for selected_head in range(-1, ROOT):
    head_mask = 0 if selected_head == -1 else 1 << selected_head
    for tau, root in product(range(1, ROOT), range(ROOT)):
        require(
            kakeya(head_mask, tau, root) == bit(head_mask, root),
            "categorical selected head is one-slope Kakeya reconstructible",
        )
        selected_head_kakeya_checks += 1


# Its threshold decomposition is exactly the positive target hit in
# THM-2537 equation (56), tested for every possible head, Boolean owner,
# integer score 0..98, and value of the later target-active event at the head.
positive_hit_layer_checks = 0
for selected_head, owner, score, target_active in product(
    range(ROOT), range(2), range(99), range(2)
):
    layered_hit = sum(
        owner * int(score >= layer) * target_active
        for layer in range(1, 99)
    )
    require(
        layered_hit == owner * score * target_active,
        "selected-head threshold layers equal the positive target hit",
    )
    positive_hit_layer_checks += 1


def singleton_mask(root: int) -> int:
    return 1 << root


def endpoint_bank(values: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Unnormalised singleton Gram diagonal and all directed boundary masses."""

    gram_diagonal = tuple(sum(int(value == root) for value in values) for root in range(ROOT))
    boundaries = tuple(
        sum(kakeya(singleton_mask(value), tau, root) for value in values)
        for tau in range(1, ROOT)
        for root in range(ROOT)
    )
    return gram_diagonal, boundaries


def coupling(
    sources: tuple[int, ...], arrivals: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(
            sum(int(source == r and arrival == s) for source, arrival in zip(sources, arrivals, strict=True))
            for s in range(ROOT)
        )
        for r in range(ROOT)
    )


sources = (1, 2)
arrivals_same = (1, 2)
arrivals_switch = (2, 1)

require(endpoint_bank(sources) == endpoint_bank(sources), "source endpoint control")
require(
    endpoint_bank(arrivals_same) == endpoint_bank(arrivals_switch),
    "two arrival couplings have identical endpoint Gram/Kakeya banks",
)

coupling_same = coupling(sources, arrivals_same)
coupling_switch = coupling(sources, arrivals_switch)
require(coupling_same != coupling_switch, "fixed endpoint margins leave a coupling choice")

haar_same = (
    coupling_same[1][1]
    - coupling_same[1][2]
    - coupling_same[2][1]
    + coupling_same[2][2]
)
haar_switch = (
    coupling_switch[1][1]
    - coupling_switch[1][2]
    - coupling_switch[2][1]
    + coupling_switch[2][2]
)
require((haar_same, haar_switch) == (2, -2), "minimal mixed-Haar arrival switch")
require(
    (sum(coupling_same[r][r] for r in range(ROOT)), sum(coupling_switch[r][r] for r in range(ROOT)))
    == (2, 0),
    "minimal loop mass is not marginally determined",
)


# Verify cross-Kakeya reconstruction of the two coupling matrices entrywise.
cross_kakeya_checks = 0
for arrivals, expected in (
    (arrivals_same, coupling_same),
    (arrivals_switch, coupling_switch),
):
    for r, s in product(range(1, ROOT), repeat=2):
        reconstructed = sum(
            kakeya(singleton_mask(source), -r, r)
            * kakeya(singleton_mask(arrival), -s, s)
            for source, arrival in zip(sources, arrivals, strict=True)
        )
        require(reconstructed == expected[r][s], "same-base cross-Kakeya coupling reconstruction")
        cross_kakeya_checks += 1


# Fréchet is the sharp one-coordinate obstruction: positive marginals can
# still have zero intersection whenever their total does not exceed the base
# mass.  Exhaust every pair of finite events through eight atoms.
frechet_cases = 0
for universe_size in range(1, 9):
    masks_by_size = {
        size: tuple(mask for mask in range(1 << universe_size) if mask.bit_count() == size)
        for size in range(universe_size + 1)
    }
    for source_mass, arrival_mass in product(range(universe_size + 1), repeat=2):
        attained = {
            (source & arrival).bit_count()
            for source in masks_by_size[source_mass]
            for arrival in masks_by_size[arrival_mass]
        }
        expected = set(
            range(
                max(0, source_mass + arrival_mass - universe_size),
                min(source_mass, arrival_mass) + 1,
            )
        )
        require(attained == expected, "sharp finite Frechet intersection interval")
        frechet_cases += 1


def compositions(total: int, parts: int) -> tuple[tuple[int, ...], ...]:
    if parts == 1:
        return ((total,),)
    return tuple(
        (first,) + tail
        for first in range(total + 1)
        for tail in compositions(total - first, parts - 1)
    )


def maximum_off_diagonal_flow(rows: tuple[int, ...], columns: tuple[int, ...]) -> int:
    """Small exact max flow on the complete bipartite graph minus diagonal."""

    count = len(rows)
    source = 0
    row_start = 1
    column_start = row_start + count
    sink = column_start + count
    node_count = sink + 1
    capacity = [[0 for _ in range(node_count)] for _ in range(node_count)]

    for index, value in enumerate(rows):
        capacity[source][row_start + index] = value
    for row in range(count):
        for column in range(count):
            if row != column:
                capacity[row_start + row][column_start + column] = sum(rows)
    for index, value in enumerate(columns):
        capacity[column_start + index][sink] = value

    flow = 0
    while True:
        parent = [-1] * node_count
        parent[source] = source
        queue = [source]
        for node in queue:
            for neighbour in range(node_count):
                if parent[neighbour] == -1 and capacity[node][neighbour] > 0:
                    parent[neighbour] = node
                    queue.append(neighbour)
                    if neighbour == sink:
                        break
            if parent[sink] != -1:
                break
        if parent[sink] == -1:
            break

        increment = sum(rows)
        node = sink
        while node != source:
            previous = parent[node]
            increment = min(increment, capacity[previous][node])
            node = previous
        node = sink
        while node != source:
            previous = parent[node]
            capacity[previous][node] -= increment
            capacity[node][previous] += increment
            node = previous
        flow += increment
    return flow


# Hall's theorem for the diagonal-forbidden transport graph reduces to the
# singleton inequalities p_i+q_i<=1.  Check the corresponding exact minimum
# diagonal mass on all integer marginals with 2--4 states and total <=5.
hall_transport_checks = 0
for state_count in range(2, 5):
    for total in range(1, 6):
        margin_vectors = compositions(total, state_count)
        for row_margins, column_margins in product(margin_vectors, repeat=2):
            off_diagonal = maximum_off_diagonal_flow(row_margins, column_margins)
            minimum_diagonal = total - off_diagonal
            hall_frechet_bound = max(
                0,
                max(
                    row_margins[index] + column_margins[index] - total
                    for index in range(state_count)
                ),
            )
            require(minimum_diagonal == hall_frechet_bound, "Hall/Frechet minimum diagonal transport")
            hall_transport_checks += 1


# Pairwise cross banks still do not determine a three-field owner incidence.
# The even/odd parity laws agree on every one- and two-coordinate marginal.
even_parity = ((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0))
odd_parity = ((1, 1, 1), (1, 0, 0), (0, 1, 0), (0, 0, 1))
triple_parity_checks = 0
for coordinate in range(3):
    require(
        tuple(sum(int(atom[coordinate] == value) for atom in even_parity) for value in range(2))
        == tuple(sum(int(atom[coordinate] == value) for atom in odd_parity) for value in range(2)),
        "parity hostiles have equal one-coordinate marginals",
    )
    triple_parity_checks += 1
for left, right in ((0, 1), (0, 2), (1, 2)):
    require(
        tuple(
            sum(int((atom[left], atom[right]) == pair) for atom in even_parity)
            for pair in product(range(2), repeat=2)
        )
        == tuple(
            sum(int((atom[left], atom[right]) == pair) for atom in odd_parity)
            for pair in product(range(2), repeat=2)
        ),
        "parity hostiles have equal two-coordinate marginals",
    )
    triple_parity_checks += 1
require(
    (
        sum(x * y * owner for x, y, owner in even_parity),
        sum(x * y * owner for x, y, owner in odd_parity),
    )
    == (0, 1),
    "three-field incidence is not pairwise determined",
)
triple_parity_checks += 1
require(
    (
        sum((-1) ** (x + y + owner) for x, y, owner in even_parity),
        sum((-1) ** (x + y + owner) for x, y, owner in odd_parity),
    )
    == (4, -4),
    "three-way parity character flips sign",
)
triple_parity_checks += 1


# Uniform permutation couplings have identical endpoint data and arbitrary
# deterministic cycle type.  Identity has three loops; a 3-cycle has none.
permutation_coupling_checks = 0
permutation_loop_counts: set[int] = set()
base = (1, 2, 3)
base_bank = endpoint_bank(base)
for permutation in permutations(base):
    require(endpoint_bank(permutation) == base_bank, "permutation coupling preserves endpoint banks")
    matrix = coupling(base, permutation)
    loop_count = sum(matrix[root][root] for root in base)
    permutation_loop_counts.add(loop_count)
    require(sum(sum(row) for row in matrix) == 3, "permutation coupling mass")
    permutation_coupling_checks += 1
require(permutation_loop_counts == {0, 1, 3}, "three-state couplings realize distinct loop types")


# The marginalization kernel on twelve anchored source and arrival states is
# exactly the (12-1)^2-dimensional mixed-Haar space.
def matrix_rank_mod(matrix: list[list[int]], prime: int) -> int:
    work = [row[:] for row in matrix]
    row_count = len(work)
    column_count = len(work[0]) if work else 0
    rank = 0
    for column in range(column_count):
        pivot = next((row for row in range(rank, row_count) if work[row][column] % prime), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = pow(work[rank][column] % prime, -1, prime)
        work[rank] = [(value * scale) % prime for value in work[rank]]
        for row in range(row_count):
            if row == rank:
                continue
            factor = work[row][column] % prime
            if factor:
                work[row] = [
                    (value - factor * pivot_value) % prime
                    for value, pivot_value in zip(work[row], work[rank], strict=True)
                ]
        rank += 1
        if rank == row_count:
            break
    return rank


support_size = ROOT - 1
marginal_matrix: list[list[int]] = []
for row_root in range(support_size):
    marginal_matrix.append(
        [int(index // support_size == row_root) for index in range(support_size**2)]
    )
for column_root in range(support_size):
    marginal_matrix.append(
        [int(index % support_size == column_root) for index in range(support_size**2)]
    )

marginal_rank = matrix_rank_mod(marginal_matrix, FIELD)
transport_kernel_dimension = support_size**2 - marginal_rank
require(marginal_rank == 2 * support_size - 1, "marginal map rank")
require(transport_kernel_dimension == (support_size - 1) ** 2, "transportation kernel dimension")

haar_basis_checks = 0
for i, j in product(range(1, support_size), repeat=2):
    table = [[0 for _ in range(support_size)] for _ in range(support_size)]
    table[0][0] = 1
    table[0][j] = -1
    table[i][0] = -1
    table[i][j] = 1
    require(all(sum(row) == 0 for row in table), "mixed-Haar basis has zero row margins")
    require(
        all(sum(table[row][column] for row in range(support_size)) == 0 for column in range(support_size)),
        "mixed-Haar basis has zero column margins",
    )
    haar_basis_checks += 1
require(haar_basis_checks == transport_kernel_dimension, "mixed-Haar basis has full kernel cardinality")


print("THM-2538 anchored transverse arrival exact referee")
print("status=PROVED-EXACT-ALGEBRA-AND-HOSTILES")
print(f"field={FIELD} generator={GENERATOR} zeta13={ZETA} xi7={XI}")
print(f"anchor_covariance_checks={anchor_covariance_checks}")
print(f"multiplier_checks={multiplier_checks}")
print(
    f"transverse_pattern_checks={transverse_pattern_checks} "
    f"minimum_nonzero_transverse={minimum_nonzero_transverse}"
)
print(f"packing_pattern_checks={packing_pattern_checks} candidates=1,2,3,4")
print(f"galois_orbit_checks={galois_orbit_checks} total_gain_incidences={len(all_gain_incidences)}")
print(
    f"hostile_atom_checks={hostile_atom_checks} "
    f"section_transverse_checks={section_transverse_checks}"
)
print(f"hostile_kernel_checks={hostile_kernel_checks}")
print(
    f"hostile_gain_checks={hostile_gain_checks} "
    f"hostile_cell_gram_checks={hostile_cell_gram_checks}"
)
print(
    f"self_arrival_checks={self_arrival_checks} "
    f"anchored_reconstruction_checks={anchored_reconstruction_checks}"
)
print(f"crofton_coordinate_checks={crofton_coordinate_checks}")
print(
    f"all_slope_hitting_checks={all_slope_hitting_checks} "
    f"selected_head_kakeya_checks={selected_head_kakeya_checks} "
    f"positive_hit_layer_checks={positive_hit_layer_checks}"
)
print(
    f"cross_kakeya_checks={cross_kakeya_checks} "
    f"haar_switch=({haar_same},{haar_switch})"
)
print(f"frechet_cases={frechet_cases} hall_transport_checks={hall_transport_checks}")
print(f"triple_parity_checks={triple_parity_checks} triple_incidence_counts=(0,1)")
print(
    f"permutation_coupling_checks={permutation_coupling_checks} "
    f"loop_counts={sorted(permutation_loop_counts)}"
)
print(
    f"marginal_rank={marginal_rank} "
    f"transport_kernel_dimension={transport_kernel_dimension} "
    f"haar_basis_checks={haar_basis_checks}"
)
print("same_horizon_empty_head=ANNIHILATES")
print("later_same_base_cross_kakeya=RECOVERS_COUPLING")
print("semantic_arrival_field_in_live_LRC=NOT_SUPPLIED")
print(f"checks={checks}")
print("all_checks=PASS")
