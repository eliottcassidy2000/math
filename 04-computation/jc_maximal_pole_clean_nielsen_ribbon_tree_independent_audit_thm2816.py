from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, factorial


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def all_pairings(points: tuple[int, ...]):
    if not points:
        yield ()
        return
    first = points[0]
    for index in range(1, len(points)):
        second = points[index]
        remainder = points[1:index] + points[index + 1 :]
        for tail in all_pairings(remainder):
            yield tuple(sorted(((first, second),) + tail))


def all_matchings(vertex_count: int, edge_count: int):
    for endpoints in combinations(range(vertex_count), 2 * edge_count):
        yield from all_pairings(endpoints)


def chords_cross(first: tuple[int, int], second: tuple[int, int]) -> bool:
    a, b = sorted(first)
    c, d = sorted(second)
    return (a < c < b < d) or (c < a < d < b)


def is_noncrossing(matching: tuple[tuple[int, int], ...]) -> bool:
    return not any(
        chords_cross(matching[left], matching[right])
        for left in range(len(matching))
        for right in range(left + 1, len(matching))
    )


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def cycle_decomposition(permutation: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    seen: set[int] = set()
    answer = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycle = []
        point = start
        while point not in seen:
            seen.add(point)
            cycle.append(point)
            point = permutation[point]
        answer.append(tuple(cycle))
    return tuple(answer)


def pole_cycles(
    vertex_count: int, matching: tuple[tuple[int, int], ...]
) -> tuple[tuple[int, ...], ...]:
    involution = list(range(vertex_count))
    for left, right in matching:
        involution[left], involution[right] = right, left
    full_cycle = tuple((index + 1) % vertex_count for index in range(vertex_count))
    return cycle_decomposition(compose(tuple(involution), full_cycle))


def rotate_matching(
    matching: tuple[tuple[int, int], ...], shift: int, vertex_count: int
) -> tuple[tuple[int, int], ...]:
    return tuple(
        sorted(
            tuple(
                sorted(((left + shift) % vertex_count, (right + shift) % vertex_count))
            )
            for left, right in matching
        )
    )


def rotation_orbit_key(
    matching: tuple[tuple[int, int], ...],
    labelled_cycles: tuple[tuple[int, ...], ...],
    vertex_count: int,
):
    candidates = []
    for shift in range(vertex_count):
        shifted_cycles = tuple(
            tuple(sorted((point + shift) % vertex_count for point in cycle))
            for cycle in labelled_cycles
        )
        candidates.append(
            (rotate_matching(matching, shift, vertex_count), shifted_cycles)
        )
    return min(candidates)


def connected_tree(vertex_count: int, edges: tuple[tuple[int, int], ...]) -> bool:
    adjacency = [set() for _ in range(vertex_count)]
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    reached = {0}
    frontier = [0]
    while frontier:
        vertex = frontier.pop()
        for neighbor in adjacency[vertex]:
            if neighbor not in reached:
                reached.add(neighbor)
                frontier.append(neighbor)
    return len(edges) == vertex_count - 1 and len(reached) == vertex_count


def canonical_rotation(entries: tuple[tuple[int, int], ...]):
    return min(entries[index:] + entries[:index] for index in range(len(entries)))


def extracted_ribbon_signature(
    matching: tuple[tuple[int, int], ...],
    labelled_cycles: tuple[tuple[int, ...], ...],
):
    face_of: dict[int, int] = {}
    for face, cycle in enumerate(labelled_cycles):
        for point in cycle:
            face_of[point] = face

    partner: dict[int, int] = {}
    dual_edges = []
    for left, right in matching:
        partner[left], partner[right] = right, left
        left_face, right_face = face_of[left], face_of[right]
        require(left_face != right_face, "planar chord became a dual loop")
        dual_edges.append(tuple(sorted((left_face, right_face))))
    edge_tuple = tuple(sorted(dual_edges))
    require(len(set(edge_tuple)) == len(edge_tuple), "dual multiedge found")
    require(connected_tree(len(labelled_cycles), edge_tuple), "dual is not a tree")

    vertex_data = []
    for face, cycle in enumerate(labelled_cycles):
        marked_positions = []
        for position, point in enumerate(cycle):
            if point in partner:
                marked_positions.append((position, face_of[partner[point]]))
        require(marked_positions, "dual-tree vertex has degree zero")
        corner_data = []
        for index, (position, neighbor) in enumerate(marked_positions):
            next_position = marked_positions[(index + 1) % len(marked_positions)][0]
            gap = (next_position - position) % len(cycle)
            if gap == 0:
                gap = len(cycle)
            require(gap > 0, "nonpositive corner arc")
            corner_data.append((neighbor, gap))
        require(
            sum(gap for _, gap in corner_data) == len(cycle),
            "corner arcs do not sum to pole length",
        )
        vertex_data.append(canonical_rotation(tuple(corner_data)))
    return edge_tuple, tuple(vertex_data)


def positive_compositions(total: int, length: int):
    if length == 1:
        yield (total,)
        return
    for first in range(1, total - length + 2):
        for tail in positive_compositions(total - first, length - 1):
            yield (first,) + tail


def prufer_tree(word: tuple[int, ...], vertex_count: int):
    degrees = [1] * vertex_count
    for vertex in word:
        degrees[vertex] += 1
    edges = []
    for vertex in word:
        leaf = next(index for index, degree in enumerate(degrees) if degree == 1)
        edges.append(tuple(sorted((leaf, vertex))))
        degrees[leaf] -= 1
        degrees[vertex] -= 1
    final = [index for index, degree in enumerate(degrees) if degree == 1]
    require(len(final) == 2, "Pruefer decoder did not end with two leaves")
    edges.append(tuple(sorted(final)))
    return tuple(sorted(edges))


def abstract_ribbon_signatures(passport: tuple[int, ...]):
    vertex_count = len(passport)
    signatures = set()
    for word in product(range(vertex_count), repeat=max(0, vertex_count - 2)):
        edges = prufer_tree(tuple(word), vertex_count)
        adjacency = [set() for _ in range(vertex_count)]
        for left, right in edges:
            adjacency[left].add(right)
            adjacency[right].add(left)

        local_options = []
        admissible = True
        for vertex, neighbors in enumerate(adjacency):
            degree = len(neighbors)
            if passport[vertex] < degree:
                admissible = False
                break
            first = min(neighbors)
            cyclic_orders = (
                (first,) + tail
                for tail in permutations(sorted(neighbors - {first}))
            )
            compositions = tuple(positive_compositions(passport[vertex], degree))
            options = []
            for cyclic_order in cyclic_orders:
                for corners in compositions:
                    options.append(
                        tuple(zip(cyclic_order, corners))
                    )
            local_options.append(tuple(options))
        if not admissible:
            continue
        for vertex_data in product(*local_options):
            signatures.add((edges, tuple(vertex_data)))
    return signatures


def ordered_compositions(total: int, length: int):
    yield from positive_compositions(total, length)


def audit_permutation_case(edge_count: int, vertex_count: int):
    genus_histogram: dict[int, int] = defaultdict(int)
    planar_records = []
    for matching in all_matchings(vertex_count, edge_count):
        faces = pole_cycles(vertex_count, matching)
        genus_numerator = edge_count + 1 - len(faces)
        require(genus_numerator >= 0, "negative chord genus")
        require(genus_numerator % 2 == 0, "nonintegral chord genus")
        genus = genus_numerator // 2
        noncrossing = is_noncrossing(matching)
        require(
            (genus == 0) == noncrossing,
            "genus-zero/noncrossing equivalence failed",
        )
        genus_histogram[genus] += 1
        if noncrossing:
            planar_records.append((matching, faces))

    catalan = comb(2 * edge_count, edge_count) // (edge_count + 1)
    require(
        len(planar_records) == comb(vertex_count, 2 * edge_count) * catalan,
        "noncrossing Catalan control failed",
    )

    orbit_sets: dict[tuple[int, ...], set] = defaultdict(set)
    signature_sets: dict[tuple[int, ...], set] = defaultdict(set)
    orbit_to_signature = {}
    signature_to_orbit = {}
    for matching, faces in planar_records:
        for labelled_cycles in permutations(faces):
            passport = tuple(len(cycle) for cycle in labelled_cycles)
            orbit = rotation_orbit_key(matching, labelled_cycles, vertex_count)
            signature = extracted_ribbon_signature(matching, labelled_cycles)
            old_signature = orbit_to_signature.setdefault(orbit, signature)
            require(old_signature == signature, "orbit changed ribbon signature")
            old_orbit = signature_to_orbit.setdefault(signature, orbit)
            require(old_orbit == orbit, "ribbon signature did not reconstruct orbit")
            orbit_sets[passport].add(orbit)
            signature_sets[passport].add(signature)

    passports = tuple(ordered_compositions(vertex_count, edge_count + 1))
    require(set(orbit_sets) == set(passports), "ordered passport universe mismatch")
    expected = factorial(edge_count - 1) * comb(
        vertex_count - edge_count - 1, edge_count - 1
    )
    reconstructed = 0
    for passport in passports:
        abstract = abstract_ribbon_signatures(passport)
        require(
            signature_sets[passport] == abstract,
            f"ribbon reconstruction mismatch for {passport}",
        )
        require(len(orbit_sets[passport]) == expected, "Nielsen count mismatch")
        reconstructed += len(abstract)
    return (
        sum(genus_histogram.values()),
        len(planar_records),
        len(passports),
        expected,
        reconstructed,
        tuple(sorted(genus_histogram.items())),
    )


def polynomial_trim(poly: list[Fraction]) -> list[Fraction]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def polynomial_add(left: list[Fraction], right: list[Fraction]):
    size = max(len(left), len(right))
    answer = [Fraction(0)] * size
    for index in range(size):
        answer[index] = (
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
        )
    return polynomial_trim(answer)


def polynomial_scale(poly: list[Fraction], scalar: Fraction):
    return polynomial_trim([scalar * coefficient for coefficient in poly])


def polynomial_multiply(left: list[Fraction], right: list[Fraction]):
    answer = [Fraction(0)] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            answer[left_degree + right_degree] += left_coefficient * right_coefficient
    return polynomial_trim(answer)


def polynomial_power(poly: list[Fraction], exponent: int):
    answer = [Fraction(1)]
    for _ in range(exponent):
        answer = polynomial_multiply(answer, poly)
    return answer


def polynomial_derivative(poly: list[Fraction]):
    if len(poly) == 1:
        return [Fraction(0)]
    return [degree * poly[degree] for degree in range(1, len(poly))]


def polynomial_evaluate(poly: list[Fraction], value: Fraction):
    answer = Fraction(0)
    for coefficient in reversed(poly):
        answer = answer * value + coefficient
    return answer


def polynomial_divmod(numerator: list[Fraction], denominator: list[Fraction]):
    remainder = numerator[:]
    quotient = [Fraction(0)] * max(1, len(numerator) - len(denominator) + 1)
    while len(remainder) >= len(denominator) and remainder != [0]:
        degree = len(remainder) - len(denominator)
        coefficient = remainder[-1] / denominator[-1]
        quotient[degree] += coefficient
        subtraction = [Fraction(0)] * degree + polynomial_scale(denominator, coefficient)
        remainder = polynomial_add(remainder, polynomial_scale(subtraction, -1))
    return polynomial_trim(quotient), polynomial_trim(remainder)


def audit_response_converse(edge_count: int, pole_parts: tuple[int, int]):
    require(edge_count == 1, "explicit response control is the one-chord chamber")
    left_part, right_part = pole_parts
    degree = left_part + right_part
    x = [Fraction(0), Fraction(1)]
    x_minus_one = [Fraction(-1), Fraction(1)]
    D = polynomial_multiply(
        polynomial_power(x, left_part),
        polynomial_power(x_minus_one, right_part),
    )
    gamma = Fraction(left_part, degree)
    A = -polynomial_evaluate(D, gamma)
    require(A != 0, "third-fibre constant vanished")
    B = polynomial_add(D, [A])
    E = [ -gamma, Fraction(1)]
    S, remainder = polynomial_divmod(B, polynomial_power(E, 2))
    require(remainder == [0], "double-zero division failed")
    require(polynomial_evaluate(S, gamma) != 0, "zero had multiplicity above two")
    T = polynomial_multiply(x, x_minus_one)
    derivative_gate = polynomial_add(
        polynomial_multiply(polynomial_derivative(D), T),
        polynomial_scale(polynomial_multiply(D, E), -degree),
    )
    require(derivative_gate == [0], "D'=NDE/T gate failed")

    C = -degree * A
    # For F=B/D, B'=D', hence the numerator of F' is -A D'.
    # The target formula C E/(D T) is equivalent after clearing D^2 T
    # to -A D' T=C E D.
    response_gate = polynomial_add(
        polynomial_scale(polynomial_multiply(polynomial_derivative(D), T), -A),
        polynomial_scale(polynomial_multiply(E, D), -C),
    )
    require(response_gate == [0], "response logarithmic derivative gate failed")
    require(len(S) - 1 == degree - 2, "simple-zero degree failed")
    return C


permutation_cases = (
    *((1, degree) for degree in range(2, 9)),
    *((2, degree) for degree in range(4, 10)),
    *((3, degree) for degree in range(6, 11)),
    *((4, degree) for degree in range(8, 11)),
    (5, 10),
)

print("THM-2816 INDEPENDENT MAXIMAL-POLE RIBBON-TREE AUDIT")
print("e,N | all matchings | planar | passports | per passport | signatures")
total_matchings = 0
total_signatures = 0
for edge_count, vertex_count in permutation_cases:
    (
        all_count,
        planar_count,
        passport_count,
        per_passport,
        signature_count,
        genus_histogram,
    ) = audit_permutation_case(edge_count, vertex_count)
    total_matchings += all_count
    total_signatures += signature_count
    print(
        f"{edge_count},{vertex_count} | {all_count} | {planar_count} | "
        f"{passport_count} | {per_passport} | {signature_count} | "
        f"genus={genus_histogram}"
    )

response_controls = []
for pole_parts in ((1, 1), (1, 4), (2, 3), (4, 5), (7, 2)):
    response_controls.append((pole_parts, audit_response_converse(1, pole_parts)))

print(f"total_all_matchings={total_matchings}")
print(f"total_reconstructed_labelled_ribbon_signatures={total_signatures}")
print(f"response_converse_controls={response_controls}")
print("source_normalization=high point infinity; beta_1=0; beta_2=1")
print("formula=(e-1)!*binom(N-e-1,e-1)")
print("ALL INDEPENDENT EXACT CHECKS PASSED")
