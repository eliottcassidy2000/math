#!/usr/bin/env python3
"""Exact independent referee for THM-983's prime common-scale H6 gate.

The proof-bearing object is a six-label subset of F_13^* together with the
directed provider/owner ratio-capacity matrix.  At a prime scale p the only
effective orders are 1 and p.  Mixed order rows die at an order-p owner, so
the exceptional calculation is just the 924 all-order-p supports.

This referee derives the cardinality formula from integer floors, rebuilds
the complete B6/B5 table, checks the formula directly throughout THM-860's
finite scale range, and scans p=23 and p=29 in two equivalent coordinates:
integer capacities and the binary high-ratio Cayley digraph.  The latter is
faithful only for these two exceptional residue rows.  A tournament completion
is reported as lossy telemetry, never used in the emptiness proof.

Only the standard library is used.  Assertions, not optimized-mode-sensitive
``assert`` statements, carry every certificate check.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations, product


MODULUS = 13
LABELS = tuple(range(1, MODULUS))
SUPPORT_SIZE = 6
SCALE_CAP = 840
EXPECTED_B6 = (1, 3, 5, 6, 6, 6, 7, 9, 11, 12, 12, 12)
EXPECTED_B5 = (1, 3, 5, 5, 5, 5, 6, 8, 10, 10, 10, 10)
# Exact spectrum of the p=23 low-ratio Cayley kernel after symmetrization,
# in the characters k=0,...,11 of Z/12.  See spectral_obstruction_23.
SPECTRUM_23 = (5, -1, -2, -1, 2, -1, 1, -1, 2, -1, -2, -1)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ratio(provider: int, owner: int) -> int:
    return provider * pow(owner, -1, MODULUS) % MODULUS


def centered(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue - modulus if 2 * residue > modulus else residue


def crt_base(provider: int, order: int, unit: int) -> int:
    """Algebraic CRT solution x=order*provider (13), x=unit (order)."""
    coefficient = (
        provider - unit * pow(order, -1, MODULUS)
    ) % MODULUS
    base = unit + order * coefficient
    require(base % MODULUS == order * provider % MODULUS,
            "mod-thirteen CRT mismatch")
    require(base % order == unit % order, "effective-order CRT mismatch")
    require(0 <= base < MODULUS * order, "CRT base outside fundamental range")
    return base


def literal_sheet_mask(
    scale: int, provider: int, owner: int, order: int, unit: int
) -> int:
    """Reconstruct the actual affine sheet mask before taking cardinality."""
    base = crt_base(provider, order, unit)
    owner_inverse = pow(owner, -1, MODULUS)
    return sum(
        1 << sheet
        for sheet in range(scale)
        if -order
        < centered(base * (owner_inverse + MODULUS * sheet),
                   MODULUS * order)
        <= order
    )


def residue_bonus(scale_residue: int, provider_owner_ratio: int) -> int:
    """The q-independent part of the order-p sheet cardinality.

    If p=13q+s and a is the representative of sr modulo 13 in {1,...,12},
    then

      #(-p < x <= p, x == pr mod 13)
        = floor((p-a)/13) - floor((-p-a)/13)
        = 2q + floor((s-a)/13) - floor((-s-a)/13).
    """
    require(1 <= scale_residue < MODULUS, "scale residue outside F_13^*")
    require(provider_owner_ratio in LABELS, "ratio outside F_13^*")
    representative = scale_residue * provider_owner_ratio % MODULUS
    require(representative in LABELS, "zero residue in prime-scale formula")
    return (
        (scale_residue - representative) // MODULUS
        - (-scale_residue - representative) // MODULUS
    )


def formula_cardinality(scale: int, provider_owner_ratio: int) -> int:
    quotient, residue = divmod(scale, MODULUS)
    require(residue != 0, "common scale must be coprime to thirteen")
    return 2 * quotient + residue_bonus(residue, provider_owner_ratio)


def direct_cardinality(scale: int, provider_owner_ratio: int) -> int:
    """Literal integer-window count, independent of the floor formula."""
    target = scale * provider_owner_ratio % MODULUS
    return sum(
        1
        for integer in range(-scale + 1, scale + 1)
        if integer % MODULUS == target
    )


def bonus_tables() -> tuple[tuple[int, ...], tuple[int, ...]]:
    six_largest = []
    five_largest = []
    for residue in LABELS:
        row = sorted(
            (residue_bonus(residue, r) for r in LABELS), reverse=True
        )
        # The length-2s window contains zero once, but provider/owner ratios
        # range only over F_13^*.  Hence the nonzero-residue bonus mass is
        # 2s-1 rather than 2s.
        require(sum(row) == 2 * residue - 1,
                "nonzero bonus mass is not 2s-1")
        require(set(row) <= {0, 1, 2}, "unexpected bonus height")
        six_largest.append(sum(row[:6]))
        five_largest.append(sum(row[:5]))
    return tuple(six_largest), tuple(five_largest)


def is_prime(integer: int) -> bool:
    if integer < 2:
        return False
    divisor = 2
    while divisor * divisor <= integer:
        if integer % divisor == 0:
            return False
        divisor += 1
    return True


def owner_capacity(scale: int, support: tuple[int, ...], owner: int) -> int:
    return sum(formula_cardinality(scale, ratio(provider, owner))
               for provider in support)


def state_cardinality(
    scale: int, provider: int, owner: int, order: int
) -> int:
    if order == 1:
        return scale if provider == owner else 0
    require(order == scale, "non-prime effective order")
    return formula_cardinality(scale, ratio(provider, owner))


def high_ratios(scale: int) -> frozenset[int]:
    _quotient, residue = divmod(scale, MODULUS)
    row = {r: residue_bonus(residue, r) for r in LABELS}
    height = max(row.values())
    return frozenset(r for r, value in row.items() if value == height)


def high_predicate(scale: int, support: tuple[int, ...], owner: int) -> bool:
    """Binary Cayley-graph threshold, valid exactly at p=23 and p=29."""
    high = high_ratios(scale)
    high_count = sum(ratio(provider, owner) in high for provider in support)
    if scale == 23:
        return high_count >= 5
    if scale == 29:
        return high_count == 5
    raise ValueError("binary exceptional predicate requested at another scale")


def multiply_support(
    support: tuple[int, ...], multiplier: int
) -> tuple[int, ...]:
    return tuple(sorted(multiplier * label % MODULUS for label in support))


def support_orbits(
    supports: tuple[tuple[int, ...], ...]
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    unseen = set(supports)
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = tuple(sorted({multiply_support(seed, a) for a in LABELS}))
        require(set(orbit) <= unseen, "multiplication orbits overlap")
        unseen.difference_update(orbit)
        orbits.append(orbit)
    return tuple(orbits)


def generated_subgroup(generators: frozenset[int]) -> frozenset[int]:
    subgroup = {1}
    changed = True
    while changed:
        changed = False
        for left in tuple(subgroup):
            for right in generators | frozenset(subgroup):
                product = left * right % MODULUS
                if product not in subgroup:
                    subgroup.add(product)
                    changed = True
    return frozenset(subgroup)


def exponent_table(generator: int) -> dict[int, int]:
    table = {}
    value = 1
    for exponent in range(len(LABELS)):
        require(value not in table, "generator cycle repeated early")
        table[value] = exponent
        value = value * generator % MODULUS
    require(value == 1 and set(table) == set(LABELS),
            "chosen label is not a generator of F_13^*")
    return table


def low_edge_count_23(support: tuple[int, ...]) -> int:
    low = set(LABELS) - set(high_ratios(23))
    return sum(
        ratio(provider, owner) in low
        for provider in support
        for owner in support
    )


def spectral_obstruction_23() -> tuple[tuple[int, ...], int]:
    """Return the exact low-kernel spectrum and its six-set Rayleigh floor.

    With generator 2, the low ratios at p=23 have exponent set
    D={2,3,6,8,9} in Z/12.  Symmetrization has Fourier eigenvalue

      cos(pi*k/3) + 2*cos(pi*k/2) + cos(2*pi*k/3) + (-1)^k,

    giving SPECTRUM_23 exactly.  Its minimum is -2.  For a six-set indicator
    x=(1/2)1+y, ||y||^2=3, so its internal directed low-edge count is
    x^T M x >= 5*36/12 - 2*3 = 9.  Six feasible owners would instead give
    at most one low provider per owner, hence at most six low edges.
    """
    exponents = exponent_table(2)
    low = set(LABELS) - set(high_ratios(23))
    low_exponents = tuple(sorted(exponents[label] for label in low))
    require(low_exponents == (2, 3, 6, 8, 9),
            "p=23 low-ratio exponent set mismatch")
    require(min(SPECTRUM_23[1:]) == -2,
            "p=23 nontrivial spectral floor mismatch")
    rayleigh_floor = 5 * SUPPORT_SIZE * SUPPORT_SIZE // len(LABELS) - 2 * 3
    require(rayleigh_floor == 9, "p=23 Rayleigh floor mismatch")
    return low_exponents, rayleigh_floor


def strongly_connected_components(
    vertices: tuple[int, ...], edges: dict[int, set[int]]
) -> tuple[tuple[int, ...], ...]:
    """Small deterministic Kosaraju implementation."""
    visited: set[int] = set()
    finish: list[int] = []

    def forward(vertex: int) -> None:
        visited.add(vertex)
        for neighbor in sorted(edges[vertex]):
            if neighbor not in visited:
                forward(neighbor)
        finish.append(vertex)

    for vertex in vertices:
        if vertex not in visited:
            forward(vertex)

    reverse = {vertex: set() for vertex in vertices}
    for source in vertices:
        for target in edges[source]:
            reverse[target].add(source)

    visited.clear()
    components = []

    def backward(vertex: int, component: list[int]) -> None:
        visited.add(vertex)
        component.append(vertex)
        for neighbor in sorted(reverse[vertex]):
            if neighbor not in visited:
                backward(neighbor, component)

    for vertex in reversed(finish):
        if vertex not in visited:
            component: list[int] = []
            backward(vertex, component)
            components.append(tuple(sorted(component)))
    return tuple(sorted(components, key=lambda c: (len(c), c)))


def hamiltonian_path_count(
    vertices: tuple[int, ...], edges: dict[int, set[int]]
) -> int:
    index = {vertex: position for position, vertex in enumerate(vertices)}
    size = 1 << len(vertices)
    paths = [[0] * len(vertices) for _ in range(size)]
    for position in range(len(vertices)):
        paths[1 << position][position] = 1
    for mask in range(size):
        for last in range(len(vertices)):
            count = paths[mask][last]
            if count == 0:
                continue
            source = vertices[last]
            for target in edges[source]:
                position = index[target]
                if mask & (1 << position) == 0:
                    paths[mask | (1 << position)][position] += count
    return sum(paths[-1])


def redei_path(
    vertices: tuple[int, ...], edges: dict[int, set[int]]
) -> tuple[int, ...]:
    path: list[int] = []
    for vertex in vertices:
        position = next(
            (i for i, old in enumerate(path) if old in edges[vertex]),
            len(path),
        )
        path.insert(position, vertex)
    require(
        all(path[i + 1] in edges[path[i]] for i in range(len(path) - 1)),
        "Redei insertion did not produce a Hamiltonian path",
    )
    return tuple(path)


def tournament_fingerprint(scale: int) -> dict[str, object]:
    """Complete bonus-comparison tournament; natural label order breaks ties."""
    _quotient, residue = divmod(scale, MODULUS)
    edges = {vertex: set() for vertex in LABELS}
    tie_edges = 0
    flips = 0
    for left, right in combinations(LABELS, 2):
        left_value = residue_bonus(residue, ratio(left, right))
        right_value = residue_bonus(residue, ratio(right, left))
        if left_value >= right_value:
            winner, loser = left, right
            tie_edges += left_value == right_value
        else:
            winner, loser = right, left
            flips += 1
        edges[winner].add(loser)
    score_histogram = Counter(len(edges[vertex]) for vertex in LABELS)
    triangles = sum(
        (
            second in edges[first]
            and third in edges[second]
            and first in edges[third]
        )
        or (
            third in edges[first]
            and second in edges[third]
            and first in edges[second]
        )
        for first, second, third in combinations(LABELS, 3)
    )
    return {
        "score_histogram": tuple(sorted(score_histogram.items())),
        "directed_triangles": triangles,
        "scc_sizes": tuple(sorted(map(len, strongly_connected_components(LABELS, edges)))),
        "tie_edges": tie_edges,
        "flips_from_natural_order": flips,
        "redei_path": redei_path(LABELS, edges),
        "hamiltonian_path_count": hamiltonian_path_count(LABELS, edges),
    }


def cayley_fingerprint(scale: int) -> dict[str, object]:
    high = high_ratios(scale)
    edges = {
        provider: {
            owner
            for owner in LABELS
            if owner != provider and ratio(provider, owner) in high
        }
        for provider in LABELS
    }
    indegrees = Counter(
        sum(vertex in edges[source] for source in LABELS) for vertex in LABELS
    )
    outdegrees = Counter(len(edges[vertex]) for vertex in LABELS)
    triangles = sum(
        (
            second in edges[first]
            and third in edges[second]
            and first in edges[third]
        )
        + (
            third in edges[first]
            and second in edges[third]
            and first in edges[second]
        )
        for first, second, third in combinations(LABELS, 3)
    )
    return {
        "high_ratios": tuple(sorted(high)),
        "outdegree_histogram": tuple(sorted(outdegrees.items())),
        "indegree_histogram": tuple(sorted(indegrees.items())),
        "directed_triangles": triangles,
        "scc_sizes": tuple(sorted(map(len, strongly_connected_components(LABELS, edges)))),
    }


def format_histogram(histogram: Counter[int]) -> str:
    return " ".join(f"{key}:{histogram[key]}" for key in sorted(histogram))


def format_fingerprint(fingerprint: dict[str, object]) -> str:
    return " ".join(f"{key}={value}" for key, value in fingerprint.items())


def main() -> None:
    require(len(LABELS) == 12, "label grammar mismatch")
    supports = tuple(combinations(LABELS, SUPPORT_SIZE))
    require(len(supports) == 924, "six-support census mismatch")

    order_words = tuple(
        word
        for word in product((1, 23), repeat=SUPPORT_SIZE)
        if all(23 in word[:omitted] + word[omitted + 1 :]
               for omitted in range(SUPPORT_SIZE))
    )
    require(len(order_words) == 57, "prime hereditary order grammar mismatch")
    require(all(sum(order == 23 for order in word) >= 2 for word in order_words),
            "hereditary prime word has fewer than two full orders")

    b6, b5 = bonus_tables()
    require(b6 == EXPECTED_B6, "B6 table mismatch")
    require(b5 == EXPECTED_B5, "B5 table mismatch")
    require(max(b6[s - 1] - s for s in LABELS) == 2,
            "six-largest symbolic tail bound mismatch")
    require(max(b5[s - 1] - s for s in LABELS) == 2,
            "five-largest symbolic tail bound mismatch")

    # The direct window implementation and derived floor formula agree on the
    # complete THM-860 range.  The displayed floor identity itself is valid
    # for arbitrary q, so the theorem's prime statement is not range-limited.
    direct_checks = 0
    for scale in range(1, SCALE_CAP + 1):
        if scale % MODULUS == 0:
            continue
        for provider_owner_ratio in LABELS:
            require(
                direct_cardinality(scale, provider_owner_ratio)
                == formula_cardinality(scale, provider_owner_ratio),
                "direct/floor cardinality mismatch",
            )
            direct_checks += 1

    # A mixed row has an order-p owner.  Every order-one provider vanishes
    # there, and the at-most-five order-p providers pay at most 10q+B5(s).
    # This is strictly less than p for every q>=1, uniformly in s.
    mixed_deficits = tuple(3 + s - b5[s - 1] for s in LABELS)
    require(min(mixed_deficits) > 0, "mixed-row deficit is not uniformly positive")

    # For all-order-p rows the scalar ceiling is 12q+B6(s).  q>=3 is
    # uniformly impossible.  Within THM-860's range only 23 and 29 survive
    # this numerical ceiling among primes p>=19.
    prime_candidates = tuple(
        scale
        for scale in range(19, SCALE_CAP + 1)
        if is_prime(scale)
        and 12 * (scale // MODULUS) + b6[scale % MODULUS - 1] >= scale
    )
    require(prime_candidates == (23, 29), "prime ceiling exceptions mismatch")

    # Reconstruct actual CRT sheet masks at the only two scalar-ceiling
    # exceptions.  This independently checks both the order-one delta law and
    # unit-independence of the order-p formula before the support scan.
    literal_mask_checks = 0
    for scale in prime_candidates:
        for provider in LABELS:
            for owner in LABELS:
                order_one_size = literal_sheet_mask(
                    scale, provider, owner, 1, 0
                ).bit_count()
                require(
                    order_one_size
                    == state_cardinality(scale, provider, owner, 1),
                    "order-one literal mask law mismatch",
                )
                literal_mask_checks += 1
                for unit in range(1, scale):
                    size = literal_sheet_mask(
                        scale, provider, owner, scale, unit
                    ).bit_count()
                    require(
                        size == state_cardinality(
                            scale, provider, owner, scale
                        ),
                        "order-p literal mask cardinality depends on unit",
                    )
                    literal_mask_checks += 1

    orbits = support_orbits(supports)
    orbit_histogram = Counter(len(orbit) for orbit in orbits)
    require(len(orbits) == 80, "multiplication-orbit count mismatch")
    require(orbit_histogram == Counter({2: 1, 4: 1, 6: 3, 12: 75}),
            "multiplication-orbit histogram mismatch")

    low_exponents_23, spectral_floor_23 = spectral_obstruction_23()
    low_edge_histogram_23 = Counter(low_edge_count_23(support)
                                    for support in supports)
    require(
        low_edge_histogram_23
        == Counter({10: 30, 12: 208, 13: 288, 14: 216,
                    16: 96, 17: 48, 18: 38}),
        "p=23 internal low-edge histogram mismatch",
    )
    require(min(low_edge_histogram_23) >= spectral_floor_23 > SUPPORT_SIZE,
            "p=23 spectral contradiction does not clear six")

    bonus_digest = sha256()
    for residue in LABELS:
        for provider_owner_ratio in LABELS:
            bonus_digest.update(bytes((residue, provider_owner_ratio,
                                       residue_bonus(residue, provider_owner_ratio))))

    print("THM-983 independent prime common-scale H6 scalar referee")
    print("scope: primitive proper AP-centred common-scale Hamming-six scalar gate only")
    print("formula: a_p(r)=2q+floor((s-a)/13)-floor((-s-a)/13), a=sr mod 13")
    print(f"direct formula checks through c={SCALE_CAP}: {direct_checks}")
    print("literal exceptional CRT-mask checks:", literal_mask_checks)
    print("B6:", " ".join(map(str, b6)))
    print("B5:", " ".join(map(str, b5)))
    print("mixed minimum deficit p-(10q+B5):", min(mixed_deficits))
    print("all-order prime ceiling exceptions p>=19:", prime_candidates)
    print("support orbits:", format_histogram(orbit_histogram))
    print("bonus digest:", bonus_digest.hexdigest())
    print("p=23 low exponents in Z/12:", low_exponents_23)
    print("p=23 symmetrized low-kernel spectrum:", SPECTRUM_23)
    print("p=23 spectral six-set floor:", spectral_floor_23)
    print("p=23 exact low-edge histogram:",
          format_histogram(low_edge_histogram_23))

    expected_ledgers = {
        23: {
            "feasible": Counter({0: 204, 1: 528, 2: 138, 3: 48, 4: 6}),
            "minimum": Counter({19: 6, 20: 204, 21: 650, 22: 64}),
            "max_feasible": 4,
            "ledger_digest": "70b8beb8154a021f9c91fd10690b6d41447b98d09a7ce087ac16159662fe5382",
            "witness_digest": "a49dd90d0a50d5815586ff42f3cb62f5a4eb302366692fe7b31bcc9c42a6cd25",
        },
        29: {
            "feasible": Counter({0: 846, 1: 72, 2: 6}),
            "minimum": Counter({25: 156, 26: 660, 27: 108}),
            "max_feasible": 2,
            "ledger_digest": "3090ed7e26e420c90813aeaacb87acadf582212350d33cdf1778f62c7f042374",
            "witness_digest": "03b35f0e6f7d9d3b706cd6f48b491a76637e7e3127383af841e47ded22c3f29f",
        },
    }

    for scale in (23, 29):
        feasible_histogram: Counter[int] = Counter()
        minimum_capacity_histogram: Counter[int] = Counter()
        ledger_digest = sha256()
        witness_digest = sha256()
        maximum_feasible = 0
        for support in supports:
            capacities = tuple(owner_capacity(scale, support, owner)
                               for owner in support)
            feasible = tuple(capacity >= scale for capacity in capacities)
            binary = tuple(high_predicate(scale, support, owner)
                           for owner in support)
            require(feasible == binary,
                    "capacity and high-ratio support scans disagree")
            feasible_count = sum(feasible)
            feasible_histogram[feasible_count] += 1
            minimum_capacity_histogram[min(capacities)] += 1
            maximum_feasible = max(maximum_feasible, feasible_count)
            failing_owner = next(owner for owner, ok in zip(support, feasible) if not ok)
            ledger_digest.update(bytes(support + capacities + feasible))
            witness_digest.update(bytes(support + (failing_owner,)))

        expected = expected_ledgers[scale]
        require(feasible_histogram == expected["feasible"],
                "feasible-owner histogram mismatch")
        require(minimum_capacity_histogram == expected["minimum"],
                "minimum-capacity histogram mismatch")
        require(maximum_feasible == expected["max_feasible"],
                "exceptional support maximum mismatch")
        require(maximum_feasible < SUPPORT_SIZE,
                "exceptional scalar support bank is nonempty")
        require(ledger_digest.hexdigest() == expected["ledger_digest"],
                "support ledger digest mismatch")
        require(witness_digest.hexdigest() == expected["witness_digest"],
                "failure-witness digest mismatch")

        print(f"p={scale} high ratios:", tuple(sorted(high_ratios(scale))))
        print(f"p={scale} feasible-owner histogram:",
              format_histogram(feasible_histogram))
        print(f"p={scale} minimum-capacity histogram:",
              format_histogram(minimum_capacity_histogram))
        print(f"p={scale} maximum feasible owners: {maximum_feasible}/6")
        print(f"p={scale} support ledger digest: {ledger_digest.hexdigest()}")
        print(f"p={scale} failure-witness digest: {witness_digest.hexdigest()}")
        print(f"p={scale} Cayley carrier:",
              format_fingerprint(cayley_fingerprint(scale)))
        print(f"p={scale} lossy tournament:",
              format_fingerprint(tournament_fingerprint(scale)))

    subgroup_29 = generated_subgroup(high_ratios(29))
    require(subgroup_29 == frozenset(LABELS),
            "p=29 high ratios do not generate F_13^*")
    print("p=23 structural certificate: low-edge Rayleigh floor 9 > 6")
    print("p=29 closure certificate: <high ratios> =", tuple(sorted(subgroup_29)))
    print("carrier audit: ratio Cayley digraph + embedded six-support is faithful;")
    print("  tournament completion forgets absolute capacities and symmetric ties;")
    print("  unlabelled sheets, units, gaps, and divisor vertices destroy owner incidence.")
    print("PASS: p=23 and p=29 exceptional 924-support banks are empty")
    print("HONEST LIMIT: composite scales, H5, non-AP/deep sheets, and global n=12 remain open")


if __name__ == "__main__":
    main()
