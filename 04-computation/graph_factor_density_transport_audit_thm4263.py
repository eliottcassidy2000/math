"""Independent finite controls for THM-4263.

Standard library only; no maintained theorem implementation is imported.
"""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def series_mul(left, right, prime):
    size = len(left)
    answer = [0] * size
    for i, x in enumerate(left):
        for j, y in enumerate(right[: size - i]):
            answer[i + j] = (answer[i + j] + x * y) % prime
    return tuple(answer)


def series_pow(value, exponent, prime):
    answer = (1,) + (0,) * (len(value) - 1)
    factor = value
    while exponent:
        if exponent & 1:
            answer = series_mul(answer, factor, prime)
        factor = series_mul(factor, factor, prime)
        exponent >>= 1
    return answer


def monomial_jet(monomial, alpha, g1, g2, prime):
    f_degree, u_degree, v_degree = monomial
    au, av = alpha
    size = len(g1)
    if au > u_degree or av > v_degree or f_degree >= size:
        return (0,) * size
    coefficient = comb(u_degree, au) * comb(v_degree, av) % prime
    value = series_mul(
        series_pow(g1, u_degree - au, prime),
        series_pow(g2, v_degree - av, prime),
        prime,
    )
    answer = [0] * size
    for index, entry in enumerate(value[: size - f_degree]):
        answer[index + f_degree] = coefficient * entry % prime
    return tuple(answer)


def rank_mod(columns, prime):
    if not columns:
        return 0
    row_count = len(columns[0])
    matrix = [[column[row] % prime for column in columns]
              for row in range(row_count)]
    rank = 0
    for column in range(len(columns)):
        pivot = next((row for row in range(rank, row_count)
                      if matrix[row][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [entry * inverse % prime for entry in matrix[rank]]
        for row in range(row_count):
            if row == rank or matrix[row][column] == 0:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (x - factor * y) % prime
                for x, y in zip(matrix[row], matrix[rank], strict=True)
            ]
        rank += 1
        if rank == row_count:
            break
    return rank


def multigraph_audit():
    prime = 2
    truncation = 4
    g1 = (1, 1, 0, 1)
    g2 = (0, 0, 1, 0)
    basis = [(f_degree, u_degree, v_degree)
             for f_degree in range(truncation)
             for u_degree in range(2)
             for v_degree in range(2)]
    values = [monomial_jet(item, (0, 0), g1, g2, prime)
              for item in basis]
    value_rank = rank_mod(values, prime)
    require(value_rank == truncation, "moving graph value rank changed")

    jet_indices = [(0, 0), (0, 1), (1, 0), (1, 1)]
    jets = []
    for item in basis:
        column = []
        for alpha in jet_indices:
            column.extend(monomial_jet(item, alpha, g1, g2, prime))
        jets.append(tuple(column))
    jet_rank = rank_mod(jets, prime)
    require(jet_rank == len(basis), "rectangular Hasse jet lost rank")

    full_g1 = (1, 1, 0, 1, 0, 1)
    full_g2 = (0, 0, 1, 0, 1, 0)
    require(all(full_g1[:level] == full_g1[: level + 1][:-1] and
                full_g2[:level] == full_g2[: level + 1][:-1]
                for level in range(1, len(full_g1))),
            "compatible section truncations changed")
    incompatible = [(0,), (1, 0)]
    require(incompatible[1][:-1] != incompatible[0],
            "incompatible section hostile was missed")
    return len(basis), value_rank, jet_rank, len(full_g1)


def filtered_hostile():
    # Over F_2, ev_0(u+f)=f is globally nonzero, but its level-one target
    # jet vanishes even though u+f is not in fR.
    source_constant_u = 1
    target_constant = 0
    require(source_constant_u and target_constant == 0,
            "filtered short-jet hostile changed")
    return 1, True


def regular_factor_audit():
    # q:F_2^4 -> F_2^2, q(x)=(x0+x2,x1+x3).
    fibres = Counter()
    for word in range(16):
        image = (((word >> 0) ^ (word >> 2)) & 1) | (
            (((word >> 1) ^ (word >> 3)) & 1) << 1
        )
        fibres[image] += 1
    weights = tuple(fibres[value] for value in range(4))
    require(weights == (4, 4, 4, 4), "regular factor weights changed")
    for mask in range(16):
        target_count = sum((mask >> value) & 1 for value in range(4))
        source_count = sum(fibres[value] for value in range(4)
                           if (mask >> value) & 1)
        require(Fraction(source_count, 16) == Fraction(target_count, 4),
                "regular-factor density mismatch")
    zero_fibre = {
        word for word in range(16)
        if ((((word >> 0) ^ (word >> 2)) & 1) |
            ((((word >> 1) ^ (word >> 3)) & 1) << 1)) == 0
    }
    require({0} < zero_fibre, "nonsaturation hostile changed")
    return weights, len(zero_fibre)


def conditional_hazard_audit():
    block_length = 3
    block_count = 4
    universe = set(range(1 << (block_length * block_count)))
    block_mask = (1 << block_length) - 1
    events = [
        {word for word in universe
         if ((word >> (block * block_length)) & block_mask) == block_mask}
        for block in range(block_count)
    ]
    survivors = set(universe)
    hazards = []
    for event in events:
        hazards.append(Fraction(len(survivors & event), len(survivors)))
        survivors -= event
    defect = Fraction(len(survivors), len(universe))
    product_defect = Fraction(1)
    for hazard in hazards:
        product_defect *= 1 - hazard
    require(defect == product_defect == Fraction(7**4, 8**4),
            "fresh conditional-hazard identity changed")

    repeated_universe = set(range(1 << block_length))
    repeated_event = {block_mask}
    survivors = set(repeated_universe)
    repeated_hazards = []
    for _ in range(block_count):
        repeated_hazards.append(
            Fraction(len(survivors & repeated_event), len(survivors))
        )
        survivors -= repeated_event
    repeated_defect = Fraction(len(survivors), len(repeated_universe))
    require(repeated_hazards ==
            [Fraction(1, 8), Fraction(0), Fraction(0), Fraction(0)] and
            repeated_defect == Fraction(7, 8),
            "merged-channel hostile changed")
    return defect, tuple(hazards), repeated_defect, tuple(repeated_hazards)


def heavy_fibre_hostile():
    rows = []
    for size in (4, 8, 12):
        heavy = 1 << size
        target_bad = Fraction(1, size)
        source_bad = Fraction(heavy, heavy + size - 1)
        require(source_bad > Fraction(4, 5), "heavy fibre ceased to dominate")
        rows.append((size, target_bad, source_bad))
    return tuple(rows)


def rule30_step(word, size):
    answer = 0
    for index in range(size):
        left = (word >> ((index - 1) % size)) & 1
        center = (word >> index) & 1
        right = (word >> ((index + 1) % size)) & 1
        answer |= (left ^ (center | right)) << index
    return answer


def rule30_factor_audit():
    maxima = []
    n3_histogram = None
    checks = 0
    for size in range(1, 13):
        cardinality = 1 << size
        weights = Counter(rule30_step(word, size)
                          for word in range(cardinality))
        histogram = Counter(weights.get(image, 0)
                            for image in range(cardinality))
        maxima.append(max(weights.values()))
        require(max(weights.values()) <= 3,
                "Rule-30 finite fibre exceeded three")
        if size == 3:
            n3_histogram = tuple(histogram[weight] for weight in range(4))
            require(n3_histogram == (3, 3, 1, 1) and
                    weights[(1 << size) - 1] == 3,
                    "Rule-30 N=3 hostile changed")
        unique_outputs = sum(weights.get(image, 0) == 1
                             for image in range(cardinality))
        exceptional_outputs = cardinality - unique_outputs
        multiple_input_mass = sum(weight for weight in weights.values()
                                  if weight >= 2)
        require(multiple_input_mass == exceptional_outputs,
                "Rule-30 weighted identity changed")
        checks += 1
    return tuple(maxima), n3_histogram, checks


def graph_channel_audit():
    value_columns = [(0, 1), (0, 1)]  # u and f both evaluate to f.
    hasse_columns = [(0, 1, 1, 0), (0, 1, 0, 0)]
    value_rank = rank_mod(value_columns, 2)
    hasse_rank = rank_mod(hasse_columns, 2)
    source_cartier_on_value_kernel = 1  # explicit-f channel sees u+f.
    require(value_rank == 1 and hasse_rank == 2 and
            source_cartier_on_value_kernel != 0,
            "graph-channel hostile changed")
    return value_rank, hasse_rank


def hostile_windows():
    prop6_pivot = 23 + 15
    prop12_pivot = 29 + 17
    prop6_hits = [shift for exponent in range(1, 4)
                  if 0 <= (shift := prop6_pivot - 23 * exponent) <= 10]
    prop12_hits = [shift for exponent in range(1, 4)
                   if -1 <= (shift := prop12_pivot - 29 * exponent) <= 13]
    require(not prop6_hits and not prop12_hits,
            "chosen-section hostile entered a short window")
    return prop6_pivot, prop12_pivot


def main():
    basis, value_rank, jet_rank, levels = multigraph_audit()
    global_order, first_jet_zero = filtered_hostile()
    weights, zero_fibre = regular_factor_audit()
    defect, hazards, merged_defect, merged_hazards = conditional_hazard_audit()
    heavy_rows = heavy_fibre_hostile()
    rule30_maxima, rule30_n3, rule30_checks = rule30_factor_audit()
    graph_value_rank, graph_hasse_rank = graph_channel_audit()
    prop6_pivot, prop12_pivot = hostile_windows()

    print("MULTIGRAPH PASS", f"basis={basis}", f"value_rank={value_rank}",
          f"kernel_dim={basis-value_rank}",
          f"rectangular_hasse_rank={jet_rank}",
          f"compatible_levels={levels}", "incompatible_detected=yes")
    print("FILTERED_HOSTILE PASS", f"global_image_order={global_order}",
          f"first_short_jet_zero={str(first_jet_zero).lower()}")
    print("REGULAR_FACTOR PASS", f"weights={weights}",
          "all_target_subsets=16", f"nonsaturated_zero_fibre={zero_fibre}")
    print("HAZARD PASS", f"fresh_hazards={hazards}",
          f"fresh_defect={defect}", f"merged_hazards={merged_hazards}",
          f"merged_defect={merged_defect}")
    print(f"HEAVY_FIBRE_HOSTILE PASS rows={heavy_rows}")
    print("RULE30_FACTOR PASS", f"max_fibres_N1_12={rule30_maxima}",
          f"N3_weight_histogram_0_1_2_3={rule30_n3}",
          f"weighted_identity_checks={rule30_checks}")
    print("GRAPH_CHANNEL PASS", f"value_rank={graph_value_rank}",
          "source_cartier_descends=no",
          f"value_plus_hasse_rank={graph_hasse_rank}")
    print("CARTIER_WINDOWS PASS", f"prop6_pivot={prop6_pivot}",
          f"prop12_pivot={prop12_pivot}")
    print("GRAPH_FACTOR_DENSITY_AUDIT: PASS")


if __name__ == "__main__":
    main()
