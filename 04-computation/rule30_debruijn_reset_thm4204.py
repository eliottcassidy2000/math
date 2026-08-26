#!/usr/bin/env python3
"""Exact audit for THM-4204: Rule 30 de Bruijn reset and saturation."""

from collections import Counter, deque
from itertools import product


StateMap = tuple[int, int, int, int]
Matrix = tuple[tuple[int, int, int, int], ...]


def check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def rule30(left: int, center: int, right: int) -> int:
    return left ^ (center | right)


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(4))
            for j in range(4)
        )
        for i in range(4)
    )


IDENTITY: Matrix = tuple(
    tuple(int(i == j) for j in range(4)) for i in range(4)
)


def debruijn_matrix(output: int) -> Matrix:
    rows = [[0] * 4 for _ in range(4)]
    for left, center, right in product((0, 1), repeat=3):
        if rule30(left, center, right) != output:
            continue
        source = 2 * left + center
        target = 2 * center + right
        rows[source][target] += 1
    return tuple(tuple(row) for row in rows)


A0 = debruijn_matrix(0)
A1 = debruijn_matrix(1)
GENERATORS = (A0, A1)

EXPECTED_A0: Matrix = (
    (1, 0, 0, 0),
    (0, 0, 0, 0),
    (0, 1, 0, 0),
    (0, 0, 1, 1),
)
EXPECTED_A1: Matrix = (
    (0, 1, 0, 0),
    (0, 0, 1, 1),
    (1, 0, 0, 0),
    (0, 0, 0, 0),
)


def matrix_word(word: str) -> Matrix:
    result = IDENTITY
    for symbol in word:
        result = matmul(result, GENERATORS[int(symbol)])
    return result


def matrix_trace(matrix: Matrix) -> int:
    return sum(matrix[i][i] for i in range(4))


def state_map(matrix: Matrix) -> StateMap:
    image = []
    for column in range(4):
        hits = [row for row in range(4) if matrix[row][column] == 1]
        check(len(hits) == 1, f"column {column} is not deterministic")
        check(
            sum(matrix[row][column] for row in range(4)) == 1,
            f"column {column} has nonunit mass",
        )
        image.append(hits[0])
    return tuple(image)  # type: ignore[return-value]


F0 = state_map(A0)
F1 = state_map(A1)
MAP_GENERATORS = (F0, F1)
MAP_IDENTITY: StateMap = (0, 1, 2, 3)


def compose(left: StateMap, right: StateMap) -> StateMap:
    return tuple(left[right[j]] for j in range(4))  # type: ignore[return-value]


def map_word(word: str) -> StateMap:
    result = MAP_IDENTITY
    for symbol in word:
        result = compose(result, MAP_GENERATORS[int(symbol)])
    return result


def transformation_rank(mapping: StateMap) -> int:
    return len(set(mapping))


def constant_matrix(state: int) -> Matrix:
    return tuple(
        tuple(int(row == state) for _ in range(4)) for row in range(4)
    )


def generated_monoid() -> tuple[list[StateMap], dict[StateMap, str]]:
    states = [MAP_IDENTITY]
    first_word = {MAP_IDENTITY: ""}
    queue = deque([MAP_IDENTITY])
    while queue:
        current = queue.popleft()
        for bit, generator in enumerate(MAP_GENERATORS):
            successor = compose(current, generator)
            if successor not in first_word:
                first_word[successor] = first_word[current] + str(bit)
                states.append(successor)
                queue.append(successor)
    return states, first_word


def fixed_points(mapping: StateMap) -> int:
    return sum(mapping[state] == state for state in range(4))


def minimize_trace_automaton(
    states: list[StateMap],
) -> tuple[list[tuple[bool, int, int]], list[list[int]]]:
    index = {state: i for i, state in enumerate(states)}
    transitions = [
        [index[compose(state, generator)] for generator in MAP_GENERATORS]
        for state in states
    ]
    unique = [fixed_points(state) == 1 for state in states]

    colors = [int(flag) for flag in unique]
    while True:
        signatures = [
            (unique[i], colors[transitions[i][0]], colors[transitions[i][1]])
            for i in range(len(states))
        ]
        palette = {signature: j for j, signature in enumerate(sorted(set(signatures)))}
        refined = [palette[signature] for signature in signatures]
        if refined == colors:
            break
        colors = refined

    blocks_by_color: dict[int, list[int]] = {}
    for i, color in enumerate(colors):
        blocks_by_color.setdefault(color, []).append(i)
    blocks = list(blocks_by_color.values())
    block_of = {
        state_index: block_index
        for block_index, block in enumerate(blocks)
        for state_index in block
    }
    raw_table: list[tuple[bool, int, int]] = []
    for block in blocks:
        representative = block[0]
        row = (
            unique[representative],
            block_of[transitions[representative][0]],
            block_of[transitions[representative][1]],
        )
        for member in block:
            check(unique[member] == row[0], "acceptance is not block-constant")
            check(
                block_of[transitions[member][0]] == row[1]
                and block_of[transitions[member][1]] == row[2],
                "transition is not block-constant",
            )
        raw_table.append(row)

    start_raw = block_of[0]
    order = []
    seen = {start_raw}
    queue = deque([start_raw])
    while queue:
        state = queue.popleft()
        order.append(state)
        for successor in raw_table[state][1:]:
            if successor not in seen:
                seen.add(successor)
                queue.append(successor)
    check(len(order) == len(blocks), "minimized quotient is not reachable")
    canonical = {raw: new for new, raw in enumerate(order)}
    table = [
        (
            raw_table[raw][0],
            canonical[raw_table[raw][1]],
            canonical[raw_table[raw][2]],
        )
        for raw in order
    ]
    canonical_blocks = [blocks[raw] for raw in order]
    return table, canonical_blocks


def matrix_add(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    return [
        [left[i][j] + right[i][j] for j in range(len(left))]
        for i in range(len(left))
    ]


def integer_matrix_product(
    left: list[list[int]], right: list[list[int]]
) -> list[list[int]]:
    size = len(left)
    return [
        [sum(left[i][k] * right[k][j] for k in range(size)) for j in range(size)]
        for i in range(size)
    ]


def integer_matrix_power(matrix: list[list[int]], exponent: int) -> list[list[int]]:
    size = len(matrix)
    result = [[int(i == j) for j in range(size)] for i in range(size)]
    base = matrix
    power = exponent
    while power:
        if power & 1:
            result = integer_matrix_product(result, base)
        base = integer_matrix_product(base, base)
        power >>= 1
    return result


def matrix_vector(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [
        sum(matrix[i][j] * vector[j] for j in range(len(vector)))
        for i in range(len(vector))
    ]


def center_direct(length: int) -> str:
    row = {0: 1}
    result = []
    for time in range(length):
        result.append(str(row.get(0, 0)))
        next_row = {}
        for site in range(-time - 1, time + 2):
            value = rule30(
                row.get(site - 1, 0),
                row.get(site, 0),
                row.get(site + 1, 0),
            )
            if value:
                next_row[site] = 1
        row = next_row
    return "".join(result)


def center_packed(length: int) -> str:
    packed = 1
    result = []
    for time in range(length):
        result.append(str((packed >> time) & 1))
        packed ^= (packed << 1) | (packed << 2)
    return "".join(result)


def isolated_support(time: int) -> str:
    row = {0: 1}
    for current_time in range(time):
        next_row = {}
        for site in range(-current_time - 1, current_time + 2):
            value = rule30(
                row.get(site - 1, 0),
                row.get(site, 0),
                row.get(site + 1, 0),
            )
            if value:
                next_row[site] = 1
        row = next_row
    return "".join(str(row.get(site, 0)) for site in range(-time, time + 1))


def embedded_isolated_row(time: int, length: int) -> str:
    check(length >= 2 * time + 1, "ring is too short for the isolated row")
    support = isolated_support(time)
    bits = ["0"] * length
    for offset, bit in enumerate(support, start=-time):
        bits[offset % length] = bit
    return "".join(bits)


def cyclic_image(word: str) -> str:
    length = len(word)
    bits = [int(bit) for bit in word]
    return "".join(
        str(rule30(bits[(i - 1) % length], bits[i], bits[(i + 1) % length]))
        for i in range(length)
    )


def direct_preimage_histogram(length: int) -> Counter[str]:
    counts: Counter[str] = Counter()
    for value in range(1 << length):
        word = "".join(str((value >> i) & 1) for i in range(length))
        counts[cyclic_image(word)] += 1
    return counts


def depth_preimage_count(target: str, depth: int) -> int:
    length = len(target)
    count = 0
    for value in range(1 << length):
        word = "".join(str((value >> i) & 1) for i in range(length))
        for _ in range(depth):
            word = cyclic_image(word)
        count += int(word == target)
    return count


def contains_cyclic_010(word: str) -> bool:
    length = len(word)
    return any(
        word[i] == "0"
        and word[(i + 1) % length] == "1"
        and word[(i + 2) % length] == "0"
        for i in range(length)
    )


def inverse_word_from_fixed_state(word: str, fixed_state: int) -> str:
    states = [0] * (len(word) + 1)
    states[-1] = fixed_state
    for i in range(len(word) - 1, -1, -1):
        states[i] = MAP_GENERATORS[int(word[i])][states[i + 1]]
    check(states[0] == fixed_state, "supplied state is not fixed")
    return "".join(str(states[i] & 1) for i in range(len(word)))


def main() -> None:
    check(A0 == EXPECTED_A0, "independently generated A0 changed")
    check(A1 == EXPECTED_A1, "independently generated A1 changed")
    check(F0 == (0, 2, 3, 3), "unexpected backward map f0")
    check(F1 == (2, 0, 1, 1), "unexpected backward map f1")
    print("PASS  local rule independently regenerates the two de Bruijn matrices")
    print(f"      f0={F0}, f1={F1}; every matrix column has mass one")

    reset_by_length: dict[int, list[str]] = {}
    for length in range(1, 5):
        words = [
            "".join(bits)
            for bits in product("01", repeat=length)
            if transformation_rank(map_word("".join(bits))) == 1
        ]
        reset_by_length[length] = words
    check(reset_by_length[1] == [], "length-one reset appeared")
    check(reset_by_length[2] == [], "length-two reset appeared")
    check(reset_by_length[3] == [], "length-three reset appeared")
    check(reset_by_length[4] == ["0010", "1010"], "reset threshold changed")
    check(matrix_word("0010") == constant_matrix(3), "0010 is not E3")
    check(matrix_word("1010") == constant_matrix(1), "1010 is not E1")
    for state in range(4):
        for generator in GENERATORS:
            check(
                matmul(constant_matrix(state), generator) == constant_matrix(state),
                "constant transformation is not suffix-absorbing",
            )
    print("PASS  shortest reset threshold is four")
    print("      length-4 reset words: 0010 -> E3, 1010 -> E1")

    monoid, first_word = generated_monoid()
    ranks = Counter(transformation_rank(state) for state in monoid)
    check(len(monoid) == 41, "transformation monoid size changed")
    check(
        sum(bool(word) for word in first_word.values()) == 40,
        "nonempty-word semigroup size changed",
    )
    check(ranks == Counter({1: 4, 2: 30, 3: 6, 4: 1}), "rank census changed")
    constants = {state for state in monoid if transformation_rank(state) == 1}
    check(
        constants == {tuple([state] * 4) for state in range(4)},
        "rank-one ideal is not the four constant maps",
    )
    print("PASS  generated inverse-boundary transformation monoid is exact")
    print("      monoid size=41, nonempty semigroup size=40;")
    print("      monoid ranks 1/2/3/4 occur 4/30/6/1 times")

    table, blocks = minimize_trace_automaton(monoid)
    expected_table = [
        (False, 1, 2),
        (False, 1, 3),
        (False, 4, 5),
        (True, 6, 7),
        (True, 4, 8),
        (False, 9, 0),
        (True, 6, 6),
        (False, 1, 1),
        (False, 6, 10),
        (False, 9, 11),
        (False, 4, 4),
        (False, 6, 12),
        (True, 9, 9),
    ]
    check(table == expected_table, "canonical minimized trace DFA changed")
    check([len(block) for block in blocks] == [2, 5, 1, 5, 2, 1, 10, 5, 2, 2, 2, 2, 2], "DFA block sizes changed")
    print("PASS  Moore minimization gives a 13-state unique-fibre DFA")
    print("      state: unique?, on-0, on-1")
    for state, row in enumerate(table):
        print(f"      {state:2d}: {int(row[0])} {row[1]:2d} {row[2]:2d}")

    size = len(table)
    transition_count = [[0] * size for _ in range(size)]
    for state, (_, zero, one) in enumerate(table):
        transition_count[state][zero] += 1
        transition_count[state][one] += 1
    nonunique = [int(not row[0]) for row in table]
    polynomial = integer_matrix_power(transition_count, 5)
    for coefficient, exponent in [(-3, 2), (-2, 1), (-2, 0)]:
        power = integer_matrix_power(transition_count, exponent)
        polynomial = matrix_add(
            polynomial,
            [[coefficient * entry for entry in row] for row in power],
        )
    check(
        matrix_vector(polynomial, nonunique) == [0] * size,
        "recurrence certificate P(T)b=0 failed",
    )

    distributions = [0] * len(monoid)
    distributions[0] = 1
    unique_counts = []
    nonunique_counts = []
    index = {state: i for i, state in enumerate(monoid)}
    for length in range(65):
        unique_count = sum(
            count
            for state, count in zip(monoid, distributions)
            if fixed_points(state) == 1
        )
        unique_counts.append(unique_count)
        nonunique_counts.append((1 << length) - unique_count)
        next_distribution = [0] * len(monoid)
        for state, count in zip(monoid, distributions):
            for generator in MAP_GENERATORS:
                next_distribution[index[compose(state, generator)]] += count
        distributions = next_distribution
    check(nonunique_counts[:5] == [1, 2, 2, 5, 10], "initial V_N changed")
    for length in range(5, len(nonunique_counts)):
        expected = (
            3 * nonunique_counts[length - 3]
            + 2 * nonunique_counts[length - 4]
            + 2 * nonunique_counts[length - 5]
        )
        check(nonunique_counts[length] == expected, f"V recurrence failed at {length}")
    print("PASS  exact nonunique-output recurrence and rational OGF certificate")
    print("      V0..V14=" + ",".join(map(str, nonunique_counts[:15])))
    print("      V_N=3V_(N-3)+2V_(N-4)+2V_(N-5), N>=5")
    print("      V(z)=(1+2z+2z^2+2z^3+2z^4)/(1-3z^3-2z^4-2z^5)")

    for length in range(1, 10):
        histogram = direct_preimage_histogram(length)
        for bits in product("01", repeat=length):
            word = "".join(bits)
            check(
                histogram[word] == matrix_trace(matrix_word(word)),
                f"direct/matrix inverse count mismatch for N={length}, word={word}",
            )
            if length >= 3 and contains_cyclic_010(word):
                check(histogram[word] == 1, f"010 reset failed for {word}")
        direct_unique = sum(count == 1 for count in histogram.values())
        check(direct_unique == unique_counts[length], f"U_{length} mismatch")
    check(
        transformation_rank(map_word("0101")) == 2
        and transformation_rank(map_word("1010")) == 1
        and fixed_points(map_word("0101")) == 1,
        "cyclic-rotation rank hostile changed",
    )
    print("PASS  direct ring enumeration agrees for every output through N=9")
    print("      hostile: M_0101 has rank two, while rotated M_1010 is constant")

    direct_center = center_direct(64)
    packed_center = center_packed(64)
    check(direct_center == packed_center, "direct and packed center prefixes differ")
    seed_prefix = "110111001100010"
    check(direct_center.startswith(seed_prefix), "seed center prefix changed")
    prefix_products = [matrix_word(direct_center[:length]) for length in range(1, 65)]
    check(
        all(transformation_rank(matrix) > 1 for matrix in prefix_products[:14]),
        "center de Bruijn product reset before length 15",
    )
    check(prefix_products[14] == constant_matrix(1), "length-15 product is not E1")
    check(
        all(matrix == constant_matrix(1) for matrix in prefix_products[14:]),
        "center prefix product did not remain E1",
    )
    print("PASS  isolated-seed center prefix first resets at length 15")
    print(f"      c_0...c_14={seed_prefix}; product=E1 for every N>=15")
    print("      dyadic product ranks/traces:")
    for length in (1, 2, 4, 8, 16, 32, 64):
        matrix = prefix_products[length - 1]
        print(
            f"      N={length:2d}: rank={transformation_rank(state_map(matrix))}, "
            f"trace={matrix_trace(matrix)}"
        )

    target = seed_prefix
    histogram_15 = direct_preimage_histogram(15)
    check(histogram_15[target] == 1, "length-15 center word is not uniquely invertible")
    predecessor = inverse_word_from_fixed_state(target, 1)
    check(predecessor == "100010000111010", "length-15 predecessor changed")
    check(cyclic_image(predecessor) == target, "direct predecessor image mismatch")
    print("PASS  independent length-15 cyclic inversion")
    print(f"      {predecessor} -> {target}")

    p_matrix = matmul(A0, A0)
    reset_zero = constant_matrix(0)
    check(
        p_matrix
        == (
            (1, 0, 0, 0),
            (0, 0, 0, 0),
            (0, 0, 0, 0),
            (0, 1, 1, 1),
        ),
        "A0^2 changed",
    )
    power = IDENTITY
    for exponent in range(1, 65):
        power = matmul(power, A0)
        if exponent >= 2:
            check(power == p_matrix, f"A0 power did not stabilize at {exponent}")
    check(matrix_word("00" + "110010") == reset_zero, "even row reset failed")
    check(matrix_word("00" + "11011110") == reset_zero, "odd row reset failed")
    check(matrix_trace(matrix_word("00" + "11001")) == 1, "time-two head failed")
    check(matrix_trace(matrix_word("00" + "1101111")) == 1, "time-three head failed")
    check(matrix_trace(matrix_word("00" + "111")) == 2, "time-one fork failed")
    for time in range(4, 65):
        support = isolated_support(time)
        if time % 2 == 0:
            check(support.startswith("110010"), f"even edge prefix failed at {time}")
        elif time >= 5:
            check(support.startswith("11011110"), f"odd edge prefix failed at {time}")
    print("PASS  physical-row edge reset is uniform in time")
    print("      even edge 110010, odd edge 11011110; two leading zeros give E0")

    for background, previous, center, following in product((0, 1), repeat=4):
        left = background ^ previous
        middle = background ^ center
        right = background ^ following
        direct = rule30(left, middle, right)
        defect_formula = (
            previous
            ^ ((1 ^ background) & (center ^ following))
            ^ (center & following)
        )
        check(direct == defect_formula, "two-background Laurent formula failed")
    for length in range(5, 12):
        sparse = ["0"] * length
        sparse[0] = "1"
        dense = ["1"] * length
        dense[-1] = "0"
        dense[0] = "0"
        check(
            cyclic_image("".join(sparse)) == cyclic_image("".join(dense)),
            f"sparse/dense coalescence failed at N={length}",
        )
    print("PASS  all 16 local background-defect cases and sparse/dense fork")
    print("      F(delta_0)=F(1+delta_-1+delta_0)")

    for length in range(2, 33):
        seed_word = "1" + "0" * (length - 1)
        dense_cyclic_word = "00" + "1" * (length - 2)
        seed_count = matrix_trace(matrix_word(seed_word))
        dense_count = matrix_trace(matrix_word(dense_cyclic_word))
        check(seed_count == 1, f"seed predecessor count failed at N={length}")
        expected_dense = (2, 1, 0)[(length - 2) % 3]
        check(dense_count == expected_dense, f"dense fork count failed at N={length}")
    for length in range(5, 12):
        histogram = direct_preimage_histogram(length)
        seed_word = "1" + "0" * (length - 1)
        dense_cyclic_word = "00" + "1" * (length - 2)
        check(histogram[seed_word] == 1, f"direct seed count failed at N={length}")
        check(
            histogram[dense_cyclic_word] == (2, 1, 0)[(length - 2) % 3],
            f"direct dense count failed at N={length}",
        )
        for time in range(1, (length - 3) // 2 + 1):
            physical = embedded_isolated_row(time, length)
            expected = 2 if time == 1 else 1
            check(
                histogram[physical] == expected,
                f"physical predecessor count failed at N={length}, t={time}",
            )
    print("PASS  direct physical-row and two-background ring counts through N=11")
    print("      #preimages(W_t)=1 for t>=2; #preimages(W_1)=2")

    depth_controls = []
    for dyadic_time, length, expected in ((2, 8, 3), (4, 16, 1)):
        target_row = embedded_isolated_row(dyadic_time, length)
        count = depth_preimage_count(target_row, dyadic_time + 1)
        check(count == expected, f"depth ancestry failed at h={dyadic_time}")
        depth_controls.append((dyadic_time, length, count))
    print("PASS  independent exhaustive dyadic depth-ancestry controls")
    for dyadic_time, length, count in depth_controls:
        print(f"      h={dyadic_time}, N={length}, depth={dyadic_time + 1}: {count}")

    print("\nSCOPE")
    print("  PROVED: shortest spatial reset; exact unique-preimage counting law;")
    print("          physical-row backward rigidity; cubic dyadic ancestry;")
    print("          center-prefix and dyadic-prefix de Bruijn saturation.")
    print("  OPEN:   Rule 30 center nonperiodicity, density, and complexity prizes.")
    print("  LOSS:   the temporal center word was retyped as one spatial cyclic output;")
    print("          after c_14 the product state contains no later center bit.")


if __name__ == "__main__":
    main()
