#!/usr/bin/env python3
"""Exact finite checks for THM-2309.

The proofs of the field-section and CRT lemmas are symbolic.  This companion
checks every finite support pattern used by the septimal construction, the
typed exact star packet, the hostile representable matroid, and all numerical
ledgers.  It deliberately uses explicit exceptions rather than ``assert`` so
that ``python -O`` runs the same checks.
"""

from itertools import combinations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def inv_mod(a, p):
    return pow(a % p, -1, p)


def rank_mod(matrix, p):
    a = [[x % p for x in row] for row in matrix]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = inv_mod(a[rank][col], p)
        a[rank] = [(scale * x) % p for x in a[rank]]
        for r in range(rows):
            if r != rank and a[r][col]:
                scale = a[r][col]
                a[r] = [
                    (a[r][c] - scale * a[rank][c]) % p
                    for c in range(cols)
                ]
        rank += 1
        if rank == rows:
            break
    return rank


def det_mod(matrix, p):
    n = len(matrix)
    require(all(len(row) == n for row in matrix), "determinant is not square")
    a = [[x % p for x in row] for row in matrix]
    det = 1
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        value = a[col][col]
        det = (det * value) % p
        scale = inv_mod(value, p)
        a[col] = [(scale * x) % p for x in a[col]]
        for r in range(col + 1, n):
            if a[r][col]:
                factor = a[r][col]
                a[r] = [
                    (a[r][c] - factor * a[col][c]) % p
                    for c in range(n)
                ]
    return det % p


def columns_to_rows(columns):
    return [list(row) for row in zip(*columns)]


def column(matrix, index):
    return [row[index] for row in matrix]


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def exact_star_packet(w, omitted):
    """Rows are ordered by P={j}+U\\{omitted}; coordinate order is j,a,b,U."""
    j, a, b = 0, 1, 2
    units = list(range(3, 9))
    included = [u for u in units if u != omitted]
    pivot = [j] + included
    rows = []
    for k in pivot:
        row = [0] * 9
        row[omitted] = w[k]
        row[k] = -w[omitted]
        rows.append(row)
    ka, kb = included[:2]
    ia = pivot.index(ka)
    ib = pivot.index(kb)
    rows[ia][omitted] += w[a]
    rows[ia][a] -= w[omitted]
    rows[ib][omitted] += w[b]
    rows[ib][b] -= w[omitted]
    return rows, pivot


def first_nonzero_vectors(p, dimension):
    for values in product(range(p), repeat=dimension):
        if any(values):
            yield list(values)


def field_section(v, p, pivot):
    """Construct R with R[pivot]=I, Rv=0, and no zero column."""
    dimension = len(pivot)
    size = len(v)
    complement = [i for i in range(size) if i not in pivot]
    require(len(complement) == 3, "expected three-column complement")

    columns = [[0] * dimension for _ in range(size)]
    for r, index in enumerate(pivot):
        columns[index][r] = 1

    target = [(-v[index]) % p for index in pivot]
    active = [index for index in complement if v[index] % p]
    inactive = [index for index in complement if not v[index] % p]
    require(active, "complement misses scalar support")

    unit_vector = [1] + [0] * (dimension - 1)
    for index in inactive:
        columns[index] = unit_vector[:]

    if len(active) == 1:
        index = active[0]
        scale = inv_mod(v[index], p)
        columns[index] = [(scale * x) % p for x in target]
        require(any(columns[index]), "unique active complement column is zero")
    else:
        residual = target[:]
        for index in active[:-2]:
            columns[index] = unit_vector[:]
            residual = [
                (residual[r] - v[index] * unit_vector[r]) % p
                for r in range(dimension)
            ]
        penultimate, last = active[-2:]
        chosen = None
        last_column = None
        for candidate in first_nonzero_vectors(p, dimension):
            remainder = [
                (residual[r] - v[penultimate] * candidate[r]) % p
                for r in range(dimension)
            ]
            solved = [
                (inv_mod(v[last], p) * remainder[r]) % p
                for r in range(dimension)
            ]
            if any(solved):
                chosen = candidate
                last_column = solved
                break
        require(chosen is not None, "failed to avoid zero final column")
        columns[penultimate] = chosen
        columns[last] = last_column

    rows = columns_to_rows(columns)
    require(rank_mod(rows, p) == dimension, "section lost row rank")
    require(
        all(any(value % p for value in column(rows, i)) for i in range(size)),
        "section has a dark column",
    )
    require(
        all(dot(row, v) % p == 0 for row in rows),
        "section row misses scalar kernel",
    )
    pivot_matrix = [[row[i] for i in pivot] for row in rows]
    require(det_mod(pivot_matrix, p) == 1, "section pivot is not identity")
    return rows


def crt_centered(a7, a13):
    match = next(x for x in range(91) if x % 7 == a7 % 7 and x % 13 == a13 % 13)
    return match if match <= 45 else match - 91


def choose_omitted_for_seven_support(v7):
    targets = [1, 2]
    units = list(range(3, 9))
    if any(v7[index] for index in targets):
        return units[0]
    return next(index for index in units if v7[index])


def check_direct_packet():
    w = [13, 13**3, 2 * 13**5, 1, 14, 27, 40, 53, 66]
    total_height = sum(abs(x) for x in w)
    units = list(range(3, 9))
    for omitted in units:
        rows, pivot = exact_star_packet(w, omitted)
        require(all(dot(row, w) == 0 for row in rows), "star row is not exact")
        pivot_matrix = [[row[i] for i in pivot] for row in rows]
        expected = pow(-w[omitted], 6, 13)
        require(det_mod(pivot_matrix, 13) == expected, "wrong star determinant")
        require(expected != 0, "star determinant is not a thirteen-unit")
        require(
            all(any(value % 13 for value in column(rows, i)) for i in range(9)),
            "grafted star has a dark column",
        )
        require(
            all(max(abs(x) for x in row) <= total_height for row in rows),
            "star row exceeds S",
        )
    return w, total_height


def check_hostile_packet():
    p = 13
    basis = [[int(i == j) for i in range(6)] for j in range(6)]
    units = basis[:5] + [[-1, -1, -1, -1, -1, 0]]
    columns = [
        [1, 2, 0, 0, 0, 0],  # selected j, inside the unit hyperplane
        [0, 0, 0, 0, 0, 1],  # target a, outside
        [0, 0, 1, 2, 0, 0],  # target b, inside
    ] + units
    rows = columns_to_rows(columns)
    scalar = [0, 0, 0] + [1] * 6
    require(rank_mod(rows, p) == 6, "hostile packet is not rank six")
    require(all(dot(row, scalar) % p == 0 for row in rows), "wrong annihilator")
    require(all(any(x % p for x in col) for col in columns), "hostile dark column")
    lam = [1] * 6
    word = [dot(lam, col) % p for col in columns]
    require(word == [3, 1, 3, 1, 1, 1, 1, 1, 8], "wrong all-unit word")
    require(all(word), "hostile word is not all-unit")

    unit_indices = list(range(3, 9))
    pivot_counts = []
    for blocker in range(3):
        count = 0
        for five_units in combinations(unit_indices, 5):
            selected = [blocker] + list(five_units)
            minor = [[row[i] for i in selected] for row in rows]
            count += det_mod(minor, p) != 0
        pivot_counts.append(count)
    require(pivot_counts == [0, 6, 0], "wrong hostile cocircuit spectrum")
    return pivot_counts, word


def check_all_septimal_support_patterns():
    checked = 0
    for mask in range(1 << 9):
        support = [i for i in range(9) if (mask >> i) & 1]
        if len(support) < 2:
            continue
        v = [0] * 9
        for i in support:
            v[i] = 1 + ((3 * i + len(support)) % 6)
        omitted = choose_omitted_for_seven_support(v)
        pivot = [0] + [i for i in range(3, 9) if i != omitted]
        complement = [i for i in range(9) if i not in pivot]
        require(any(v[i] for i in complement), "chosen complement misses support")
        field_section(v, 7, pivot)
        checked += 1
    require(checked == 502, "wrong septimal support-pattern count")
    return checked


def check_crt_lift(w, total_height):
    v7 = [x % 7 for x in w]
    v13 = [x % 13 for x in w]
    require(sum(bool(x) for x in v7) >= 2, "sample lacks septimal support")
    omitted = choose_omitted_for_seven_support(v7)
    pivot = [0] + [i for i in range(3, 9) if i != omitted]
    rows7 = field_section(v7, 7, pivot)
    rows13 = field_section(v13, 13, pivot)
    centered = [
        [crt_centered(rows7[r][c], rows13[r][c]) for c in range(9)]
        for r in range(6)
    ]

    # In the typed sample, coordinate u_1 has scalar value one.
    bezout = [0] * 9
    bezout[3] = 1
    require(dot(bezout, w) == 1, "wrong Bezout vector")
    exact = []
    for row in centered:
        value = dot(row, w)
        require(value % 91 == 0, "CRT row is not in the kernel modulo 91")
        lifted = [row[i] - value * bezout[i] for i in range(9)]
        require(dot(lifted, w) == 0, "Bezout lift is not exact")
        require(
            all((lifted[i] - row[i]) % 91 == 0 for i in range(9)),
            "Bezout lift changed its CRT residue",
        )
        require(
            max(abs(x) for x in lifted) <= 45 * (1 + total_height),
            "Bezout height bound failed",
        )
        exact.append(lifted)

    for p in (7, 13):
        pivot_matrix = [[row[i] for i in pivot] for row in exact]
        require(det_mod(pivot_matrix, p) == 1, "CRT owner pivot failed")
        require(
            all(any(value % p for value in column(exact, i)) for i in range(9)),
            "CRT packet has a dark column",
        )
    return omitted


def main():
    w, total_height = check_direct_packet()
    pivot_counts, hostile_word = check_hostile_packet()
    support_patterns = check_all_septimal_support_patterns()
    crt_omitted = check_crt_lift(w, total_height)

    unanchored_13 = 9 * 12**5
    normalized_13 = 9 * 12**4
    unanchored_7 = 3 * 6**5
    normalized_7 = 3 * 6**4
    unanchored_91 = unanchored_13 * unanchored_7
    normalized_91 = normalized_13 * normalized_7
    require(unanchored_91 == 52242776064, "wrong unanchored CRT constant")
    require(normalized_91 == 725594112, "wrong normalized CRT constant")
    require(6 * 45 // 2 == 135, "wrong CRT height coefficient")

    print("THM-2309 exact companion")
    print("direct_star_omissions=6")
    print(f"typed_scalar_l1_height={total_height}")
    print(f"hostile_owner_pivot_counts={pivot_counts}")
    print(f"hostile_all_unit_word={hostile_word}")
    print(f"septimal_support_patterns_checked={support_patterns}")
    print(f"sample_crt_omitted_coordinate={crt_omitted}")
    print(f"unanchored_13={unanchored_13}")
    print(f"normalized_13={normalized_13}")
    print(f"unanchored_7={unanchored_7}")
    print(f"normalized_7={normalized_7}")
    print(f"unanchored_91={unanchored_91}")
    print(f"normalized_91={normalized_91}")
    print("direct_height=3*S*(13^n-1)")
    print("crt_height=135*(1+S*B(w))*(91^n-1)")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
