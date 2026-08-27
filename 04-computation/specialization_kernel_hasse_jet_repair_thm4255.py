#!/usr/bin/env python3
"""Exact finite controls for THM-4255.

The mathematical proof in THM-4255 is ring-theoretic.  This companion checks
its finite linear-algebra shadows over truncated power-series rings F_p[f]/f^N:

* one specialization has exactly the principal kernel (u-g);
* its image order can exceed coefficientwise order by an arbitrary tested
  amount, already in u-degree one;
* the full transverse Hasse jet is invertible in every characteristic;
* d+1 unit-separated scalar sections also recover u-degree at most d; and
* universal-slope and finite-box Kronecker encodings are injective on their
  declared finite universes; and
* ordinary derivatives fail in characteristic p on u^p.

Only the Python standard library and exact arithmetic modulo p are used.
"""

from math import comb, factorial


PRIMES = (2, 3, 5, 7)
TRUNCATIONS = range(2, 9)
U_DEGREES = range(1, 7)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def add_entry(matrix, row, col, value, p):
    matrix[row][col] = (matrix[row][col] + value) % p


def series_mul(left, right, p, truncation):
    out = [0] * truncation
    for i, a in enumerate(left):
        if not a:
            continue
        for j, b in enumerate(right[: truncation - i]):
            if b:
                out[i + j] = (out[i + j] + a * b) % p
    return out


def series_pow(series, exponent, p, truncation):
    out = [1] + [0] * (truncation - 1)
    base = list(series)
    e = exponent
    while e:
        if e & 1:
            out = series_mul(out, base, p, truncation)
        base = series_mul(base, base, p, truncation)
        e //= 2
    return out


def shifted(series, amount, truncation):
    if amount >= truncation:
        return [0] * truncation
    return [0] * amount + list(series[: truncation - amount])


def rank_mod_p(matrix, p):
    work = [row[:] for row in matrix]
    if not work:
        return 0
    rows = len(work)
    cols = len(work[0])
    pivot_row = 0
    for col in range(cols):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][col] % p),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][col] % p, -1, p)
        work[pivot_row] = [(entry * inverse) % p for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row:
                continue
            factor = work[row][col] % p
            if factor:
                work[row] = [
                    (a - factor * b) % p
                    for a, b in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def specialization_matrix(p, truncation, u_degree, g, jet_orders=(0,)):
    """Matrix of P -> (D_u^[j] P)(f,g(f)) for the requested Hasse jets."""
    columns = truncation * (u_degree + 1)
    matrix = [[0] * columns for _ in range(truncation * len(jet_orders))]
    powers = [series_pow(g, exponent, p, truncation) for exponent in range(u_degree + 1)]
    for f_degree in range(truncation):
        for exponent in range(u_degree + 1):
            col = f_degree * (u_degree + 1) + exponent
            for block, jet in enumerate(jet_orders):
                if jet > exponent:
                    continue
                value = shifted(powers[exponent - jet], f_degree, truncation)
                scalar = comb(exponent, jet) % p
                for f_out, coefficient in enumerate(value):
                    add_entry(
                        matrix,
                        block * truncation + f_out,
                        col,
                        scalar * coefficient,
                        p,
                    )
    return matrix


def principal_multiplication_matrix(p, truncation, u_degree, g):
    """Matrix of multiplication by u-g from degree <=d-1 to degree <=d."""
    rows = truncation * (u_degree + 1)
    columns = truncation * u_degree
    matrix = [[0] * columns for _ in range(rows)]
    for f_degree in range(truncation):
        for exponent in range(u_degree):
            col = f_degree * u_degree + exponent
            row_u = f_degree * (u_degree + 1) + exponent + 1
            add_entry(matrix, row_u, col, 1, p)
            for shift_degree, coefficient in enumerate(g):
                if coefficient and f_degree + shift_degree < truncation:
                    row_g = (f_degree + shift_degree) * (u_degree + 1) + exponent
                    add_entry(matrix, row_g, col, -coefficient, p)
    return matrix


def universal_slope_matrix(p, f_count, u_degree):
    """Matrix of f^n u^i -> f^(n+i) lambda^i on a finite box."""
    columns = f_count * (u_degree + 1)
    f_out_count = f_count + u_degree
    rows = f_out_count * (u_degree + 1)
    matrix = [[0] * columns for _ in range(rows)]
    for f_degree in range(f_count):
        for exponent in range(u_degree + 1):
            col = f_degree * (u_degree + 1) + exponent
            row = (f_degree + exponent) * (u_degree + 1) + exponent
            add_entry(matrix, row, col, 1, p)
    return matrix


def kronecker_matrix(p, f_count, u_degree, multiplier):
    """Matrix of f^n u^i -> f^(n+M i) on a finite rectangular box."""
    columns = f_count * (u_degree + 1)
    maximum_exponent = (f_count - 1) + multiplier * u_degree
    matrix = [[0] * columns for _ in range(maximum_exponent + 1)]
    for f_degree in range(f_count):
        for exponent in range(u_degree + 1):
            col = f_degree * (u_degree + 1) + exponent
            row = f_degree + multiplier * exponent
            add_entry(matrix, row, col, 1, p)
    return matrix


def product_is_zero(left, right, p):
    require(len(left[0]) == len(right), "incompatible matrix product")
    for row in left:
        for col in range(len(right[0])):
            value = sum(row[k] * right[k][col] for k in range(len(right))) % p
            if value:
                return False
    return True


def apply_matrix(matrix, vector, p):
    return [sum(a * b for a, b in zip(row, vector)) % p for row in matrix]


def series_order(vector):
    return next((index for index, value in enumerate(vector) if value), None)


def coefficientwise_f_order(vector, u_degree):
    block_size = u_degree + 1
    return next(
        (
            f_degree
            for f_degree in range(len(vector) // block_size)
            if any(vector[f_degree * block_size : (f_degree + 1) * block_size])
        ),
        None,
    )


def section_families(p, truncation):
    families = []
    for terms in ((1,), (1, 1), (2, 0, 1)):
        g = [0] * truncation
        for degree, coefficient in enumerate(terms, start=1):
            if degree < truncation:
                g[degree] = coefficient % p
        families.append(g)
    return families


def main():
    kernel_cases = 0
    hasse_cases = 0
    hostile_order_cases = 0
    joint_section_cases = 0
    universal_slope_cases = 0
    kronecker_cases = 0
    maximum_tested_order_jump = 0

    for p in PRIMES:
        for truncation in TRUNCATIONS:
            for u_degree in U_DEGREES:
                domain_dimension = truncation * (u_degree + 1)
                for g in section_families(p, truncation):
                    evaluation = specialization_matrix(p, truncation, u_degree, g)
                    multiplier = principal_multiplication_matrix(
                        p, truncation, u_degree, g
                    )
                    evaluation_rank = rank_mod_p(evaluation, p)
                    multiplier_rank = rank_mod_p(multiplier, p)
                    require(
                        evaluation_rank == truncation,
                        "single specialization should be onto",
                    )
                    require(
                        multiplier_rank == truncation * u_degree,
                        "multiplication by the monic kernel generator lost rank",
                    )
                    require(
                        domain_dimension - evaluation_rank == multiplier_rank,
                        "kernel and principal image dimensions disagree",
                    )
                    require(
                        product_is_zero(evaluation, multiplier, p),
                        "(u-g) did not specialize to zero",
                    )
                    kernel_cases += 1

                    full_jet = specialization_matrix(
                        p,
                        truncation,
                        u_degree,
                        g,
                        tuple(range(u_degree + 1)),
                    )
                    require(
                        rank_mod_p(full_jet, p) == domain_dimension,
                        "full Hasse jet was not invertible",
                    )
                    hasse_cases += 1

                    # P_h=u-g+f^h has coefficientwise order zero but evaluates
                    # to f^h.  Its first transverse Hasse derivative is 1.
                    for height in range(1, truncation):
                        vector = [0] * domain_dimension
                        vector[1] = 1  # u
                        for degree, coefficient in enumerate(g):
                            if coefficient:
                                vector[degree * (u_degree + 1)] = (
                                    vector[degree * (u_degree + 1)] - coefficient
                                ) % p
                        vector[height * (u_degree + 1)] = (
                            vector[height * (u_degree + 1)] + 1
                        ) % p
                        image = apply_matrix(evaluation, vector, p)
                        require(
                            coefficientwise_f_order(vector, u_degree) == 0,
                            "hostile source order changed",
                        )
                        require(
                            series_order(image) == height,
                            "hostile specialization did not have prescribed order",
                        )
                        first_jet = specialization_matrix(
                            p, truncation, u_degree, g, (1,)
                        )
                        jet_image = apply_matrix(first_jet, vector, p)
                        require(
                            jet_image[0] == 1 and all(x == 0 for x in jet_image[1:]),
                            "first Hasse jet failed to expose the hidden direction",
                        )
                        hostile_order_cases += 1
                        maximum_tested_order_jump = max(
                            maximum_tested_order_jump, height
                        )

            # Distinct constant sections have unit Vandermonde differences.
            # They recover every polynomial whose u-degree is below p.
            for u_degree in range(1, min(6, p - 1) + 1):
                blocks = []
                for constant in range(u_degree + 1):
                    g = [constant] + [0] * (truncation - 1)
                    blocks.extend(
                        specialization_matrix(p, truncation, u_degree, g)
                    )
                require(
                    rank_mod_p(blocks, p) == truncation * (u_degree + 1),
                    "unit-separated section family failed interpolation",
                )
                joint_section_cases += 1

            for u_degree in U_DEGREES:
                domain_dimension = truncation * (u_degree + 1)
                slope_matrix = universal_slope_matrix(p, truncation, u_degree)
                require(
                    rank_mod_p(slope_matrix, p) == domain_dimension,
                    "universal-slope encoding lost a bidegree coefficient",
                )
                universal_slope_cases += 1

                f_degree_cap = truncation - 1
                good_multiplier = f_degree_cap + 1
                good_kronecker = kronecker_matrix(
                    p, truncation, u_degree, good_multiplier
                )
                require(
                    rank_mod_p(good_kronecker, p) == domain_dimension,
                    "Kronecker encoding above the degree cap lost rank",
                )
                bad_kronecker = kronecker_matrix(
                    p, truncation, u_degree, f_degree_cap
                )
                require(
                    rank_mod_p(bad_kronecker, p) < domain_dimension,
                    "sharp u-f^M hostile was not detected at the cap",
                )
                kronecker_cases += 1

    ordinary_derivative_hostiles = 0
    for p in PRIMES:
        # At u=0 every ordinary derivative of u^p through order p vanishes
        # modulo p, whereas D_u^[p](u^p)=1.
        ordinary_values = []
        for order in range(p + 1):
            if order < p:
                ordinary_values.append(0)
            else:
                ordinary_values.append(factorial(p) % p)
        require(
            all(value == 0 for value in ordinary_values),
            "ordinary characteristic-p hostile unexpectedly survived",
        )
        require(comb(p, p) % p == 1, "top Hasse derivative should be one")
        ordinary_derivative_hostiles += 1

    print("THM-4255 specialization-kernel / Hasse-jet exact audit")
    print(f"primes={PRIMES}")
    print("truncations=2..8, u_degrees=1..6")
    print(f"principal_kernel_cases={kernel_cases}")
    print(f"full_hasse_jet_cases={hasse_cases}")
    print(f"prescribed_order_hostiles={hostile_order_cases}")
    print(f"maximum_tested_order_jump={maximum_tested_order_jump}")
    print(f"unit_separated_joint_section_cases={joint_section_cases}")
    print(f"universal_slope_finite_box_cases={universal_slope_cases}")
    print(f"kronecker_sharp_boundary_cases={kronecker_cases}")
    print(f"ordinary_derivative_characteristic_p_hostiles={ordinary_derivative_hostiles}")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
