#!/usr/bin/env python3
"""Exact companion for THM-3353 split-prime parabolic transplants.

Universe
--------
* every rational prime p == 1 (mod 4) below 5,000;
* both determinant-matched cusp gates A and D;
* the first 24 lawful rows on each compiled arithmetic branch;
* direct 2x2 parameter words and an independent 3x3 Berggren evaluation;
* every prime p == 3 (mod 4) below 500 as a no-root hostile;
* five CRT rank controls for each of p=5,13,17.

All checks are exact integer arithmetic.  Explicit exceptions make normal and
optimized Python execute the same validation path.
"""

from math import gcd, isqrt


U2 = ((2, -1), (1, 0))
A2 = ((2, 1), (1, 0))
D2 = ((1, 2), (0, 1))
PARAMETER_MATRICES = {"U": U2, "A": A2, "D": D2}

U3 = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
A3 = ((1, 2, 2), (2, 1, 2), (2, 2, 3))
D3 = ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3))
BERGGREN_MATRICES = {"U": U3, "A": A3, "D": D3}

I2 = ((1, 0), (0, 1))
SWAP = ((0, 1), (1, 0))
ROOT_PARAMETERS = (2, 1)
ROOT_TRIPLE = (3, 4, 5)
PRIME_LIMIT = 5000
ROW_COUNT = 24


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def matrix_vector(matrix, vector):
    return tuple(
        sum(matrix[i][j] * vector[j] for j in range(len(vector)))
        for i in range(len(matrix))
    )


def matrix_multiply(left, right):
    rows = len(left)
    middle = len(right)
    cols = len(right[0])
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(middle))
              for j in range(cols))
        for i in range(rows)
    )


def identity_matrix(size):
    return tuple(tuple(int(i == j) for j in range(size)) for i in range(size))


def matrix_power(matrix, exponent):
    require(exponent >= 0, "negative matrix exponent")
    result = identity_matrix(len(matrix))
    base = matrix
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = matrix_multiply(result, base)
        base = matrix_multiply(base, base)
        remaining //= 2
    return result


def transpose(matrix):
    return tuple(tuple(matrix[j][i] for j in range(len(matrix)))
                 for i in range(len(matrix[0])))


def determinant2(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def word_matrix(word):
    matrix = I2
    for letter in word:
        matrix = matrix_multiply(PARAMETER_MATRICES[letter], matrix)
    return matrix


def word_matrix3(word):
    matrix = identity_matrix(3)
    for letter in word:
        matrix = matrix_multiply(BERGGREN_MATRICES[letter], matrix)
    return matrix


def apply_parameter_word(word, vector=ROOT_PARAMETERS):
    for letter in word:
        vector = matrix_vector(PARAMETER_MATRICES[letter], vector)
    return vector


def apply_berggren_word(word, vector=ROOT_TRIPLE):
    for letter in word:
        vector = matrix_vector(BERGGREN_MATRICES[letter], vector)
    return vector


def euclid_triple(parameters):
    m, n = parameters
    return (m * m - n * n, 2 * m * n, m * m + n * n)


def norm2(vector):
    return vector[0] * vector[0] + vector[1] * vector[1]


def is_prime(number):
    if number < 2:
        return False
    if number % 2 == 0:
        return number == 2
    divisor = 3
    while divisor * divisor <= number:
        if number % divisor == 0:
            return False
        divisor += 2
    return True


def split_representation(prime):
    for b in range(1, isqrt(prime) + 1):
        a2 = prime - b * b
        a = isqrt(a2)
        if a > b and a * a == a2:
            require(gcd(a, b) == 1 and (a - b) % 2 == 1,
                    f"nonprimitive split representation at p={prime}")
            return a, b
    raise RuntimeError(f"missing split representation at p={prime}")


def odd_cusp_tail(target):
    """Unique word R with R(3,1)=target for a coprime odd target."""
    x, y = target
    require(x > y > 0 and gcd(x, y) == 1 and x % 2 == y % 2 == 1,
            f"bad odd cusp target {target}")
    inverse_letters = []
    while (x, y) != (3, 1):
        old_sum = x + y
        if x < 2 * y:
            inverse_letters.append("U")
            x, y = y, 2 * y - x
        elif x < 3 * y:
            inverse_letters.append("A")
            x, y = y, x - 2 * y
        elif x > 3 * y:
            inverse_letters.append("D")
            x, y = x - 2 * y, y
        else:
            require(y == 1, f"nonprimitive ratio-three boundary {(x, y)}")
            x, y = (3, 1)
            continue
        require(x > y > 0 and x + y < old_sum,
                f"cusp descent failed at {(x, y)}")
        require(gcd(x, y) == 1 and x % 2 == y % 2 == 1,
                f"cusp descent lost type at {(x, y)}")
    return "".join(reversed(inverse_letters))


def ceil_div(numerator, denominator):
    return -((-numerator) // denominator)


def compile_prime(prime):
    a, b = split_representation(prime)
    d = a * a - b * b
    e = 2 * a * b
    orientation = 1 if d > e else -1
    base_h = ((d, e), (-orientation * e, orientation * d))
    target = matrix_vector(base_h, (1, 1))
    require(target == (d + e, abs(d - e)),
            f"positive cusp mismatch at p={prime}")
    require(matrix_multiply(transpose(base_h), base_h)
            == ((prime * prime, 0), (0, prime * prime)),
            f"orthogonality mismatch at p={prime}")
    tail = odd_cusp_tail(target)
    require(apply_parameter_word(tail, (3, 1)) == target,
            f"tail reconstruction failed at p={prime}")

    branches = []
    for gate in "AD":
        suffix = gate + tail
        suffix_matrix = word_matrix(suffix)
        determinant = determinant2(suffix_matrix)
        require(abs(determinant) == 1, f"nonunimodular suffix at p={prime}")
        require(matrix_vector(suffix_matrix, (1, 1)) == target,
                f"boundary collision failed at p={prime}, gate={gate}")

        toggle = base_h
        if determinant2(toggle) != prime * prime * determinant:
            toggle = matrix_multiply(toggle, SWAP)
        require(determinant2(toggle) == prime * prime * determinant,
                f"determinant gauge failed at p={prime}, gate={gate}")
        require(matrix_multiply(transpose(toggle), toggle)
                == ((prime * prime, 0), (0, prime * prime)),
                f"gauged orthogonality failed at p={prime}, gate={gate}")
        require(matrix_vector(toggle, (1, 1)) == target,
                f"gauged cusp mismatch at p={prime}, gate={gate}")

        second_column = (suffix_matrix[0][1], suffix_matrix[1][1])
        toggle_transpose_column = matrix_vector(transpose(toggle), second_column)
        u0 = tuple(prime - value for value in toggle_transpose_column)
        require(u0[0] == u0[1] + 1,
                f"intercept is not consecutive at p={prime}, gate={gate}: {u0}")
        residue_lift = u0[1]
        residue = residue_lift % (prime * prime)
        start = max(1, ceil_div(2 - residue_lift, prime * prime))
        branches.append({
            "gate": gate,
            "suffix": suffix,
            "suffix_matrix": suffix_matrix,
            "toggle": toggle,
            "u0": u0,
            "residue_lift": residue_lift,
            "residue": residue,
            "start": start,
        })

    left, right = branches
    require(left["residue_lift"] + right["residue_lift"] + 1
            == 2 * prime - prime * prime,
            f"complementary lift identity failed at p={prime}")
    require((left["residue"] + right["residue"] + 1) % prime == 0,
            f"local conjugation identity failed at p={prime}")
    require(left["residue"] % prime != right["residue"] % prime,
            f"two gates hit one local root at p={prime}")
    return {
        "prime": prime,
        "a": a,
        "b": b,
        "target": target,
        "tail": tail,
        "branches": branches,
    }


def p_valuation(number, prime):
    valuation = 0
    while number % prime == 0:
        number //= prime
        valuation += 1
    return valuation


def audit_branch(compilation, branch):
    prime = compilation["prime"]
    q = (prime - 1) // 2
    suffix = branch["suffix"]
    start = branch["start"]
    toggle = branch["toggle"]
    residue_lift = branch["residue_lift"]
    residue = branch["residue"]
    suffix_matrix2 = branch["suffix_matrix"]
    suffix_matrix3 = word_matrix3(suffix)
    prefix_matrix2 = matrix_power(D2, q)
    prefix_matrix3 = matrix_power(D3, q)

    require((2 * residue * residue + 2 * residue + 1) % prime == 0,
            f"compiled residue is not a root at p={prime}")
    require((2 * residue * residue + 2 * residue + 1) % (prime * prime) != 0,
            f"compiled residue accidentally Hensel-lifted at p={prime}")

    for s in range(start, start + ROW_COUNT):
        t = prime * prime * s + residue_lift
        require(t >= 2, f"unlawful U-spine index at p={prime}, s={s}")
        input_parameters = (t + 1, t)
        intermediate2 = matrix_vector(prefix_matrix2, ROOT_PARAMETERS)
        intermediate2 = matrix_vector(matrix_power(U2, s - 1), intermediate2)
        output_parameters = matrix_vector(suffix_matrix2, intermediate2)
        numerator = matrix_vector(toggle, input_parameters)
        require(numerator[0] % prime == numerator[1] % prime == 0,
                f"nonintegral rational toggle at p={prime}, s={s}")
        toggled_parameters = (numerator[0] // prime, numerator[1] // prime)
        require(output_parameters == toggled_parameters,
                f"word/toggle mismatch at p={prime}, s={s}")
        require(output_parameters[0] > output_parameters[1] > 0,
                f"output left Euclid chamber at p={prime}, s={s}")
        require(gcd(*output_parameters) == 1
                and (output_parameters[0] - output_parameters[1]) % 2 == 1,
                f"output lost primitivity at p={prime}, s={s}")
        require(norm2(input_parameters) == norm2(output_parameters),
                f"hypotenuse changed at p={prime}, s={s}")
        require(p_valuation(norm2(input_parameters), prime) == 1,
                f"toggle prime valuation is not one at p={prime}, s={s}")

        parameter_triple = euclid_triple(output_parameters)
        intermediate3 = matrix_vector(prefix_matrix3, ROOT_TRIPLE)
        intermediate3 = matrix_vector(matrix_power(U3, s - 1), intermediate3)
        direct_triple = matrix_vector(suffix_matrix3, intermediate3)
        require(parameter_triple == direct_triple,
                f"2x2/3x3 route mismatch at p={prime}, s={s}")

        input_depth = t - 1
        output_depth = q + (s - 1) + len(suffix)
        require(q > 0 and input_depth > 0,
                f"empty-prefix hostile failed at p={prime}, s={s}")
        expected_cost = ((prime * prime + 1) * s + residue_lift
                         + q + len(suffix) - 2)
        expected_jump = abs((prime * prime - 1) * s + residue_lift
                            - q - len(suffix))
        require(input_depth + output_depth == expected_cost,
                f"path-cost formula failed at p={prime}, s={s}")
        require(abs(input_depth - output_depth) == expected_jump,
                f"depth-jump formula failed at p={prime}, s={s}")


def roots_mod_prime(prime):
    return tuple(t for t in range(prime)
                 if (2 * t * t + 2 * t + 1) % prime == 0)


def crt_pair(residue, modulus, new_residue, new_modulus):
    step = ((new_residue - residue) * pow(modulus, -1, new_modulus)) % new_modulus
    return ((residue + modulus * step) % (modulus * new_modulus),
            modulus * new_modulus)


def crt_rank_control(compilation, branch, extra_primes):
    prime = compilation["prime"]
    residue = branch["residue"]
    modulus = prime * prime
    value = residue
    for extra in extra_primes:
        roots = roots_mod_prime(extra)
        require(len(roots) == 2, f"bad CRT split prime {extra}")
        value, modulus = crt_pair(value, modulus, roots[0], extra)
    if value < 2:
        value += modulus
    require(value % (prime * prime) == residue,
            f"CRT lost compiled branch at p={prime}")
    hypotenuse = 2 * value * (value + 1) + 1
    require(p_valuation(hypotenuse, prime) == 1,
            f"CRT changed toggle valuation at p={prime}")
    for extra in extra_primes:
        require(hypotenuse % extra == 0,
                f"CRT lost extra prime {extra} at p={prime}")
    return value, hypotenuse


def main():
    split_primes = [p for p in range(5, PRIME_LIMIT)
                    if p % 4 == 1 and is_prime(p)]
    inert_primes = [p for p in range(3, 500)
                    if p % 4 == 3 and is_prime(p)]
    compilations = []
    branch_rows = 0
    max_tail = (0, None)
    all_start_one = True

    for prime in split_primes:
        compilation = compile_prime(prime)
        compilations.append(compilation)
        if len(compilation["tail"]) > max_tail[0]:
            max_tail = (len(compilation["tail"]), prime)
        for branch in compilation["branches"]:
            audit_branch(compilation, branch)
            branch_rows += ROW_COUNT
            all_start_one = all_start_one and branch["start"] == 1

    for prime in inert_primes:
        require(roots_mod_prime(prime) == (),
                f"inert-prime hostile has a root at p={prime}")

    controls = []
    extras = (13, 17, 29, 37, 41, 53, 61)
    for prime in (5, 13, 17):
        compilation = next(item for item in compilations if item["prime"] == prime)
        available = tuple(q for q in extras if q != prime)
        branch = compilation["branches"][0]
        for rank in range(1, 6):
            t, hypotenuse = crt_rank_control(compilation, branch, available[:rank])
            controls.append((prime, rank, t, hypotenuse.bit_length()))

    p5 = next(item for item in compilations if item["prime"] == 5)
    p13 = next(item for item in compilations if item["prime"] == 13)
    require(p5["tail"] == "DD", "p=5 cusp tail changed")
    require(p5["branches"][0]["suffix"] == "ADD"
            and p5["branches"][0]["residue_lift"] == 1,
            "THM-3345 branch was not recovered")
    require(p13["tail"] == "AA", "p=13 cusp tail control changed")

    print("THM-3353 split-prime parabolic branch-transplant exact audit")
    print(f"split_prime_universe=p<5000,count={len(split_primes)}")
    print(f"branches={2 * len(split_primes)},rows_per_branch={ROW_COUNT},"
          f"exact_rows={branch_rows}")
    print(f"inert_hostiles=p<500,count={len(inert_primes)},all_rootless=True")
    print(f"all_sampled_branch_starts_equal_one={all_start_one}")
    print(f"max_cusp_tail_length={max_tail[0]},at_p={max_tail[1]}")
    for compilation in compilations[:6]:
        pieces = []
        for branch in compilation["branches"]:
            pieces.append(
                f"{branch['gate']}:r0={branch['residue_lift']},"
                f"r={branch['residue']},S={branch['suffix']}"
            )
        print(f"p={compilation['prime']}={compilation['a']}^2+"
              f"{compilation['b']}^2,cusp={compilation['target']},"
              f"tail={compilation['tail'] or '-'} | " + "; ".join(pieces))
    print("p5_known_branch=U^(25s)->DD U^(s-1) ADD")
    print("p5_companion_branch=U^(25s-18)->DD U^(s-1) DDD")
    print("p13_branches=D^6 U^(s-1) AAA ; D^6 U^(s-1) DAA")
    print("crt_controls(prime,extra_rank,t,hypotenuse_bits):")
    for row in controls:
        print(f"  {row}")
    print("normalization=all checks exact; compare stdout after LF normalization")
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
