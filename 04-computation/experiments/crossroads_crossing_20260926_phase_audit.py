"""Independent exact audit; imports none of the implementations under review."""

from fractions import Fraction as F
from math import lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def floor_log2(x):
    k = x.numerator.bit_length() - x.denominator.bit_length()
    if x < F(2)**k:
        k -= 1
    require(F(2)**k <= x < F(2)**(k+1), "floor logarithm")
    return k


def cuts(gains):
    return sorted([F(1), F(2)] + [F(2)**(floor_log2(g)+1)/g for g in gains])


def language(gains):
    words = set()
    boundaries = cuts(gains)
    for left, right in zip(boundaries, boundaries[1:]):
        t = (left+right)/2
        prefixes = [0] + [floor_log2(t*g) for g in gains]
        words.add(tuple(b-a for a, b in zip(prefixes, prefixes[1:])))
    return words


def phase_controls():
    thresholds = []
    perturbed_checks = saturation_checks = largest_bits = 0
    for length in range(1, 25):
        ideal = [F(3, 2)**j for j in range(1, length+1)]
        boundaries = cuts(ideal)
        rho = min(b/a for a, b in zip(boundaries, boundaries[1:]))
        minimum = 1
        while (1+F(1, 3*minimum))**length >= rho:
            minimum += 1
        thresholds.append(minimum)
        expected_language = language(ideal)
        require(len(expected_language) == length+1, "rotation complexity")
        for mode in range(3):
            gain = F(1)
            perturbed = []
            for j in range(length):
                n = minimum if mode == 0 else minimum+j if mode == 1 else 2*minimum+3*j
                gain *= F(3, 2)*(1+F(1, 3*n))
                perturbed.append(gain)
            require(language(perturbed) == expected_language, "perturbed language changed")
            perturbed_checks += 1

        # Independent source selection: round DOWN to -1 modulo 2^(L+1).
        found_words = set()
        for left, right in zip(boundaries, boundaries[1:]):
            t = (left+right)/2
            prefixes = [0] + [floor_log2(t*g) for g in ideal]
            expected = tuple(b-a for a, b in zip(prefixes, prefixes[1:]))
            exponent = 100
            modulus = 2**(length+1)
            while True:
                center = (t*2**exponent).__floor__()
                source = center-(center+1) % modulus
                current = source
                actual = []
                valuations = []
                for _ in range(length):
                    old_height = current.bit_length()-1
                    odd_half_step = (3*current+1)//2
                    actual.append(odd_half_step.bit_length()-1-old_height)
                    numerator = 3*current+1
                    a = (numerator & -numerator).bit_length()-1
                    valuations.append(a)
                    current = numerator >> a
                if tuple(actual) == expected and all(a == 1 for a in valuations):
                    break
                exponent *= 2
                require(exponent <= 1000, "saturation source search")
            require(source > 0 and current > source, "positive growing saturation")
            found_words.add(tuple(actual))
            largest_bits = max(largest_bits, source.bit_length())
            saturation_checks += 1
        require(found_words == expected_language, "unsaturated rotation language")
    print("independent_cut_stability_checks", perturbed_checks)
    print("height_thresholds_L1_to_L24", thresholds)
    print("actual_all_a1_saturation_witnesses", saturation_checks)
    print("largest_saturation_source_bits", largest_bits)


def multiplicative_order(a, p):
    residue, order = a % p, 1
    while residue != 1:
        residue = residue*a % p
        order += 1
    return order


def valuation(x, p):
    x = F(x)
    require(x != 0, "valuation at zero")

    def integer_valuation(z):
        z = abs(z)
        count = 0
        while z % p == 0:
            z //= p
            count += 1
        return count

    return integer_valuation(x.numerator)-integer_valuation(x.denominator)


def evaluate(coefficients, x):
    result = F(0)
    for coefficient in coefficients:
        result = result*x+coefficient
    return result


def crt(congruences):
    representative, modulus = 0, 1
    for residue, next_modulus in congruences:
        representative += modulus*((residue-representative)*pow(modulus, -1, next_modulus) % next_modulus)
        modulus *= next_modulus
    return representative, modulus


def polynomial_controls():
    # Descending coefficients; n+5 intentionally kills the first B=1
    # expanding negative fixed point r=-5, exercising root avoidance.
    polynomials = [(1, 0), (1, 1), (1, -1), (1, 5), (1, 0, 1),
                   (1, 0, -1), (1, 0, 0, 2), (3, 1), (9,)]
    total = 0
    for primes in ((2,), (2, 3), (2, 3, 5, 11), (2, 3, 7, 13)):
        zero_run = lcm(1, *[multiplicative_order(2, p) for p in primes if p >= 5])
        progression = lcm(1, *[multiplicative_order(base, p) for p in primes if p >= 5 for base in (2, 3)])
        odd_run = 1
        excluded_roots = 0
        while True:
            if 3**odd_run > 2**(odd_run+zero_run):
                fixed = -F(3**odd_run-2**odd_run, 3**odd_run-2**(odd_run+zero_run))
                if all(evaluate(coefficients, fixed) for coefficients in polynomials):
                    break
                excluded_roots += 1
            odd_run += progression
        if zero_run == 1:
            require(excluded_roots == 1 and odd_run == 3, "deliberate root exclusion")
        length = odd_run+zero_run
        multiplier = F(3**odd_run, 2**length)
        precision = {p: 1+max(valuation(evaluate(coefficients, fixed), p) for coefficients in polynomials) for p in primes}
        target_vector = tuple(valuation(evaluate(coefficients, fixed), p) for p in primes for coefficients in polynomials)
        largest_bits = 0
        for repetitions in (1, 2, 4, 8):
            congruences = []
            for p, e in precision.items():
                exponent = e+(repetitions*length if p == 2 else 0)
                modulus = p**exponent
                residue = fixed.numerator*pow(fixed.denominator, -1, modulus) % modulus
                congruences.append((residue, modulus))
            source, modulus = crt(congruences)
            required_source = 2**(repetitions*length)*10**6
            source += max(0, (required_source-source)//modulus+1)*modulus
            current = source
            minimum = source
            word = "1"*odd_run+"0"*zero_run
            for j in range(repetitions*length):
                require(current & 1 == int(word[j % length]), "actual parity shadow")
                current = (3*current+1)//2 if current & 1 else current//2
                minimum = min(minimum, current)
                require(all(evaluate(coefficients, current) for coefficients in polynomials), "intermediate polynomial root")
            require(F(current) == fixed+multiplier**repetitions*(source-fixed), "affine endpoint")
            require(F(current, source) > multiplier**repetitions, "expanding endpoint ratio")
            require(minimum > 10**6, "finite core avoidance")
            for value in (source, current):
                actual_vector = tuple(valuation(evaluate(coefficients, value), p) for p in primes for coefficients in polynomials)
                require(actual_vector == target_vector, "polynomial valuation vector changed")
            largest_bits = max(largest_bits, source.bit_length(), current.bit_length())
            total += 1
        print("polynomial_prime_set", primes, "B", zero_run, "L", progression,
              "a", odd_run, "precision", precision, "repetitions", [1, 2, 4, 8],
              "largest_endpoint_bits", largest_bits, "excluded_candidate_roots", excluded_roots, "PASS")
    print("polynomial_shadow_blocks", total, "polynomial_features", len(polynomials))


def main():
    phase_controls()
    polynomial_controls()
    print("PASS: independent integer/Fraction audit; no implementation under review imported")


if __name__ == "__main__":
    main()
