#!/usr/bin/env python3
"""Independent exact controls; no imports from the producer packet.

The all-order proof is audited in REVIEW.md. These finite checks do not
substitute for it. Failures remain enabled under python -O.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from math import isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRODUCER = HERE.parent / "session13-nonlinear-degree"
CHECKS = 0


def require(condition, *context):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(context)


def bits(number, length):
    return tuple((number >> index) & 1 for index in range(length))


def integer(word):
    return sum(value << index for index, value in enumerate(word))


def rotate(word):
    return word[-1:] + word[:-1]


def twisted(word):
    return (1 - word[-1],) + word[:-1]


def orbit(word, action):
    result = [word]
    new = action(word)
    while new != word:
        require(new not in result, "non-permutation action")
        result.append(new)
        new = action(new)
    return result


def partition(length, action):
    unseen = set(range(1 << length))
    classes = []
    while unseen:
        root = min(unseen)
        numbers = [integer(word) for word in orbit(bits(root, length), action)]
        unseen.difference_update(numbers)
        classes.append(numbers)
    return classes


def primitive_root(word):
    for length in range(1, len(word) + 1):
        if len(word) % length == 0 and word == word[:length] * (len(word) // length):
            return word[:length]
    raise RuntimeError("primitive root missing")


def allowed(word):
    return sum(primitive_root(word)) % 2 == 0


def quotient(word):
    return tuple(value ^ word[-1] for value in word[:-1])


def polynomial(word, terms):
    q = quotient(word)
    return sum(all(q[index] for index in support) for support in terms) % 2


def vector(word, terms):
    # Recover inverse shifts by indexing a complete FORWARD orbit.
    forward = orbit(word, twisted)
    return tuple(polynomial(forward[(-index) % len(forward)], terms)
                 for index in range(len(word)))


def anf_via_subsets(truth):
    """Coefficient S is the XOR over all truth values at subsets of S."""
    result = []
    for support in range(len(truth)):
        value = 0
        subset = support
        while True:
            value ^= truth[subset]
            if subset == 0:
                break
            subset = (subset - 1) & support
        result.append(value)
    return result


def canonical(word, action):
    return min(map(integer, orbit(word, action)))


def witnesses():
    examples = {
        2: (),
        4: ((1,), (0, 1)),
        6: ((0,), (1,), (0, 2), (0, 3)),
        8: ((0,), (0, 4), (1, 5), (3, 4), (4, 5),
            (1, 2, 3), (3, 5, 6)),
    }
    report = []
    for length, terms in examples.items():
        table = [vector(bits(value, length), terms) for value in range(1 << length)]
        fibers = Counter(map(integer, table))
        image_types = {}
        for source, target in enumerate(table):
            word = bits(source, length)
            require(allowed(target), "target landing", length, source)
            require(table[integer(twisted(word))] == rotate(target),
                    "equivariance", length, source)
            source_type = canonical(word, twisted)
            target_type = canonical(target, rotate)
            require(image_types.get(source_type, target_type) == target_type,
                    "well-defined type map")
            image_types[source_type] = target_type
        target_words = {value for value in range(1 << length)
                        if allowed(bits(value, length))}
        target_types = {canonical(bits(value, length), rotate) for value in target_words}
        require(set(fibers) == target_words, "labelled surjectivity", length)
        require(len(image_types) == len(set(image_types.values())) == len(target_types),
                "type bijection", length)
        full_scalar = [word[0] for word in table]
        small_scalar = full_scalar[:1 << (length - 1)]
        full_anf = anf_via_subsets(full_scalar)
        small_anf = anf_via_subsets(small_scalar)
        full_degree = max([mask.bit_count() for mask, value in enumerate(full_anf) if value] or [0])
        small_degree = max([mask.bit_count() for mask, value in enumerate(small_anf) if value] or [0])
        require(full_degree == small_degree == max(map(len, terms), default=0), "degree")
        require({mask for mask, value in enumerate(small_anf) if value}
                == {sum(1 << index for index in term) for term in terms}, "ANF")
        reps = sorted(image_types)
        report.append({
            "m": length, "degree": full_degree,
            "source_representatives": reps,
            "images": [integer(table[value]) for value in reps],
            "target_type_representatives": sorted(target_types),
            "target_labelled": len(target_words),
            "fiber_size_histogram": dict(sorted(Counter(fibers.values()).items())),
            "full_scalar_weight": sum(full_scalar),
            "full_scalar_Walsh": sum(1 - 2 * bit for bit in full_scalar),
            "quotient_scalar_weight": sum(small_scalar),
            "quotient_scalar_Walsh": sum(1 - 2 * bit for bit in small_scalar),
        })
    return report


def dyadic_orbits():
    report = []
    for length in (1, 2, 4, 8, 16):
        sources = partition(length, twisted)
        targets = [row for row in partition(length, rotate)
                   if allowed(bits(row[0], length))]
        expected = (1 << length) // (2 * length)
        require(all(len(row) == 2 * length for row in sources), "source free", length)
        require(len(sources) == len(targets) == expected, "type count", length)
        mass = sum(row[0].bit_count() for row in targets)
        require(2 * mass == length * (expected - 1), "target weight sum", length)
        target_ids = {min(row) for row in targets}
        for row in targets:
            root = row[0]
            if root:
                complement = canonical(tuple(1 - bit for bit in bits(root, length)), rotate)
                require(complement in target_ids, "complement closure", length, root)
                require(root.bit_count() + complement.bit_count() == length,
                        "complement weight", length, root)
        report.append({"m": length, "source_types": len(sources),
                       "target_types": len(targets), "source_orbit_length": 2 * length,
                       "target_weight_mass": mass})
    return report


def fixed_subsets():
    report = []
    # Some non-dyadic e also exercise the stated all-length restriction lemma.
    for small, repeat in ((1, 3), (2, 3), (3, 3), (4, 3)):
        large = small * repeat
        embedded_source = {}
        embedded_target = {}
        for value in range(1 << small):
            word = bits(value, small)
            lift_source = tuple(word[index % small] ^ ((index // small) % 2)
                                for index in range(large))
            lift_target = word * repeat
            embedded_source[value] = integer(lift_source)
            embedded_target[value] = integer(lift_target)
            next_word = twisted(word)
            next_lift = tuple(next_word[index % small] ^ ((index // small) % 2)
                              for index in range(large))
            require(twisted(lift_source) == next_lift, "source intertwining")
            require(rotate(lift_target) == rotate(word) * repeat, "target intertwining")
            require(allowed(word) == allowed(lift_target), "primitive predicate")
        fixed = set()
        for value in range(1 << large):
            word = bits(value, large)
            full = word + tuple(1 - bit for bit in word)
            rotated = full[-2 * small:] + full[:-2 * small]
            if rotated == full:
                fixed.add(value)
        require(fixed == set(embedded_source.values()), "exact source fixed set", small, repeat)
        small_sources = partition(small, twisted)
        small_targets = [row for row in partition(small, rotate)
                         if allowed(bits(row[0], small))]
        source_image_types = {canonical(bits(embedded_source[row[0]], large), twisted)
                              for row in small_sources}
        target_image_types = {canonical(bits(embedded_target[row[0]], large), rotate)
                              for row in small_targets}
        require(len(source_image_types) == len(small_sources), "source orbit reflection")
        require(len(target_image_types) == len(small_targets), "target orbit reflection")
        require(len(small_sources) == len(small_targets), "finite equal sizes")
        report.append({"e": small, "q": repeat, "m": large,
                       "fixed_source_words": len(fixed), "types": len(small_sources)})
    # Oddness is essential: at m=2e, T^(2e) is complementation, hence has no fixed point.
    for small in (1, 2, 3, 4):
        large = 2 * small
        require(all(tuple(1 - bit for bit in bits(value, large)) != bits(value, large)
                    for value in range(1 << large)), "even-repeat hostile")
    return report


def walsh_all_functions():
    report = []
    for variables in range(1, 5):
        # Convert ANF coefficient to truth-set parity through explicit subsets.
        parity_masks = []
        for support in range(1 << variables):
            subset = support
            mask = 0
            while True:
                mask |= 1 << subset
                if subset == 0:
                    break
                subset = (subset - 1) & support
            parity_masks.append(mask)
        positive_degree = 0
        for truth in range(1 << (1 << variables)):
            degree = max([support.bit_count()
                          for support, mask in enumerate(parity_masks)
                          if (truth & mask).bit_count() % 2] or [0])
            if degree:
                positive_degree += 1
                walsh = (1 << variables) - 2 * truth.bit_count()
                require(walsh % (1 << ((variables + degree - 1) // degree)) == 0,
                        "divisibility", variables, truth, degree)
        report.append({"variables": variables,
                       "all_functions": 1 << (1 << variables),
                       "positive_degree_functions": positive_degree})
    return report


def quadratic_input():
    length = 8
    source_reps = [min(row) for row in partition(length, twisted)]
    target_reps = [min(row) for row in partition(length, rotate)
                   if allowed(bits(row[0], length))]
    supports = [support for degree in range(3)
                for support in combinations(range(length - 1), degree)]
    columns = []
    for support in supports:
        words = [integer(vector(bits(value, length), (support,))) for value in source_reps]
        columns.append(sum(word << (length * index) for index, word in enumerate(words)))
    labels = []
    for value in range(1 << length):
        word = bits(value, length)
        labels.append(target_reps.index(canonical(word, rotate)) if allowed(word) else -1)
    text = f"{len(columns)} {len(source_reps)}\n"
    text += "".join(f"{column % (1 << 64)} {column >> 64}\n" for column in columns)
    text += " ".join(map(str, labels)) + "\n"
    original = (PRODUCER / "input8.txt").read_bytes()
    require(text.encode() == original, "complete independent basis byte match")
    (HERE / "input8.independent.txt").write_text(text)
    return {"dimension": len(columns), "all_coefficient_vectors": 1 << len(columns),
            "independent_input_sha256": sha256(text.encode()).hexdigest(),
            "producer_input_byte_identical": True}


def elementary_examples():
    report = []
    for prime in (17, 97, 193, 641, 257, 65537):
        require(all(prime % factor for factor in range(2, 1 + isqrt(prime))),
                "prime example", prime)
        length = (prime - 1) // 2
        dyadic = length & -length
        exponent = dyadic.bit_length() - 1
        report.append({"prime": prime, "m": length, "e": dyadic,
                       "a": exponent, "lower_bound": (dyadic + exponent - 2) // exponent})
    # Certified small divisors of F_n, including composite F_5 and F_6.
    fermat = []
    for n, prime in ((2, 17), (3, 257), (4, 65537), (5, 641), (6, 274177)):
        require(all(prime % factor for factor in range(2, 1 + isqrt(prime))),
                "Fermat divisor prime")
        require(pow(2, 1 << n, prime) == prime - 1, "Fermat divisibility", n, prime)
        require(pow(2, 1 << (n + 1), prime) == 1, "order upper", n, prime)
        require(all(pow(2, 1 << k, prime) != 1 for k in range(n + 1)),
                "order exact", n, prime)
        exponent = ((prime - 1) // 2 & -((prime - 1) // 2)).bit_length() - 1
        require(exponent >= n, "dyadic growth", n, prime)
        fermat.append({"n": n, "prime_divisor": prime,
                       "order_of_2": 1 << (n + 1), "v2_m": exponent})
    return {"table": report, "Fermat_divisor_controls": fermat}


def small_census():
    report = []
    for length in (2, 3, 4, 6):
        source_reps = [min(row) for row in partition(length, twisted)]
        target_reps = [min(row) for row in partition(length, rotate)
                       if allowed(bits(row[0], length))]
        target_lookup = [
            target_reps.index(canonical(bits(value, length), rotate))
            if allowed(bits(value, length)) else -1
            for value in range(1 << length)
        ]
        for degree in (1, 2):
            supports = [support for order in range(degree + 1)
                        for support in combinations(range(length - 1), order)]
            basis = []
            for support in supports:
                words = [integer(vector(bits(value, length), (support,)))
                         for value in source_reps]
                basis.append(sum(word << (length * index)
                                 for index, word in enumerate(words)))
            target_valued = 0
            bijections = 0
            for coefficients in range(1 << len(supports)):
                packed = 0
                for index, column in enumerate(basis):
                    if (coefficients >> index) & 1:
                        packed ^= column
                labels = [target_lookup[(packed >> (length * index)) & ((1 << length) - 1)]
                          for index in range(len(source_reps))]
                if -1 not in labels:
                    target_valued += 1
                    bijections += len(set(labels)) == len(target_reps)
            report.append({"m": length, "degree_bound": degree,
                           "maps": 1 << len(supports), "target_valued": target_valued,
                           "type_bijections": bijections})
    original = json.loads((PRODUCER / "small_enumeration.audit.json").read_text())
    require(report == original, "independent complete small census")
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=str(HERE / "independent_checks.json"))
    args = parser.parse_args()
    result = {
        "status": "FINITE-EXACT independent controls, not the all-order proof",
        "witnesses": witnesses(),
        "dyadic_orbits": dyadic_orbits(),
        "fixed_subsets": fixed_subsets(),
        "walsh_all_functions": walsh_all_functions(),
        "quadratic_basis": quadratic_input(),
        "examples": elementary_examples(),
        "small_census": small_census(),
    }
    result["explicit_checks"] = CHECKS
    result["script_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    Path(args.output).write_text(serialized)
    print(serialized, end="")


if __name__ == "__main__":
    main()
