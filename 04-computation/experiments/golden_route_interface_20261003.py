"""Cross-check literal phi digits, exact residue guards, and timed route words.

The universe is explicit in main. Independent ordinary iteration provides
finite complete suffixes; no claim of complete infinite-family coverage.
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from itertools import product


def load(name):
    path = Path(__file__).with_name(name + ".py")
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def need(condition, reason):
    if not condition:
        raise ValueError(reason)


def ordinary(n):
    return 3*n+1 if n % 2 else n//2


def suffix_to_one(n, bound=10000):
    bits = []
    while n != 1:
        need(len(bits) < bound, "declared finite suffix search exhausted")
        bits.append(n % 2)
        n = ordinary(n)
    return bits


def main():
    digit = load("golden_digit_carry_20261003")
    clock = load("golden_prime_clocks_20261003")
    route = load("collatz_mixed_return_compiler_20261003")
    pair_add, pair_mul = digit.pair_add, digit.pair_mul

    def parity_value(bits):
        return digit.evaluate_laurent({-j-1: bit for j, bit in enumerate(bits) if bit})

    def twice_theta_home(bits):
        # Theta(1)=phi/2 for the ordinary 1,4,2 cycle.
        return pair_add(pair_mul((2, 0), parity_value(bits)),
                        pair_mul((0, 1), digit.phi_power(-len(bits))))

    need(pair_mul((0, 1), pair_add((1, 0), pair_mul((-1, 0), digit.phi_power(-3))))
         == pair_mul((2, 0), digit.phi_power(-1)), "Theta(1)=phi/2")

    instances = complete_suffixes = numerator_checks = 0
    max_source = max_suffix_time = max_radix = 0
    first_examples = {}
    stored_edges = {}
    raw_odd_edges = max_tournament_vertices = 0
    for head, m in product(((1,), route.HEAD), range(1, 7)):
        row = route.compile_family(head, route.CYCLE, m)
        for lift in (0, 1, 2):
            _, source, _ = route.instantiate(row, lift)
            verified = route.verify_instance(row, lift)
            source_word = digit.canonical_phi_integer(source)
            source_poly, radix = digit.parse_phi_word(source_word)
            modulus = 2**row["K"] * 3**row["p"]
            period = clock.guard_clock(row["K"], row["p"])
            need(digit.read_phi_word(source_word) == (source, 0), "literal source exact value")
            need(digit.evaluate_laurent(source_poly, modulus) ==
                 digit.read_phi_word(source_word, modulus) == (source % modulus, 0),
                 "literal source modular reader")
            encoded = digit.read_phi_word(source_word, 2**row["K"])
            need(encoded == (row["residue"], 0), "source word passes exact dyadic guard")
            hostile_word = digit.canonical_phi_integer(source+2)
            need(digit.read_phi_word(hostile_word, 2**row["K"]) !=
                 (row["residue"], 0), "neighbor guard rejection")

            # The residue reader needs only the radix phase, but exact source
            # arithmetic still needs the full radix and the complete word.
            unpointed = digit.read_phi_word(source_word.replace(".", ""), modulus)
            phase = pair_mul(unpointed, clock.phi_power(-(radix % period), modulus), modulus)
            need(phase == (source % modulus, 0), "finite radix clock")

            # Independently replay ordinary steps and evaluate their positions.
            n = source
            bits = []
            for _ in range(verified["ordinary_time"]):
                bits.append(n % 2)
                n = ordinary(n)
            need(n == verified["endpoint"], "independent ordinary endpoint")
            need(parity_value(bits) == tuple(verified["golden_constant"]), "time polynomial")
            need(digit.phi_power(-len(bits)) == tuple(verified["golden_tail_multiplier"]),
                 "retained tail time")
            need("11" not in "".join(map(str, bits)), "ordinary parity admissibility")

            # A full route is only asserted for these finite endpoint tests.
            suffix = suffix_to_one(n)
            joined = suffix_to_one(source)
            need(joined == bits+suffix, "first-hit suffix splice")
            composed = pair_add(pair_mul((2, 0), parity_value(bits)),
                                pair_mul(digit.phi_power(-len(bits)), twice_theta_home(suffix)))
            need(composed == twice_theta_home(joined), "exact golden suffix composition")
            complete_suffixes += 1

            # A complete first-hit route also has an intrinsic SCC-word encoding:
            # reverse sibling index j becomes a cyclic tournament block R_(2j+3).
            # Keep only this established block-size representation, not an
            # invented orientation of symmetric product or color data.
            current, indices = source, []
            while current != 1:
                target, a = route.odd_step(current)
                base_a = 2 if target % 3 == 1 else 1
                need(target % 3 != 0 and a >= base_a and (a-base_a) % 2 == 0,
                     "compressed route phase legality")
                indices.append((a-base_a)//2)
                need(current not in stored_edges or stored_edges[current] == (target, a),
                     "shared suffix edge consistency")
                stored_edges[current] = target, a
                raw_odd_edges += 1
                current = target
            block_sizes = [2*j+3 for j in indices]+[1]
            recovered = 1
            for block in reversed(block_sizes[:-1]):
                j = (block-3)//2
                a = (2 if recovered % 3 == 1 else 1)+2*j
                numerator = 2**a*recovered-1
                need(numerator % 3 == 0, "inverse block integrality")
                recovered = numerator//3
                need(recovered > 1, "first-hit root padding excluded")
            need(recovered == source, "tournament SCC-size route recovery")
            max_tournament_vertices = max(max_tournament_vertices, sum(block_sizes))

            # Check the two-shift representation of 3n+1 at the same sources.
            raw = digit.polynomial_add(digit.convolve(source_poly, {2: 1, -2: 1}), {0: 1})
            normal, _ = digit.parse_phi_word(digit.canonical_phi_integer(3*source+1))
            R, Q = digit.carry_certificate(raw, normal)
            need(digit.recover_raw(normal, R, Q) == raw, "numerator carry recovery")
            need(digit.evaluate_laurent(normal) == (3*source+1, 0), "numerator value")
            numerator_checks += 1
            if head not in first_examples:
                first_examples[head] = dict(source=source, source_word=source_word, guard=(row["residue"], row["K"]),
                             endpoint=n, prefix_time=len(bits), suffix_time=len(suffix),
                             total_time=len(joined), radix=radix, modulus=modulus, period=period)
            instances += 1
            max_source = max(max_source, source)
            max_suffix_time = max(max_suffix_time, len(suffix))
            max_radix = max(max_radix, radix)

    # The evaluated prefix polynomial alone loses timing and endpoint data.
    need(parity_value([1]) == parity_value([1, 0]) == digit.phi_power(-1), "equal prefix values")
    need(ordinary(3) == 10 and ordinary(ordinary(3)) == 5, "different timed endpoints")
    need(digit.phi_power(-1) != digit.phi_power(-2), "tail multipliers differ")
    # Modulo 2, three nonzero phi states really cycle; zero is a fourth state.
    cycle = [digit.to_zeck_registers(clock.phi_power(j, 2)) for j in range(3)]
    cycle = [tuple(x % 2 for x in state) for state in cycle]
    need(cycle == [(1, 0), (0, 1), (1, 1)], "three nonzero auxiliary colors")
    need(digit.to_zeck_registers((0, 0)) == (0, 0), "neutral color retained")

    # Group the actual positional exponents by their phi-mod-2 clock phase.
    signature_checks = 0
    signatures = {}
    for n in range(1001):
        poly, _ = digit.parse_phi_word(digit.canonical_phi_integer(n))
        counts = tuple(sum(c for k, c in poly.items() if k % 3 == j) % 2 for j in range(3))
        need(counts[1] == counts[2] and (counts[0]+counts[2]) % 2 == n % 2,
             "integer support phase constraint")
        signatures[n] = counts
        signature_checks += 1
    need(signatures[3] == signatures[9] == (0, 1, 1), "phase signature primality hostile")
    need(signatures[5] == signatures[105] == (1, 0, 0), "second phase signature hostile")

    # Fourier characters are checked by their exact exponents, never by floats.
    # D_d(v)=M v+d e1 pulls frequency k back to M^T k plus phase k1*d.
    fourier_checks = 0
    for modulus in range(2, 12):
        for A, B, r, s in product(range(modulus), repeat=4):
            for d in (0, 1):
                direct = (r*(B+d)+s*(A+B)) % modulus
                dual = (s*A+(r+s)*B+r*d) % modulus
                need(direct == dual, "Fourier dual digit rule")
                fourier_checks += 1
    # Reflected states have identical cosine data at every mode, but digit1
    # destroys that ambiguity. Mode (1,1), modulus3, distinguishes the outputs.
    v, negv = (1, 0), (2, 0)
    need(all((r*v[0]+s*v[1]+r*negv[0]+s*negv[1]) % 3 == 0
             for r, s in product(range(3), repeat=2)), "all initial cosines agree")
    after_v, after_negv = (v[1]+1, v[0]+v[1]), (negv[1]+1, negv[0]+negv[1])
    phase_v, phase_negv = sum(after_v) % 3, sum(after_negv) % 3
    need(phase_v != phase_negv and phase_v != (-phase_negv) % 3,
         "cosine-only digit update fails")

    # Decorated arithmetic retains the aggregate raw Laurent polynomial.
    # A state (n,q) decodes to sigma(n)+f*q, f=t^2-t-1.
    f = {2: 1, 1: -1, 0: -1}
    sigma = lambda n: digit.parse_phi_word(digit.canonical_phi_integer(n))[0]
    add, mul = digit.polynomial_add, digit.convolve
    decorated_checks = 0
    for a, b in product(range(1, 9), repeat=2):
        sa, sb = sigma(a), sigma(b)
        qa, qb = {-(a % 3): (-1)**a}, {-(b % 4): (-1)**b, 5: 2}
        raw_a, raw_b = add(sa, mul(f, qa)), add(sb, mul(f, qb))
        acarry = digit.carry_laurent(add(sa, sb), sigma(a+b))
        ccarry = digit.carry_laurent(mul(sa, sb), sigma(a*b))
        summed_q = add(acarry, add(qa, qb))
        product_q = add(ccarry, add(mul(qa, sb), add(mul(sa, qb), mul(f, mul(qa, qb)))))
        need(add(sigma(a+b), mul(f, summed_q)) == add(raw_a, raw_b), "decorated addition")
        need(add(sigma(a*b), mul(f, product_q)) == mul(raw_a, raw_b), "decorated multiplication")
        decorated_checks += 2

    # Concurrent denominator-module bridge, checked by cleared denominators.
    phi5_minus1 = pair_add(digit.phi_power(5), (-1, 0))
    need(phi5_minus1 == pair_mul(pair_mul((-1, 0), digit.phi_power(3)), (-4, 1)),
         "monodromy ideal generator")
    need(pair_mul((-4, 1), (-1, 7)) == (11, -22), "minus5 denominator ideal")
    need(pair_mul((-4, 1), (1, 4)) == (0, -11), "rational1/13 shared denominator ideal")
    need(clock.pair_norm((-1, 7)) == -55 and clock.pair_norm((1, 4)) == -11,
         "primitive numerator is not a unit modulo11")

    print("FINITE-EXACT cross-interface universe: heads(1),(1,2); m=1..6; lifts=0,1,2")
    print(f"Sources / neighbor guard rejections / clock phases: {instances} each")
    print(f"Exact two-shift numerator carry certificates: {numerator_checks}")
    print(f"Complete suffixes independently found and spliced: {complete_suffixes}; search bound=10000")
    print(f"Maximum tested source={max_source}; maximum radix={max_radix}; longest suffix={max_suffix_time}")
    print("First interfaces:", first_examples)
    print(f"Complete SCC-size route encodings: {instances}; largest tournament={max_tournament_vertices} vertices")
    print(f"Exact shared-suffix storage: {raw_odd_edges} route-edge occurrences -> {len(stored_edges)} distinct odd edges")
    print("Auxiliary mod2 cycle:", cycle, "; zero remains separate; no identification with ordered R/K/B grammar")
    print(f"Exponent-phase signature constraints: {signature_checks} integers; primes3,5 share signatures with composites9,105.")
    print(f"Exact Fourier-dual digit checks: {fourier_checks}; all states/modes moduli2..11, digits0,1")
    print("Cosine hostile: mod3 states(1,0),(-1,0) share every cosine; after digit1 mode(1,1) has phases2 and0.")
    print(f"Decorated carry-memory arithmetic: {decorated_checks} exact sum/product controls, a,b=1..8")
    print("Concurrent ideal bridge: (Phi5(phi))=(phi^5-1)=(phi-4); -5 and rational1/13 coordinates share it.")
    print("Timing hostile: source3 prefixes [1] and [1,0] both evaluate to phi^-1, but endpoints are10 and5.")
    print("PROVED interface is value/guard/timed-certificate preservation; FINITE-EXACT complete routes only for tested sources.")
    print("PASS: all checks active with and without python -O")


if __name__ == "__main__":
    main()
