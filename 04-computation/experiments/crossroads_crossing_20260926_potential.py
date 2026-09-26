"""Exact controls for dyadic crossings; no sampled law is an orbit theorem.

Run: python 04-computation/experiments/crossroads_crossing_20260926_potential.py
All pass/fail controls use integer or Fraction arithmetic, including the
Sturmian language and its finite-height sufficient condition. Logarithms
are used only to print diagnostic crossing/carry identities.
"""
from fractions import Fraction
from math import lcm, log2


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def v2(n):
    return (n & -n).bit_length() - 1


def vp(n, prime):
    exponent = 0
    while n % prime == 0:
        n //= prime
        exponent += 1
    return exponent


def multiplicative_order(base, prime):
    require(prime >= 5, "order only requested for primes>=5")
    value = base % prime
    order = 1
    while value != 1:
        value = value*base % prime
        order += 1
        require(order <= prime-1, "non-prime or invalid multiplicative order")
    return order


def height(n):
    return n.bit_length() - 1


def odd_step(n, multiplier=3, shift=1):
    value = multiplier*n + shift
    a = v2(value)
    return value >> a, a


def upward(n):
    return height((3*n+1)//2) - height(n)


def floor_log2_fraction(q):
    k = q.numerator.bit_length() - q.denominator.bit_length()
    if k >= 0:
        return k - (q.numerator < (q.denominator << k))
    return k - ((q.numerator << -k) < q.denominator)


def sturmian_data(length):
    # t=2**phase. The cut at {-j*alpha}, alpha=log2(3/2), is
    # t=2**ceil(log2(3**j))/3**j. All comparisons are rational.
    cuts = [Fraction(1), Fraction(2)]
    power = 1
    for j in range(1, length+1):
        power *= 3
        cuts.append(Fraction(1 << power.bit_length(), power))
    cuts.sort()
    require(len(set(cuts)) == length+2, "cut collision")
    rho = min(b/a for a, b in zip(cuts, cuts[1:]))
    language = set()
    for left, right in zip(cuts, cuts[1:]):
        t = (left+right)/2
        totals = [floor_log2_fraction(t*Fraction(3, 2)**j)
                  for j in range(length+1)]
        language.add(tuple(b-a for a, b in zip(totals, totals[1:])))
    require(len(language) == length+1, "Sturmian complexity")
    require(all(set(word) <= {0, 1} for word in language), "alphabet")
    return language, rho


def all_odd_realizations(length):
    """One large exact source for EVERY length-L rotation word."""
    cuts = [Fraction(1), Fraction(2)]
    for j in range(1, length+1):
        power = 3**j
        cuts.append(Fraction(1 << power.bit_length(), power))
    cuts.sort()
    modulus = 1 << (length+1)
    found = set()
    for left, right in zip(cuts, cuts[1:]):
        target = (left+right)/2
        exponent = 5*length+20
        while True:
            center = (target.numerator << exponent)//target.denominator
            source = center + ((-1-center) % modulus)
            current = source
            word, valuations = [], []
            for _ in range(length):
                word.append(upward(current))
                current, a = odd_step(current)
                valuations.append(a)
            expected_floors = [floor_log2_fraction(target*Fraction(3,2)**j)
                               for j in range(length+1)]
            expected = tuple(b-a for a,b in zip(expected_floors,expected_floors[1:]))
            if tuple(word) == expected:
                require(set(valuations) == {1}, "all-odd saturation valuations")
                require(current > source, "all-odd saturation grows")
                found.add(tuple(word))
                break
            exponent += 1
    return found


def safe_height(length, rho):
    # A sufficient odd-source lower bound; not asserted optimal.
    def safe(m):
        return Fraction(3*m+1, 3*m)**length < rho
    high = 1
    while not safe(high):
        high *= 2
    low = 1
    while low < high:
        middle = (low+high)//2
        if safe(middle):
            high = middle
        else:
            low = middle+1
    return low


def laminar_excursions(values):
    stack, pairs, unmatched_down = [], [], []
    for i, (a, b) in enumerate(zip(values, values[1:])):
        delta = height(b)-height(a)
        require(delta in (-1, 0, 1), "nearest-neighbor height")
        if delta == 1:
            stack.append(i)
        elif delta == -1:
            if stack:
                pairs.append((stack.pop(), i+1))
            else:
                unmatched_down.append(i+1)
    for a, b in pairs:
        require(height(values[a]) == height(values[b]), "return height")
        require(all(height(x) > height(values[a]) for x in values[a+1:b]),
                "strict excursion interior")
    for a, b in pairs:
        for c, d in pairs:
            require(not (a < c < b < d), "nonlaminar pairs")
    return pairs, stack, unmatched_down


def main():
    limit = 100_000
    exact_edges = 0
    for n in range(1, limit+1, 2):
        nxt, a = odd_step(n)
        u = upward(n)
        require(u in (0, 1), "up-crossing alphabet")
        require(height(nxt)-height(n) == u-(a-1), "crossing flux")
        require(not (u == 0 and upward(nxt) == 0), "00 impossible")
        nxt2, _ = odd_step(nxt)
        if min(n, nxt, nxt2) >= 7:
            require((u, upward(nxt), upward(nxt2)) != (1, 1, 1),
                    "111 above small core")
        exact_edges += 1
    print(f"exact_odd_sources={exact_edges}, range=1..{limit}, flux/no00/no111 PASS")

    tables = []
    controls = 0
    for length in range(1, 17):
        language, rho = sturmian_data(length)
        safe = safe_height(length, rho)
        eligible = 0
        for n in range(1, limit+1, 2):
            current = n
            states, word = [], []
            for _ in range(length):
                states.append(current)
                word.append(upward(current))
                current, _ = odd_step(current)
            if min(states) >= safe:
                require(tuple(word) in language, "finite-height language control")
                eligible += 1
        controls += eligible
        tables.append((length, len(language), safe, eligible))
    print("length, language_size, sufficient_minimum, eligible_exact_controls")
    for row in tables:
        print(*row, sep=", ")
    print(f"finite_height_language_controls={controls} PASS")

    saturated = 0
    for length in range(1, 17):
        language, _ = sturmian_data(length)
        require(all_odd_realizations(length) == language,
                "every finite Sturmian word has all-a=1 realization")
        saturated += len(language)
    print(f"all_a1_saturation: lengths1..16, words={saturated} PASS")

    reset_count = 0
    for k in range(3, 202, 2):
        source = ((1 << (k+2))-5)//3
        target, a = odd_step(source)
        require(target == (1 << k)-1 and a == 2, "shallow reset family")
        require(v2(source+1) == 1 and v2(target+1) == k, "resource reset")
        require((source % 3 != 0) == (k % 6 in (1, 5)), "internal-image residue")
        endpoint = target
        internal_allowed = k % 6 in (1, 5)
        for _ in range(k-1):
            if internal_allowed:
                require(endpoint % 3 != 0 and endpoint % 5 != 0,
                        "full odd macroblock survives mod30")
            endpoint, a = odd_step(endpoint)
            require(a == 1, "reset-followed-by-growing-run")
        require(endpoint == 2*3**(k-1)-1 and v2(endpoint+1) == 1,
                "same-resource expanding macroblock")
        require(source % 5 != 0 and target % 5 != 0, "mod30 does not remove reset")
        reset_count += 1
    for n in range(3, limit+1, 2):
        r = v2(n+1)
        nxt, a = odd_step(n)
        if a == 1:
            require(r >= 2 and v2(nxt+1) == r-1, "resource consumption")
            q_before = 3**r*((n+1) >> r)
            s = v2(nxt+1)
            q_after = 3**s*((nxt+1) >> s)
            require(q_before == q_after, "integer run invariant")
    print(f"resource_controls: shallow_resets={reset_count}, run_invariant alloddn<=100000 PASS")

    prime_sets = [(2,), (3,), (11,), (2, 3, 11), (2, 3, 5, 11),
                  (2, 3, 5, 7, 11, 13, 17, 19)]
    for primes in prime_sets:
        period = 6
        for prime in primes:
            if prime >= 5:
                period = lcm(period, multiplicative_order(2, prime),
                             multiplicative_order(3, prime))
        if primes == (2, 3, 11):
            require(period == 30, "2,3,11 period")
        if primes == (2, 3, 5, 11):
            require(period == 60, "2,3,5,11 period")
        for parameter in range(1, 4):
            k = 1+period*parameter
            source = ((1 << (k+2))-5)//3
            z = 2*3**(k-1)-1
            final, a = odd_step(z)
            require(a == 2 and final == (3**k-1)//2, "finite-prime final reset")
            require(tuple(vp(source+1,p) for p in primes) == tuple(vp(2,p) for p in primes),
                    "finite-prime source resource")
            require(tuple(vp(final+1,p) for p in primes) == tuple(vp(2,p) for p in primes),
                    "finite-prime endpoint resource")
            require(final > source and final % 3 != 0 and final % 5 != 0,
                    "finite-prime expansion and mod30 hostile")
        print(f"finite_prime_resources={primes}, period={period}, parameters=1..3 PASS")

    values = [27]
    while values[-1] != 1:
        x = values[-1]
        values.append(x//2 if x % 2 == 0 else (3*x+1)//2)
    pairs, residual, downs = laminar_excursions(values)
    growing = [(values[a], values[b], b-a) for a, b in pairs if values[b] > values[a]]
    require((27, 31, 3) in growing, "growing completed excursion hostile")
    require(len(residual) == 0, "27 final stack closes")
    require(len(downs) == height(27), "endpoint flux")
    print(f"orbit27: steps={len(values)-1}, matched_excursions={len(pairs)}, "
          f"unmatched_down={len(downs)}, growing_returns={growing}")

    # An independent exact multiplicative version of the logarithmic identity.
    for n0 in (1, 3, 27, 703, 10087, 35655, 626331, 2**80-1):
        n = n0
        correction, valuations, ups = Fraction(1), 0, 0
        samples = []
        for j in range(80):
            samples.append(n)
            correction *= Fraction(3*n+1, 3*n)
            ups += upward(n)
            n, a = odd_step(n)
            valuations += a
            require(Fraction(n, n0) == Fraction(3**(j+1), 2**valuations)*correction,
                    "independent multiplicative carry identity")
        C = log2(correction.numerator)-log2(correction.denominator)
        # Use exact integer heights: taking floating log2(n) modulo one
        # spuriously wraps a value just below a large power of two.
        diagnostic = ups-80*log2(1.5)-C-(log2(n0)-height(n0))+(log2(n)-height(n))
        print(f"source={n0}, odd_steps=80, up_crossings={ups}, even_halvings={valuations-80}, "
              f"carry_bits={C:.12f}, float_identity_error={diagnostic:.3g}")

    # Positive cycles show why finite carry cannot be silently assumed.
    require(odd_step(5, 3, -1) == (7, 1), "3n-1 hostile1")
    require(odd_step(7, 3, -1) == (5, 2), "3n-1 hostile2")
    require([odd_step(x, 5, 1)[0] for x in (13, 33, 83)] == [33, 83, 13],
            "5n+1 cycle hostile")
    require(upward(1) == 1, "ordinary 1-cycle crossing hostile")
    print("positive_cycle_controls: 3n-1 odd cycle(5,7); 5n+1(13,33,83); "
          "3n+1 odd fixed point1, all reciprocal sums divergent PASS")
    print("PASS: exact arithmetic controls; no numerical proof of a hypothetical orbit")


if __name__ == "__main__":
    main()
