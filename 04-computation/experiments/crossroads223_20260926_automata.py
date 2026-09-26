"""Exact finite residue/return-language probes; Python standard library only.

Run from repository root:
  python 04-computation/experiments/crossroads223_20260926_automata.py
No simulation is used to infer termination of an unbounded integer orbit.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb, isqrt


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def prime(p):
    return p >= 2 and all(p % q for q in range(2, isqrt(p) + 1))


def order(a, p):
    x = 1
    for k in range(1, p):
        x = x * a % p
        if x == 1:
            return k
    raise RuntimeError("not a multiplicative unit")


def affine(word, q=3, sigma=1):
    """Return integers A,B,C with word map x -> (A*x+C)/B."""
    a, b, c = 1, 1, 0
    for bit in word:
        a, c = q**bit * a, q**bit * c + sigma * bit * b
        b *= 2
    return a, b, c


def carry_dp(length, ones, p):
    """Normalized affine translation C/2^length, by chronological bits."""
    inv2 = pow(2, -1, p)
    states = {(0, 0): 1}
    for j in range(length):
        nxt = Counter()
        for (e, c), count in states.items():
            if e + length - j - 1 >= ones:
                nxt[e, c * inv2 % p] += count
            if e < ones:
                nxt[e + 1, (3 * c + 1) * inv2 % p] += count
        states = nxt
    return Counter({c: count for (e, c), count in states.items() if e == ones})


def carry_direct(length, ones, p):
    """Independent denominator-cleared closed formula, no affine recurrence."""
    weights = [[pow(2, i, p) * pow(3, ones - j - 1, p) % p
                for i in range(length)] for j in range(ones)]
    scale = pow(pow(2, length, p), -1, p)
    result = Counter()
    for positions in combinations(range(length), ones):
        carry = sum(weights[j][i] for j, i in enumerate(positions))
        result[carry * scale % p] += 1
    return result


def orbit_lengths(a, c, p):
    unseen, lengths = set(range(p)), []
    while unseen:
        start = min(unseen)
        x, length = start, 0
        while x in unseen:
            unseen.remove(x)
            x, length = (a * x + c) % p, length + 1
        check(x == start, "affine map was not a permutation")
        lengths.append(length)
    return Counter(lengths)


def expected_orbits(a, c, p):
    if a == 1:
        return Counter({1: p}) if c == 0 else Counter({p: 1})
    d = order(a, p)
    return Counter({1: 1, d: (p - 1) // d})


def valuation(n, p):
    check(n != 0, "valuation of zero needs infinity")
    exponent = 0
    while n % p == 0:
        n //= p
        exponent += 1
    return exponent


def sink_components(targets):
    """Kosaraju decomposition, returning only closed components."""
    size = len(targets)
    reverse = [[] for _ in range(size)]
    for v, edges in enumerate(targets):
        for w in edges:
            reverse[w].append(v)
    seen, finishing = set(), []
    for root in range(size):
        stack = [(root, False)]
        while stack:
            v, done = stack.pop()
            if done:
                finishing.append(v)
            elif v not in seen:
                seen.add(v)
                stack.append((v, True))
                stack.extend((w, False) for w in targets[v] if w not in seen)
    component, groups = {}, []
    for root in reversed(finishing):
        if root in component:
            continue
        group, stack = [], [root]
        while stack:
            v = stack.pop()
            if v in component:
                continue
            component[v] = len(groups)
            group.append(v)
            stack.extend(reverse[v])
        groups.append(group)
    return [group for j, group in enumerate(groups)
            if all(component[w] == j for v in group for w in targets[v])]


def tight_cycle(targets, weights, potential):
    size = len(targets)
    tight = [[w for w in targets[v] if potential[v] == weights[v] + potential[w]]
             for v in range(size)]
    color, path, position = [0]*size, [], {}
    for root in range(size):
        if color[root]:
            continue
        stack = [(root, 0)]
        color[root], position[root] = 1, 0
        path.append(root)
        while stack:
            v, j = stack[-1]
            if j == len(tight[v]):
                stack.pop()
                path.pop()
                position.pop(v)
                color[v] = 2
                continue
            w = tight[v][j]
            stack[-1] = (v, j+1)
            if color[w] == 1:
                return path[position[w]:]
            if color[w] == 0:
                color[w], position[w] = 1, len(path)
                path.append(w)
                stack.append((w, 0))
    raise RuntimeError("no zero-weight cycle certificate")


def minimum_cycle_mean(targets, parity):
    """Karp candidate, certified by exact integer potential + tight cycle."""
    size = len(targets)
    distance = [0]*size
    for _ in range(size):
        distance = [parity[v] + min(distance[w] for w in targets[v])
                    for v in range(size)]
    terminal, distance, maxima = distance, [0]*size, [-1.0]*size
    for step in range(size):
        for v in range(size):
            maxima[v] = max(maxima[v], (terminal[v]-distance[v])/(size-step))
        distance = [parity[v] + min(distance[w] for w in targets[v])
                    for v in range(size)]
    value = Fraction(min(maxima)).limit_denominator(size)
    weights = [value.denominator*b-value.numerator for b in parity]
    potential = [0]*size
    for _ in range(size):
        newer = [min(0, weights[v]+min(potential[w] for w in targets[v]))
                 for v in range(size)]
        if newer == potential:
            break
        potential = newer
    check(all(potential[v] <= weights[v]+potential[w]
              for v in range(size) for w in targets[v]), "minimum-mean lower certificate")
    cycle = tight_cycle(targets, weights, potential)
    check(sum(weights[v] for v in cycle) == 0, "minimum-mean upper cycle certificate")
    return value, cycle


def q7_fixed_adversary(level, escape):
    """Top lift except the negative states -1,...,-escape are replaced."""
    half = 2**(level-1)
    representatives = [p + (0 if half-escape <= p < half else half)
                       for p in range(half)]
    targets = [tuple({((7*s+sign)//2)%half for sign in (-1,1)}) if s%2
               else (s//2%half,) for s in representatives]
    best, witness, component_sizes = Fraction(-1), None, []
    for group in sink_components(targets):
        index = {v: j for j, v in enumerate(group)}
        local = [tuple(index[w] for w in targets[v]) for v in group]
        value, cycle = minimum_cycle_mean(local, [v%2 for v in group])
        component_sizes.append(len(group))
        if value > best:
            best = value
            witness = [representatives[group[v]] for v in cycle]
    return best, witness, sorted(component_sizes)


def q7_probe():
    print("Q7 fixed-adversary exact certificates; level escape value sink_sizes cycle_length")
    for level in (5, 7, 9, 11):
        for escape in (0, 1, 3, 7):
            value, cycle, sizes = q7_fixed_adversary(level, escape)
            if escape == 0:
                check(value == Fraction(1,3), "negative-adversary inherited control")
            print(level, escape, str(value), sizes, len(cycle))
    value, cycle, sizes = q7_fixed_adversary(13, 1)
    check(value == Fraction(12,41), "larger one-point escape control")
    print(13, 1, str(value), sizes, len(cycle))
    value, cycle, sizes = q7_fixed_adversary(7, 1)
    check(value == Fraction(2,9), "one-point escape hostile")
    signed = [v if v < 64 else v-128 for v in cycle]
    print("Level-7 escape witness in signed representatives:", signed)


def main():
    p = 223
    check(prime(p), "223 not prime")
    check(2**37 - 1 == p * 616318177, "Mersenne factor")
    check(order(2, p) == 37 and order(3, p) == 222, "orders")
    check(pow(3, 180, p) == 2, "discrete logarithm")
    print("223 prime; 2^37-1=223*616318177; ord(2)=37; ord(3)=222; 2=3^180")

    # Paley orientation is a relation on all distinct pairs, not a trajectory.
    chi = lambda z: 1 if pow(z % p, (p - 1) // 2, p) == 1 else -1
    paley_checks = 0
    for sigma, bit in product((-1, 1), (0, 1)):
        a = pow(3, bit, p) * pow(2, -1, p) % p
        c = sigma * bit * pow(2, -1, p) % p
        for x in range(p):
            for y in range(p):
                if x == y:
                    continue
                check(chi((a * y + c) - (a * x + c)) == (-1)**bit * chi(y-x),
                      "Paley covariance")
                paley_checks += 1
    print("Paley covariance checks:", paley_checks, "; even preserves, odd reverses, either sign")

    affine_checks = 0
    for modulus in (7, 31, 223):
        for length in range(1, 9):
            for word in product((0, 1), repeat=length):
                a, b, c = affine(word)
                aa, cc = a * pow(b, -1, modulus) % modulus, c * pow(b, -1, modulus) % modulus
                check(orbit_lengths(aa, cc, modulus) == expected_orbits(aa, cc, modulus),
                      "finite affine orbit classification")
                affine_checks += 1
    print("Affine permutation classifications:", affine_checks, "(all words lengths 1..8; p=7,31,223)")

    resonances = [(r, e) for r in range(1, 50) for e in range(r+1)
                  if (e - 180*r) % 222 == 0]
    expanding = [(r, e) for r, e in resonances if 3**e > 2**r]
    check(resonances[0] == (21, 6) and expanding == [(26, 18), (31, 30)], "resonance list")
    print("First resonance:", resonances[0], "; expanding resonances through length 49:", expanding)
    for length, ones in ((21, 6), (26, 18)):
        counts = carry_dp(length, ones, p)
        independent = carry_direct(length, ones, p)
        check(counts == independent, "DP/direct carry histograms differ")
        check(sum(counts.values()) == comb(length, ones), "word universe")
        print("Carry histogram", (length, ones), "total", sum(counts.values()),
              "zero", counts[0], "nonzero range", (min(counts[c] for c in range(1,p)),
                                                      max(counts[c] for c in range(1,p))))
        if (length, ones) == (26, 18):
            check(counts[0] == 7735, "zero-carry count")
            primitive_zero = counts[0] - comb(13, 9)
            primitive_total = comb(26, 18) - comb(13, 9)
            check((primitive_zero, primitive_total) == (7020, 1561560), "primitive words")
            check(primitive_zero % 26 == primitive_total % 26 == 0, "necklace integrality")
            print("Primitive necklaces: zero", primitive_zero//26,
                  "; nonzero", (primitive_total-primitive_zero)//26)

    doubled = 0
    for positions in combinations(range(13), 9):
        word = tuple(int(i in positions) for i in range(13))
        a, b, c = affine(word)
        aa, cc = a*pow(b,-1,p)%p, c*pow(b,-1,p)%p
        check(aa == p-1, "13-block multiplier")
        check((aa*aa, (aa+1)*cc%p) == ((p-1)**2, 0), "reflection squares")
        a2, b2, c2 = affine(word+word)
        check(a2%p == b2%p and c2%p == 0, "doubled word identity")
        doubled += 1
    check(doubled == 715, "doubled universe")
    print("All 715 imprimitive (26,18) words are doubled reflections, hence residue identities")

    # Powers do not alter the characteristic-zero rational fixed-point quotient.
    power_checks = 0
    for q, sigma in product((3, 5), (-1, 1)):
        for length in range(1, 8):
            for word in product((0, 1), repeat=length):
                a, b, c = affine(word, q, sigma)
                for repeats in (2, 3, 7):
                    ad, bd, cd = affine(word*repeats, q, sigma)
                    check(Fraction(cd, bd-ad) == Fraction(c,b-a), "rational fixed point changed")
                    power_checks += 1
    print("Fixed-point invariance under word powers:", power_checks, "exact rational checks")

    word = (1,)*18 + (0,)*8
    a, b, c = affine(word)
    ad, bd, cd = affine(word*p)
    check(c%p == 17 and valuation(b-a,p) == 1, "translation obstruction control")
    check((valuation(cd,p), valuation(bd-ad,p)) == (1,2), "power valuation cancellation")
    check(Fraction(c,b-a) == Fraction(cd,bd-ad) == Fraction(-77431669,64062325),
          "223 repeats cannot repair denominator")
    print("Explicit 1^18 0^8 control: source -77431669/64062325; (v223(C),v223(B-A))")
    print("  before (0,1), after 223 repeats (1,2); finite residue return, same nonintegral source")

    # Every finite all-ones language fragment occurs at positive integers,
    # while its mod-223 coordinate can remain the apparent fixed state -1.
    for horizon in range(1, 129):
        n = p * 2**horizon - 1
        for j in range(horizon+1):
            check(n == p * 3**j * 2**(horizon-j) - 1, "positive all-ones family")
            check(n%p == p-1, "fixed auxiliary residue")
            if j < horizon:
                check(n%2 == 1, "odd prefix")
                n = (3*n+1)//2
    n = 7*2**37
    check(n%p == 7 and n//2**37 == 7 and n != 7, "37-zero false integer cycle control")
    for q, sigma, word, start in ((3,1,(1,0),1), (3,-1,(1,1,0),5), (5,1,(1,1,0,0,0),1)):
        n = start
        for bit in word:
            check(n%2 == bit, "genuine cycle parity control")
            n = (q*n+sigma)//2 if bit else n//2
        check(n == start, "genuine signed cycle control")
    print("Positive controls: horizons 1..128; 37-zero false closure; genuine +3, -3, +5 cycles")
    q7_probe()
    print("PASS: finite exact arithmetic only; no Collatz convergence claim")


if __name__ == "__main__":
    main()
