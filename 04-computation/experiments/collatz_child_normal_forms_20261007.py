"""Exact guarded normal forms for three negative-cycle inverse branches.

No ROOT oracle is used by the production APIs.  A terminal router state is
an obligation, not a convergence certificate.  Run with Python 3, also -O.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product


# Label: (binary exponent, ternary exponent, carry, actual forward word).
GENERATORS = {
    1: (1, 1, 1, (1,)),
    5: (3, 2, 5, (1, 2)),
    17: (11, 7, 2363, (1, 1, 1, 2, 1, 1, 4)),
}


def natural(n):
    if type(n) is not int or n < 0:
        raise ValueError("expected an exact nonnegative integer")
    return n


def odd(n):
    if type(n) is not int or n <= 0 or n % 2 != 1:
        raise ValueError("expected an exact positive odd integer")
    return n


def letters(word):
    if type(word) is not tuple or any(type(g) is not int or g not in GENERATORS for g in word):
        raise ValueError("expected a tuple of generator labels 1,5,17")
    return word


def valuation(n, p):
    if type(n) is not int or n == 0 or type(p) is not int or p < 2:
        raise ValueError("nonzero integer and integer base >=2 required")
    n = abs(n)
    a = 0
    while n % p == 0:
        a, n = a + 1, n // p
    return a


@dataclass(frozen=True)
class Carrier:
    A: int
    B: int
    C: int

    def __post_init__(self):
        for x in (self.A, self.B, self.C):
            natural(x)


IDENTITY = Carrier(0, 0, 0)


def encode(word):
    letters(word)
    A = B = C = 0
    for g in word:
        a, b, c, _ = GENERATORS[g]
        C = (1 << a) * C + c * 3**B
        A, B = A + a, B + b
    return Carrier(A, B, C)


def decode(record):
    if type(record) is not Carrier:
        raise ValueError("exact Carrier required")
    A, B, C = record.A, record.B, record.C
    result = []
    while (A, B, C) != (0, 0, 0):
        if A <= 0 or B <= 0 or C <= 0 or C % 2 != 1 or C % 3 == 0:
            raise ValueError("not a generated carrier")
        candidates = []
        for g, (a, b, c, _) in GENERATORS.items():
            if A >= a and B >= b:
                modulus = 3**b
                if (C * pow(2, -A, modulus) - c * pow(2, -a, modulus)) % modulus == 0:
                    candidates.append(g)
        if len(candidates) != 1:
            raise ValueError("no unique first native letter")
        g = candidates[0]
        a, b, c, _ = GENERATORS[g]
        numerator = C - (1 << (A - a)) * c
        if numerator < 0 or numerator % (3**b):
            raise ValueError("invalid first-letter carry")
        A, B, C = A - a, B - b, numerator // (3**b)
        result.append(g)
    word = tuple(result)
    if encode(word) != record:
        raise ValueError("failed carrier authentication")
    return word


def compose(first, second):
    """Apply first, then second; authenticate both records."""
    decode(first)
    decode(second)
    return Carrier(first.A + second.A, first.B + second.B,
                   (1 << second.A) * first.C + second.C * 3**first.B)


def native_cell(record):
    decode(record)
    modulus = 3**record.B
    return (record.C * pow(2, -record.A, modulus)) % modulus, modulus


def apply(record, source):
    odd(source)
    residue, modulus = native_cell(record)
    if source % modulus != residue:
        raise ValueError("source is outside the exact native cylinder")
    child = ((1 << record.A) * source - record.C) // modulus
    odd(child)
    return child


def forward_word(word):
    letters(word)
    return tuple(a for g in reversed(word) for a in GENERATORS[g][3])


def step(n):
    odd(n)
    a = valuation(3 * n + 1, 2)
    return (3 * n + 1) >> a, a


def replay(n, word, root=False):
    odd(n)
    if type(word) is not tuple or any(type(a) is not int or a <= 0 for a in word):
        raise ValueError("exact positive valuation word required")
    if type(root) is not bool:
        raise ValueError("root flag must be bool")
    for a in word:
        if n == 1:
            raise ValueError("word crosses its first ROOT hit")
        n, actual = step(n)
        if a != actual:
            raise ValueError("incorrect actual valuation")
    if root and n != 1:
        raise ValueError("word does not end at ROOT")
    return n


def next_letter(n):
    odd(n)
    choices = [g for g in GENERATORS if n % (3**GENERATORS[g][1]) ==
               GENERATORS[g][2] * pow(2, -GENERATORS[g][0], 3**GENERATORS[g][1]) % (3**GENERATORS[g][1])]
    if len(choices) > 1:
        raise ArithmeticError("disjoint native guards violated")
    return choices[0] if choices else None


@dataclass(frozen=True)
class Reduction:
    source: int
    terminal: int
    word: tuple


def route(source):
    odd(source)
    n, word = source, []
    cap = 11 * source.bit_length()
    while (g := next_letter(n)) is not None:
        if len(word) >= cap:
            raise ArithmeticError("proved height bound violated")
        child = apply(encode((g,)), n)
        if child >= n:
            raise ArithmeticError("strict child inequality violated")
        word.append(g)
        n = child
    return Reduction(source, n, tuple(word))


def discharge(reduction, terminal_root_word):
    """Extract a source ROOT suffix from an explicitly supplied terminal proof."""
    if type(reduction) is not Reduction:
        raise ValueError("exact Reduction required")
    odd(reduction.source)
    odd(reduction.terminal)
    letters(reduction.word)
    if route(reduction.source) != reduction:
        raise ValueError("forged reduction")
    replay(reduction.terminal, terminal_root_word, root=True)
    prefix = forward_word(reduction.word)
    if terminal_root_word[:len(prefix)] != prefix:
        raise ValueError("terminal proof lacks the required source occurrence")
    suffix = terminal_root_word[len(prefix):]
    replay(reduction.source, suffix, root=True)
    return suffix


def two_log_three_power(unit, depth):
    """Exact discrete log modulo 3**depth; returns the canonical exponent."""
    natural(depth)
    if type(unit) is not int or unit % 3 == 0 or depth < 1:
        raise ValueError("unit and positive ternary precision required")
    exponent = 0 if unit % 3 == 1 else 1
    period, modulus = 2, 3
    for _ in range(1, depth):
        modulus *= 3
        candidates = [exponent + j * period for j in range(3)
                      if pow(2, exponent + j * period, modulus) == unit % modulus]
        if len(candidates) != 1:
            raise ArithmeticError("primitive base-two lifting failed")
        exponent = candidates[0]
        period *= 3
    return exponent, period


def mersenne_phase(record):
    """K>=1 in returned class iff 2**K-1 is native; None means empty."""
    decode(record)
    if record == IDENTITY:
        return 0, 1
    D = (1 << record.A) + record.C
    if D % 3 == 0:
        return None
    modulus = 3**record.B
    exponent, period = two_log_three_power(D * pow(2, -record.A, modulus), record.B)
    return exponent, period


def exponential_state(record, K):
    """Exact child value. This explicit API expands 2**K; phase() need not."""
    natural(K)
    if K < 1:
        raise ValueError("positive Mersenne exponent required")
    phase = mersenne_phase(record)
    if phase is None or (K - phase[0]) % phase[1]:
        raise ValueError("exponent is outside the native phase")
    return apply(record, (1 << K) - 1)


def cylinder_mass(word):
    return F(1, 3**encode(word).B)


def literal_root(n, cap=10000):
    """Experiment-only discovery, never called by any production API above."""
    result = []
    for _ in range(cap):
        if n == 1:
            return tuple(result)
        n, a = step(n)
        result.append(a)
    raise RuntimeError("finite experimental discovery cap exceeded")


def main():
    checks = 0

    def check(predicate):
        nonlocal checks
        checks += 1
        if not predicate:
            raise ArithmeticError("exact control failed")

    def rejects(call):
        try:
            call()
        except ValueError:
            check(True)
        else:
            check(False)

    # Independently recover the inherited inverse branch carriers.
    for g, (a, b, c, w) in GENERATORS.items():
        P, Q, B = 1, 1, 0
        for exponent in w:
            P, Q, B = 3 * P, (1 << exponent) * Q, 3 * B + Q
        check((P, Q, B) == (3**b, 1 << a, c))
        check(c == g * (3**b - (1 << a)))
    check(2 * 2048**11 < 2187**11)

    words = [w for depth in range(7) for w in product(GENERATORS, repeat=depth)]
    phase_expansions = 0
    for word in words:
        record = encode(word)
        check(decode(record) == word)
        residue, modulus = native_cell(record)
        least = residue if residue % 2 else residue + modulus
        if least <= 0:
            least += 2 * modulus
        for lift in (0, 1, 7):
            source = least + 2 * modulus * lift
            n = source
            for g in word:
                child = apply(encode((g,)), n)
                check(0 < child < n)
                check(F(child + 1, n + 1) <= F(2048, 2187))
                n = child
            check(apply(record, source) == n)
            check(replay(n, forward_word(word)) == source)
        for cut in range(len(word) + 1):
            check(compose(encode(word[:cut]), encode(word[cut:])) == record)
        phase = mersenne_phase(record)
        check((phase is None) == bool(word and word[0] == 1))
        if phase is not None and word:
            exponent, period = phase
            check(exponent % 2 == 1)
            D = (1 << record.A) + record.C
            check(pow(2, exponent + record.A, 3**record.B) == D % (3**record.B))
            check(period == 2 * 3**(record.B - 1))
            if exponent <= 4096:
                child = exponential_state(record, exponent)
                check(replay(child, forward_word(word)) == (1 << exponent) - 1)
                phase_expansions += 1

    survival = F(973, 2187)
    for depth in range(1, 7):
        level = [w for w in words if len(w) == depth]
        check(sum(map(cylinder_mass, level), F(0)) == survival**depth)
        check(sum((3 * cylinder_mass(w) for w in level if w[0] != 1), F(0)) ==
              F(244, 729) * survival**(depth - 1))
    counts = {g: 0 for g in GENERATORS}
    for j in range(2187):
        g = next_letter(2 * j + 1)
        if g is not None:
            counts[g] += 1
    check(counts == {1: 729, 5: 243, 17: 1})

    max_depth, deepest_source, terminals = 0, 1, set()
    for n in range(1, 16384, 2):
        reduced = route(n)
        check(apply(encode(reduced.word), n) == reduced.terminal)
        check(next_letter(reduced.terminal) is None)
        check(len(reduced.word) <= 11 * n.bit_length())
        check(replay(reduced.terminal, forward_word(reduced.word)) == n)
        terminals.add(reduced.terminal)
        if len(reduced.word) > max_depth:
            max_depth, deepest_source = len(reduced.word), n
        if n <= 1023:
            proof = literal_root(reduced.terminal)
            check(discharge(reduced, proof) == literal_root(n))

    # Concrete arbitrary-depth closure, terminal counterfamily, and order loss.
    check(mersenne_phase(encode((5,))) == (5, 6))
    check(mersenne_phase(encode((17,))) == (733, 1458))
    check(route(31) == Reduction(31, 27, (5,)))
    for t in range(32):
        K = 5 + 18 * t
        child = exponential_state(encode((5,)), K)
        check(child % 3 == 0)
        check(next_letter(child) is None)
        check(replay(child, (1, 2)) == (1 << K) - 1)
    for r in range(1, 5):
        K = 2 + 3**(2*r-1)
        n = (1 << K) - 1
        check(valuation(n + 5, 3) == 2*r)
        check(route(n).word[:r] == (5,) * r)
        child = apply(encode((5,) * r), n)
        check((child + 5) % 9 != 0)
    check(encode((1, 17)) == Carrier(12, 8, 9137))
    check(encode((17, 1)) == Carrier(12, 8, 6913))
    check(encode((5,) * 4) == Carrier(12, 8, 12325))
    check(encode((1, 5)) == Carrier(4, 3, 23))
    check(encode((5, 1)) == Carrier(4, 3, 19))
    r15, m15 = native_cell(encode((1, 5)))
    r51, m51 = native_cell(encode((5, 1)))
    check(m15 == m51 == 27 and r15 != r51)
    check(F(23 - 19, 27) == F(4, 27))

    for bad in (True, 3.0, 0, 2, -1):
        rejects(lambda bad=bad: route(bad))
    for bad in ((True,), (1.0,), (2,), [1]):
        rejects(lambda bad=bad: encode(bad))
    for args in ((True, 0, 0), (0, 0.0, 0), (-1, 1, 1)):
        rejects(lambda args=args: Carrier(*args))
    for bad in (Carrier(1, 1, 5), Carrier(0, 1, 0), Carrier(1, 0, 1)):
        rejects(lambda bad=bad: decode(bad))
    rejects(lambda: apply(encode((5,)), 27))
    rejects(lambda: exponential_state(encode((1,)), 5))
    rejects(lambda: exponential_state(encode((5,)), True))
    rejects(lambda: two_log_three_power(3, 2))
    rejects(lambda: replay(True, (), root=True))
    rejects(lambda: replay(1.0, (), root=True))
    rejects(lambda: discharge(Reduction(True, 1, ()), ()))
    rejects(lambda: discharge(Reduction(31, 3, (5,)), (1, 4)))
    rejects(lambda: discharge(route(31), (1, 2)))
    rejects(lambda: replay(1, (2,), root=True))
    check(discharge(route(1), ()) == ())

    print("Guarded child normal forms: exact deterministic router, not universal ROOT closure.")
    print("Generators: G1=(2n-1)/3; G5=(8n-5)/9; G17=(2048n-2363)/2187.")
    print("Native guards: 2 mod3; 4 mod9; 2170 mod2187 (pairwise disjoint).")
    print(f"Generated words: {len(words)} (length0..6); explicit Mersenne phase replays: {phase_expansions}.")
    print(f"Odd-source census1..16383: 8192; distinct terminals: {len(terminals)}; maximum router depth: {max_depth} at {deepest_source}.")
    print("Supplied terminal ROOT proof transfers: all512 odd sources1..1023; experimental discovery separate.")
    print("Finite-depth source survival=(973/2187)^k; odd-exponent survival=(244/729)(973/2187)^(k-1).")
    print("Slope alias: (1,17),(17,1),(5,5,5,5) share(A,B)=(12,8), carries9137,6913,12325.")
    print("Mersenne phases: G1 empty; G5 K5mod6; G17 K733mod1458.")
    print("Hostile: every K5mod18 gives a3-divisible immediate terminal, beginning31->27.")
    print("Native child is smaller; its retained forward word returns to the original source, so ROOT is still an obligation.")
    print(f"Exact checks: {checks}")


if __name__ == "__main__":
    main()
