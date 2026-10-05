"""Exact finite affine groups and guarded lifts of every group element.

No assertion that an arbitrary supplied source reaches a supplied hub.
All certificates construct an ancestor of that hub; large exponents stay
symbolic. Checks deliberately survive python -O. Stdlib only.
"""

from collections import deque
from dataclasses import dataclass
from functools import lru_cache
from fractions import Fraction
from itertools import product
from math import gcd, lcm


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def integer(value, minimum, name):
    if type(value) is not int or value < minimum:
        raise ValueError(name)


def modulus_ok(m):
    integer(m, 5, "modulus")
    if gcd(m, 6) != 1:
        raise ValueError("modulus must be coprime to six")


def compose(f, g, m):
    """f after g, for affine pairs (slope, intercept)."""
    return (f[0] * g[0] % m, (f[0] * g[1] + f[1]) % m)


def inverse(f, m):
    a = pow(f[0], -1, m)
    return (a, -a * f[1] % m)


def generators(m):
    return ((2 % m, 0), (2 * pow(3, -1, m) % m, -pow(3, -1, m) % m))


def order(a, m):
    if gcd(a, m) != 1:
        raise ValueError("nonunit")
    k, x = 1, a % m
    while x != 1:
        x = x * a % m
        k += 1
    return k


@lru_cache(None)
def slope_group(m):
    modulus_ok(m)
    found, queue = {1}, deque([1])
    while queue:
        x = queue.popleft()
        for a in (2, 3):
            y = a * x % m
            if y not in found:
                found.add(y)
                queue.append(y)
    return frozenset(found)


@lru_cache(None)
def affine_group(m):
    """Independent BFS using inverse Collatz generators, not the formula."""
    modulus_ok(m)
    gens = generators(m)
    found, queue = {(1, 0)}, deque([(1, 0)])
    while queue:
        g = queue.popleft()
        for f in gens:
            h = compose(f, g, m)
            if h not in found:
                found.add(h)
                queue.append(h)
    return frozenset(found)


def prime_powers(n):
    factors, p = [], 2
    while p * p <= n:
        if n % p == 0:
            q = 1
            while n % p == 0:
                q *= p
                n //= p
            factors.append(q)
        p += 1
    if n > 1:
        factors.append(n)
    return factors


def word_data(word):
    """Chronological D/E inverse word: (2^length*x-B)/3^reds."""
    e, b = 0, 0
    for letter in word:
        if letter == "D":
            b *= 2
        elif letter == "E":
            b = 2 * b + 3**e
            e += 1
        else:
            raise ValueError("unknown inverse letter")
    return len(word), e, b


def word_mod(word, m):
    d, e = generators(m)
    g = (1, 0)
    for letter in word:
        g = compose(d if letter == "D" else e, g, m)
    return g


@dataclass(frozen=True)
class Template:
    modulus: int
    word: str
    marks: tuple  # (affine pair, prefix length), immediately after E
    reds: int
    guard: int


@lru_cache(None)
def cover_template(m):
    """Euler circuit in the positive two-generator Cayley graph, then E."""
    group, gens = affine_group(m), generators(m)
    used = {g: 0 for g in group}
    stack, incoming, reverse_letters = [(1, 0)], [], []
    while stack:
        g = stack[-1]
        i = used[g]
        if i < 2:
            used[g] += 1
            stack.append(compose(gens[i], g, m))
            incoming.append("DE"[i])
        else:
            stack.pop()
            if incoming:
                reverse_letters.append(incoming.pop())
    circuit = "".join(reversed(reverse_letters))
    check(len(circuit) == 2 * len(group), "Euler edge count")
    g, marks = (1, 0), []
    for j, letter in enumerate(circuit, 1):
        g = compose(gens[letter == "E"], g, m)
        if letter == "E":
            marks.append((g, j))
    check(g == (1, 0), "Euler circuit closes")
    check(len(marks) == len(group) and {g for g, _ in marks} == group,
          "E-marked prefixes cover entire group exactly once")
    word = circuit + "E"
    length, reds, carry = word_data(word)
    ternary = 3**reds
    guard = carry * pow(2**length, -1, ternary) % ternary
    check(guard % 3 != 0 and reds == len(group) + 1, "unit master guard")
    return Template(m, word, tuple(marks), reds, guard)


def ternary_log_two(value, exponent):
    """Digit lifting; returns k in [0, 2*3^(exponent-1))."""
    integer(exponent, 1, "ternary exponent")
    if gcd(value, 3) != 1:
        raise ValueError("ternary nonunit")
    k, period, modulus = (0 if value % 3 == 1 else 1), 2, 3
    for _ in range(1, exponent):
        modulus *= 3
        choices = [k + digit * period for digit in range(3)
                   if pow(2, k + digit * period, modulus) == value % modulus]
        check(len(choices) == 1, "unique next ternary discrete-log digit")
        k = choices[0]
        period *= 3
    return k


@dataclass(frozen=True)
class InverseCertificate:
    hub: int
    modulus: int
    target: tuple
    doublings: int
    word: str

    def endpoint_mod(self, q):
        integer(q, 1, "output modulus")
        length, reds, carry = word_data(self.word)
        denominator = 3**reds
        numerator = (pow(2, self.doublings + length, q * denominator)
                     * self.hub - carry) % (q * denominator)
        check(numerator % denominator == 0, "endpoint division precision")
        return numerator // denominator

    def verify(self):
        integer(self.hub, 1, "hub")
        if self.hub % 2 == 0 or self.hub % 3 == 0:
            raise ValueError("hub must be a positive odd ternary unit")
        modulus_ok(self.modulus)
        integer(self.doublings, 1, "doublings")
        if (type(self.word) is not str or not self.word or self.word[-1] != "E"
                or any(x not in "DE" for x in self.word)):
            raise ValueError("inverse word must end in E")
        if (type(self.target) is not tuple or len(self.target) != 2
                or any(type(x) is not int for x in self.target)
                or any(not 0 <= x < self.modulus for x in self.target)):
            raise ValueError("affine target")
        length, reds, carry = word_data(self.word)
        if self.doublings < 2 * length:
            raise ValueError("strict-height guard")
        if (pow(2, self.doublings + length, 3**reds) * self.hub - carry) % 3**reds:
            raise ValueError("ternary legality guard")
        actual = compose(word_mod(self.word, self.modulus),
                         (pow(2, self.doublings, self.modulus), 0), self.modulus)
        if actual != self.target:
            raise ValueError("affine target mismatch")
        if self.endpoint_mod(3) == 0:
            raise ValueError("output is a ternary leaf")
        check(self.endpoint_mod(2) == 1, "last E gives odd endpoint")
        return True

    def odd_inverse_exponents(self):
        pending, result = self.doublings, []
        for letter in self.word:
            if letter == "D":
                pending += 1
            else:
                result.append(pending + 1)
                pending = 0
        return tuple(result)


def lift_all(m, hub):
    """One source-aware guard calculation lifts every affine target."""
    integer(hub, 1, "hub")
    if hub % 2 == 0 or hub % 3 == 0:
        raise ValueError("hub must be a positive odd ternary unit")
    template = cover_template(m)
    ternary = 3**template.reds
    wanted = template.guard * pow(hub, -1, ternary) % ternary
    k = ternary_log_two(wanted, template.reds)
    joint = lcm(2 * 3 ** (template.reds - 1), order(2, m))
    threshold = 2 * len(template.word)
    if k < threshold:
        k += ((threshold - k + joint - 1) // joint) * joint
    scaling = (pow(2, k, m), 0)
    result = {}
    for g, j in template.marks:
        target = compose(g, scaling, m)
        result[target] = InverseCertificate(hub, m, target, k, template.word[:j])
    check(set(result) == affine_group(m), "right coset reindexing")
    return result


def stepwise_mod(cert, q):
    """Independent prefix reader consumes one ternary digit per E."""
    integer(q, 1, "output modulus")
    precision = q * 3 ** cert.word.count("E")
    x = pow(2, cert.doublings, precision) * cert.hub % precision
    for letter in cert.word:
        if letter == "D":
            x = 2 * x % precision
        else:
            check(x % 3 == 2, "independent intermediate integer guard")
            x = (2 * x - 1) // 3
            precision //= 3
            x %= precision
    check(precision == q, "precision budget consumed exactly")
    return x


def repeat_refuel(word, copies, hub, m=5):
    """Refuel a fixed repeated inverse block; target is computed, not imposed."""
    integer(copies, 1, "copies")
    integer(hub, 1, "hub")
    modulus_ok(m)
    if hub % 2 == 0 or hub % 3 == 0:
        raise ValueError("hub must be a positive odd ternary unit")
    if type(word) is not str or not word or word[-1] != "E":
        raise ValueError("repeated block must end in E")
    length, reds, carry = word_data(word)
    power = reds * copies + 1  # sufficient extra digit for a ternary-unit endpoint
    ternary = 3**power
    axis = carry * pow(2**length - 3**reds, -1, ternary) % ternary
    wanted = axis * pow(hub, -1, ternary) % ternary
    a = ternary_log_two(wanted, power)
    period = lcm(2 * 3 ** (power - 1), order(2, m))
    threshold = 2 * length * copies
    if a < threshold:
        a += ((threshold - a + period - 1) // period) * period
    repeated = word * copies
    target = compose(word_mod(repeated, m), (pow(2, a, m), 0), m)
    return InverseCertificate(hub, m, target, a, repeated)


def main():
    print("FINITE AFFINE GROUPS AND GUARDED LIFTS; arbitrary supplied-source convergence OPEN")
    moduli = [m for m in range(5, 152) if gcd(m, 6) == 1]
    mismatches, group_sum = [], 0
    for m in moduli:
        group, slopes = affine_group(m), slope_group(m)
        check(group == {(a, b) for a in slopes for b in range(m)}, "semidirect formula")
        d, e = generators(m)
        commutator = compose(compose(compose(d, e, m), inverse(d, m), m), inverse(e, m), m)
        check(commutator == (1, -pow(3, -1, m) % m), "unit commutator translation")
        group_sum += len(group)
        local_product = 1
        for q in prime_powers(m):
            local_product *= len(slope_group(q))
        if local_product != len(slopes):
            mismatches.append((m, len(slopes), local_product))
    check(mismatches[0] == (95, 36, 72), "least composite-product hostile")
    check(pow(2, 31, 95) == 3 and order(2, 95) == 36, "exact cyclic H95")
    def character(a, p):
        return 1 if pow(a, (p - 1) // 2, p) == 1 else -1
    correlated = {a for a in range(95) if gcd(a, 95) == 1 and character(a, 5) == character(a, 19)}
    check(slope_group(95) == correlated and 7 not in correlated, "quadratic correlation equals H95")
    print("Group universe: all", len(moduli), "moduli5..151 coprime6;", group_sum, "affine elements")
    print("Full local-product mismatches (modulus, true slope size, local product):", mismatches)
    print("Mod95: H=<2>, order36; 3=2^31; G order3420, not6840; missing slope7")

    # Explicit group-derived-series checks in a separate small universe.
    commutators = 0
    for m in (5, 7, 11):
        group = affine_group(m)
        offsets = set()
        for f, g in product(group, repeat=2):
            h = compose(compose(compose(f, g, m), inverse(f, m), m), inverse(g, m), m)
            check(h[0] == 1, "derived group contains translations only")
            offsets.add(h[1])
            commutators += 1
        check(offsets == set(range(m)), "all translations already commutators")
    print("Independent derived-series pair controls:", commutators)

    lift_checks = reader_checks = literal_marks = literal_edges = 0
    widths = []
    for m in (5, 7, 11, 19):
        template = cover_template(m)
        for hub in (1, 5, 7):
            lifts = lift_all(m, hub)
            widths.append((m, hub, next(iter(lifts.values())).doublings.bit_length()))
            for target, cert in lifts.items():
                check(cert.verify(), "compressed certificate")
                check(cert.endpoint_mod(m) == (target[0] * hub + target[1]) % m, "target action")
                for q in (2, 9, 19, 95, 4096):
                    check(cert.endpoint_mod(q) == stepwise_mod(cert, q), "independent modular reader")
                    reader_checks += 1
                exponents = cert.odd_inverse_exponents()
                check(len(exponents) == cert.word.count("E") and
                      sum(exponents) == cert.doublings + len(cert.word), "odd rank/cost")
                lift_checks += 1

        # Choose a legal hub for the SMALL exponent2L; unlike arbitrary hubs,
        # this permits complete expansion and independent literal replay.
        k = 2 * len(template.word)
        ternary = 3**template.reds
        hub = template.guard * pow(pow(2, k, ternary), -1, ternary) % ternary
        if hub % 2 == 0:
            hub += ternary
        lifts = lift_all(m, hub)
        check(next(iter(lifts.values())).doublings == k, "small literal lifting exponent")
        x, action, marked = hub << k, (pow(2, k, m), 0), {}
        gens = generators(m)
        for j, letter in enumerate(template.word, 1):
            old = x
            if letter == "D":
                x *= 2
                check(x % 2 == 0 and x // 2 == old, "literal D reverse edge")
            else:
                check(x % 3 == 2, "literal E guard")
                x = (2 * x - 1) // 3
                check(x % 2 == 1 and (3 * x + 1) // 2 == old, "literal E reverse edge")
            check(x > hub, "strict first-hit hub height")
            action = compose(gens[letter == "E"], action, m)
            literal_edges += 1
            if letter == "E" and j < len(template.word):
                check(x % 3 != 0, "marked state remains extendable")
                cert = lifts[action]
                check(cert.endpoint_mod(x + 1) == x, "full literal endpoint")
                marked[action] = x
                literal_marks += 1
        check(set(marked) == affine_group(m), "literal all-map coverage")
    print("Compressed all-map lifts:", lift_checks, "; independent modular readers:", reader_checks)
    print("Initial exponent BIT LENGTHS (modulus, hub, bits):", widths)
    print("Literal high-hub controls:", literal_edges, "D/E edges;", literal_marks, "marked maps")

    blocks = ("E", "DE", "DDE", "EE", "EDE", "DEDE", "DDEE")
    fuel_checks = refill_checks = refill_literal_edges = 0
    for block, copies in product(blocks, range(1, 5)):
        length, reds, carry = word_data(block)
        q, d = 2**length, 3**reds
        for start in range(1, 101):
            expected = ((q - d) * start - carry) % 3 ** (reds * copies) == 0
            x, legal = start, True
            for letter in block * copies:
                if letter == "D":
                    x *= 2
                elif x % 3 != 2:
                    legal = False
                    break
                else:
                    x = (2 * x - 1) // 3
                check(x > 0, "independent repetition positivity")
            check(legal == expected, "exact repeated-block fuel criterion")
            fuel_checks += 1
        for hub in (1, 5, 7):
            cert = repeat_refuel(block, copies, hub)
            check(cert.verify(), "compressed repeated-block certificate")
            for output_modulus in (9, 95, 4096):
                check(cert.endpoint_mod(output_modulus) == stepwise_mod(cert, output_modulus),
                      "independent repeated-block reader")
            refill_checks += 1
        # Independent literal positivity test: prescribe the small exponent
        # first, then choose a positive odd hub in its exact fuel class.
        a = 2 * length * copies
        ternary = 3 ** (reds * copies + 1)
        axis = carry * pow(q - d, -1, ternary) % ternary
        hub = axis * pow(pow(2, a, ternary), -1, ternary) % ternary
        if hub % 2 == 0:
            hub += ternary
        cert = repeat_refuel(block, copies, hub)
        x = hub << cert.doublings
        for letter in cert.word:
            if letter == "D":
                x *= 2
            else:
                check(x % 3 == 2, "literal repetition guard")
                x = (2 * x - 1) // 3
            check(x > hub, "literal repetition first-hit height")
            refill_literal_edges += 1
        check(x % 3 != 0 and x % 2 == 1 and cert.endpoint_mod(x + 1) == x,
              "literal refuel endpoint")
    check((2 * 5 - 1) // 3 == 3, "exact one-digit fuel can end at ternary leaf")
    check((2 * 11 - 1) // 3 == 7 and (-11 - 1) % 9 != 0,
          "extra unit digit is sufficient, not necessary")
    print("Repetition fuel:", fuel_checks, "exact guards;", refill_checks,
          "refuel certificates;", refill_literal_edges, "literal positivity edges")

    # Solvability and group identities do not supply partial-map legality.
    check(word_mod("E", 5)[0] != 0 and (2 * 1 - 1) % 3 != 0, "E(1) modular/formal hostile")
    for n in range(1, 101):
        x = Fraction(n)
        for letter in ("T0", "E", "T0", "T1", "D", "D"):
            x = {"T0": lambda y: y / 2, "E": lambda y: (2 * y - 1) / 3,
                 "T1": lambda y: (3 * y + 1) / 2, "D": lambda y: 2 * y}[letter](x)
        check(x == n + 1, "formal translation identity")
        # If the first two letters are legal, E returns an odd integer,
        # which makes the following T0 illegal. This is an empty domain.
        first_two_legal = n % 2 == 0 and (n // 2) % 3 == 2
        check(not first_two_legal or ((n - 1) // 3) % 2 == 1,
              "formal translation has no guarded integer input")
    cert = next(iter(lift_all(5, 1).values()))
    bad = [InverseCertificate(True, 5, cert.target, cert.doublings, cert.word),
           InverseCertificate(3, 5, cert.target, cert.doublings, cert.word),
           InverseCertificate(1, 6, cert.target, cert.doublings, cert.word),
           InverseCertificate(1, 5, cert.target, 1, cert.word),
           InverseCertificate(1, 5, (True, 0), cert.doublings, cert.word),
           InverseCertificate(1, 5, cert.target, cert.doublings, "D")]
    for invalid in bad:
        try:
            invalid.verify()
        except ValueError:
            pass
        else:
            raise RuntimeError("malformed certificate accepted")
    print("Hostiles: E(1); formal +1 word has empty integer domain (100 controls); six bad certificates")
    print("Scope: constructed ancestors of supplied hubs; hub-home proof remains a premise")


if __name__ == "__main__":
    main()
