"""Exact affine ports, adaptive cylinder cover, and unbounded root schemas.

No external packages, floating point, or convergence assumption.
"""

from dataclasses import dataclass


def check(test, label):
    if not test:
        raise RuntimeError(label)


def step(n, sign=1):
    return (3 * n + sign) // 2 if n % 2 else n // 2


@dataclass(frozen=True)
class Port:
    word: tuple
    sign: int = 1

    @property
    def length(self):
        return len(self.word)

    @property
    def slope_numerator(self):
        return 3 ** sum(self.word)

    @property
    def modulus(self):
        return 2 ** self.length

    @property
    def carry(self):
        carry = 0
        for j, bit in enumerate(self.word):
            carry = 3 ** bit * carry + bit * 2 ** j
        return carry

    @property
    def residue(self):
        q = self.modulus
        return (-self.sign * self.carry * pow(self.slope_numerator, -1, q)) % q

    @property
    def intercept(self):
        return (self.slope_numerator * self.residue + self.sign * self.carry) // self.modulus

    @property
    def first_descending_lift(self):
        gap = self.modulus - self.slope_numerator
        check(gap > 0, "port is not contractive")
        return max(0, (self.intercept - self.residue) // gap + 1)

    def input(self, k):
        return self.residue + self.modulus * k

    def output(self, k):
        return self.intercept + self.slope_numerator * k


def adaptive_cover(depth, sign):
    active = [Port((), sign)]
    leaves = []
    for _ in range(depth):
        next_active = []
        for node in active:
            for bit in (0, 1):
                child = Port(node.word + (bit,), sign)
                if child.slope_numerator < child.modulus:
                    leaves.append(child)
                else:
                    next_active.append(child)
        active = next_active
    return leaves, active


def direct_word(n, length, sign=1):
    word = []
    for _ in range(length):
        word.append(n % 2)
        n = step(n, sign)
    return tuple(word), n


def inverse_rise_port(target, a, b, sign=1):
    check(target > 0 and target % 2 and target % 3, "target must be odd and prime to3")
    check(a >= 2 and b >= a, "strict-rank inverse constructor")
    check(sign in (1, -1), "sign domain")
    numerator = 2 ** b * target + sign
    check(numerator % (3 ** a) == 0, "3-adic port matching")
    u = numerator // (3 ** a)
    n = 2 ** a * u - sign
    check(u % 2 == 1 and n > target, "oddness and integer rank")
    return n


def main():
    for sign in (1, -1):
        leaves, residual = adaptive_cover(14, sign)
        positive_fringe = set()
        lift_tests = 0
        for port in leaves + residual:
            p, q = port.slope_numerator, port.modulus
            independent = sum(bit * 2 ** j * 3 ** sum(port.word[j + 1:])
                              for j, bit in enumerate(port.word))
            check(port.carry == independent, "carry independent sum")
            for k in (0, 1, 2, 17):
                n = port.input(k)
                word, out = direct_word(n, port.length, sign)
                check(word == port.word and out == port.output(k), "port direct orbit")
                if q > p:
                    check((out < n) == (k >= port.first_descending_lift), "exact threshold")
                lift_tests += 1
            if q > p:
                for k in range(port.first_descending_lift):
                    if port.input(k) > 0:
                        positive_fringe.add(port.input(k))
        # Independent native-residue membership, not the word-tree recurrence.
        native_residual = []
        for r in range(2 ** 14):
            n, power = r, 1
            has_slope_exit = False
            for j in range(1, 15):
                if n % 2:
                    power *= 3
                n = step(n, sign)
                if power < 2 ** j:
                    has_slope_exit = True
                    break
            if not has_slope_exit:
                native_residual.append(r)
        check(sorted(p.residue for p in residual) == native_residual, "independent residual classes")
        cover_mass = sum(2 ** (14 - p.length) for p in leaves)
        check(cover_mass + len(residual) == 2 ** 14, "exact disjoint cover")
        print("ADAPTIVE", sign, "depth14", "ports", len(leaves),
              "residual_classes", len(residual), "positive_fringe", sorted(positive_fringe),
              "lift_controls", lift_tests)
        print("residual_first16", native_residual[:16])

    # The counter grammar works at arbitrary a; these are finite controls.
    diagonal = 0
    for a in range(1, 65):
        port = Port((1,) * a + (0,) * a)
        r = pow(3 ** a, -1, 2 ** a)
        check(port.residue == 2 ** a * r - 1, "balanced-cylinder residue")
        expected_threshold = 1 if a == 1 else 0
        check(port.first_descending_lift == expected_threshold, "balanced threshold")
        for k in (0, 1, 3, 100):
            n, out = port.input(k), port.output(k)
            _, actual = direct_word(n, 2 * a)
            check(out == actual, "balanced actual output")
            check(out < n or n == out == 1, "balanced descent")
            diagonal += 1
    print("balanced_counter_grammar_controls", diagonal, "a1..64_lifts0,1,3,100")

    # A parameterized complete root proof, without searching an orbit.
    for sign in (1, -1):
        roots = []
        for a in range(2, 10):
            for t in range(4):
                b = 3 ** (a - 1) * ((2 * t + 1) if sign == 1 else 2 * (t + 1))
                n = inverse_rise_port(1, a, b, sign)
                # Only a odd steps are directly iterated; all b halvings are
                # independently checked as a power-of-two identity.
                word, peak = direct_word(n, a, sign)
                check(word == (1,) * a and peak == 2 ** b, "complete root family")
                check(peak >> (b - 2) == 4, "compressed halving to root4")
                roots.append((a, t, b, n.bit_length()))
        print("complete_root_schema_controls", sign, len(roots), "a2..9_t0..3")
        for a in range(2, 5):
            b = 3 ** (a - 1) * (1 if sign == 1 else 2)
            n = inverse_rise_port(1, a, b, sign)
            print("root_example", sign, "a", a, "b", b, "n", n, "steps_to4", a + b - 2)
        # A root certificate beyond this finite cover's lookahead, represented
        # by counters rather than millions of individual halving steps.
        a = 15
        b = 3 ** (a - 1) * (1 if sign == 1 else 2)
        n = inverse_rise_port(1, a, b, sign)
        word, after_cut = direct_word(n, 14, sign)
        check(word == (1,) * 14 and after_cut > n, "root source survives depth14 cutoff")
        _, peak = direct_word(n, a, sign)
        check(peak == 2 ** b and peak >> (b - 2) == 4, "root beyond cutoff")
        print("root_beyond_cutoff", sign, "a", a, "b", b,
              "source_bits", n.bit_length(), "residue_mod16384", n % 16384)

    # Exhaust all exponent phases for each modulus, including failures.
    phases = 0
    for a in range(1, 9):
        modulus, period = 3 ** a, 2 * 3 ** (a - 1)
        check(len({pow(2, b, modulus) for b in range(period)}) == period,
              "2 generates every unit modulo3^a")
        for b in range(1, 2 * period + 1):
            check(((pow(2, b, modulus) + 1) % modulus == 0)
                  == (b % period == period // 2), "plus integrality iff phase")
            check(((pow(2, b, modulus) - 1) % modulus == 0)
                  == (b % period == 0), "minus integrality iff phase")
            phases += 2
    print("integrality_iff_controls", phases, "a1..8_b1..two_periods_both_signs")

    # Port matching to any odd target coprime to3, followed by its supplied proof.
    matches = 0
    for sign in (1, -1):
        for target in range(1, 102, 2):
            if target % 3 == 0:
                continue
            for a in range(2, 7):
                modulus, period = 3 ** a, 2 * 3 ** (a - 1)
                values = [b for b in range(period)
                          if (pow(2, b, modulus) * target + sign) % modulus == 0]
                check(len(values) == 1, "unique port phase")
                b = values[0]
                while b < a:
                    b += period
                n = inverse_rise_port(target, a, b, sign)
                _, peak = direct_word(n, a, sign)
                check(peak == 2 ** b * target, "generic matched port")
                matches += 1
    print("generic_inverse_port_phase_controls", matches, "targets<=101_a2..6_both_signs")

    # Exact sign-sensitive controls and bounded-lookahead hostility.
    for a in range(1, 65):
        n = 2 ** (a + 1) - 1
        for j in range(1, a + 1):
            check(direct_word(n, j)[1] == 3 ** j * 2 ** (a + 1 - j) - 1 > n,
                  "unbounded growing prefix")
    # Derive the second cycle's shortcut word independently if convention differs.
    for start in (5, 17):
        n, word = start, []
        while not word or n != start:
            check(len(word) < 100, "minus control cycle length")
            word.append(n % 2)
            n = step(n, -1)
        port = Port(tuple(word), -1)
        check(port.slope_numerator > port.modulus, "minus cycles have expanding slope")
        check(port.slope_numerator * start - port.carry == port.modulus * start,
              "negative carry exactly compensates")
        print("minus_cycle_hostile", start, "word", ''.join(map(str, word)),
              "p", port.slope_numerator, "q", port.modulus, "carry", port.carry)
    print("ALL CHECKS PASS; residual remains explicit and OPEN")


if __name__ == "__main__":
    main()
