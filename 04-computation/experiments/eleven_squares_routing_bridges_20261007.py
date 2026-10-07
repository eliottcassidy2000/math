#!/usr/bin/env python3
"""Exact ordered guard routing; no orbit-convergence oracle or imported audits.

Run: python -B -O -X utf8 this_file.py
The companion note proves the general statements and specifies all scopes.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import permutations, product
from math import gcd


CHECKS = Counter()


def check(ok, label):
    CHECKS[label] += 1
    if not ok:
        raise RuntimeError(f"{label}: check {CHECKS[label]}")


def integer(x, *, positive=False, natural=False):
    if type(x) is not int or (positive and x <= 0) or (natural and x < 0):
        raise ValueError("exact integer outside the declared domain")


def maximum(a, b):
    return b if a is None else a if b is None else max(a, b)


def crt(a, m, b, n):
    """Intersection of two integer residue classes, including noncoprime moduli."""
    for x in (a, m, b, n):
        integer(x)
    if m <= 0 or n <= 0:
        raise ValueError("positive moduli required")
    g = gcd(m, n)
    if (b - a) % g:
        return None
    v = n // g
    t = 0 if v == 1 else ((b - a) // g * pow(m // g, -1, v)) % v
    modulus = m * v
    return (a + m * t) % modulus, modulus


@dataclass(frozen=True)
class Packet:
    """A partial increasing affine map plus a native guard and token summary.

    A packet is not a proof that the map is a Collatz operation. Such local
    semantics must accompany its leaf in the retained labelled proof tree.
    cut=None means -infinity; all finite cuts are strict input lower bounds.
    """
    p: int = 1
    q: int = 1
    b: int = 0
    residue: int = 0
    modulus: int = 1
    cut: object = None
    delta: int = 0
    need: int = 0
    empty: bool = False

    def __post_init__(self):
        for x in (self.p, self.q, self.b, self.residue, self.modulus,
                  self.delta, self.need):
            integer(x)
        if self.p <= 0 or self.q <= 0 or self.modulus <= 0 or self.need < 0:
            raise ValueError("positive slope/modulus and nonnegative credit requirement")
        if not 0 <= self.residue < self.modulus or type(self.empty) is not bool:
            raise ValueError("canonical guard and Boolean empty marker required")
        if self.cut is not None and type(self.cut) is not F:
            raise ValueError("cut must be an exact Fraction or None")
        if gcd(gcd(self.p, self.q), abs(self.b)) != 1:
            raise ValueError("canonical affine numerator/denominator required")
        if (self.p*self.residue+self.b) % self.q or (self.p*self.modulus) % self.q:
            raise ValueError("native progression must have integral outputs")
        if self.need < max(0, -self.delta):
            raise ValueError("credit requirement must fund the terminal balance")

    def at(self, n):
        integer(n)
        return F(self.p * n + self.b, self.q)

    def allows(self, n, credit=0):
        integer(n)
        integer(credit, natural=True)
        return (not self.empty and credit >= self.need
                and (n - self.residue) % self.modulus == 0
                and (self.cut is None or n > self.cut))

    def then(self, other):
        if type(other) is not Packet:
            raise ValueError("Packet required")
        if self.empty or other.empty:
            return EMPTY
        # Pull back the second native guard without cancelling its modulus.
        modulus = self.q * other.modulus
        rhs = self.q * other.residue - self.b
        g = gcd(self.p, modulus)
        if rhs % g:
            return EMPTY
        modulus //= g
        r = 0 if modulus == 1 else (rhs // g * pow(self.p // g, -1, modulus)) % modulus
        guard = crt(self.residue, self.modulus, r, modulus)
        if guard is None:
            return EMPTY
        pulled_cut = (None if other.cut is None
                      else (self.q * other.cut - self.b) / self.p)
        p = other.p * self.p
        q = other.q * self.q
        b = other.p * self.b + other.b * self.q
        div = gcd(gcd(p, q), abs(b))
        return Packet(p // div, q // div, b // div, *guard,
                      maximum(self.cut, pulled_cut),
                      self.delta + other.delta,
                      max(self.need, other.need - self.delta))


IDENTITY = Packet()
EMPTY = Packet(empty=True)


def valuation(a):
    """A strict pre-root ordinary odd Collatz step, with exact valuation a."""
    integer(a, positive=True)
    q = 1 << a
    modulus = 2 * q
    return Packet(3, q, 1, ((q - 1) * pow(3, -1, modulus)) % modulus,
                  modulus, F(1))


def compile_word(word):
    out = IDENTITY
    for a in word:
        out = out.then(valuation(a))
    return out


class PacketTree:
    """Balanced ordered storage; each edit costs <=ceil(log2(d)) combinations.

    No bit-cost bound is hidden here. Exact coefficients and moduli can grow.
    Leaves retain their distinct types; EMPTY does not discard those records.
    """
    def __init__(self, leaves):
        self.leaves = list(leaves)
        if any(type(x) is not Packet for x in self.leaves):
            raise ValueError("Packet leaves required")
        self.width = 1
        while self.width < len(self.leaves):
            self.width *= 2
        self.nodes = [IDENTITY] * (2 * self.width)
        self.nodes[self.width:self.width + len(self.leaves)] = self.leaves
        for i in range(self.width - 1, 0, -1):
            self.nodes[i] = self.nodes[2*i].then(self.nodes[2*i+1])

    @property
    def root(self):
        return self.nodes[1]

    def replace(self, index, packet):
        integer(index, natural=True)
        if index >= len(self.leaves) or type(packet) is not Packet:
            raise ValueError("valid leaf index and Packet required")
        self.leaves[index] = packet
        i = self.width + index
        self.nodes[i] = packet
        work = 0
        while i > 1:
            i //= 2
            self.nodes[i] = self.nodes[2*i].then(self.nodes[2*i+1])
            work += 1
        return work


NATIVE = {
    "H": Packet(729, 1024, 669, 155, 2048, F(0), 2, 0),
    "G": Packet(9, 8, 5, 11, 16, F(0), -1, 1),
    "A": Packet(81, 128, 85, 187, 256, F(0), 1, 0),
    "B": Packet(81, 128, 73, 7, 256, F(0), 1, 0),
    "L": Packet(9, 16, -3, 219, 256, F(0), 4, 0),
}


def direct_word(word, n):
    """Independent strict-root replay; no precomputed affine data."""
    x = n
    for a in word:
        if x <= 1 or x % 2 == 0:
            return None
        z = 3*x + 1
        count = 0
        while z % 2 == 0:
            z //= 2
            count += 1
        if count != a:
            return None
        x = z
    return x


def direct_native(word, n, credit):
    x = n
    for letter in word:
        item = NATIVE[letter]
        if not item.allows(x, credit):
            return None
        y = item.at(x)
        if y.denominator != 1:
            raise RuntimeError("native leaf failed its local integrality premise")
        x = int(y)
        credit += item.delta
    return x, credit


def reversal(d):
    names = list(range(d))
    swaps = [(i, d - 1 - i) for i in range(d // 2)]
    for i, j in swaps:
        names[i], names[j] = names[j], names[i]
    return names, swaps


def main():
    for d in range(33):
        names, swaps = reversal(d)
        check(names == list(reversed(range(d))), "labelled_axis_reversal")
        check(len(swaps) == d // 2, "labelled_axis_reversal")
        # Expose/process/restore, retaining named full payloads and masks.
        axes = [(j, j + 2, ("payload", j), ("mask", j)) for j in range(d)]
        original = axes[:]
        moves = 0
        visited = []
        for i in range(d):
            if i != d - 1:
                axes[i], axes[-1] = axes[-1], axes[i]
                moves += 1
            visited.append(axes[-1][0])
            if i != d - 1:
                axes[i], axes[-1] = axes[-1], axes[i]
                moves += 1
            check(axes == original, "expose_restore_payloads")
        check(visited == list(range(d)) and moves == 2*max(0, d-1), "expose_restore_payloads")

    for m in range(1, 10):
        for n in range(1, 10):
            for a in range(m):
                for b in range(n):
                    out = crt(a, m, b, n)
                    actual = [x for x in range(m*n) if x % m == a and x % n == b]
                    predicted = [] if out is None else list(range(out[0], m*n, out[1]))
                    check(actual == predicted, "generalized_CRT")

    generic = [Packet(p, 1, b, r, m, F(b, 2), delta, max(0, -delta))
               for p, b, r, m, delta in
               [(1, -2, 0, 2, 1), (2, 1, 1, 3, -1), (3, -1, 0, 3, 0),
                (2, -2, 1, 2, 2), (3, 2, 2, 4, -2), (1, 0, 2, 5, 0)]]
    generic += [valuation(1), valuation(3), IDENTITY, EMPTY]
    for a, b, c in product(generic, repeat=3):
        check(a.then(b).then(c) == a.then(b.then(c)), "generic_associativity")
    for a, b in product(generic, repeat=2):
        combined = a.then(b)
        for x in range(-5, 13):
            for credit in (0, 2, 4):
                direct = False
                if a.allows(x, credit):
                    y = a.at(x)
                    check(y.denominator == 1, "generic_integrality")
                    direct = b.allows(int(y), credit+a.delta)
                check(combined.allows(x, credit) == direct, "generic_pullback_semantics")

    words = [()]
    for length in range(1, 6):
        words.extend(product(range(1, 5), repeat=length))
    for word in words:
        packet = compile_word(word)
        tree = PacketTree([valuation(a) for a in word])
        check(packet == tree.root, "ordered_reassociation")
        if word:
            q = 1 << sum(word)
            check(packet.modulus == 2*q, "full_native_word_guard")
            for t in (0, 1, 3):
                n = packet.residue + packet.modulus*t
                replay = direct_word(word, n)
                check(packet.allows(n) == (replay is not None), "literal_strict_root_replay")
                if replay is not None:
                    check(packet.at(n) == replay, "literal_strict_root_replay")
            # A consecutive wrong odd residue is not repaired by the same slope.
            n = packet.residue + 2
            check(packet.allows(n) == (direct_word(word, n) is not None), "off_guard")
            for i in (0, len(word)-1):
                changed = list(word)
                changed[i] += 1
                work = tree.replace(i, valuation(changed[i]))
                check(tree.root == compile_word(changed), "source_guard_updates")
                check(work == tree.width.bit_length()-1, "source_guard_updates")
                tree.replace(i, valuation(word[i]))

    for length in range(5):
        for word in product(NATIVE, repeat=length):
            packet = PacketTree([NATIVE[x] for x in word]).root
            seq = IDENTITY
            for x in word:
                seq = seq.then(NATIVE[x])
            check(packet == seq, "native_tree_reassociation")
            # Includes empty-guard and unfunded controls, and unbounded-height lifts.
            seeds = [1, 27, 155, 219]
            if not packet.empty:
                seeds += [packet.residue + packet.modulus*t for t in (0, 1, 3)]
            for n in seeds:
                for credit in (0, packet.need, packet.need+2):
                    out = direct_native(word, n, credit)
                    check(packet.allows(n, credit) == (out is not None), "native_guard_and_budget")
                    if out is not None:
                        check(out == (packet.at(n), credit+packet.delta), "native_guard_and_budget")

    p12, p21 = compile_word((1, 2)), compile_word((2, 1))
    check((p12.p, p12.q, p12.b, p12.residue, p12.modulus) == (9, 8, 5, 11, 16), "order_hostile")
    check((p21.p, p21.q, p21.b, p21.residue, p21.modulus) == (9, 8, 7, 9, 16), "order_hostile")
    padded = compile_word((4, 2))
    check(padded.residue == 5 and padded.cut == 5 and not padded.allows(5), "root_padding_hostile")
    check(direct_word((4,), 5) == 1 and direct_word((4, 2), 133) == 19, "root_padding_hostile")
    check(NATIVE["L"].at(27) == 15 and not NATIVE["L"].allows(27), "native_guard_cancellation")
    check(NATIVE["L"].then(NATIVE["B"]).empty, "native_guard_cancellation")
    lg = NATIVE["L"].then(NATIVE["G"])
    check((lg.p, lg.q, lg.b, lg.residue, lg.modulus, lg.delta, lg.need)
          == (81, 128, 53, 219, 256, 3, 0), "native_guard_cancellation")

    for a in range(7):
        center = -F(1, 3-(1 << a))
        for x in range(-5, 6):
            check(F(3*x+1, 1 << a)-center == F(3, 1 << a)*(x-center),
                  "changing_fixed_point_charts")
    w11 = compile_word((1, 1))
    check((w11.p, w11.q, w11.b, w11.residue, w11.modulus) == (9, 4, 5, 7, 8),
          "changing_fixed_point_charts")

    # Parent's separately discovered block replacement, not an unguarded swap.
    w = tuple(map(int, "111112223113"))
    v = tuple(map(int, "121313121121"))
    pw, pv = compile_word(w), compile_word(v)
    check((pw.p, pw.q) == (pv.p, pv.q) and pv.b-pw.b == 4*pw.p, "whole_block_replacement")
    check(pw.residue == 257727 and pw.modulus == 1048576, "whole_block_replacement")
    for t in range(16):
        n = 257727 + 1048576*t
        check(pw.allows(n) and pv.allows(n-4), "whole_block_replacement")
        check(direct_word(w, n) == direct_word(v, n-4), "whole_block_replacement")
        x = n
        for a in w:
            x = (3*x+1) // (1 << a)
            check(x > n, "whole_block_growth")

    # Nontrivial congruence pullback with a nonunit coefficient.
    f = Packet(3, 1, 0, 0, 1)
    check(f.then(Packet(1, 1, 0, 1, 3)).empty, "nonunit_pullback")
    check(f.then(Packet(1, 1, 0, 0, 9)).modulus == 3, "nonunit_pullback")

    malformed = [lambda: valuation(True), lambda: valuation(1.0), lambda: valuation(0),
                 lambda: Packet(p=True), lambda: Packet(cut=1.0),
                 lambda: Packet(residue=1), lambda: Packet(p=2, q=2, b=0),
                 lambda: PacketTree([1]), lambda: PacketTree([IDENTITY]).replace(True, IDENTITY),
                 lambda: IDENTITY.allows(True), lambda: IDENTITY.allows(1, -1),
                 lambda: Packet(p=1, q=2), lambda: Packet(delta=-1)]
    for operation in malformed:
        try:
            operation()
        except ValueError:
            check(True, "malformed_exact_inputs")
        else:
            check(False, "malformed_exact_inputs")

    print("Exact labelled routing and guarded ordered packet compiler")
    print("Imported multiplication and square-packing theorems: accepted, not re-audited")
    print("Axis dimensions: 0..32; strict valuation words: letters 1..4, lengths 0..5")
    print("Native H/G/A/B/L programs: lengths 0..4; native guards and token budgets retained")
    print("12 versus 21: same (P,Q)=(9,8), carries 5/7, guards 11/9 mod16")
    print("42: root-padding source5 rejected by exact strict cut5; source133 ->19 accepted")
    print("LB empty; LG: (81x+53)/128 on219 mod256, delta3, need0")
    print("Whole-block replacement: 257727+1048576t joins its smaller child n-4; t0..15 checked")
    print("Tree updates: logarithmically many packet combinations; no unit-bit-cost or coverage claim")
    for key, value in sorted(CHECKS.items()):
        print(f"{key}: {value}")
    print(f"TOTAL: {sum(CHECKS.values())}")


if __name__ == "__main__":
    main()
