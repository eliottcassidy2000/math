"""Exact clock-pair composition and periodic-tail certificates.

No search in this file assumes universal Collatz convergence. Universes are
all 288 maps on labelled sets of sizes 1..4, clocks 0..6, all 15,625 triples
of clock pairs in [0,4]^2, and complete Terras period words through 12.
"""
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import gcd


def need(condition, message):
    if not condition:
        raise ValueError(message)


def integer(n):
    need(type(n) is int, "exact integer required")
    return n


def clock(n):
    integer(n)
    need(n >= 0, "nonnegative clock required")
    return n


def iterate(step, x, n):
    integer(x)
    clock(n)
    for _ in range(n):
        x = step(x)
    return x


def clock_product(left, right):
    need(type(left) is tuple and len(left) == 2, "clock pair required")
    need(type(right) is tuple and len(right) == 2, "clock pair required")
    a, b = map(clock, left)
    c, d = map(clock, right)
    return a + max(c-b, 0), d + max(b-c, 0)


@dataclass(frozen=True)
class Receipt:
    x: int
    y: int
    a: int
    b: int


def verify(step, receipt):
    need(type(receipt) is Receipt, "typed receipt required")
    integer(receipt.x)
    integer(receipt.y)
    clock(receipt.a)
    clock(receipt.b)
    need(iterate(step, receipt.x, receipt.a) ==
         iterate(step, receipt.y, receipt.b), "false endpoint equality")
    return receipt


def reverse(step, receipt):
    verify(step, receipt)
    return Receipt(receipt.y, receipt.x, receipt.b, receipt.a)


def compose(step, left, right):
    verify(step, left)
    verify(step, right)
    need(left.y == right.x, "source seam mismatch")
    a, b = clock_product((left.a, left.b), (right.a, right.b))
    result = Receipt(left.x, right.y, a, b)
    return verify(step, result)


def close_walk(step, receipts):
    need(type(receipts) is tuple and bool(receipts), "nonempty receipt tuple required")
    total = verify(step, receipts[0])
    for receipt in receipts[1:]:
        total = compose(step, total, receipt)
    need(total.x == total.y, "walk is not closed")
    return total


def tail_period(step, loops):
    need(type(loops) is tuple and bool(loops), "nonempty loop tuple required")
    source = loops[0].x if type(loops[0]) is Receipt else None
    integer(source)
    periods, entries = [], []
    for loop in loops:
        verify(step, loop)
        need(loop.x == source == loop.y, "loops require one base source")
        if loop.a != loop.b:
            periods.append(abs(loop.a-loop.b))
            entries.append(min(loop.a, loop.b))
    need(bool(periods), "zero-charge loops do not certify a period")
    period, entry = 0, max(entries)
    for p in periods:
        period = gcd(period, p)
    z = iterate(step, source, entry)
    need(iterate(step, z, period) == z, "period consequence failed")
    return entry, period, z


def terras(n):
    integer(n)
    return (3*n+1)//2 if n % 2 else n//2


def root_from_loops(loops):
    need(type(loops) is tuple and bool(loops), "nonempty loops required")
    for loop in loops:
        need(type(loop) is Receipt, "typed receipt required")
        integer(loop.x)
        integer(loop.y)
        need(loop.x > 0 and loop.y > 0, "positive Terras sources required")
    entry, period, z = tail_period(terras, loops)
    need(period == 2, "this root rule requires gcd exactly two")
    need(z in (1, 2), "positive two-period classification failed")
    path = [loops[0].x]
    while path[-1] != 1 and len(path) <= entry + 1:
        path.append(terras(path[-1]))
    need(path[-1] == 1, "root extraction exceeded proved deadline")
    return tuple(path)


def periodic_points(length):
    """All positive T-points fixed by T^length, using all parity words."""
    clock(length)
    need(length > 0, "positive period required")
    points = set()
    for bits in product((0, 1), repeat=length):
        p, c, q = 1, 0, 1
        for bit in bits:
            if bit:
                p, c = 3*p, 3*c+q
            q *= 2
        if q <= p or c == 0:
            continue
        candidate = Fraction(c, q-p)
        if candidate.denominator != 1:
            continue
        n = candidate.numerator
        z = n
        for bit in bits:
            if z % 2 != bit:
                break
            z = terras(z)
        else:
            need(z == n, "affine periodic candidate failed replay")
            points.add(n)
    return tuple(sorted(points))


def main():
    checks = 0

    def check(ok, message):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(message)

    # Associativity and an independent partial-translation representation.
    pairs = list(product(range(5), repeat=2))
    for r, s, t in product(pairs, repeat=3):
        check(clock_product(clock_product(r, s), t) ==
              clock_product(r, clock_product(s, t)), "associativity")
    for r, s in product(pairs, repeat=2):
        a, b = clock_product(r, s)
        check(a-b == r[0]-r[1]+s[0]-s[1], "additive charge")
        for n in range(17):
            first = None if n < s[1] else n+s[0]-s[1]
            literal = None if first is None or first < r[1] else first+r[0]-r[1]
            normalized = None if n < b else n+a-b
            check(literal == normalized, "independent tail-domain composition")

    maps, loop_pairs = 0, 0
    for size in range(1, 5):
        for values in product(range(size), repeat=size):
            maps += 1
            step = lambda x, values=values: values[x]
            for x in range(size):
                loops = [Receipt(x, x, a, b) for a in range(7) for b in range(a)
                         if iterate(step, x, a) == iterate(step, x, b)]
                for left, right in product(loops, repeat=2):
                    entry, p, z = tail_period(step, (left, right))
                    # Independently find the first return of the extracted z.
                    zz, first = step(z), 1
                    while zz != z and first <= size:
                        zz, first = step(zz), first+1
                    check(zz == z and p % first == 0, "gcd is a certified tail period")
                    loop_pairs += 1
                # An actual edge followed by a delayed equality, for all y.
                y = step(x)
                r = Receipt(x, y, 1, 0)
                s = Receipt(y, step(y), 1, 0)
                check(compose(step, r, s) == Receipt(x, step(y), 2, 0), "literal chain")

    # Zero lag cannot erase the entry barrier of a genuine merge.
    delayed = verify(terras, Receipt(5, 4, 3, 3))
    check(delayed.x != delayed.y, "zero charge does not mean equal sources")
    neutral = close_walk(terras, (delayed, reverse(terras, delayed)))
    check((neutral.a, neutral.b) == (3, 3), "inverse keeps entry idempotent")
    # A zero-charge closed proof also exists for the nonperiodic map n -> n+1.
    shift = lambda n: n+1
    check(verify(shift, Receipt(0, 0, 4, 4)).a == 4, "zero loop has no periodic consequence")

    # Root-free syntax: two verified closed clock words, no root label supplied.
    root_rows = []
    for n in (1, 2, 3, 27, 223, 233, 322):
        x, first = n, 0
        while x != 1 and first < 10000:
            x, first = terras(x), first+1
        check(x == 1, "declared finite source has a root word")
        loops = (Receipt(n, n, first+4, first), Receipt(n, n, first+6, first))
        path = root_from_loops(loops)
        check(len(path)-1 == first and path[-1] == 1, "root readout from gcd(4,6)")
        root_rows.append((n, first))

    # Exact common-future joins for the owner's three integers.
    join = compose(terras, Receipt(223, 233, 7, 15), Receipt(233, 322, 0, 10))
    check(join == Receipt(223, 322, 7, 25), "triplet clock composition")
    check(iterate(terras, 223, 7) == iterate(terras, 322, 25) == 425,
          "triplet actual endpoint")

    census = []
    for length in range(1, 13):
        points = periodic_points(length)
        check(points == ((1, 2) if length % 2 == 0 else ()), "complete short-period census")
        census.append((length, points))

    hostiles = [
        lambda: clock_product((True, 0), (0, 0)),
        lambda: clock_product((1.0, 0), (0, 0)),
        lambda: clock_product((-1, 0), (0, 0)),
        lambda: verify(terras, Receipt(5, 4, 0, 0)),
        lambda: compose(terras, Receipt(5, 4, 3, 3), Receipt(2, 1, 1, 0)),
        lambda: tail_period(shift, (Receipt(0, 0, 4, 4),)),
        lambda: root_from_loops((Receipt(-1, -1, 4, 0), Receipt(-1, -1, 6, 0))),
        lambda: root_from_loops((Receipt(0, 0, 4, 0), Receipt(0, 0, 6, 0))),
        lambda: root_from_loops((Receipt(1, 1, 4, 0),)),
        lambda: root_from_loops(()),
    ]
    for action in hostiles:
        try:
            action()
        except ValueError:
            check(True, "hostile rejected")
        else:
            check(False, "hostile accepted")

    print("PROVED: exact common-future clock composition and nonzero-charge period extraction")
    print("clock_triples=15625; tail_compositions=625; test_tail_points_each=17")
    print(f"complete_finite_maps={maps}; checked_loop_pairs={loop_pairs}")
    print("triplet: T^7(223)=T^15(233)=T^25(322)=425; closed triangle has charge zero")
    print("finite_root_rows=" + repr(root_rows))
    print("complete_period_words_lengths1..12=8190; positive fixed points are(1,2) exactly at even lengths")
    print("gcd(4,6)=2 certifies ROOT from verified closed receipts; no universal receipt existence")
    print(f"typed_guard_hostiles={len(hostiles)}; checks={checks}")


if __name__ == "__main__":
    main()
