"""Search both binary words of the physical six-vertex core at a prescribed lower count."""

import argparse
from fractions import Fraction
from check import p, perfect_shape, two_slots, word_count
from replay import inverse_one_sided


def addresses(target, max_width=None):
    """Yield every (v,y) producing the specified lower endpoint within the width limit."""
    if target <= 0 or target % 2 == 0:
        return
    cap = (target + 10) // 14

    def padded(u, slots, tag_complement):
        h = 25 * (1 + p(u))
        t = int(u == 0)
        while (slots << t) <= cap:
            v = ((u + 1) << t) - 1
            slope = (25 * (u + 1) + 28 * slots) << t
            for parity in (0, 1):
                k = h + slope - 28 * tag_complement - 10 - 25 * (v % 2 == parity == 0)
                if parity == 0 and target == k:
                    yield v, 0
                for y in inverse_one_sided(target, k, h, slope, parity):
                    assert word_count(v, y) == target
                    yield v, y
            t += 1

    yield from padded(0, 2, 1)
    width, least_slots = 2, 2
    while least_slots <= cap and (max_width is None or width <= max_width):
        for u in range(1 << (width - 1), 1 << width, 2):
            slots = two_slots(u)
            if slots <= cap:
                complement = (1 << width) - int(format(u, 'b')[::-1], 2)
                yield from padded(u, slots, complement)
        width += 1
        least_slots *= 3


def density_bound(width=10):
    grid = 10**9
    rounded = 0
    for b in range(2, width + 1):
        for u in range(1 << (b-1), 1 << b, 2):
            c = 25*(u+1) + 28*two_slots(u)
            rounded += (2*grid + c-1) // c
    coefficients = Fraction(1, 81) + Fraction(rounded, grid) + Fraction(3, 28) * Fraction(2, 3)**(width-1)
    return 4 * coefficients


def checks():
    for b in range(2, 9):
        for u in range(1 << (b-1), 1 << b, 2):
            slots = two_slots(u)
            assert 2*3**(b-2) <= slots <= 2*3**b
            assert two_slots(2*u+1) == 2*slots
    for v in range(1, 9):
        for y in range(12):
            N = word_count(v, y)
            assert (v, y) in addresses(N)
    # Entire small inverse fibers compared with the original two-word formula.
    brute = {}
    for v in range(1, 32):
        for y in range(32):
            value = word_count(v, y)
            if value <= 1500:
                brute.setdefault(value, set()).add((v, y))
    # L+10 >= 53(v+1)(y+1) rules out any further v or y here.
    for N in range(1, 1501, 2):
        assert set(addresses(N)) == brute.get(N, set())
    bound = density_bound()
    assert bound < Fraction(133, 1000)
    print('Upper density among odd integers using either endpoint is below 0.133.')
    print('Checked all inverse fibers through 1500 and two-slot bounds through width 8.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--lower', type=int, help='prescribed lower (odd) endpoint')
    group.add_argument('--exponent', type=int, help='search for both neighbors of the Euclid center')
    parser.add_argument('--max-width', type=int, help='cap the normalized even prefix, leaving all tail lengths free')
    args = parser.parse_args()
    if args.lower is None and args.exponent is None:
        checks()
    else:
        if args.exponent is not None and args.exponent < 2:
            parser.error('the exponent must be at least two')
        N = args.lower if args.lower is not None else perfect_shape(args.exponent)-1
        if N < 1 or N % 2 == 0:
            parser.error('the lower count must be positive and odd')
        matches = 0
        for v, y in addresses(N, args.max_width):
            print(f'lower {N}, upper {N+2}, receiver {v}, right address {y}', flush=True)
            matches += 1
        scope = 'full two-word family' if args.max_width is None else f'normalized width through {args.max_width}'
        print(f'{matches} solutions in the {scope}')
