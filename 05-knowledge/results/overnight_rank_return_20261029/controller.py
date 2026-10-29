"""Exact certificates for the (d,r)=(3,3) chronological Collatz parser.

Run from any directory. No asserts: checks remain live under python -O.
The finite sweep supplements the all-integer inequalities in CONTROLLER.md.
"""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def step(n):
    t = 3 * n + 1
    a = (t & -t).bit_length() - 1
    return t >> a, a


def f(n):
    return n ** 3 + 5 * (n == 3)


def block(n):
    source, q, word, values = n, 1, [], [n]
    while q:
        target, a = step(n)
        nq = q - 1 + 3 * (a % 3 == 0)
        require(128 * 4 ** nq * f(target) <= 125 * 4 ** q * f(n),
                ('potential', n, q, a))
        word.append(a)
        values.append(target)
        n, q = target, nq
    return n, tuple(word), values


def affine(word):
    c = d = 0
    for a in word:
        c = 3 * c + 2 ** d
        d += a
    j = len(word)
    mod = 2 ** (d + 1)
    residue = ((2 ** d - c) * pow(3 ** j, -1, mod)) % mod
    return j, d, c, residue, mod


WORDS = [(1,), (3, 1, 1, 1), (3, 3, 1, 1, 1, 1, 1),
         (3, 1, 3, 1, 1, 1, 1), (3, 1, 1, 3, 1, 1, 1)]


def failure(n):
    return any(n % affine(w)[4] == affine(w)[3] for w in WORDS)


def first_reference_descent(n):
    """Finite-height experiment only; no claim of universal termination."""
    source, word, cuts = n, [], []
    for _ in range(10000):
        n, w, _ = block(n)
        word.extend(w)
        cuts.append(len(word))
        if n < source:
            return tuple(word), tuple(cuts)
    raise ValueError(('controller budget exhausted', source))


def guarded_cell(word, cuts):
    """Derive exact integer interval/cylinder without simulating a source."""
    low, high = 3, None
    q = 1
    actual_cuts = []
    for k, a in enumerate(word, 1):
        require(a >= 1, 'positive letters')
        q += 3 * (a % 3 == 0) - 1
        if q == 0:
            actual_cuts.append(k)
            q = 1
    require(tuple(actual_cuts) == cuts, 'parser cuts')
    require(cuts[-1] == len(word), 'terminal cut')
    for cut in cuts:
        j, d, c, _, _ = affine(word[:cut])
        gap = 2 ** d - 3 ** j
        if cut == cuts[-1]:
            require(gap > 0, 'positive terminal descent gap')
            low = max(low, c // gap + 1)
        elif gap > 0:
            bound = c // gap
            high = bound if high is None else min(high, bound)
    j, d, c, residue, mod = affine(word)
    return dict(word=list(word), cuts=list(cuts), lower=low, upper=high,
                residue=residue, modulus=mod, odd_clock=j, binary_clock=d,
                carry=c)


def receipt(height):
    require(height >= 3, 'height at least three')
    # Infinite-range inequalities; the remaining boundary is explicit.
    require(Fraction(1331, 1372) < Fraction(125, 128), 'n>=7 bound')
    rho = Fraction(31250, 35937)
    require(4 * rho ** 10 < 1, 'ten-step descent above 31')
    require(Fraction(358, 47) < 13, 'one-branch nonminimal bound')
    require(Fraction(32354, 1909) < 19, 'two-branch nonminimal bound')
    for n in range(1, 19, 2):
        target, w, _ = block(n)
        require(sum(a % 3 == 0 for a in w) < 2, 'two-branch minimum')
    small = set()
    for n in range(1, 32, 2):
        while n not in small:
            small.add(n)
            n, _ = step(n)
    require(max(small) == 3077 and len(small) == 53, 'small invariant set')
    require(all(step(n)[0] in small for n in small), 'small-set closure')
    critical = []
    for n in range(1, 3078, 2):
        target, w, _ = block(n)
        m = sum(a % 3 == 0 for a in w)
        require(m < 3 or target < n, ('small boundary', n))
        critical.append([n, target, len(w), m])
    # Independent word arithmetic verifies the claimed cylinders, including
    # their exact valuations rather than just their endpoint parity.
    cylinders = []
    for w in WORDS:
        j, d, c, r, mod = affine(w)
        require(3 ** j > 2 ** d, 'expanding failure branch')
        for offset in range(5):
            n = r + offset * mod
            target, got, _ = block(n)
            require(got == w and 2 ** d * target == 3 ** j * n + c,
                    ('cylinder', n))
        cylinders.append(dict(word=list(w), residue=r, modulus=mod,
                              numerator=3 ** j, carry=c, denominator=2 ** d))
    mass = sum(Fraction(2, row['modulus']) for row in cylinders)
    require(mass == Fraction(1059, 2048), 'failure mass')
    failures, max_steps, max_source = {}, 0, None
    for n in range(3, height + 1, 2):
        target, w, _ = block(n)
        require((target >= n) == failure(n), ('failure classifier', n))
        if len(w) > max_steps:
            max_steps, max_source = len(w), n
        if target >= n:
            m = sum(a % 3 == 0 for a in w)
            failures[m] = failures.get(m, 0) + 1
    # Build exact source-preserving first-reference-descent cells and
    # independently cover the finite universe by their arithmetic guards.
    bank = {}
    control_height = min(height, 10001)
    for n in range(3, control_height + 1, 2):
        word, cuts = first_reference_descent(n)
        bank[word] = guarded_cell(word, cuts)
    cells = list(bank.values())
    owners = {}
    for index, cell in enumerate(cells):
        end = min(control_height, cell['upper'] or control_height)
        low, mod, r = cell['lower'], cell['modulus'], cell['residue']
        start = low + (r - low) % mod
        for n in range(start, end + 1, mod):
            require(n not in owners, ('overlapping source guards', n))
            owners[n] = index
    require(set(owners) == set(range(3, control_height + 1, 2)),
            'finite source cover')
    return dict(status='PROVED scoped inequalities / FINITE-EXACT certificates',
                height=height, failure_cylinders=cylinders,
                failure_odd_haar_mass=str(mass), small_closed_set=sorted(small),
                finite_boundary=critical, failure_histogram=failures,
                longest_block=dict(steps=max_steps, source=max_source),
                controller_height=control_height, controller_cells=cells,
                scope='Single-block termination is universal; repeated-reference '
                      'descent is certified only through controller_height.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--height', type=int, default=100001)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    data = receipt(args.height)
    encoded = json.dumps(data, indent=2, sort_keys=True) + '\n'
    path = Path(__file__).with_suffix('.json')
    if args.write:
        path.write_text(encoded, encoding='utf-8')
    elif path.exists():
        require(json.loads(path.read_text(encoding='utf-8')) == json.loads(encoded),
                'stored receipt differs; use same height')
    print(json.dumps(dict(height=args.height, sha256=hashlib.sha256(
        encoded.encode()).hexdigest(), failure_histogram=data['failure_histogram'],
        longest_block=data['longest_block'], controller_cells=len(data['controller_cells']),
        controller_height=data['controller_height']), sort_keys=True))


if __name__ == '__main__':
    main()
