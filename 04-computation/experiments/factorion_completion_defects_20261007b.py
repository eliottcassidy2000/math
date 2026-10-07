"""Finite global certificates for digit-factorial maps and honest completion defects.

The digit histogram is sufficient for the next integer, not for the source.
An independent full integer core checks the compressed cycle portrait.
The modified all-to-1 dynamics mark exactly which edges were changed.
"""
from array import array
from dataclasses import dataclass, replace
from itertools import combinations_with_replacement
from math import factorial, comb

CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def natural(n, lower=0):
    need(type(n) is int and n >= lower, 'exact integer with stated lower bound')


def digits(n, base):
    natural(n, 1)
    natural(base, 2)
    out = []
    while n:
        n, r = divmod(n, base)
        out.append(r)
    return tuple(reversed(out))


def histogram(n, base):
    counts = [0]*base
    for d in digits(n, base):
        counts[d] += 1
    return tuple(counts)


def value(counts):
    need(type(counts) is tuple and len(counts) >= 2, 'typed base histogram')
    need(all(type(c) is int and c >= 0 for c in counts), 'nonnegative exact counts')
    need(sum(counts[1:]) > 0, 'canonical positive numeral has a nonzero digit')
    return sum(c*factorial(d) for d, c in enumerate(counts))


def digit_map(n, base):
    return value(histogram(n, base))


def canonical_cycle(cycle):
    need(bool(cycle) and len(set(cycle)) == len(cycle), 'nonempty simple marked cycle')
    index = cycle.index(min(cycle))
    return tuple(cycle[index:]+cycle[:index])


@dataclass(frozen=True)
class Portrait:
    base: int
    digit_cut: int
    integer_core: int
    histograms: tuple
    image: tuple
    cycles: tuple
    histogram_distance: dict


def finite_parameters(base):
    natural(base, 2)
    maximum = factorial(base-1)
    d = 2
    while base**(d-1) <= d*maximum:
        d += 1
    return d, (d-1)*maximum


def compile_portrait(base):
    d, bound = finite_parameters(base)
    weights = tuple(factorial(x) for x in range(base))
    evaluation = {}
    for length in range(1, d):
        for word in combinations_with_replacement(range(base), length):
            if word[-1] == 0:
                continue
            counts = [0]*base
            for x in word:
                counts[x] += 1
            evaluation[tuple(counts)] = sum(weights[x] for x in word)
    need(len(evaluation) == comb(d-1+base, base)-d,
         'all positive numeral histograms through the digit cutoff')
    successor = {h: histogram(n, base) for h, n in evaluation.items()}
    need(all(h in evaluation for h in successor.values()), 'closed finite histogram carrier')
    distance, cycles = {}, set()
    for start in evaluation:
        path, positions, h = [], {}, start
        while h not in distance and h not in positions:
            positions[h] = len(path)
            path.append(h)
            h = successor[h]
        if h in positions:
            cut = positions[h]
            cycle = path[cut:]
            integers = [evaluation[x] for x in cycle]
            need(len(set(integers)) == len(integers), 'cycle evaluation loses no phase')
            cycles.add(canonical_cycle(integers))
            for x in cycle:
                distance[x] = 0
            path = path[:cut]
        for x in reversed(path):
            distance[x] = distance[successor[x]]+1
    image = tuple(sorted(set(evaluation.values())))
    need(all(n <= bound for n in image), 'all image integers stay in the forward-invariant core')
    need(all(digit_map(n, base) in image for n in image), 'smaller image core is invariant')
    for cycle in cycles:
        need(all(digit_map(n, base) == cycle[(i+1) % len(cycle)]
                 for i, n in enumerate(cycle)), 'reconstructed integer cycle in actual order')
    return Portrait(base, d, bound, tuple(evaluation), image, tuple(sorted(cycles)), distance)


def rank(portrait, n):
    """Explicit global integer rank for arrival at the complete cycle set."""
    ds = digits(n, portrait.base)
    excess = max(0, len(ds)-(portrait.digit_cut-1))
    if excess:
        return (len(portrait.histograms)+1)*excess
    if any(n in cycle for cycle in portrait.cycles):
        return 0
    return 1+portrait.histogram_distance[histogram(n, portrait.base)]


def audit_portrait(portrait):
    need(type(portrait) is Portrait, 'typed finite portrait')
    for field in ('base', 'digit_cut', 'integer_core'):
        natural(getattr(portrait, field), 1)
    need(type(portrait.histograms) is tuple and type(portrait.image) is tuple
         and type(portrait.cycles) is tuple and type(portrait.histogram_distance) is dict,
         'exact containers for the compiled carrier')
    for h in portrait.histograms:
        need(type(h) is tuple and all(type(x) is int for x in h), 'exact histogram coordinates')
    need(all(type(x) is int for x in portrait.image), 'exact image integers')
    need(all(type(c) is tuple and all(type(x) is int for x in c) for c in portrait.cycles),
         'exact cycle integers')
    need(all(type(h) is tuple and all(type(x) is int for x in h) and type(v) is int
             for h, v in portrait.histogram_distance.items()), 'exact distance table')
    need(compile_portrait(portrait.base) == portrait, 'recompute complete finite portrait')
    return portrait


def _classify(portrait, source):
    """Finite certificate from the proved rank; no arbitrary search cap."""
    natural(source, 1)
    cycles = {n: cycle for cycle in portrait.cycles for n in cycle}
    path, n = [source], source
    while n not in cycles:
        previous = rank(portrait, n)
        n = digit_map(n, portrait.base)
        need(rank(portrait, n) < previous, 'global rank decreases until an enumerated cycle')
        path.append(n)
    return tuple(path), cycles[n]


def classify(portrait, source):
    return _classify(audit_portrait(portrait), source)


def direct_core(base, bound):
    """Independent arithmetic path: integer digit recurrence, no histograms."""
    weights = tuple(factorial(d) for d in range(base))
    successor = array('I', [0])*(bound+1)
    for n in range(1, bound+1):
        successor[n] = weights[n] if n < base else successor[n//base]+weights[n % base]
    need(max(successor) <= bound, 'independent full core is invariant')
    done = bytearray(bound+1)
    cycles = set()
    for source in range(1, bound+1):
        if done[source]:
            continue
        path, n = [], source
        while not done[n]:
            done[n] = 1
            path.append(n)
            n = successor[n]
        if done[n] == 1:
            cycle = path[path.index(n):]
            cycles.add(canonical_cycle(cycle))
        for n in path:
            done[n] = 2
    return successor, tuple(sorted(cycles))


def completed_potential(successor, cycles):
    """Modify one marked edge per nonroot cycle, retaining the defect support."""
    bound = len(successor)-1
    cuts = {min(cycle): len(cycle) for cycle in cycles if 1 not in cycle}
    potential = array('I', [0])*(bound+1)
    for source in range(2, bound+1):
        if potential[source]:
            continue
        path, seen, n = [], set(), source
        while n != 1 and potential[n] == 0:
            need(n not in seen, 'modified completion has no unresolved cycle')
            seen.add(n)
            path.append(n)
            n = 1 if n in cuts else successor[n]
        height = potential[n]
        for n in reversed(path):
            height += 1
            potential[n] = height
    defects = {}
    for n in range(1, bound+1):
        residual = int(n != 1)+potential[successor[n]]-potential[n]
        if residual:
            defects[n] = residual
    need(defects == cuts, 'original dynamics retain exactly one cycle-length defect per artificial edge')
    for cycle in cycles:
        need(sum(int(n != 1)+potential[successor[n]]-potential[n] for n in cycle)
             == sum(n != 1 for n in cycle), 'cycle sum exposes the unchanged telescoping obstruction')
    return potential, defects


def rejects(callback):
    try:
        callback()
    except (ValueError, TypeError):
        return True
    return False


def main():
    expected = {
        6: ((1,), (2,), (25,), (26,)),
        10: ((1,), (2,), (145,), (169, 363601, 1454),
             (871, 45361), (872, 45362), (40585,)),
    }
    rows = []
    for base in (6, 10):
        portrait = compile_portrait(base)
        successor, cycles = direct_core(base, portrait.integer_core)
        need(cycles == portrait.cycles == expected[base],
             'independent complete integer and compressed histogram cycle classifications')
        potential, defects = completed_potential(successor, cycles)
        # Complete check on the compressed image, plus small and enormous sources.
        for n in portrait.image:
            if n not in {v for c in cycles for v in c}:
                need(rank(portrait, successor[n]) < rank(portrait, n), 'rank on every image-core integer')
        for source in tuple(range(1, 300))+tuple(base**k-1 for k in (1, 2, 7, 8, 31, 1000)):
            path, cycle = _classify(portrait, source)
            need(path[-1] in cycle and all(digit_map(x, base) == y for x, y in zip(path, path[1:])),
                 'all-height emitted certificate replays')
        for k in range(portrait.digit_cut, portrait.digit_cut+20):
            need(base**(k-1) > k*factorial(base-1), 'strict digit-count descent boundary')
        rows.append((base, portrait.digit_cut, portrait.integer_core,
                     len(portrait.histograms), len(portrait.image), len(cycles), max(potential)))
        print('BASE', base, 'CYCLES', cycles)
        print('MODIFIED_ROOT_DEFECTS', sorted(defects.items()))
        need(classify(portrait, 169)[1] in cycles, 'public certificate API authenticates the portrait')
        need(rejects(lambda: audit_portrait(replace(portrait, integer_core=True))), 'forged numeric metadata')
        need(rejects(lambda: audit_portrait(replace(portrait, cycles=((1,),)))), 'missing terminal cycles rejected')

    need(digits(25, 6) == (4, 1) and digits(26, 6) == (4, 2), 'base-6 values use decimal labels')
    need(digit_map(169, 10) == 363601 and digit_map(196, 10) == 363601,
         'permutation quotient preserves the immediate future')
    need(169 != 196 and histogram(169, 10) == histogram(196, 10),
         'same future does not recover original source')
    need(digit_map(145, 10) == 145 and digit_map(145, 10) != 1,
         'modified all-to-1 completion is not a certificate for the original map')
    for bad in (True, 3.0, 0, -1):
        need(rejects(lambda bad=bad: digits(bad, 10)), 'invalid positive source type/domain')
    need(rejects(lambda: value((0,)*10)), 'empty histogram rejected')
    print('GLOBAL_PORTRAIT_ROWS(base,digit_cut,integer_core,histograms,image,cycles,max_modified_rank)', rows)
    print('PROVED: every positive integer enters the completely classified finite cycle portrait')
    print('COMPLETION: every artificial root edge remains marked; its defect equals its original cycle length')
    print('COLLATZ_TRANSFER_OPEN: an authenticated all-source absorbing core or pointwise paid rank is still needed')
    print('CHECKS', CHECKS)


if __name__ == '__main__':
    main()
