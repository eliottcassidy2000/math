"""Exact palette-action and marked-seam controls, not a six-colouring search."""
from fractions import Fraction as Q
from itertools import combinations_with_replacement, permutations, product
from math import floor

CHECKS = 0


def need(ok, why):
    if not ok:
        raise ValueError(why)


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def audit_perm(p):
    need(type(p) is tuple and all(type(x) is int for x in p)
         and set(p) == set(range(len(p))), 'exact permutation')


def compose(p, q):
    audit_perm(p)
    audit_perm(q)
    need(len(p) == len(q), 'same palette')
    return tuple(p[q[i]] for i in range(len(p)))


def phase(m):
    need(type(m) is int and m >= 1, 'positive number of phase orbits')
    return tuple(3*b+(j+1) % 3 for b in range(m) for j in range(3))


def centralizer(m):
    """All commuting permutations, by orbit permutation and phase offset."""
    phase(m)
    return tuple(tuple(3*s[b]+(j+d[b]) % 3 for b in range(m) for j in range(3))
                 for s in permutations(range(m)) for d in product(range(3), repeat=m))


def seam_solution(n, seams):
    """Solve x_u xor x_v = label, retaining parallel edges and loops."""
    need(type(n) is int and n >= 0 and type(seams) is tuple, 'exact seam source')
    adj = [[] for _ in range(n)]
    for e in seams:
        need(type(e) is tuple and len(e) == 3
             and all(type(v) is int for v in e), 'exact seam row')
        u, v, bit = e
        need(0 <= u < n and 0 <= v < n and bit in (0, 1), 'seam domain')
        adj[u].append((v, bit))
        adj[v].append((u, bit))
    out = [None]*n
    for root in range(n):
        if out[root] is not None:
            continue
        out[root] = 0
        stack = [root]
        while stack:
            u = stack.pop()
            for v, bit in adj[u]:
                proposed = out[u] ^ bit
                if out[v] is None:
                    out[v] = proposed
                    stack.append(v)
                elif out[v] != proposed:
                    return None
    return tuple(out)


def colour9(x, y):
    need(type(x) is Q and type(y) is Q, 'exact rational test point')
    return floor(5*x/3) % 3, floor(5*y/3) % 3


def rejecting(fn):
    try:
        fn()
    except ValueError:
        return True
    return False


def main():
    print('Scope: fixed palette actions, a proved nine-colour control, finite seam systems.')
    r = phase(2)
    ident = tuple(range(6))
    direct = tuple(p for p in permutations(range(6)) if compose(p, r) == compose(r, p))
    check(set(direct) == set(centralizer(2)), 'independent full S6 centralizer')
    involutions = tuple(p for p in direct if compose(p, p) == ident)
    fixed = sorted(sum(p[i] == i for i in range(6)) for p in involutions)
    check(len(direct) == 18 and fixed == [0, 0, 0, 6], 'six-palette reflection obstruction')
    print('Two phase-three orbits: centralizer18; involution fixed-colour counts', fixed)
    bank9 = centralizer(3)
    id9 = tuple(range(9))
    refs9 = tuple(p for p in bank9 if compose(p, p) == id9)
    check(len(bank9) == 162 and len(refs9) == 10, 'three-orbit centralizer and involutions')
    check(sorted(sum(p[i] == i for i in range(9)) for p in refs9) == [3]*9+[9],
          'first possible fixed and moving phase orbits')
    print('Three phase-three orbits: centralizer162; nine nonidentity reflections fix3 colours.')

    # Opposite phase on the two sides repairs commutation, but not the axis.
    same = tuple(3*b+(j+1) % 3 for b in range(2) for j in range(3))
    opposite = tuple(3*b+(j+(1 if b == 0 else -1)) % 3
                     for b in range(2) for j in range(3))
    mirror = tuple(3*(1-b)+(-j) % 3 for b in range(2) for j in range(3))
    check(compose(mirror, same) != compose(same, mirror), 'reversing same phase does not commute')
    check(compose(mirror, opposite) == compose(opposite, mirror), 'retained side repairs phase rule')
    check(all(mirror[i] != i for i in range(6)), 'axis has no colour under paired reflection')
    # Changing the translation cycle type admits the formal group relation.
    alternative = (1, 2, 0, 3, 4, 5)
    reflection = (0, 1, 2, 4, 3, 5)
    check(compose(alternative, reflection) == compose(reflection, alternative)
          and compose(reflection, reflection) == ident, '3+1+1+1 formal relaxation')
    print('Opposite half-plane phases commute with paired reflection; their axis still has no fixed colour.')
    print('Cycle type3+1+1+1 admits the symmetry algebra; no metric colouring follows.')

    # Exact global proof in the note; finite checks include negative boundaries.
    side = Q(3, 5)
    points = tuple(Q(i, 10) for i in range(-24, 25))
    for x, y in product(points, repeat=2):
        c = colour9(x, y)
        check(colour9(y, x) == (c[1], c[0]), 'reflection including half-open boundaries')
        check(colour9(x+side, y+side) == ((c[0]+1) % 3, (c[1]+1) % 3), 'diagonal phase')
    check(2*side*side < 1 and 2*side > 1, 'same-cell diameter and different-cell gap')
    # Same-colour cells differ by multiples of three; this is a bounded
    # independent endpoint check, not the all-plane proof by itself.
    for i, j in product(range(-6, 7), repeat=2):
        if i % 3 == j % 3 == 0 and (i, j) != (0, 0):
            gapx, gapy = max(0, abs(i)-1)*side, max(0, abs(j)-1)*side
            check(gapx*gapx+gapy*gapy > 1, 'same-colour square separation')
    print('Sharp positive control: square side3/5, colour(floor(5x/3),floor(5y/3)) mod3, palette9.')

    cases = 0
    for n in range(4):
        pairs = tuple(combinations_with_replacement(range(n), 2))
        for choices in product(range(4), repeat=len(pairs)):
            seams = tuple((u, v, bit) for (u, v), mask in zip(pairs, choices)
                          for bit in (0, 1) if mask >> bit & 1)
            got = seam_solution(n, seams)
            literal = tuple(x for x in product((0, 1), repeat=n)
                            if all(x[u] ^ x[v] == b for u, v, b in seams))
            check((got is not None) == bool(literal), 'all small signed seam systems')
            if got is not None:
                check(got in literal, 'actual exported frame assignment')
            cases += 1
    # Lower half-plane contains a unit edge, upper half-plane its equilateral
    # third point. The two crossing edges impose contradictory parallel seams.
    check(seam_solution(2, ((0, 1, 0), (0, 1, 1))) is None, 'cross-seam triangle hostile')
    check(seam_solution(3, ((0, 1, 1), (1, 2, 1), (0, 2, 0))) is not None,
          'zero-holonomy positive control')
    # In orthonormal (u,v) coordinates the squared vertical differences are
    # (sqrt3/2)^2, while horizontal differences have square1/4.
    check(Q(1, 4)+Q(3, 4) == 1, 'exact equilateral cross-edge distances')
    print('Literal seam systems including loops and parallel constraints:', cases)
    for fn in (lambda: phase(True), lambda: compose((False, 1), (0, 1)),
               lambda: seam_solution(True, ()),
               lambda: seam_solution(2, ((False, 1, 0),)),
               lambda: seam_solution(2, ((0, 1, 0.0),)),
               lambda: seam_solution(2, [(0, 1, 0)]),
               lambda: seam_solution(2, ((0, 2, 0),)),
               lambda: colour9(0.0, Q(1))):
        check(rejecting(fn), 'exact-domain hostile')
    print(f'PASS: {CHECKS} exact checks. Six-colour plane construction remains OPEN.')


if __name__ == '__main__':
    main()
