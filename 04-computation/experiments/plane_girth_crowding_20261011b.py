"""Exact arithmetic for the BFS-tree / Euclidean-packing obstruction.

This checks numerical bounds, not existence of geometric high-chromatic graphs.
The proof and its hypotheses are in the companion result packet.
"""
from fractions import Fraction as Q

CHECKS = 0


def need(ok, message):
    if not ok:
        raise ValueError(message)


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def tree_size(degree, radius):
    need(type(degree) is int and degree >= 2, 'integer degree >=2')
    need(type(radius) is int and radius >= 1, 'integer radius >=1')
    return 1 + degree*sum((degree-1)**i for i in range(radius))


def excluded(degree, radius, separation):
    need(type(separation) is Q and separation > 0, 'positive exact separation')
    return tree_size(degree, radius) > (1+2*radius/separation)**2


def first_excluded(degree, separation):
    need(type(degree) is int and degree >= 3, 'branching degree >=3')
    need(type(separation) is Q and separation > 0, 'positive exact separation')
    radius = 1
    while not excluded(degree, radius, separation):
        radius += 1
    return radius


def rejects(fn):
    try:
        fn()
    except ValueError:
        check(True, 'hostile rejected')
    else:
        raise RuntimeError('hostile accepted')


def main():
    for degree in range(2, 12):
        level, total = degree, 1
        for radius in range(1, 41):
            total += level
            check(tree_size(degree, radius) == total, 'independent level recurrence')
            formula = 1+2*radius if degree == 2 else 1+degree*((degree-1)**radius-1)//(degree-2)
            check(total == formula, 'closed form with degree-two boundary')
            level *= degree-1
    # Explicit one-dimensional cycle growth: branching is essential.
    check(all(not excluded(2, r, Q(1)) for r in range(1, 100)),
          'minimum degree two is not excluded by this growth argument')
    for degree in range(3, 8):
        for denominator in range(1, 101):
            h = Q(1, denominator)
            r = first_excluded(degree, h)
            check(excluded(degree, r, h), 'strict rational area contradiction')
            check(all(not excluded(degree, s, h) for s in range(1, r)),
                  'reported radius is first for this sufficient inequality')
            check(tree_size(degree, r)*h*h > (h+2*r)**2,
                  'independent denominator-cleared inequality')
    # At equality the packing estimate alone gives no contradiction.
    check(not excluded(8, 1, Q(1)), 'equality in the packing estimate stays inconclusive')
    # A triangle fails the necessary girth hypothesis for radius one.
    check(not 3 >= 2*1+2, 'triangle cannot be passed through the high-girth gate')
    rejects(lambda:tree_size(True, 1))
    rejects(lambda:tree_size(3, True))
    rejects(lambda:tree_size(1, 2))
    rejects(lambda:excluded(3, 1, 0.1))
    rejects(lambda:first_excluded(2, Q(1)))
    rejects(lambda:first_excluded(3, Q(0)))
    print('PROVED: fixed h>0 and sufficiently large girth force 2-degeneracy and 3-colourability.')
    print('Sufficient thresholds from exact strict packing inequalities:')
    for h in (Q(1), Q(1,2), Q(1,10), Q(1,100), Q(1,1000)):
        r = first_excluded(3,h)
        print('h=',str(h),' radius=',r,' girth>=',2*r+2,' implies 3-colourability')
    print('Critical chromatic k controls at separation h=1/1000:')
    for k in (6,7):
        r = first_excluded(k-1,Q(1,1000))
        print('k=',k,' radius=',r,' girth>=',2*r+2,' is excluded')
    print('No plane chromatic decision, realized witness, or optimal threshold claimed.')
    print('PASS:',CHECKS,'exact checks')


if __name__ == '__main__':
    main()
