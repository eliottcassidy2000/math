"""Integer quadratic components and an obstruction to finite monotone rank charts.

The proofs are in the companion note. Exact finite controls are independent
of the Collatz convergence conjecture and survive Python optimization.
"""
from fractions import Fraction
from itertools import product
from math import gcd, isqrt

import collatz_branch_toll_rank_20261004 as inherited


def need(ok, message):
    if not ok:
        raise ValueError(message)


def parameter(c):
    need(type(c) is int and c in (0, 1, 2), 'one of the three exact quadratic parameters')


def quadratic(x, c):
    parameter(c)
    return x*x-c


def core(c):
    parameter(c)
    return tuple(range(-2, 3)) if c == 2 else (-1, 0, 1)


def escape_address(x, c):
    """Return (component root, forward depth, sign), only outside the finite core."""
    parameter(c)
    need(type(x) is int and x not in core(c), 'escaping exact integer')
    sign = 1 if x > 0 else -1
    value, depth = abs(x), 0
    boundary = 3 if c == 2 else 2
    while True:
        root = isqrt(value+c)
        if root*root != value+c:
            break
        need(boundary <= root < value, 'inverse stripping strictly lowers magnitude')
        value, depth = root, depth+1
    return value, depth, sign


def decode_address(address, c):
    parameter(c)
    base, depth, sign = address
    boundary = 3 if c == 2 else 2
    need(type(base) is int and base >= boundary and isqrt(base+c)**2 != base+c,
         'canonical component root')
    need(type(depth) is int and depth >= 0 and type(sign) is int and sign in (-1, 1),
         'exact depth and sign')
    value = base
    for _ in range(depth):
        value = quadratic(value, c)
    return sign*value


def odd_step(n):
    need(type(n) is int and n > 0 and n % 2, 'positive odd source')
    z = 3*n+1
    a = (z & -z).bit_length()-1
    return z >> a, a


def increasing_run(length, unit):
    need(type(length) is int and length >= 1, 'positive exact run length')
    need(type(unit) is int and unit > 0 and unit % 2, 'positive odd cofactor')
    nodes = tuple(3**j*2**(length+1-j)*unit-1 for j in range(length+1))
    for x, y in zip(nodes, nodes[1:]):
        need(odd_step(x) == (y, 1) and y > x, 'actual strictly increasing one-run')
    return nodes


def one_run_rank(n):
    """Proper negative-anchor rank; decreases on valuation-one edges only."""
    need(type(n) is int and n > 0 and n % 2, 'positive odd rank input')
    k = ((n+1) & -(n+1)).bit_length()-1
    t = (n+1) >> k
    return 3**k*t, k


def main():
    print('QUADRATIC COMPONENTS AND RANK ATLASES: exact mechanisms, no global Collatz closure')
    controls, roots = 0, []
    for c in (0, 1, 2):
        K = core(c)
        need(all(quadratic(x,c) in K for x in K), 'finite forward-invariant core')
        critical, x = [], 0
        while x not in critical:
            critical.append(x)
            x = quadratic(x,c)
        cycle = critical[critical.index(x):]
        bases = set()
        for x in range(-5000, 5001):
            if x in K:
                continue
            address = escape_address(x,c)
            need(decode_address(address,c) == x, 'lossless signed component codec')
            b,j,s = address
            need(escape_address(quadratic(x,c),c) == (b,j+1,1),
                 'forward transition retains component and consumes sign')
            if x > 0:
                bases.add(b)
                need(quadratic(x,c) > x, 'integer escape outside core')
            controls += 1
        expected = 5000-isqrt(5000+c)
        need(len(bases) == expected, 'exact root counting formula')
        roots.append((c, len(bases)))
        print('x^2-'+str(c)+': core',K,'; critical orbit',critical,'; eventual cycle',cycle,
              '; escaping roots b<=5000:',len(bases))
    print('Exact signed codec/transition checks:',controls,'; roots count X-floor(sqrt(X+c))')
    for c,b in ((0,2),(1,2),(2,3)):
        chain = [b]
        for _ in range(4):
            chain.append(quadratic(chain[-1],c))
        need(escape_address(-chain[-1],c) == (b,4,-1), 'explicit escaping component')
        print('Escaping component c='+str(c)+': +/-',chain)

    # Bounded real dynamics with exponentially increasing arithmetic denominator.
    a,b,C = 3,4,5
    rational_checks = []
    for j in range(9):
        x = Fraction(2*a,C)
        need(a*a+b*b == C*C and gcd(gcd(a,b),C) == 1, 'primitive signed Gaussian lift')
        need(-2 <= x <= 2 and x.denominator == 5**(2**j), 'bounded real / unbounded denominator')
        rational_checks.append((j,x.denominator.bit_length()))
        aa,bb,CC = a*a-b*b,2*a*b,C*C
        need(Fraction(2*aa,CC) == quadratic(x,2), 'exact trace semiconjugacy')
        a,b,C = aa,bb,CC
    print('Chebyshev rational orbit6/5: denominator=5^(2^j); denominator bit lengths:',rational_checks)
    for z in (Fraction(2),Fraction(3,2),Fraction(-2),Fraction(1,3)):
        need(z*z+1/(z*z) == quadratic(z+1/z,2), 'J(z^2)=J(z)^2-2')

    # Exhaustive arbitrary chart assignments in a small finite proof control.
    # Each chart n^d is proper and increasing; selecting history-dependent
    # labels still cannot evade the repeated-chart obstruction.
    atlas_assignments = 0
    for q in range(1,5):
        nodes = increasing_run(q,101)
        for labels in product(range(1,q+1),repeat=q+1):
            values = [n**degree for n,degree in zip(nodes,labels)]
            need(not all(y<x for x,y in zip(values,values[1:])),
                 'no arbitrary q-chart selection ranks this q-edge increasing block')
            witness = next((i,j) for i in range(q+1) for j in range(i+1,q+1)
                           if labels[i] == labels[j])
            i,j = witness
            need(values[j] > values[i], 'repeated monotone chart contradiction')
            atlas_assignments += 1
        # A finite unrolled atlas works before a chart has to repeat.
        values = [nodes[i]**(q-i) for i in range(q)]
        need(all(y<x for x,y in zip(values,values[1:])), 'finite unrolled positive control')
    need(atlas_assignments == 1114, 'complete small chart-assignment universe')
    print('Arbitrary polynomial-chart assignments rejected:',atlas_assignments,
          '; q1..4, q+1 increasing states; finite unrolled controls pass')
    macro_assignments = 0
    for q in range(1,4):
        nodes = increasing_run(3*q,101)
        for lengths in product(range(1,4),repeat=q):
            indices = [0]
            for length in lengths:
                indices.append(indices[-1]+length)
            selected = [nodes[i] for i in indices]
            for labels in product(range(1,q+1),repeat=q+1):
                values = [n**degree for n,degree in zip(selected,labels)]
                need(not all(y<x for x,y in zip(values,values[1:])),
                     'bounded forward macros do not remove repeated-chart obstruction')
                macro_assignments += 1
    need(macro_assignments == 2262, 'complete bounded-macro assignment universe')
    print('Bounded macro controls:',macro_assignments,'; q1..3, all macro lengths1..3 and all chart assignments')

    # The incoming unbounded valuation coordinate escapes the chart hypothesis.
    one_controls = 0
    for length in range(1,65):
        for unit in (1,3,101):
            nodes = increasing_run(length,unit)
            ranks = [one_run_rank(n) for n in nodes]
            need(all(E>=n+1 for n,(E,_) in zip(nodes,ranks)), 'proper one-run rank bound')
            need(all(b<a for a,b in zip(ranks,ranks[1:])), 'unbounded precision pays every one edge')
            one_controls += length
    need(odd_step(9) == (7,2) and one_run_rank(7) > one_run_rank(9),
         'one-run rank fails at a reset even when the integer decreases')
    need([one_run_rank(n)[0] for n in (27,41,31)] == [63,63,243],
         'reset-two energy debt remains visible')
    print('Valuation-coordinate edge controls:',one_controls,
          '; reset hostile9->7 raises energy15->27;27->41->31 gives63,63,243')

    # A quadratic's sign fold is an actual common-future relation. The incoming
    # Collatz root-centered energy fold is only a rank equality.
    for n in range(-100,101):
        for c in (0,1,2):
            need(quadratic(n,c) == quadratic(-n,c), 'quadratic sign fold preserves the next value')
    need(inherited.rank(3) == inherited.rank(-1), 'incoming equal-energy fold')
    need(inherited.step(3)[0] == 5 and inherited.step(-1)[0] == -1,
         'equal-rank fold is not a Collatz common-future edge')
    need(inherited.step(5)[0] == 1, 'positive3 and negative-1 have different exact basins')
    print('Fold hostile: rank(3)=rank(-1), but3->5->1 and-1->-1; sign cannot be dropped')

    for thunk in (lambda: escape_address(1,0),lambda: escape_address(3,True),
                  lambda: decode_address((4,0,1),0),lambda: increasing_run(True,1),
                  lambda: one_run_rank(2)):
        try:
            thunk()
        except ValueError:
            pass
        else:
            raise ValueError('invalid type or domain accepted')
    print('PASS: finite critical orbit, proper graph reduction, and global attraction are different predicates')


if __name__ == '__main__':
    main()
