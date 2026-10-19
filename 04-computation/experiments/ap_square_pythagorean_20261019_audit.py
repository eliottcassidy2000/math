#!/usr/bin/env python3
"""Independent standard-library audit of AP square-sum triangles.

Run with ordinary Python or python -O: checks use explicit exceptions, not assert.
The exhaustive triangle census uses two independent routes:
(1) primitive Pythagorean triples plus all integer-square scales;
(2) the middle square x+z=2y, enumerating the smaller pair-sum root.
Finite checks support, but do not replace, the accompanying exact proofs.
"""
from itertools import combinations
from math import gcd, isqrt, lcm, pi, sin, sqrt, log

MAX_VERTEX = 1_000_000
TREE_DEPTH = 8
HYPO_MAX = 100_000


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def square(n):
    return n >= 0 and isqrt(n) ** 2 == n


def euclid(m, n):
    return n*n-m*m, 2*m*n, n*n+m*m


def admissible(m, n):
    a, b, c = euclid(m, n)
    return c*c > 4*a*b


def cone(m, n):
    if not admissible(m, n):
        return None
    if n < 2*m:
        return 'L'
    if n > 3*m:
        return 'R'
    raise RuntimeError(f'Admissible point outside exact coarse cone separation: {(m,n)}')


def child(m, n, letter):
    if letter == 'L':
        return n, 2*n-m
    if letter == 'M':
        return n, 2*n+m
    if letter == 'R':
        return m, n+2*m
    raise ValueError(letter)


def primitive_parameters(c_max):
    for n in range(2, isqrt(c_max)+1):
        for m in range(1, n):
            if m*m+n*n <= c_max and (m+n) % 2 == 1 and gcd(m,n) == 1:
                yield m, n


def verify_triangle(triangle):
    x, y, z = triangle
    require(0 < x < y < z and x+z == 2*y, f'Bad AP: {triangle}')
    sums = (x+y, x+z, y+z)
    require(all(square(v) for v in sums), f'Bad pair-sum squares: {triangle}')
    roots = tuple(isqrt(v) for v in sums)
    a, b, c = roots
    require(a % 2 == b % 2 == c % 2 == 0, 'Odd pair-sum root')
    require((a+c) % 4 == (c-a) % 4 == 0, 'Pythagorean inverse not integral')
    p, q = (a+c)//4, (c-a)//4
    g = gcd(p, q)
    A, B, C = p//g, q//g, b//(2*g)
    require(A > B > 0 and A*A+B*B == C*C, 'Inverse Pythagorean failure')
    require(gcd(A,B) == 1 and (A+B) % 2 == 1, 'Nonprimitive inverse')
    expected = tuple(2*g*g*(C*C+k*4*A*B) for k in (-1,0,1))
    require(expected == triangle, 'Inverse reconstruction failure')
    require(all(gcd(v,w) == 2*g*g for v,w in combinations(triangle,2)), 'Vertex gcd failure')
    require(all(gcd(v,w) == 2*g for v,w in combinations(roots,2)), 'Root gcd failure')
    require((y-x) % (96*g*g) == 0, 'Difference divisibility failure')
    require(square(y//2) and not square(x//2) and not square(z//2), 'Twice-square position failure')
    return g, A, B, C


def ppt_triangle_census(limit):
    triangles = set()
    primitive = set()
    for m, n in primitive_parameters(isqrt(limit//2)):
        a, b, c = euclid(m,n)
        if c*c <= 4*a*b:
            continue
        base = tuple(2*(c*c+k*4*a*b) for k in (-1,0,1))
        if base[2] > limit:
            continue
        require(base not in primitive, 'Duplicate primitive triangle')
        primitive.add(base)
        for g in range(1, isqrt(limit//base[2])+1):
            triangle = tuple(g*g*v for v in base)
            require(triangle not in triangles, 'Duplicate square-scaled triangle')
            require(verify_triangle(triangle)[0] == g, 'Incorrect recovered scale')
            triangles.add(triangle)
    return triangles, primitive


def middle_square_census(limit):
    triangles = set()
    # Since x+z=2y is a square, y=2s^2 for a positive integer s.
    for s in range(1, isqrt(limit//2)+1):
        y = 2*s*s
        # Let a^2=x+y. Positivity gives a^2>y, and x<y gives a<2s.
        # The upper-vertex cutoff requires z=3y-a^2<=limit.
        lower_square = max(y+1, 3*y-limit)
        lower = isqrt(lower_square)
        if lower*lower < lower_square:
            lower += 1
        for a in range(lower, 2*s):
            if square(4*y-a*a):
                triangle = (a*a-y, y, 3*y-a*a)
                require(triangle not in triangles, 'Middle-square duplicate')
                verify_triangle(triangle)
                triangles.add(triangle)
    return triangles


def audit_tree():
    nodes = [('',1,2)]
    seen = set()
    admissible_count = 0
    entry_counts = {'L3':0,'L4':0,'R3':0,'R4':0}
    for word, m, n in nodes:
        require((m,n) not in seen, 'Duplicate Euclid tree node')
        seen.add((m,n))
        require(0 < m < n and gcd(m,n) == 1 and (m+n)%2 == 1, 'Invalid tree parameter')
        kind = cone(m,n)
        children = {letter:child(m,n,letter) for letter in 'LMR'}
        require(not admissible(*children['M']), 'M child entered admissible cone')
        if kind is not None:
            admissible_count += 1
            allowed = [letter for letter,pair in children.items() if admissible(*pair)]
            require(allowed == [kind], 'Incorrect admissible child set')
            require(word.endswith(kind*3), 'Admissible word lacks terminal triple letter')
        for letter in 'LR':
            pair = (m,n)
            first = None
            for step in range(1,5):
                pair = child(*pair,letter)
                if admissible(*pair) and first is None:
                    first = step
            require(cone(*pair) == letter, 'Four repeated steps missed the required cone')
            if not word.endswith(letter):
                require(first in (3,4), 'A maximal run first entered outside steps 3 and 4')
                entry_counts[letter+str(first)] += 1
        if word.endswith('LLLL') or word.endswith('RRRR'):
            require(kind is not None, 'Terminal four-letter suffix is inadmissible')
        if len(word) < TREE_DEPTH:
            nodes.extend((word+letter,*pair) for letter,pair in children.items())
    return len(seen), admissible_count, entry_counts


def audit_hypotenuse_counts():
    thresholds = (100, 1000, 10000, HYPO_MAX)
    counts = {x:0 for x in thresholds}
    total = {x:0 for x in thresholds}
    all_admissible = []
    for m,n in primitive_parameters(HYPO_MAX):
        a,b,c = euclid(m,n)
        status = admissible(m,n)
        for x in thresholds:
            if c <= x:
                total[x] += 1
                counts[x] += status
        if status:
            all_admissible.append((c,m,n))
        # Independent exact polynomial checks, including opposite-parity boundaries.
        quartic = n**4-8*m*n**3+2*m*m*n*n+8*m**3*n+m**4
        require((quartic > 0) == status and quartic != 0, 'Cone quartic/boundary failure')
        require(not admissible(*child(m,n,'M')), 'Direct census M-cone failure')
        if status:
            expected = cone(m,n)
            actual = ''.join(q for q in 'LMR' if admissible(*child(m,n,q)))
            require(actual == expected, 'Direct census exact child failure')
    require(min(all_admissible) == (41,4,5), 'Minimum admissible hypotenuse failure')
    return thresholds,total,counts


def audit_windmill(size=8):
    pairs = [(k,k+1) for k in range(4,4+size)]
    triples = [euclid(*pair) for pair in pairs]
    H = lcm(*(c for _,_,c in triples))
    center = 2*H*H
    vertices = {center}
    expected_edges = set()
    for a,b,c in triples:
        require(c*c > 4*a*b, 'Inadmissible L-spine point')
        g = H//c
        triangle = tuple(2*g*g*(c*c+k*4*a*b) for k in (-1,0,1))
        require(triangle[1] == center, 'Windmill centers differ')
        require(triangle[0] not in vertices and triangle[2] not in vertices, 'Windmill endpoint collision')
        verify_triangle(triangle)
        vertices.update(triangle)
        expected_edges.update(tuple(sorted(edge)) for edge in combinations(triangle,2))
    actual_edges = {tuple(sorted(edge)) for edge in combinations(vertices,2) if square(sum(edge))}
    require(len(vertices) == 2*size+1 and expected_edges <= actual_edges, 'Windmill subgraph failure')
    return size,len(vertices),len(actual_edges)-len(expected_edges)


def simpson(f,a,b,n=20000):
    require(n % 2 == 0, 'Simpson interval count is odd')
    h=(b-a)/n
    return h/3*(f(a)+f(b)+4*sum(f(a+h*k) for k in range(1,n,2))+2*sum(f(a+h*k) for k in range(2,n,2)))


def main():
    ppt,primitive = ppt_triangle_census(MAX_VERTEX)
    brute = middle_square_census(MAX_VERTEX)
    require(ppt == brute, f'Census disagreement: PPT-only={ppt-brute}, brute-only={brute-ppt}')
    minimum = min(ppt,key=lambda t:t[2])
    require(minimum == (482,3362,6242), 'Minimum largest-vertex failure')
    require((386,8450,16514) in ppt, 'R^3 example missing')
    require(verify_triangle(minimum) == (1,40,9,41), 'Low-spine example inverse failure')
    print(f'PASS independent exhaustive triangle census, z <= {MAX_VERTEX}: {len(ppt)} total, {len(primitive)} primitive (g=1).')
    print(f'PASS smallest largest vertex: {minimum}; pair-sum roots (62,82,98).')
    print('PASS R^3 example: (386,8450,16514); smallest largest vertex does not mean smallest first vertex.')
    n_nodes,n_admissible,entries = audit_tree()
    print(f'PASS Euclid tree depth {TREE_DEPTH}: {n_nodes} distinct nodes, {n_admissible} admissible; exact child rules and suffix tests.')
    print(f'PASS first entry in maximal terminal runs: {entries}; uniform L^4/R^4 tests.')
    thresholds,total,counts = audit_hypotenuse_counts()
    print('PASS direct primitive-Pythagorean census, cone quartic, exact children, no equality boundary:')
    for x in thresholds:
        print(f'  C <= {x}: all={total[x]}, admissible={counts[x]}, X/(6*pi)={x/(6*pi):.9f}')
    size,vertices,extra = audit_windmill()
    print(f'PASS common-center {size}-triangle windmill: {vertices} distinct vertices, {extra} extra edges in this sample; theorem claims a subgraph only.')
    J=simpson(lambda u:1/sqrt(1+2*sin(u)),0,pi/6)
    kappa=J/(sqrt(2)*pi*pi)
    require(abs(kappa-0.0309567443988252174) < 1e-13, 'Counting constant quadrature mismatch')
    print(f'PASS numerical counting constant: kappa={kappa:.17g}; kappa/2={kappa/2:.17g}.')
    print('ANALYTIC: P(N)=kappa*sqrt(N)+O(N^(1/4)*log N); all scales T(N)=(kappa/2)*sqrt(N)*log N+O(sqrt N).')
    print('ANALYTIC: coprime opposite-parity lattice density is 4/pi^2; admissible Euclid angular width is pi/12.')
    print('NOTE: finite endpoint checks agree with the four-squares-in-AP theorem; this script does not reprove that theorem.')
    print('ALL AUDIT CHECKS PASSED (no Python assert statements used).')


if __name__ == '__main__':
    main()
