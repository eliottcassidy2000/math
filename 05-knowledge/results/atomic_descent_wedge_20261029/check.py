"""Exact scoped controls; standard library, no removable proof assertions."""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from hashlib import sha256
from pathlib import Path
import json

CHECKS = 0


def require(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(x) for x in zip(*a)]


def mm(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def add(a, b, scale=1):
    return [[x+scale*y for x, y in zip(r, s)] for r, s in zip(a, b)]


def power(a, n, modulus=None):
    result = eye(len(a))
    while n:
        if n & 1:
            result = mm(result, a)
        a = mm(a, a)
        if modulus:
            result = [[x % modulus for x in row] for row in result]
            a = [[x % modulus for x in row] for row in a]
        n //= 2
    return result


def det(a):
    return sum((-1)**sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p)))
               * prod(a[i][p[i]] for i in range(len(p)))
               for p in permutations(range(len(a))))


def prod(xs):
    result = 1
    for x in xs:
        result *= x
    return result


def wedge(a):
    pairs = list(combinations(range(len(a)), 2))
    return [[a[i][k]*a[j][l]-a[i][l]*a[j][k] for k, l in pairs] for i, j in pairs]


def bivector(u, v):
    return [u[i]*v[j]-u[j]*v[i] for i, j in combinations(range(4), 2)]


def pf(p):
    return p[0]*p[5]-p[1]*p[4]+p[2]*p[3]


def mv(a, v):
    return [sum(x*y for x, y in zip(row, v)) for row in a]


def spectral_checks():
    c = [[1, 0, 1], [1, 1, 1], [1, 1, 0]]
    w = wedge(c)
    s = [[-1, 0, -1], [0, 0, -1], [0, 1, -1]]
    require(det(c) == -1 and det(s) == -1, 'unimodular source and intertwiner')
    require(mm(w, s) == mm(s, add(c, eye(3), -1)), 'exterior shift conjugacy')
    require(add(add(add(power(c, 3), power(c, 2), -2), c, -1), eye(3)) == [[0]*3 for _ in range(3)], 'source cubic')
    require(add(add(add(power(w, 3), power(w, 2)), w, -2), eye(3), -1) == [[0]*3 for _ in range(3)], 'exterior cubic')
    v, counts = [0, 1, 0], []
    for m in range(101):
        counts.append(sum(v))
        v = mv(transpose(c), v)
    require(counts[9] == 2065 and counts[6] == 182, 'order13 and order10 counts')
    for m in range(3, len(counts)):
        require(counts[m] == 2*counts[m-1]+counts[m-2]-counts[m-3], 'count recurrence')
    n = add(c, eye(3), -3)
    require(all(x % 7 == 0 for row in power(n, 3) for x in row), 'nilpotence mod7')
    require(power(c, 7, 7) == [[3*x for x in row] for row in eye(3)], 'seventh power mod7')
    for m, value in enumerate(counts):
        require(value % 7 == (pow(3, m, 7)*(1-m*(m-1)//2)) % 7, 'mod7 closed law')
        require((value % 7 == 0) == (m % 7 in (2, 6)), 'exact zero classes mod7')
    require(power(c, 12, 13) == eye(3), 'mod13 period divides12')
    g = [[1, 0, 0, 0]] + [[0]+row for row in c]
    b = wedge(g)
    j = [[0]*6 for _ in range(6)]
    for i, k, sign in ((0, 5, 1), (1, 4, -1), (2, 3, 1)):
        j[i][k] = j[k][i] = sign
    require(mm(mm(transpose(b), j), b) == [[-x for x in row] for row in j], 'wedge pairing anti-isometry')
    require(mm(b, j) != mm(j, b), 'not a Hodge-invariant transfer')
    for u in product((-1, 0, 1), repeat=4):
        for v in product((-1, 0, 1), repeat=4):
            p = bivector(u, v)
            require(pf(p) == 0 and mv(b, p) == bivector(mv(g, u), mv(g, v)), 'literal decomposable transport')
    require(pf([1, 0, 0, 0, 0, 1]) == 1, 'sum of decomposables need not decompose')
    perm = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
    require(sum(wedge(perm)[i][i] for i in range(6)) == 0, 'signed transposition trace')
    require(sum(abs(wedge(perm)[i][i]) for i in range(6)) == 2, 'unsigned transposition trace')
    # Letter e: four-state acceptance and six-state legality are different functors.
    ae = [[int(j == target) for j in range(4)] for target in (3, 2, 3, 3)]
    le = [[int(j == target) for j in range(6)] for target in (5, 2, 4, 5, 4, 5)]
    require(len({tuple(row) for row in le}) == 3, 'three independent legality images')
    we = wedge(ae)
    nonzero = [row for row in we if any(row)]
    require(nonzero and all(row == nonzero[0] or row == [-x for x in nonzero[0]] for row in nonzero), 'rank-one exterior acceptance e')
    # X=d/dx, Y=d/dy+x d/dz: {x,y}=1, {x,z}=x.
    for x in range(-3,4):
        pp=[[0,1,x,0],[-1,0,0,0],[-x,0,0,0],[0,0,0,0]]
        dp=lambda i,j: int((i,j)==(0,2))-int((i,j)==(2,0))
        jacobi=sum(pp[i][0]*dp(j,k) for i,j,k in ((0,1,2),(1,2,0),(2,0,1)))
        require(pf([1,x,0,0,0,0]) == 0 and jacobi == 1, 'decomposable field fails Jacobi')
    return dict(C=c, exterior=w, intertwiner=s, counts_order4_to17=counts[:14],
                mod7_zero_orders=[3, 6], mod13_period_divisor=12,
                source_discriminant=49, field_jacobi_defect=1,
                letter_e_ranks=dict(exterior_acceptance=1, legality=3))


def mul27(a, b):
    z = [0]*5
    for i in range(3):
        for j in range(3):
            z[i+j] += a[i]*b[j]
    # x^3=x+1 over F3.
    for k in (4, 3):
        z[k-2] += z[k]
        z[k-3] += z[k]
    return tuple(x % 3 for x in z[:3])


def projective_checks():
    vectors = list(product(range(3), repeat=3))
    nonzero = [v for v in vectors if any(v)]
    canon = lambda v: min(v, tuple((-x) % 3 for x in v))
    points = sorted({canon(v) for v in nonzero})
    require(len(vectors) == 27 and len(nonzero) == 26 and len(points) == 13, '27/26/13 source quotient')
    lines = [{i for i, v in enumerate(points) if sum(a*b for a, b in zip(v, normal)) % 3 == 0} for normal in points]
    require(all(len(line) == 4 for line in lines), 'four points per line')
    require(all(len(a & b) == 1 for a, b in combinations(lines, 2)), 'unique line intersection')
    generator, x, orbit = (0, 2, 0), (1, 0, 0), []
    for _ in range(26):
        orbit.append(x)
        x = mul27(x, generator)
    require(len(set(orbit)) == 26 and x == (1, 0, 0), 'primitive field element')
    squares = {mul27(v, v) for v in nonzero}
    require(len(squares) == 13 and len({canon(v) for v in squares}) == 13, 'square identifies signed pairs with norm-one subgroup')
    return dict(vector_count=27, nonzero_count=26, projective_points=13,
                projective_line_F27=28, projective_line_F9=10, line_sizes=4)


def arithmetic_checks():
    s = {2, 3}; stages = [sorted(s)]
    for _ in range(6):
        s = s | {a+b for a, b in combinations(sorted(s), 2)}
        stages.append(sorted(s))
    for t in range(4, 7):
        m = 27*2**(t-4)+1
        require(stages[t] == [2, 3, 5]+list(range(7, m+1)), 'summand frontier')
    require([(a, 13-a) for a in sorted(s) if a < 13-a and 13-a in s] == [(2, 11), (3, 10), (5, 8)], 'thirteen parent fibre')
    square_rows = []
    for m in range(1, 51):
        for a in (1, m*m, m*(2*m+1), 3*m*m+1):
            d1 = sum(a+i for i in range(m+1))-sum(a+i for i in range(m+1, 2*m+1))
            d2 = sum((a+i)**2 for i in range(m+1))-sum((a+i)**2 for i in range(m+1, 2*m+1))
            require(d1 == a-m*m, 'linear balance')
            require(d2 == (a-m*(2*m+1))*(a+m), 'square balance factorization')
        if m <= 3:
            a=m*(2*m+1)
            square_rows.append(dict(left=list(range(a,a+m+1)), right=list(range(a+m+1,a+2*m+1)), total=sum((a+i)**2 for i in range(m+1))))
    require(sum((181-i)**15 for i in range(4))-sum((181+i)**15 for i in range(1, 4)) == -61355600116348100287252240367, 'higher-power unit-width hostile retained')
    # Elementary normalized square lift of coprime additive triples.
    from math import gcd
    lift_examples = []
    for a in range(1, 70):
        for b in range(a+1, 70):
            if gcd(a, b) != 1:
                continue
            h = 4 if (a+b) % 2 == 0 else 1
            aa, bb, cc = (a-b)**2//h, 4*a*b//h, (a+b)**2//h
            require(aa+bb == cc and gcd(aa, bb) == 1, 'primitive additive-to-square lift')
            if (a,b) in ((1,2),(1,4),(2,3),(1,8)):
                lift_examples.append([a,b,a+b,aa,bb,cc])
    return dict(frontiers=[max(row) for row in stages], square_rows=square_rows, square_lifts=lift_examples)


def odd_step(n, coefficient=3):
    z = coefficient*n+1
    a = (z & -z).bit_length()-1
    return z >> a, a


def weight(m):
    return F(1, 2*4**(m.bit_length()-1))


def descent_checks():
    orbit, n = [], 13
    for _ in range(3):
        orbit.append(n)
        n, _ = odd_step(n, 5)
    require(orbit == [13,33,83] and n == 13, '5x+1 cycle')
    require(sum(weight((n+1)//2) for n in orbit) == F(69,2048), 'persistent source atoms')
    for modulus in (1,2,3,10,13,26,27,28):
        for length in range(1,9):
            n0 = modulus*2**(length+1)*17-1
            n = n0
            for j in range(length):
                n, a = odd_step(n)
                require(a == 1 and n % modulus == (-1) % modulus, 'fixed-residue hostile')
                require(n+1 == 3**(j+1)*2**(length-j)*modulus*17, 'marked-source clock')
            require(n > n0, 'growth across finite-phase return')
    bank = list(range(1,1025))
    current = {m:2*m-1 for m in bank}
    alive = {m for m in bank if m != 1}
    lower = []
    first_violation = None
    for t in range(65):
        mass = sum((weight(m) for m in alive), F(0))
        lower.append(str(mass))
        if mass > F(1, 2**t) and first_violation is None:
            first_violation = dict(odd_clock=t, exact_lower_mass=str(mass), proposed_bound=str(F(1,2**t)))
        for m in list(alive):
            current[m], _ = odd_step(current[m])
            if current[m] == 1:
                alive.remove(m)
    require(first_violation is not None, 'geometric rate is not free')
    # Independent native trace for a named slow atom.
    n, j, divisions = 27, 0, 0
    while n != 1:
        z = 3*n+1; a=0
        while z % 2 == 0:
            z//=2; a+=1
        n=z; j+=1; divisions+=a
    require(j == 41 and divisions == 70 and weight(14) == F(1,128), '27 retains its atom through odd clock40')
    return dict(cycle5=orbit, cycle_mass=str(F(69,2048)), geometric_half_mass_hostile=first_violation,
                source27=dict(odd_clock=j, division_clock=divisions, source_mass=str(weight(14))),
                survivor_lower_bounds=lower,
                scope='Finite mass lower bounds only; universal 3x+1 decay remains open')


def atomic_failure_checks():
    blocks = [(1,1,1),(4,6,179),(7,11,16177),(7,11,12289),(7,11,9697)]
    q = sum((F(1,4**d) for _,d,_ in blocks), F(0))
    require(q == F(1049603,4194304), 'baseline transfer coefficient')
    rows, states = [], [(0,0,0)]
    for k in range(7):
        base, spikes, least = F(0), F(0), None
        for j,d,c in states:
            rn = ((2**d-c)*pow(3**j,-1,2**(d+1))) % 2**(d+1)
            r = ((rn+1)//2) % 2**d
            base += F(1,4**d)
            if r:
                spikes += weight(r)
            least = rn if least is None else min(least,rn)
        require(base == q**k, 'complete-word baseline recurrence')
        rows.append(dict(blocks=k, baseline=str(base), spikes=str(spikes),
                         atomic_mass=str(base+spikes), least_source=least))
        states = [(j+jj,d+dd,3**jj*c+2**d*cc)
                  for j,d,c in states for jj,dd,cc in blocks] if k < 6 else states
    # Independent literal source simulation, with exact omitted-source tail.
    partial = [F(0)]*7
    for m in range(1,1024):
        n=2*m-1
        partial[0] += weight(m)
        for k in range(1,7):
            source=n; queue=1
            while queue:
                n,a=odd_step(n)
                queue += -1+3*int(a % 3 == 0)
            if n <= source:
                break
            partial[k] += weight(m)
    for k,row in enumerate(rows):
        exact=F(row['atomic_mass'])
        require(partial[k] <= exact <= partial[k]+F(1,1024), 'independent atom-bank enclosure')
        if k:
            require(exact <= F(rows[k-1]['atomic_mass']), 'nested survivor mass')
    for length in range(9):
        for r in range(2**length):
            mass=F(1,4**length)+(weight(r) if r else 0)
            children=2*F(1,4**(length+1))+(weight(r) if r else 0)+weight(r+2**length)
            require(mass == children, 'exact binary refinement conservation')
            require(weight(r+2**length) == F(1,2*4**length), 'new atom born from half tail')
    return dict(q=str(q), rows=rows, independent_source_indices=[1,1023], tail='1/1024')


def main():
    result = dict(spectral=spectral_checks(), projective=projective_checks(),
                  arithmetic=arithmetic_checks(), descent=descent_checks(),
                  atomic_failure=atomic_failure_checks())
    result.update(status='PASS', explicit_checks=CHECKS,
                  source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
