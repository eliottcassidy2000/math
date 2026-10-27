"""Independent exact audit of all-period golden fibres. Python standard library only.
No import/use of producer code or floating point. The completeness of the seven
carries is proved in audit.md; this verifies their entire finite transition graph.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import product
import json
CHECKS = 0

def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)
CARRIES = ((-2, 1), (-1, 1), (-1, 0), (0, 0), (1, 0), (1, -1), (2, -1))
vertices = tuple(((a, b, p, q) for (a, b), (p, q) in product(CARRIES, product(range(2), repeat=2))))
V = set(vertices)
edges = {}
for a, b, p, q in vertices:
    edges[a, b, p, q] = tuple(((b - u + v, a + b, u, v) for u, v in product(range(2), repeat=2) if not (p and u) and (not (q and v)) and ((b - u + v, a + b, u, v) in V)))
serial = 0
ids = {}
lows = {}
stack = []
active = set()
components = []

def visit(x):
    global serial
    ids[x] = lows[x] = serial
    serial += 1
    stack.append(x)
    active.add(x)
    for y in edges[x]:
        if y not in ids:
            visit(y)
            lows[x] = min(lows[x], lows[y])
        elif y in active:
            lows[x] = min(lows[x], ids[y])
    if ids[x] == lows[x]:
        component = []
        while True:
            y = stack.pop()
            active.remove(y)
            component.append(y)
            if x == y:
                break
        components.append(tuple(sorted(component)))
for x in vertices:
    if x not in ids:
        visit(x)
cyclic = sorted((c for c in components if len(c) > 1 or c[0] in edges[c[0]]))
expected = sorted((tuple(sorted(c)) for c in (((0, 0, 0, 0), (0, 0, 1, 1)), ((1, 0, 0, 0), (-1, 1, 1, 0)), ((-1, 0, 0, 0), (1, -1, 0, 1)), ((-2, 1, 1, 0), (2, -1, 0, 1)))))
check(cyclic == expected, 'cyclic == expected')
internal = []
for c in cyclic:
    E = [(x, y) for x in c for y in edges[x] if y in c]
    if c != tuple(sorted(((0, 0, 0, 0), (0, 0, 1, 1)))):
        check(len(E) == 2 and all((x != y for x, y in E)), 'len(E)==2 and all(x!=y for x,y in E)')
    else:
        check(len(E) == 3 and all((x[:2] == y[:2] == (0, 0) for x, y in E)), 'len(E)==3 and all(x[:2]==y[:2]==(0,0) for x,y in E)')
    internal.append(E)
owner = {x: i for i, c in enumerate(components) for x in c}
cedges = {i: set() for i in range(len(components))}
for x in V:
    for y in edges[x]:
        if owner[x] != owner[y]:
            cedges[owner[x]].add(owner[y])
ranks = {}

def rank(i):
    if i not in ranks:
        ranks[i] = 0 if not cedges[i] else 1 + max((rank(j) for j in cedges[i]))
    return ranks[i]
for i in cedges:
    rank(i)
check(all((ranks[i] > ranks[j] for i in cedges for j in cedges[i])), 'all(ranks[i]>ranks[j] for i in cedges for j in cedges[i])')

def plus(x, y):
    return (x[0] + y[0], x[1] + y[1])

def times(x, y):
    a, b = x
    c, d = y
    return (a * c + b * d, a * d + b * c + b * d)

def pp(n):
    ans = (1, 0)
    for _ in range(n):
        ans = times(ans, (0, 1))
    return ans

def read(w):
    ans = (0, 0)
    for b in w:
        ans = plus(times(ans, (0, 1)), (b, 0))
    return ans

def denom(n):
    a, b = pp(n)
    return (a - 1, b)

def norm(x):
    a, b = x
    return a * a + a * b - b * b

def quotient(x, y):
    a, b = times(x, (y[0] + y[1], -y[1]))
    d = norm(y)
    return (Fraction(a, d), Fraction(b, d))

def eqmod(x, y, n):
    q = quotient((x[0] - y[0], x[1] - y[1]), denom(n))
    return all((t.denominator == 1 for t in q))

def rep(x, n):
    return tuple((t % 1 for t in quotient(x, denom(n))))

def legal(w):
    return all((not (w[i] and w[(i + 1) % len(w)]) for i in range(len(w))))

def cycle(w):
    k = h = B = 0
    for b in w:
        if b:
            k += 1
            B = 3 * B + 2 ** h
        else:
            h += 1
    x = Fraction(B, 2 ** h - 3 ** k)
    y = x
    for b in w:
        check(y.denominator % 2 == 1 and y.numerator % 2 == b, 'y.denominator%2==1 and y.numerator%2==b')
        y = 3 * y + 1 if b else y / 2
    check(x == y, 'x==y')
    return x
Lprev, L = (2, 1)
census = []
for n in range(1, 17):
    if n > 1:
        Lprev, L = (L, Lprev + L)
    groups = defaultdict(list)
    sources = {}
    for w in product(range(2), repeat=n):
        if not legal(w):
            continue
        groups[rep(read(w), n)].append(w)
        x = cycle(w)
        check(x not in sources, 'x not in sources')
        sources[x] = w
        check(eqmod(read(w[1:] + w[:1]), times((0, 1), read(w)), n), 'eqmod(read(w[1:]+w[:1]),times((0,1),read(w)),n)')
    check(sum(map(len, groups.values())) == L, 'sum(map(len,groups.values()))==L')
    check(len(groups) == abs(norm(denom(n))) == L - 1 - (-1) ** n, 'len(groups)==abs(norm(denom(n)))==L-1-(-1)**n')
    bad = [set(g) for g in groups.values() if len(g) > 1]
    if n % 2:
        check(not bad, 'not bad')
    else:
        special = {tuple([0] * n), tuple([0, 1] * (n // 2)), tuple([1, 0] * (n // 2))}
        check(bad == [special], 'bad==[special]')
        check(all((eqmod(read(w), (0, 0), n) for w in special)), 'all(eqmod(read(w),(0,0),n) for w in special)')
        check({cycle(w) for w in special} == {Fraction(0), Fraction(-1), Fraction(-2)}, '{cycle(w) for w in special}=={Fraction(0),Fraction(-1),Fraction(-2)}')
    census.append({'n': n, 'words': len(sources), 'quotient': len(groups)})
lifts = []
for m in range(1, 11):
    for r in range(1, 9):
        Q = (0, 0)
        for j in range(r):
            Q = plus(Q, pp(m * j))
        check(times(denom(m), Q) == denom(m * r), 'times(denom(m),Q)==denom(m*r)')
        check(eqmod(Q, (r, 0), m), 'eqmod(Q,(r,0),m)')
        check(eqmod(times(Q, Q), times((r, 0), Q), m * r), 'eqmod(times(Q,Q),times((r,0),Q),m*r)')
        for x in product(range(-3, 4), repeat=2):
            check(quotient(times(Q, x), denom(m * r)) == quotient(x, denom(m)), 'quotient(times(Q,x),denom(m*r))==quotient(x,denom(m))')
        lifts.append([m, r])
m, r = (5, 11)
Q = (0, 0)
for j in range(r):
    Q = plus(Q, pp(m * j))
check(abs(norm(denom(55))) == 312119004989, 'abs(norm(denom(55)))==312119004989')
eleven = [times(Q, (a, 0)) for a in range(11)]
check(len({rep(x, 55) for x in eleven}) == 11, 'len({rep(x,55) for x in eleven})==11')
check(all((eqmod(x, (0, 0), 5) for x in eleven)), 'all(eqmod(x,(0,0),5) for x in eleven)')
check(all((eqmod(times(x, y), (0, 0), 55) for x in eleven for y in eleven)), 'all(eqmod(times(x,y),(0,0),55) for x in eleven for y in eleven)')
check(read((1, 0, 0)) == read((0, 1, 1)) and (not legal((0, 1, 1))), 'read((1,0,0))==read((0,1,1)) and not legal((0,1,1))')
check(eqmod(read((1, 0, 1)), read((0, 1, 0)), 3) and (not legal((1, 0, 1))), 'eqmod(read((1,0,1)),read((0,1,0)),3) and not legal((1,0,1))')
# Independent Smith control: enumerate quotient elements using field division,
# then count the [r] kernel directly rather than deriving it from Smith factors.
from math import gcd
smith_rows = []
for m in range(1, 13):
    a, b = denom(m)
    N, d = abs(norm((a,b))), gcd(a,b)
    check(N % (d*d) == 0, 'first Smith factor divides the second')
    representatives = {rep((0,0),m):(0,0)}
    queue = [(0,0)]
    for x in queue:
        for step in ((1,0),(0,1)):
            y=plus(x,step); key=rep(y,m)
            if key not in representatives:
                representatives[key]=y;queue.append(y)
    check(len(representatives)==N, 'direct quotient group order')
    rows=[]
    for r in range(1,23):
        actual=sum(eqmod(times((r,0),x),(0,0),m) for x in representatives.values())
        predicted=gcd(r,d)*gcd(r,N//d)
        check(actual==predicted, 'direct kernel equals Smith product')
        check((actual==1)==(gcd(r,N)==1), 'lossless iff coprime')
        rows.append([r,actual,N//actual])
    smith_rows.append({'m':m,'smith':[d,N//d],'return_kernel_image':rows})
# The 2/3/11 operator palette, including Paley orientation and two 5-cycles.
H={pow(x,2,11) for x in range(1,11)}
check(H=={1,3,4,5,9}, 'quadratic residues mod 11')
check({2*x%11 for x in H}=={(-x)%11 for x in H}, 'two swaps Paley colours')
check(all(((3*x)%11 in H)==(x in H) for x in range(1,11)), 'three preserves Paley arcs')
check(pow(3,5,11)==1 and all(pow(3,j,11)!=1 for j in range(1,5)), 'three order five')
check(pow(4,4,11)==3, 'three is four standard Collatz shifts on R5')

receipt = {'status': 'PASS', 'graph_vertices': len(V), 'graph_edges': sum(map(len, edges.values())), 'cyclic_sccs': cyclic, 'internal_recurrent_edges': internal, 'components_and_ranks': [{'vertices': c, 'rank': ranks[i]} for i, c in enumerate(components)], 'all_edges': [[x, edges[x]] for x in sorted(V)], 'independent_census': census, 'lift_pairs_checked': lifts, 'floating_point': False, 'smith_rows': smith_rows, 'smith_cases': 264, 'checks': CHECKS, 'elevenfold_ring_order': abs(norm(denom(55))), 'elevenfold_distinct_lifts': len({rep(x, 55) for x in eleven}), 'scope': 'All-n conclusion uses carry bound and graph proof; finite census is only a control.'}
print(json.dumps(receipt, indent=2))
