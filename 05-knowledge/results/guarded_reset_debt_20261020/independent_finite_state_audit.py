#!/usr/bin/env python3
"""Independent exact graph/LP/guard controls; standard library only."""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from random import Random
import json

checks = 0

def need(test, label):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(label)

def cycle_data(n, edges):
    # Enumerate vertex permutations and all parallel-edge choices. This is
    # independent of the producer's DFS simple-cycle enumeration.
    lookup = {}
    for edge in edges:
        lookup.setdefault(edge[:2], []).append(edge)
    output = []
    for length in range(1,n+1):
        for vertices in permutations(range(n), length):
            if vertices[0] != min(vertices):
                continue
            lists = [lookup.get((vertices[j],vertices[(j+1)%length]), [])
                     for j in range(length)]
            for cycle in product(*lists):
                multiplier = F(1)
                paid = 0
                for _,_,l,d in cycle:
                    multiplier *= l
                    paid += d > 0
                output.append((multiplier,paid))
    return output

def criterion(data):
    return all(l >= 1 and (l > 1 or p == 0) for l,p in data)

def gaussian(matrix, rhs):
    a = [list(row)+[c] for row,c in zip(matrix,rhs)]
    n = len(a)
    for column in range(n):
        pivot = next((i for i in range(column,n) if a[i][column]), None)
        if pivot is None:
            return None
        a[column], a[pivot] = a[pivot], a[column]
        scale = a[column][column]
        a[column] = [x/scale for x in a[column]]
        for i in range(n):
            if i != column:
                scale = a[i][column]
                a[i] = [x-scale*y for x,y in zip(a[i],a[column])]
    return tuple(row[-1] for row in a)

def vertex_feasible(n, edges):
    # Independent exact LP oracle: enumerate all vertices of b>=1 plus the
    # edge inequalities. A nonempty polyhedron contained in b>=1 has a
    # vertex: minimize sum(b), whose sublevel sets here are compact.
    rows = []
    for i in range(n):
        a = [F(0)]*n
        a[i] = F(1)
        rows.append((tuple(a),F(1)))
    for u,v,l,d in edges:
        a = [F(0)]*n
        a[u] -= 1
        a[v] += l
        rows.append((tuple(a),d))
    for indices in combinations(range(len(rows)),n):
        b = gaussian([rows[i][0] for i in indices],[rows[i][1] for i in indices])
        if b is not None and all(sum(x*y for x,y in zip(a,b)) >= c for a,c in rows):
            return b
    return None

def valid(edges,b):
    return all(x > 0 for x in b) and all(l*b[v]-b[u] >= d for u,v,l,d in edges)

def rational_construct(n, edges, data):
    N = 1
    for l,p in data:
        if p:
            ratio = F(2**p-1)/(l-1)
            N = max(N,(ratio.numerator+ratio.denominator-1)//ratio.denominator)
    kappa = F(N+1,N)
    for l,p in data:
        need(kappa**p <= l, 'rational perturbation preserves every cycle')
    # Multiplicative Floyd-Warshall, distinct from producer Bellman-Ford.
    dist = [[None]*n for _ in range(n)]
    for i in range(n):
        dist[i][i] = F(1)
    for u,v,l,d in edges:
        mu = l/(kappa if d > 0 else 1)
        dist[u][v] = mu if dist[u][v] is None else min(dist[u][v],mu)
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if dist[i][k] is not None and dist[k][j] is not None:
                    candidate = dist[i][k]*dist[k][j]
                    if dist[i][j] is None or candidate < dist[i][j]:
                        dist[i][j] = candidate
    b = [1/min(dist[i][j] for i in range(n) if dist[i][j] is not None) for j in range(n)]
    scale = F(1)
    for u,v,l,d in edges:
        slack = l*b[v]-b[u]
        need(slack >= 0 and (not d or slack > 0), 'constructive original edge slack')
        if d:
            scale = max(scale,d/slack)
    b = [scale*x for x in b]
    need(valid(edges,b), 'fully rational barrier certificate')
    return b

def compare(n,edges):
    data = cycle_data(n,edges)
    good = criterion(data)
    feasible = vertex_feasible(n,edges)
    need(good == (feasible is not None), 'cycle criterion iff exact vertex LP')
    if good:
        need(valid(edges,feasible), 'independent LP certificate')
        rational_construct(n,edges,data)
    return good

options = [None]+list(product([F(2,3),F(1),F(3,2)],[F(0),F(2)]))
exhaustive = feasible_count = 0
for selected in product(options,repeat=4):
    edges = [(u,v,*item) for (u,v),item in zip(product(range(2),repeat=2),selected)
             if item is not None]
    feasible_count += compare(2,edges)
    exhaustive += 1

rng = Random(904712)
random_counts = {3:400,4:50}
for n,count in random_counts.items():
    for _ in range(count):
        edges = [(rng.randrange(n),rng.randrange(n),
                  rng.choice([F(1,3),F(1,2),F(2,3),F(1),F(3,2),F(2),F(3)]),
                  rng.choice([F(0),F(1),F(7,3)]))
                 for _ in range(rng.randrange(7))]
        compare(n,edges)

for n in range(1,6):
    need(compare(n,[]), 'empty graph has positive rational barriers')
for edges,want in [([(0,0,F(1,2),F(0))],False),
                   ([(0,0,F(1),F(1))],False),
                   ([(0,0,F(1),F(0))],True),
                   ([(0,1,F(1,2),F(0)),(1,0,F(2),F(0))],True),
                   ([(0,1,F(1,2),F(1)),(1,0,F(2),F(0))],False),
                   ([(0,1,F(1,2),F(0)),(1,0,F(2),F(1))],False),
                   ([(0,1,F(1,2),F(1)),(1,0,F(3),F(2))],True),
                   ([(0,0,F(1),F(0)),(0,1,F(1),F(5))],True)]:
    need(compare(2,edges) == want, 'boundary hostile verdict')

# General telescoping, without assuming barriers are feasible.
for _ in range(1200):
    length = rng.randrange(1,16)
    b = [F(rng.randrange(1,25),rng.randrange(1,7)) for _ in range(length+1)]
    edges = [(F(rng.randrange(1,12),rng.randrange(1,12)),F(rng.randrange(6),3))
             for _ in range(length)]
    multiplier = F(1)
    debt = slack = F(0)
    x0 = x = F(rng.randrange(1,30),3)
    for j,(l,d) in enumerate(edges):
        slack += multiplier*(l*b[j+1]-b[j]-d)
        debt += multiplier*d
        multiplier *= l
        x = (x+d)/l
    need(multiplier*b[-1]-b[0]-debt == slack, 'state-dependent telescoping equality')
    need(x-b[-1] == (x0-b[0]-slack)/multiplier, 'exact barrier-excess equality')

def valuation(n):
    a = 0
    while n % 2 == 0:
        a += 1
        n //= 2
    return a

def execute(n,word):
    for a in word:
        if n <= 1:
            return None
        value = 3*n+1
        if valuation(value) != a:
            return None
        n = value//2**a
    return n

def physical(word):
    slope,intercept = F(1),F(0)
    for a in word:
        slope,intercept = F(3,2**a)*slope, F(3,2**a)*intercept+F(1,2**a)
    return 1/slope,intercept/slope

prefixes = 0
for K in range(1,129):
    for t in range(1,6):
        initial = n = 2**(K+1)*t-1
        for j in range(K):
            need(n == 3**j*2**(K+1-j)*t-1, 'all-ones exact endpoint')
            need(n > 1 and valuation(3*n+1) == 1, 'all-ones ROOT-free exact valuation')
            nxt = (3*n+1)//2
            need(nxt > n >= initial, 'all-ones strict growth')
            n = nxt
        prefixes += 1

reset_words = []
for k in range(1,6):
    for word in product(range(1,5),repeat=k):
        l,d = physical(word)
        if all(a == 1 for a in word):
            need(l < 1 and d > 0 and d > l-1, 'all-ones cannot satisfy unit budget')
        else:
            reset_words.append(word)
for t in range(1,7):
    n = 2**6*t-1
    need(all(execute(n,word) is None for word in reset_words), 'finite reset bank uncovered family')

physical_words = [(1,),(2,),(3,),(1,1),(1,2),(2,1),(2,2),(1,1,3)]
physical_graphs = physical_feasible = 0

def path_exists(state,n,depth,labels):
    if depth == 0:
        return True
    for u,v,word in labels:
        if u == state:
            output = execute(n,word)
            if output is not None and path_exists(v,output,depth-1,labels):
                return True
    return False

for selected in product([None]+physical_words,repeat=4):
    labels = [(u,v,word) for (u,v),word in zip(product(range(2),repeat=2),selected)
              if word is not None]
    edges = [(u,v,*physical(word)) for u,v,word in labels]
    physical_graphs += 1
    if criterion(cycle_data(2,edges)):
        physical_feasible += 1
        rational_construct(2,edges,cycle_data(2,edges))
        L = max([1]+[len(word) for _,_,word in labels])
        for t in range(1,4):
            n = 2**(2*L+1)*t-1
            need(all(not path_exists(state,n,2,labels) for state in range(2)),
                 'finite-state no-go from every initial state')

for V in range(1,6):
    n = 2**(V+1)-1
    need(n > 1 and not path_exists(0,n,V,[]), 'empty-edge L=1 boundary')
need(valid([(0,1,*physical((1,))),(1,0,*physical((2,2)))],[F(5),F(8)]),
     'physical state-dependent barrier example')
need(execute(31,(1,)) == 47 and execute(47,(2,2)) is None, 'next guard fails despite barriers')
need(F(3+3,2) == 3, 'strict-cycle abstract finite-head non-target')

print(json.dumps({
    'status':'PASS', 'checks':checks,
    'exact_LP_method':'enumeration of rational basic feasible vertices, not Fourier-Motzkin',
    'rational_construction_method':'cycle perturbation plus multiplicative Floyd-Warshall',
    'exhaustive_two_state_graphs':exhaustive,
    'exhaustive_two_state_feasible':feasible_count,
    'sampled_multigraphs_by_vertices':random_counts,
    'random_path_identities':1200,
    'all_ones_source_prefixes':prefixes,
    'maximum_all_ones_prefix':128,
    'reset_bank_words':len(reset_words),
    'exhaustive_physical_two_state_graphs':physical_graphs,
    'physical_two_state_graphs_with_barriers':physical_feasible,
    'global_Collatz_coverage_or_convergence_claim':False,
},indent=2,sort_keys=True))
