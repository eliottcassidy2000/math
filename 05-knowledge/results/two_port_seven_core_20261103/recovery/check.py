"""Exact controls for recovered sum/difference and Petersen predicates."""
from argparse import ArgumentParser
from collections import deque
from hashlib import sha256
from itertools import combinations, product
from math import isqrt
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
SOURCE = Path(__file__).read_bytes()
CHECKS = 0
P_EDGES = ((0,1),(0,4),(0,5),(1,2),(1,6),(2,3),(2,7),(3,4),
           (3,8),(4,9),(5,7),(5,8),(6,8),(6,9),(7,9))


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def adjacent(n, edges):
    a = [set() for _ in range(n)]
    for u, v in edges:
        a[u].add(v); a[v].add(u)
    return a


def bipartition(a):
    colors, parent = {}, {}
    for root in range(len(a)):
        if root in colors:
            continue
        colors[root] = 0; parent[root] = None; queue = deque([root])
        while queue:
            u = queue.popleft()
            for v in sorted(a[u]):
                if v not in colors:
                    colors[v] = 1-colors[u]; parent[v] = u; queue.append(v)
                elif colors[v] == colors[u]:
                    left = []; x = u
                    while x is not None:
                        left.append(x); x = parent[x]
                    right = []; x = v
                    while x not in left:
                        right.append(x); x = parent[x]
                    return None, left[:left.index(x)+1] + right[::-1]
    return colors, None


def gauge_controls():
    count = 0
    for n in range(1,6):
        pairs = list(combinations(range(n),2))
        for mask in range(1 << len(pairs)):
            edges = [e for k,e in enumerate(pairs) if mask >> k & 1]
            colors, _ = bipartition(adjacent(n,edges))
            gauges = [s for s in product((-1,1), repeat=n)
                      if all(s[u] == -s[v] for u,v in edges)]
            require(bool(gauges) == (colors is not None), 'gauge iff bipartite')
            for signs in gauges:
                for u,v in edges:
                    left = [signs[w] * int(w in (u,v)) for w in range(n)]
                    right = [signs[u] * (int(w==u)-int(w==v)) for w in range(n)]
                    require(left == right, 'unsigned/signed incidence identity')
            count += 1
    first = None
    for n in range(1,47):
        edges = [(u-1,v-1) for u,v in combinations(range(1,n+1),2)
                 if isqrt(u+v)**2 == u+v]
        _, cycle = bipartition(adjacent(n,edges))
        if cycle is not None:
            require(len(cycle)%2==1 and all(tuple(sorted((u,v))) in edges
                    for u,v in zip(cycle,cycle[1:]+cycle[:1])), 'odd square-sum cycle')
            first = {'n':n,'cycle':[u+1 for u in cycle]}; break
    require(first is not None, 'square-sum odd-cycle hostile found')
    require(((-1)-1)**2 != 0 and ((-1)**2-1**2)**2 == 0, 'mirror-factor hostile')
    return {'all_graphs_orders1_to5': count, 'first_nonbipartite_square_sum':first}


# Integer coefficient arithmetic in Z[zeta_5], Phi_5=1+x+x^2+x^3+x^4.
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def neg(a): return tuple(-x for x in a)
def mul(a,b):
    c = [0]*7
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j] += x*y
    for k in range(6,3,-1):
        value = c[k]
        for j in range(1,5): c[k-j] -= value
    return tuple(c[:4])


ONE = (1,0,0,0)
ROOTS = [ONE]
for _ in range(4): ROOTS.append(mul(ROOTS[-1],(0,1,0,0)))


def conjugate(a):
    result = (0,0,0,0)
    for k,x in enumerate(a):
        result = add(result,tuple(x*y for y in ROOTS[-k % 5]))
    return result


def norm(a): return mul(a,conjugate(a))


def hamilton_paths(a, cycles=False, first_only=False):
    n = len(a); answers = []
    def walk(word, used):
        if len(word)==n:
            if not cycles or word[0] in a[word[-1]]:
                answers.append(word[:])
            return bool(answers) and first_only
        for v in sorted(a[word[-1]]):
            if not (used >> v & 1) and walk(word+[v],used | 1<<v):
                return True
        return False
    for root in ([0] if cycles else range(n)):
        if walk([root],1<<root): break
    return answers


def petersen_geometry():
    pairs = list(combinations(range(5),2))
    values = [add(ROOTS[u],ROOTS[v]) for u,v in pairs]
    for i,j in combinations(range(10),2):
        require((norm(add(values[i],values[j])) == ONE) == set(pairs[i]).isdisjoint(pairs[j]),
                'unit-sum Kneser predicate')
    points = values + [neg(x) for x in values]
    require(len(set(points))==20, 'signed two-sums are distinct')
    edges = []
    for u,v in combinations(range(20),2):
        actual = norm(add(points[u],neg(points[v]))) == ONE
        expected = (u//10 != v//10 and set(pairs[u%10]).isdisjoint(pairs[v%10]))
        require(actual == expected, 'faithful unit-distance double cover')
        if actual: edges.append((u,v))
    a = adjacent(20,edges)
    require(len(edges)==30 and all(len(row)==3 for row in a), 'cubic double cover')
    midpoint_lengths = {norm(add(points[v],neg(points[w])))
                        for neighbors in a for v,w in combinations(sorted(neighbors),2)}
    require(len(midpoint_lengths)>1, 'midpoint line-graph map has unequal edge lengths')
    cycle = hamilton_paths(a,cycles=True,first_only=True)[0]
    require(len(set(cycle))==20, 'double-cover Hamilton cycle')
    return {'cyclotomic_coefficients':points,'base_pairs':pairs,'unit_edges':edges,
            'midpoint_squared_lengths_times4':sorted(midpoint_lengths),
            'hamilton_cycle':cycle,'coordinate_convention':'first ten s_ab=zeta^a+zeta^b; last ten -s_ab; coefficient basis 1,zeta,zeta^2,zeta^3'}


def perfect_matchings(a):
    answers = []
    def recur(unused, edges):
        if not unused:
            answers.append(edges); return
        u = min(unused)
        for v in sorted(a[u]&unused):
            recur(unused-{u,v},edges+[(u,v)])
    recur(set(range(len(a))),[])
    return answers


def color_petersen(lists):
    a = adjacent(10,P_EDGES); coloring = {}
    for u in (9,4,6,7):
        coloring[u] = min(set(lists[u])-{coloring[v] for v in a[u] if v in coloring})
    cycle = [0,1,2,3,8,5]
    reduced = {u:set(sorted(set(lists[u])-{coloring[v] for v in a[u] if v in coloring})[:2]) for u in cycle}
    require(all(len(v)==2 for v in reduced.values()), 'two colors remain on even cycle')
    if all(reduced[u]==reduced[cycle[0]] for u in cycle):
        pair = sorted(reduced[cycle[0]])
        coloring.update({u:pair[k%2] for k,u in enumerate(cycle)})
    else:
        i = next(i for i,u in enumerate(cycle) if reduced[u]-reduced[cycle[i-1]])
        cycle = cycle[i:]+cycle[:i]
        coloring[cycle[0]] = min(reduced[cycle[0]]-reduced[cycle[-1]])
        for before,u in zip(cycle,cycle[1:]):
            coloring[u] = min(reduced[u]-{coloring[before]})
    require(all(coloring[u] in lists[u] for u in range(10)) and
            all(coloring[u]!=coloring[v] for u,v in P_EDGES), 'actual list coloring')
    return [coloring[u] for u in range(10)]


def petersen_controls():
    a = adjacent(10,P_EDGES)
    paths = hamilton_paths(a); cycles = hamilton_paths(a,cycles=True)
    require(paths and not cycles, 'Petersen has paths but no Hamilton cycle')
    matchings = perfect_matchings(a); lengths = []
    for matching in matchings:
        remaining = [set(row) for row in a]
        for u,v in matching: remaining[u].remove(v); remaining[v].remove(u)
        unused = set(range(10)); sizes=[]
        while unused:
            stack=[min(unused)]; component=set(stack)
            while stack:
                u=stack.pop()
                for v in remaining[u]-component: component.add(v);stack.append(v)
            sizes.append(len(component)); unused-=component
        lengths.append(sorted(sizes))
    require(len(matchings)==6 and all(x==[5,5] for x in lengths), 'all two-factors are two odd pentagons')
    branchsets = [{i,i+5} for i in range(5)]
    require(all(any(v in a[u] for u in left for v in right) for left,right in combinations(branchsets,2)), 'explicit K5 minor')
    suppressed = {u:set(a[u])-{0} for u in range(1,10)}
    for u in (1,4,5):
        v,w=sorted(suppressed.pop(u));suppressed[v].remove(u);suppressed[w].remove(u)
        suppressed[v].add(w);suppressed[w].add(v)
    left,right={2,8,9},{3,6,7}
    require(all(suppressed[u]==(right if u in left else left) for u in left|right), 'explicit K33 subdivision')
    choices=list(combinations(range(5),3));count=0
    for groups in product(choices,repeat=3):
        color_petersen([groups[u%3] for u in range(10)]);count+=1
    return {'directed_hamilton_paths':len(paths),'hamilton_cycles':0,'first_path':paths[0],
            'perfect_matchings':matchings,'complement_cycle_lengths':lengths,
            'k5_branchsets':[sorted(s) for s in branchsets], 'k33_sides':[sorted(left),sorted(right)],
            'list_controls':count,'list_control_rule':'every three-element subset of {0,...,4} assigned independently to the three vertex congruence classes modulo3',
            'vertex_three_coloring':color_petersen([tuple(range(3))]*10)}


def main():
    parser=ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    result={'gauge':gauge_controls(),'signed_pentagon':petersen_geometry(),'petersen':petersen_controls()}
    require(Path(__file__).read_bytes()==SOURCE,'source stable')
    result.update(status='FINITE-EXACT PASS',checks=CHECKS,source_sha256=sha256(SOURCE).hexdigest())
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(text,encoding='utf8',newline='\n')
    else:print(text,end='')


if __name__=='__main__':main()
