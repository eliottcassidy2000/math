#!/usr/bin/env python3
"""Independent six-port Hamiltonian gadget audit; Python standard library only.
No imports from the producer or inherited verifiers. All checks survive python -O.
"""
import hashlib
import itertools
import json
from pathlib import Path

PAIRS = (
    ((0, 1), (2, 3), (4, 5)),
    ((0, 1), (2, 5), (3, 4)),
    ((0, 3), (1, 2), (4, 5)),
    ((0, 5), (1, 2), (3, 4)),
    ((0, 5), (1, 4), (2, 3)),
)
BOUNDARY = tuple('r%d' % i for i in range(6))

def require(value, text):
    if not value:
        raise RuntimeError(text)

def canon_edge(u, v):
    require(u != v, 'loop')
    return tuple(sorted((u, v)))

def graph(vertices, edges):
    vertices = tuple(vertices)
    edges = tuple(sorted(canon_edge(*e) for e in edges))
    require(len(vertices) == len(set(vertices)), 'duplicate vertex')
    require(len(edges) == len(set(edges)), 'parallel edge')
    require(all(u in vertices and v in vertices for u, v in edges), 'missing endpoint')
    return vertices, edges

def patch(hub):
    v = BOUNDARY + tuple('a%d' % i for i in range(6)) + tuple('v%d' % i for i in range(6))
    e = []
    for i in range(6):
        e += [('r%d' % i, 'a%d' % i), ('a%d' % i, 'v%d' % i), ('v%d' % i, 'v%d' % ((i + 1) % 6))]
    if hub:
        v += ('c',)
        e += [('c', 'v%d' % i) for i in range(6)]
    return graph(v, e)

def degree_subsets(g, boundary=BOUNDARY):
    """Enumerate EVERY edge subset satisfying degree bounds, with safe pruning.
    Required interior degrees are two. Ports range from zero through two.
    This does not assume spokes forced and does not assume acyclicity.
    """
    vertices, edges = g
    ix = {v: i for i, v in enumerate(vertices)}
    ee = [(ix[u], ix[v]) for u, v in edges]
    mandatory = [v not in boundary for v in vertices]
    deg = [0] * len(vertices)
    remaining = [0] * len(vertices)
    for u, v in ee:
        remaining[u] += 1
        remaining[v] += 1
    def visit(j, selected):
        if j == len(ee):
            if all(not mandatory[i] or deg[i] == 2 for i in range(len(vertices))):
                yield tuple(selected)
            return
        u, v = ee[j]
        remaining[u] -= 1
        remaining[v] -= 1
        for take in (False, True):
            if take:
                deg[u] += 1
                deg[v] += 1
            if all(deg[t] <= 2 and (not mandatory[t] or deg[t] + remaining[t] >= 2) for t in (u, v)):
                if take:
                    selected.append(edges[j])
                yield from visit(j + 1, selected)
                if take:
                    selected.pop()
            if take:
                deg[u] -= 1
                deg[v] -= 1
        remaining[u] += 1
        remaining[v] += 1
    yield from visit(0, [])

def components(vertices, edges):
    adj = {v: [] for v in vertices}
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    groups = []
    seen = set()
    for v in vertices:
        if v in seen:
            continue
        stack = [v]
        seen.add(v)
        group = []
        while stack:
            u = stack.pop()
            group.append(u)
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        groups.append(group)
    return adj, groups

def profile(g, boundary=BOUNDARY):
    out = {}
    all_subsets = []
    cycle_rejections = 0
    for selected in degree_subsets(g, boundary):
        adj, groups = components(g[0], selected)
        # A connected max-degree-two component is a cycle exactly when every
        # one of its vertices has degree two.
        if any(all(len(adj[v]) == 2 for v in group) for group in groups):
            cycle_rejections += 1
            continue
        d = tuple(len(adj[v]) for v in boundary)
        matching = []
        for group in groups:
            if len(group) == 1 and not adj[group[0]]:
                require(group[0] in boundary, 'isolated mandatory interior')
                continue
            ends = [boundary.index(v) for v in group if len(adj[v]) == 1]
            require(len(ends) == 2, 'invalid path component')
            matching.append(tuple(sorted(ends)))
        m = tuple(sorted(matching))
        out[d, m] = out.get((d, m), 0) + 1
        all_subsets.append({'degree': d, 'matching': m, 'edges': selected})
    return out, all_subsets, cycle_rejections

def connectivity(a, b):
    vertices = tuple(sorted({x for pair in a for x in pair}))
    return int(len(components(vertices, a + b)[1]) == 1)

K = [[connectivity(a, b) for b in PAIRS] for a in PAIRS]

def matrix_times(m, v):
    return [sum(a * b for a, b in zip(row, v)) for row in m]

def parity_rank(m):
    a = [[x & 1 for x in row] for row in m]
    row = 0
    for col in range(len(a[0])):
        pivot = next((j for j in range(row, len(a)) if a[j][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        for j in range(row + 1, len(a)):
            if a[j][col]:
                a[j] = [x ^ y for x, y in zip(a[j], a[row])]
        row += 1
    return row

def determinant(m):
    if not m:
        return 1
    if len(m) == 1:
        return m[0][0]
    return sum((-1)**j * x * determinant([r[:j] + r[j+1:] for r in m[1:]]) for j, x in enumerate(m[0]) if x)

def probe(m):
    vertices = list(BOUNDARY)
    edges = []
    for i, (u, v) in enumerate(m):
        w = 'probe%d' % i
        vertices.append(w)
        edges += [(BOUNDARY[u], w), (w, BOUNDARY[v])]
    return graph(vertices, edges)

def glue(g, h):
    require(set(g[0]) & set(h[0]) == set(BOUNDARY), 'wrong overlap')
    require(not set(g[1]) & set(h[1]), 'edge overlap')
    return graph(g[0] + tuple(v for v in h[0] if v not in BOUNDARY), g[1] + h[1])

def count_hamiltonian_orders(g):
    """Independent vertex-order DFS, not forest pairing or edge-subset gluing."""
    vertices, edges = g
    adj, _ = components(vertices, edges)
    start = min(vertices, key=lambda v: (len(adj[v]), v))
    visited = {start}
    path = [start]
    oriented = []
    def walk(v):
        if len(visited) == len(vertices):
            if start in adj[v]:
                oriented.append(tuple(path))
            return
        for w in adj[v]:
            if w not in visited:
                visited.add(w)
                path.append(w)
                walk(w)
                path.pop()
                visited.remove(w)
    walk(start)
    require(len(oriented) % 2 == 0, 'orientation imbalance')
    unique = sorted(set(min(p, (p[0],) + tuple(reversed(p[1:]))) for p in oriented))
    require(2 * len(unique) == len(oriented), 'reversal quotient')
    return len(unique), unique

def nc_matchings(v):
    if not v:
        yield ()
        return
    for j in range(1, len(v), 2):
        for a in nc_matchings(v[1:j]):
            for b in nc_matchings(v[j+1:]):
                yield tuple(sorted(((v[0], v[j]),) + a + b))

def arbitrary_matchings(v):
    if not v:
        yield ()
        return
    for j in range(1, len(v)):
        for rest in arbitrary_matchings(v[1:j] + v[j+1:]):
            yield ((v[0], v[j]),) + rest

def all_live_states(k):
    for degree in itertools.product(range(3), repeat=k):
        ends = tuple(i for i, x in enumerate(degree) if x == 1)
        if not ends or len(ends) % 2:
            continue
        for m in nc_matchings(ends):
            yield degree, m

def external_samples():
    # All are planar: outer boundary hexagon, an inner rim, radial spokes,
    # and optional center spokes. Removing any edges preserves that embedding.
    for index in range(40):
        vertices = BOUNDARY + tuple('x%d' % i for i in range(6))
        edges = [(BOUNDARY[i], 'x%d' % i) for i in range(6)]
        edges += [('x%d' % i, 'x%d' % ((i + 1) % 6)) for i in range(6)]
        if index & 1:
            vertices += ('z',)
            edges += [('z', 'x%d' % i) for i in range(6)]
        edges += [(BOUNDARY[i], BOUNDARY[(i + 1) % 6]) for i in range(6) if ((index * 37 + i * 11) >> (i % 4)) & 1]
        if index >= 4:
            edges = [e for j, e in enumerate(edges) if (j * 29 + index * 13) % 23 != 0]
        yield index, graph(vertices, edges)

def drawing_check(g):
    rim = [(2,0),(1,2),(-1,2),(-2,0),(-1,-2),(1,-2)]
    xy = {}
    for i, (x,y) in enumerate(rim):
        xy['v%d'%i] = (x,y)
        xy['a%d'%i] = (2*x,2*y)
        xy['r%d'%i] = (3*x,3*y)
    if 'c' in g[0]:
        xy['c'] = (0,0)
    def sub(x,y):
        return x[0]-y[0], x[1]-y[1]
    def cross(x,y):
        return x[0]*y[1]-x[1]*y[0]
    def dot(x,y):
        return x[0]*y[0]+x[1]*y[1]
    require(len(set(xy.values()))==len(g[0]), 'drawing duplicate coordinates')
    for v in g[0]:
        if v not in BOUNDARY:
            for i in range(6):
                p,q=xy[BOUNDARY[i]],xy[BOUNDARY[(i+1)%6]]
                require(cross(sub(q,p),sub(xy[v],p))>0, 'vertex outside disk')
    comparisons=0
    for e,f in itertools.combinations(g[1],2):
        common=set(e)&set(f)
        if common:
            z=next(iter(common))
            u=next(v for v in e if v!=z)
            v=next(v for v in f if v!=z)
            a,b=sub(xy[u],xy[z]),sub(xy[v],xy[z])
            require(cross(a,b)!=0 or dot(a,b)<0, 'incident edge overlap')
        else:
            p,q=map(xy.get,e)
            r,t=map(xy.get,f)
            a,b,w=sub(q,p),sub(t,r),sub(r,p)
            det=cross(a,b)
            if det:
                numer1,numer2=cross(w,b),cross(w,a)
                # Segment parameter tests, performed as exact integer inequalities.
                lo,hi=sorted((0,det))
                require(not (lo<=numer1<=hi and lo<=numer2<=hi), 'nonincident edge crossing')
            elif cross(w,a)==0:
                aa=dot(a,a)
                vals=sorted((dot(w,a),dot(sub(t,p),a)))
                require(vals[1]<0 or vals[0]>aa, 'collinear nonincident overlap')
        comparisons+=1
    return {'coordinates':xy,'edge_pairs_checked':comparisons,'status':'PASS'}

def walsh_check():
    W=[[1,1],[1,-3]]
    require(determinant(W)==-4, 'Walsh determinant')
    values=[]
    for x,y,z in itertools.product((-1,1),repeat=3):
        values.append(x*y-x*z+y*z)
    require(values.count(1)==6 and values.count(-3)==2, 'Walsh tournament values')
    seen={}
    for a in itertools.product(range(8),repeat=2):
        response=tuple(x%8 for x in matrix_times(W,a))
        parity=tuple(x%2 for x in a)
        require(response not in seen or seen[response]==parity,'Walsh mod-eight recovery')
        seen[response]=parity
    require(tuple(x%4 for x in matrix_times(W,(1,-1)))==(0,0), 'Walsh mod-four hostile')
    return {'matrix':W,'determinant':-4,'smith_diagonal':[1,4],
            'labelled_values':values,'mod_eight_vectors_exhausted':64,
            'mod_four_collision_coefficients':[[0,0],[1,-1]]}

def main():
    A, B = patch(True), patch(False)
    computed = {}
    for name, g, wanted in [('A', A, [3,1,1,3,1]), ('B', B, [1,0,0,1,0])]:
        p, subsets, rejects = profile(g)
        require(all(d == (1,)*6 for d, m in p), 'unexpected degree sector')
        require(all(m in PAIRS for d, m in p), 'crossing matching')
        vector = [p.get(((1,)*6, m), 0) for m in PAIRS]
        require(vector == wanted, 'wrong gadget counts')
        require(rejects == 0, 'cycle-containing degree-valid subset')
        computed[name] = {'vertices': g[0], 'edges': g[1], 'profile': vector, 'subsets': subsets, 'rejected_cycles': rejects}
    difference = [a-b for a,b in zip(computed['A']['profile'], computed['B']['profile'])]
    require(matrix_times(K, difference) == [2]*5, 'integer identity')
    require(determinant(K) == -2 and parity_rank(K) == 4, 'matrix invariants')
    # A unit determinant 4x4 minor plus determinant -2 proves the Smith form.
    unit_minor = None
    for rows in itertools.combinations(range(5), 4):
        for cols in itertools.combinations(range(5), 4):
            sub = [[K[i][j] for j in cols] for i in rows]
            if abs(determinant(sub)) == 1:
                unit_minor = {'rows': rows, 'cols': cols, 'determinant': determinant(sub)}
                break
        if unit_minor:
            break
    require(unit_minor is not None, 'no unit 4-minor')
    probes = []
    for i, m in enumerate(PAIRS):
        p, _, _ = profile(probe(m))
        require(p == {((1,)*6,m):1}, 'probe not a unit')
        ac, aw = count_hamiltonian_orders(glue(A, probe(m)))
        bc, bw = count_hamiltonian_orders(glue(B, probe(m)))
        require((ac,bc) == [(3,1),(2,0),(2,0),(3,1),(2,0)][i], 'direct cycle counts')
        probes.append({'matching':m,'A_count':ac,'B_count':bc,'A_cycles':aw,'B_cycles':bw})
    states = list(all_live_states(6))
    require(len(states) == 365, 'live states')
    full = [[int(all(x+y==2 for x,y in zip(d,e))) * connectivity(m,n) if tuple(i for i,x in enumerate(d) if x==1) == tuple(i for i,x in enumerate(e) if x==1) else 0 for e,n in states] for d,m in states]
    degrees = [sum(row) for row in full]
    require(degrees.count(1)==362 and degrees.count(2)==3, 'full block degrees')
    require(parity_rank(full)==364, 'full mod-two rank')
    noncentral = [i for i,(d,m) in enumerate(states) if d!=(1,)*6]
    require(all(sum(full[i][j] for j in noncentral)==1 for i in noncentral), 'remaining block not permutation')
    # Exhaust every residue vector modulo four: responses recover all parity bits.
    seen = {}
    for f in itertools.product(range(4), repeat=5):
        r = tuple(x%4 for x in matrix_times(K,f))
        parity = tuple(x%2 for x in f)
        require(r not in seen or seen[r] == parity, 'mod-four parity ambiguity')
        seen[r] = parity
        recovery = (r[3]%2, ((r[2]+r[4]-r[1])%4)//2,
                    ((r[1]+r[4]-r[2])%4)//2, r[0]%2,
                    ((r[1]+r[2]-r[4])%4)//2)
        require(recovery == parity, 'mod-four recovery formula')
    unrestricted = []
    for m in arbitrary_matchings(tuple(range(6))):
        response = sum(connectivity(m, n)*x for n, x in zip(PAIRS,difference))
        require(response in (2,4,6), 'arbitrary exterior parity response')
        ac, _ = count_hamiltonian_orders(glue(A,probe(m)))
        bc, _ = count_hamiltonian_orders(glue(B,probe(m)))
        require(ac-bc==response, 'arbitrary physical probe')
        unrestricted.append({'matching':m,'difference':response,'A_count':ac,'B_count':bc,'noncrossing':m in PAIRS})
    require(len(unrestricted)==15, 'arbitrary matching count')
    samples = []
    for i, x in external_samples():
        xp, _, _ = profile(x)
        mass = sum(c for (d,m),c in xp.items() if d==(1,)*6)
        ac, _ = count_hamiltonian_orders(glue(A,x))
        bc, _ = count_hamiltonian_orders(glue(B,x))
        require(ac-bc==2*mass, 'external identity sample')
        samples.append({'index':i,'vertices':len(x[0]),'edges':len(x[1]),'all_one_mass':mass,'A_count':ac,'B_count':bc})
    return {'status':'PASS','method':'Independent exhaustive degree-pruned edge subsets and independent Hamiltonian vertex-order DFS.',
            'gadget_results':computed,'drawings':{'A':drawing_check(A),'B':drawing_check(B)},'order_three_walsh':walsh_check(),'K':K,'integer_difference':difference,'K_times_difference':matrix_times(K,difference),
            'determinant':determinant(K),'rank_mod_two':parity_rank(K),'unit_four_minor':unit_minor,
            'smith_diagonal':[1,1,1,1,2],'full_live_states':len(states),'full_rank_mod_two':parity_rank(full),
            'full_smith_diagonal_multiplicities':{'1':364,'2':1},'probes':probes,
            'mod_four_vectors_exhausted':4**5,'mod_four_response_count':len(seen),'additional_planar_exteriors':samples,'arbitrary_matching_responses':unrestricted}

if __name__=='__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
