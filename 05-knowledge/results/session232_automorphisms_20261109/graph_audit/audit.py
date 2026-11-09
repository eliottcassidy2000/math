#!/usr/bin/env python3
"""Independent session232 audit. Standard library only; no producer imports."""
import ast
from collections import Counter, defaultdict, deque
from itertools import permutations, product
import hashlib
import json
from pathlib import Path

SOURCE = Path(__file__).resolve().parent.parent / 'boundary_graphs'
if not SOURCE.exists():
    SOURCE = Path(__file__).resolve().parent.parent / 'session232-action-graphs'

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

def matrix(n, mask):
    edges = iter((mask >> k) & 1 for k in range(n*(n-1)//2))
    a = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            a[i][j] = bool(next(edges))
            a[j][i] = not a[i][j]
    return a

def cayley(n, steps):
    return [[(j-i) % n in steps for j in range(n)] for i in range(n)]

def stack(a, b):
    n = len(a)
    return [row + [True]*len(b) for row in a] + [[False]*n + row for row in b]

def brute(a):
    return [p for p in permutations(range(len(a))) if all(a[u][v] for u,v in zip(p, p[1:]))]

def dfs(a):
    n = len(a)
    out = []
    def grow(p, unseen):
        if not unseen:
            out.append(p)
        else:
            for v in range(n):
                if unseen & (1 << v) and (not p or a[p[-1]][v]):
                    grow(p+(v,), unseen ^ (1 << v))
    grow((), (1 << n)-1)
    return out

def autos(a):
    n = len(a)
    if n <= 7:
        candidates = permutations(range(n))
    else:
        blocks = defaultdict(list)
        for v in range(n):
            blocks[sum(a[v])].append(v)
        blocks = list(blocks.values())
        def candidate_iter():
            for images in product(*(permutations(b) for b in blocks)):
                p = list(range(n))
                for block,image in zip(blocks,images):
                    for v,w in zip(block,image):
                        p[v] = w
                yield tuple(p)
        candidates = candidate_iter()
    return [p for p in candidates if all(a[u][v] == a[p[u]][p[v]] for u in range(n) for v in range(n))]

def encode(a):
    # All five gadget vertices are created together; unlike producer ordering.
    g = [set() for _ in a]
    for u in range(len(a)):
        for v in range(len(a)):
            if not a[u][v]:
                continue
            l,r,x,y,z = range(len(g),len(g)+5)
            g.extend(set() for _ in range(5))
            for s,t in ((u,l),(l,r),(r,v),(l,x),(r,y),(y,z)):
                g[s].add(t); g[t].add(s)
    return g

def relabel(g):
    # Reverse IDs. Decoder and generic automorphism test receive only adjacency.
    n = len(g)
    return [{n-1-v for v in g[n-1-u]} for u in range(n)]

def decode(g):
    degree = [len(adj) for adj in g]
    core = set(range(len(g)))
    queue = deque(v for v in core if degree[v] < 2)
    while queue:
        v = queue.popleft()
        if v not in core:
            continue
        core.remove(v)
        for w in g[v] & core:
            degree[w] -= 1
            if degree[w] < 2:
                queue.append(w)
    roles = {}
    for v in core:
        branches = g[v] - core
        require(len(branches) <= 1, 'at most one pendant branch')
        chain,previous = list(branches),v
        if chain:
            current = chain[0]
            while True:
                nxt = g[current] - {previous}
                require(len(nxt) <= 1 and not nxt & core, 'pendant chain has no branch or extra core attachment')
                if not nxt:
                    break
                previous,current = current,next(iter(nxt))
                chain.append(current)
        require(len(chain) in (0,1,2), 'valid pendant lengths')
        roles[v] = len(chain)
    originals = sorted(v for v in core if roles[v] == 0)
    index = {v:i for i,v in enumerate(originals)}
    recovered = [[False]*len(originals) for _ in originals]
    for l in core:
        if roles[l] != 1:
            continue
        ns = g[l] & core
        require(len(ns) == 2, 'L core degree')
        us = [v for v in ns if roles[v] == 0]
        rs = [v for v in ns if roles[v] == 2]
        require(len(us) == len(rs) == 1, 'L roles')
        u,r = us[0],rs[0]
        rn = g[r] & core
        require(len(rn) == 2 and l in rn, 'R core degree')
        v = next(iter(rn-{l}))
        require(roles[v] == 0, 'R original neighbor')
        recovered[index[u]][index[v]] = True
    return recovered, originals, roles

def generic_automorphisms(g):
    """Paired individualization/refinement on uncoloured adjacency only."""
    n = len(g)
    degree = tuple(map(len,g))
    def recurse(left,right):
        while True:
            sl = [(left[v],tuple(sorted(left[w] for w in g[v]))) for v in range(n)]
            sr = [(right[v],tuple(sorted(right[w] for w in g[v]))) for v in range(n)]
            labels = {s:k for k,s in enumerate(sorted(set(sl+sr)))}
            nl = tuple(labels[s] for s in sl)
            nr = tuple(labels[s] for s in sr)
            if Counter(nl) != Counter(nr):
                return
            stable = len(set(nl)) == len(set(left))
            left,right = nl,nr
            if stable:
                break
        lc,rc = defaultdict(list),defaultdict(list)
        for v in range(n):
            lc[left[v]].append(v); rc[right[v]].append(v)
        mult = [c for c in lc if len(lc[c]) > 1]
        if not mult:
            p = tuple(rc[left[v]][0] for v in range(n))
            require(all({p[w] for w in g[v]} == g[p[v]] for v in range(n)), 'terminal adjacency map')
            yield p
            return
        c = min(mult,key=lambda t:(len(lc[t]),t))
        v = lc[c][0]
        fresh = max(left)+1
        ll = list(left); ll[v] = fresh
        for w in rc[c]:
            rr = list(right); rr[w] = fresh
            yield from recurse(tuple(ll),tuple(rr))
    return list(recurse(degree,degree))

def cycle_lengths(p):
    todo = set(range(len(p)))
    lengths = []
    while todo:
        v = min(todo); w = v; length = 0
        while w in todo:
            todo.remove(w); length += 1; w = p[w]
        require(w == v, 'permutation cycle closure')
        lengths.append(length)
    return lengths

def check(a, expected=None, generic=True):
    n = len(a)
    paths = brute(a) if n <= 7 else dfs(a)
    ends = [sum(p[-1] == v for p in paths) for v in range(n)]
    if expected is not None:
        require(ends == expected, 'endpoint counts')
    group = autos(a)
    vtodo = set(range(n)); boundary = []
    while vtodo:
        v = min(vtodo)
        vorbit = {p[v] for p in group}; vtodo -= vorbit
        st = [p for p in group if p[v] == v]
        endpoint_orbits = {min(tuple(p[x] for x in path) for p in st) for path in paths if path[-1] == v}
        whole_over_vorbit = {min(tuple(p[x] for x in path) for p in group) for path in paths if path[-1] in vorbit}
        require(len(endpoint_orbits) == len(whole_over_vorbit), 'endpoint fibre bijection count')
        require(len(endpoint_orbits)*len(st) == ends[v], 'free endpoint action')
        boundary.append((v,len(st),ends[v],len(endpoint_orbits)))
    require(sum(b[3] for b in boundary)*len(group) == len(paths), 'global decomposition')
    g = relabel(encode(a))
    recovered, originals, roles = decode(g)
    require(recovered == [list(reversed(row)) for row in reversed(a)], 'decode reverse-labelled tournament')
    require(Counter(roles.values()) == {0:n,1:n*(n-1)//2,2:n*(n-1)//2}, 'intrinsic roles')
    require(len(g) == n+5*n*(n-1)//2, 'graph order')
    require(sum(map(len,g)) == 6*n*(n-1), 'twice graph edges')
    require(sum(len(s)==1 for s in g) == n*(n-1), 'graph leaves')
    result = {'n':n,'H':len(paths),'aut':len(group),'endpoint':ends,'boundary':boundary,'graph_order':len(g)}
    if generic:
        ga = generic_automorphisms(g)
        require(len(ga) == len(group), 'independent graph aut count')
        original_set = set(originals)
        require(all({p[v] for v in originals} == original_set for p in ga), 'graph aut preserves original role')
        restrictions = {tuple(originals.index(p[v]) for v in originals) for p in ga}
        require(len(restrictions) == len(ga), 'trivial graph restriction kernel')
        require(all(all(k % 2 for k in cycle_lengths(p)) for p in ga), 'all graph automorphism cycles odd')
        result['graph_aut'] = len(ga)
        result['graph_restriction_kernel'] = 1
        result['graph_all_signs_even'] = True
    return result

def explicit_action(a,root):
    paths = dfs(a); ids = {p:i for i,p in enumerate(paths)}
    # Generator rotates only the initial C3.
    generator = (1,2,0)+tuple(range(3,len(a)))
    adjacency = [set() for _ in paths]
    for p,i in ids.items():
        j = ids[tuple(generator[v] for v in p)]
        adjacency[i].add(j); adjacency[j].add(i)
    def components(subset):
        sizes = Counter(); todo = set(subset)
        while todo:
            component = set(); pending = [min(todo)]
            while pending:
                u = pending.pop()
                if u in component:
                    continue
                component.add(u); pending.extend(adjacency[u]-component)
            require(component <= todo, 'whole components stay in restricted action graph')
            todo -= component; sizes[len(component)] += 1
        return dict(sizes)
    return {'whole':components(range(len(paths))), 'endpoint':components(i for i,p in enumerate(paths) if p[-1]==root)}

def main():
    require(not any(isinstance(node,ast.Assert) for node in ast.walk(ast.parse((SOURCE/'verify.py').read_text()))), 'producer has no optimization-removable asserts')
    maxima = []
    for n in range(1,6):
        maxima.append(max(len(brute(matrix(n,m))) for m in range(1 << (n*(n-1)//2))))
    require(maxima == [1,1,3,5,15], 'independent exhaustive maxima')
    b7,b21,c3 = matrix(6,26),matrix(7,290),cayley(3,{1})
    output = {'exhaustive_brute_maxima':maxima,'witness7':check(b7,[1,5,7,3,2,1]),'witness21':check(b21,[3,9,21,9,7,3,3])}
    for name,b in [('7',b7),('21',b21)]:
        a = stack(c3,b)
        output['tower'+name] = check(a)
        output['action'+name] = explicit_action(a,5)
        require(output['action'+name] == {'whole':{3:19 if name=='7' else 55},'endpoint':{3:int(name)}}, 'explicit triangle action graph')
    output['R7'] = check(cayley(7,{1,2,3}),[25]*7)
    output['P7'] = check(cayley(7,{1,2,4}),[27]*7)
    total = 0
    for n in (3,4):
        for mask in range(1 << (n*(n-1)//2)):
            check(matrix(n,mask)); total += 1
    output['complete_small_graph_checks'] = total
    output['producer_sha256'] = {name:hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() for name in ('PROOF.md','verify.py','results.json')}
    print(json.dumps(output,indent=2,sort_keys=True))

if __name__ == '__main__':
    main()
