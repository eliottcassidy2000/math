#!/usr/bin/env python3
"""Exact one-bit compiler; full-record graphs and native controls, not a census."""
from collections import deque
from itertools import combinations, combinations_with_replacement, permutations, product
from pathlib import Path
import argparse
import hashlib
import json
import random

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
CHECKS = 0

def need(value, label):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(label)

def shortest(adj, a, b):
    seen = {a: None}
    queue = deque([a])
    while queue:
        x = queue.popleft()
        if x == b:
            out = [x]
            while seen[x] is not None:
                x = seen[x]
                out.append(x)
            return out[::-1]
        for y in sorted(adj[x]):
            if y not in seen:
                seen[y] = x
                queue.append(y)
    raise ValueError('endpoint records disconnected')

def compile_route(paths, qmaps):
    """qmaps retain the full record->resource map; no kernel quotient of states."""
    s = len(paths)
    records = [p[0] for p in paths]
    layers = [0] * s
    done = set()
    events = []
    def vertical(i, kind):
        layers[i] ^= 1
        events.append((kind, i, records[i], layers[i]))
    while len(done) < s:
        i = next(j for j in range(s) if j not in done)
        while True:
            need(records[i] == paths[i][0], 'new mover retains its full start record')
            for y in paths[i][1:]:
                oldq, newq = qmaps[i][records[i]], qmaps[i][y]
                passive = [j for j in range(s) if j != i]
                need(len({qmaps[j][records[j]] for j in passive}) == len(passive),
                     'passive base resource injection')
                if oldq != newq:
                    occupant = next((j for j in passive if qmaps[j][records[j]] == newq), None)
                    if occupant is not None and layers[occupant] == layers[i]:
                        vertical(occupant, 'passive')
                records[i] = y
                events.append(('horizontal', i, y, layers[i]))
            done.add(i)
            occupant = next((j for j in range(s) if j != i and
                             qmaps[j][records[j]] == qmaps[i][records[i]]), None)
            if occupant is None:
                break
            need(occupant not in done, 'injective targets exclude completed handoff')
            need(layers[occupant] != layers[i], 'handoff shares base in opposite layers')
            i = occupant
    for i in range(s):
        if layers[i]:
            vertical(i, 'reset')
    return events

def verify(paths, qmaps, adjs, events):
    """Independent literal replay: validates state changes, capacities and full labels."""
    records = [p[0] for p in paths]
    layers = [0] * len(paths)
    visit = [[p[0]] for p in paths]
    counts = dict(horizontal=0, passive=0, reset=0)
    need(len({qmaps[i][x] for i, x in enumerate(records)}) == len(paths), 'start injection')
    need(len({qmaps[i][p[-1]] for i, p in enumerate(paths)}) == len(paths), 'target injection')
    for kind, i, x, layer in events:
        oldx, oldlayer = records[i], layers[i]
        need(x in adjs[i] and layer in (0, 1), 'authenticated full record and layer')
        if kind == 'horizontal':
            need(x in adjs[i][oldx] and layer == oldlayer, 'literal local horizontal edge')
            visit[i].append(x)
        else:
            need(kind in ('passive', 'reset') and x == oldx and layer == 1-oldlayer,
                 'vertical changes only real workspace bit')
        records[i], layers[i] = x, layer
        counts[kind] += 1
        need(len({(qmaps[j][y], layers[j]) for j, y in enumerate(records)}) == len(paths),
             'every intermediate normalized resource injects')
    need(visit == paths, 'every full-record route traversed exactly once and in order')
    need(records == [p[-1] for p in paths] and layers == [0]*len(paths), 'full endpoint and layer return')
    ell = sum(len(p)-1 for p in paths)
    need(counts['horizontal'] == ell and counts['passive'] <= ell and counts['reset'] <= len(paths),
         'exact horizontal and bounded vertical clocks')
    need(len(events) <= 2*ell+len(paths), 'event bound')
    return counts

def assignments(qmaps):
    out = []
    def visit(cur, used):
        i = len(cur)
        if i == len(qmaps):
            out.append(tuple(cur)); return
        for x, q in sorted(qmaps[i].items()):
            if q not in used:
                visit(cur+[x], used | {q})
    visit([], set())
    return out

def graph(vertices, edges):
    out = {v: set() for v in vertices}
    for x, y in edges:
        out[x].add(y); out[y].add(x)
    return out

def connected(adj):
    if not adj: return False
    seen = {next(iter(adj))}
    while True:
        new = seen | {y for x in seen for y in adj[x]}
        if new == seen: return len(seen) == len(adj)
        seen = new

def exhaustive_graphs():
    # Every connected labelled graph on a nonempty subset of 3 resources,
    # independently selected per slot. Slot permutation is the only quotient.
    local = []
    for size in (1, 2, 3):
        for vertices in combinations(range(3), size):
            pairs = list(combinations(vertices, 2))
            for mask in range(1 << len(pairs)):
                adj = graph(vertices, [e for j, e in enumerate(pairs) if mask >> j & 1])
                if connected(adj): local.append(adj)
    cases = routes = 0
    for s in range(1, 4):
        for indices in combinations_with_replacement(range(len(local)), s):
            adjs = [local[i] for i in indices]
            qmaps = [{x: x for x in adj} for adj in adjs]
            ends = assignments(qmaps)
            if not ends: continue
            cases += 1
            for a, b in product(ends, repeat=2):
                paths = [shortest(adj, x, y) for adj, x, y in zip(adjs, a, b)]
                verify(paths, qmaps, adjs, compile_route(paths, qmaps))
                routes += 1
    return dict(local_connected_graphs=len(local), feasible_menu_multisets=cases,
                ordered_endpoint_routes=routes, universe='all nonempty connected local graphs on subsets of 3 labelled resources; 1..3 slots')

def seeded_kernels():
    rng = random.Random(2026102301)
    cases = routes = noninjective = same_resource_edges = 0
    for case in range(600):
        n = rng.randrange(2, 10)
        s = rng.randrange(1, n+1)
        adjs, qmaps = [], []
        starts, targets = rng.sample(range(n), s), rng.sample(range(n), s)
        chosen_a, chosen_b = [], []
        for i in range(s):
            # Every resource has two distinct kernel records. They remain
            # different vertices, while q is deliberately 2-to-1.
            vertices = list(range(2*n))
            order = rng.sample(vertices, len(vertices))
            edges = set(tuple(sorted(e)) for e in zip(order, order[1:]))
            edges.update((2*j, 2*j+1) for j in range(n))
            for u, v in combinations(vertices, 2):
                if rng.randrange(8) == 0: edges.add((u, v))
            adjs.append(graph(vertices, edges))
            qmaps.append({x: x//2 for x in vertices})
            chosen_a.append(2*starts[i]+rng.randrange(2))
            chosen_b.append(2*targets[i]+rng.randrange(2))
        paths = [shortest(adj, a, b) for adj, a, b in zip(adjs, chosen_a, chosen_b)]
        verify(paths, qmaps, adjs, compile_route(paths, qmaps))
        same_resource_edges += sum(qmaps[i][a] == qmaps[i][b] for i, p in enumerate(paths) for a, b in zip(p, p[1:]))
        noninjective += len(qmaps)
        cases += 1; routes += 1
    # Deliberately revisit a consumed base key in a non-simple resource walk.
    adjs = [graph(range(5), [(0,1),(1,2),(2,3),(3,4)]), graph(range(2),[(0,1)])]
    qmaps = [{0:0,1:1,2:0,3:2,4:3}, {0:2,1:1}]
    paths = [[0,1,2,3,4], [0,1]]
    verify(paths, qmaps, adjs, compile_route(paths, qmaps))
    return dict(seed=2026102301, random_cases=cases, noninjective_record_maps=noninjective,
                same_resource_kernel_edges=same_resource_edges, extra_repeated_resource_walk=1)

P = ((3,0,6,7,1,5,4,2), (3,1,6,7,0,5,4,2))
G = (4,0,2,3,1,5,6,7)

def source(r):
    pairs = list(combinations(range(8),2))
    core = [[0]*8 for _ in range(8)]
    for bit, (u,v) in enumerate(pairs):
        core[u][v] = (29558785 >> (27-bit)) & 1
        core[v][u] = 1-core[u][v]
    return tuple(tuple(0 if u == v else core[u%8][v%8] if u//8 == v//8 else int(u//8 < v//8)
                       for v in range(8*r)) for u in range(8*r))

def path(r, x):
    return tuple(8*j+u for j in range(r) for u in P[(x>>j)&1])

def partner(A, p):
    clock = {u:j for j,u in enumerate(p)}
    return tuple(tuple(0 if u==v else A[u][v] if abs(clock[u]-clock[v])==1 else A[v][u]
                       for v in range(len(A))) for u in range(len(A)))

def word(A,p):
    return tuple(A[p[i]][p[j]] for i in range(len(p)) for j in range(i+2,len(p)))

def native():
    rows=[]
    for r in range(1,6):
        A=source(r); n=len(A)
        pp=[path(r,x) for x in range(1<<r)]
        kk=[partner(A,p) for p in pp]
        need(all(all(A[u][v] for u,v in zip(p,p[1:])) for p in pp), 'native actual path legality')
        ww=[word(A,p) for p in pp]
        rw=[word(tuple(zip(*A)),p[::-1]) for p in pp]
        allwords=ww+rw+[tuple(1-b for b in w) for w in ww+rw]
        need(len(set(allwords)) == 4*(1<<r), 'full normalized K/converse reservoir injects')
        for x in range(1<<r):
            need(partner(tuple(zip(*A)),pp[x][::-1])==tuple(zip(*kk[x])), 'literal converse compatibility')
            for j in range(r):
                y=x^(1<<j)
                perm=G if not (x>>j)&1 else tuple(G.index(u) for u in range(8))
                gg=tuple(8*j+perm[u%8] if u//8==j else u for u in range(n))
                need(all(kk[x][u][v]==kk[y][gg[u]][gg[v]] for u in range(n) for v in range(n)),
                     'actual directed cube edge partner isomorphism')
                need([i for i in range(n) if pp[x][i]!=pp[y][i]]==[8*j+1,8*j+4], 'exact two owner clocks')
                need(gg[pp[x][0]]==pp[y][0] and gg[pp[x][-1]]==pp[y][-1], 'marked target endpoint transport')
        rows.append(dict(blocks=r,owners=n,paths=1<<r,reservoir_words=len(allwords)))
    # Fully occupied cube endpoints; exact generic routing is then replayed
    # on native paths at every emitted event, including all pins and clocks.
    controls=[]
    for r in range(1,5):
        m=1<<r
        adj={x:{x^(1<<j) for j in range(r)} for x in range(m)}
        qmaps=[{x:x for x in adj} for _ in range(m)]
        paths=[shortest(adj,x,m-1-x) for x in range(m)]
        events=compile_route(paths,qmaps)
        verify(paths,qmaps,[adj]*m,events)
        states=list(range(m)); layers=[0]*m
        for kind,i,x,layer in events:
            old=path(r+1,states[i]+m*layers[i]); new=path(r+1,x+m*layer)
            need(sum(a!=b for a,b in zip(old,new))==2,'native event costs two actual owner clocks')
            states[i],layers[i]=x,layer
            need(len(set(states[j]+m*layers[j] for j in range(m)))==m,'native exact tile injection')
        controls.append(dict(base_blocks=r,slots=m,events=len(events),horizontal=r*m,clock_work=2*len(events)))
    return dict(cube_reconstruction=rows,full_cube_routes=controls)

def kernel_source(r):
    A=source(r); p=len(A); n=p+3
    return tuple(tuple(0 if u==v else A[u][v] if u<p and v<p else
                       int(u<p) if (u<p)!=(v<p) else int((v-u)%3==1)
                       for v in range(n)) for u in range(n))

def kernel_path(r,b,k):
    return path(r,b)+tuple(8*r+(k+j)%3 for j in range(3))

def native_kernels():
    rows=[]
    for r in range(1,5):
        A=kernel_source(r); n=len(A)
        records=list(product(range(1<<r),range(3)))
        paths={z:kernel_path(r,*z) for z in records}
        partners={z:partner(A,p) for z,p in paths.items()}
        words={z:word(A,p) for z,p in paths.items()}
        need(len(set(words.values()))==1<<r,'native kernel fibre has exactly one normalized word per cube key')
        for b in range(1<<r):
            need(len({words[b,k] for k in range(3)})==1,'three full kernel records do not multiply capacity')
        family=set(words.values())
        converse={word(tuple(zip(*A)),p[::-1]) for p in paths.values()}
        reservoir=family|converse|{tuple(1-v for v in w) for w in family|converse}
        need(len(reservoir)==4*(1<<r),'kernel normalized K/converse reservoir has correct unit count')
        for (b,k),p in paths.items():
            need(all(A[u][v] for u,v in zip(p,p[1:])),'native kernel full Hamiltonian path')
            need(partner(tuple(zip(*A)),p[::-1])==tuple(zip(*partners[b,k])), 'kernel converse compatibility')
            # Actual source/kernel automorphism and its induced partner map.
            for delta in (1,2):
                kk=(k+delta)%3
                gg=tuple(u if u<8*r else 8*r+(u-8*r+delta)%3 for u in range(n))
                need(all(A[u][v]==A[gg[u]][gg[v]] for u in range(n) for v in range(n)), 'literal C3 source automorphism')
                need(tuple(gg[u] for u in p)==paths[b,kk],'full kernel owner lift transport')
                need(all(partners[b,k][u][v]==partners[b,kk][gg[u]][gg[v]] for u in range(n) for v in range(n)), 'native kernel target isomorphism')
                need(sum(x!=y for x,y in zip(p,paths[b,kk]))==3,'kernel event costs three owner clocks')
            for j in range(r):
                bb=b^(1<<j)
                perm=G if not (b>>j)&1 else tuple(G.index(u) for u in range(8))
                gg=tuple(8*j+perm[u%8] if u//8==j else u for u in range(n))
                need(all(partners[b,k][u][v]==partners[bb,k][gg[u]][gg[v]] for u in range(n) for v in range(n)), 'prefix cube target isomorphism retains kernel')
                need(sum(x!=y for x,y in zip(p,paths[bb,k]))==2,'prefix cube move costs two owner clocks')
            # Independently establish the partner strong connectivity used
            # to separate inverse tiles from decomposable source tiles.
            for start in range(n):
                seen={start}; todo=[start]
                while todo:
                    u=todo.pop()
                    for v in range(n):
                        if partners[b,k][u][v] and v not in seen:
                            seen.add(v);todo.append(v)
                need(len(seen)==n,'native kernel partner strong')
        rows.append(dict(prefix_blocks=r,owners=n,full_records=len(records),normalized_tiles=1<<r,
                         normalized_reservoir_words=len(reservoir),kernel_fibre_size=3))
    routes=[]
    for r in range(1,4):
        m=1<<r; adj={}
        for b in range(m):
            for k in range(3):
                adj[3*b+k]={3*(b^(1<<j))+k for j in range(r)}|{3*b+(k+1)%3,3*b+(k+2)%3}
        qmaps=[{z:z//3 for z in adj} for _ in range(m)]
        pp=[shortest(adj,3*i+i%3,3*(m-1-i)+(i+1)%3) for i in range(m)]
        events=compile_route(pp,qmaps);verify(pp,qmaps,[adj]*m,events)
        current=[p[0] for p in pp]; layers=[0]*m; clockwork=rotations=0
        for kind,i,z,layer in events:
            old=kernel_path(r+1,current[i]//3+m*layers[i],current[i]%3)
            new=kernel_path(r+1,z//3+m*layer,z%3)
            isrotation=kind=='horizontal' and current[i]//3==z//3
            cost=sum(a!=b for a,b in zip(old,new))
            need(cost==(3 if isrotation else 2),'full native kernel route event clock cost')
            need(all(old[8*j+t]==new[8*j+t] for j in range(r+1) for t in (0,2,3,5,6,7)), 'all prefix shared pins retained')
            current[i],layers[i]=z,layer
            need(len({(current[j]//3,layers[j]) for j in range(m)})==m,'native q-layer capacity counts kernel only once')
            clockwork+=cost;rotations+=int(isrotation)
        need(clockwork==2*len(events)+rotations,'exact 2-clock plus kernel surcharge')
        routes.append(dict(base_prefix_blocks=r,source_owners=8*r+3,stabilized_owners=8*r+11,
                           slots=m,events=len(events),kernel_rotations=rotations,owner_clock_work=clockwork))
    return dict(exact_native_fibres=rows,full_record_routes=routes)

def shared_owner_hostile():
    # All states at workspace dimensions0..3; legal one-slot cube edges.
    rows=[]
    for w in range(4):
        menus=([x for x in range(1<<(2+w)) if not x&2],
               [x for x in range(1<<(2+w)) if x&2])
        states={x for x in product(*menus) if (x[0]&1)==(x[1]&1)}
        need(len(states)==2*(1<<(2*w)),'complete shared-owner states retain all workspace choices')
        for state in states:
            need(path(2+w,state[0])[1]==path(2+w,state[1])[1],'literal dynamic owner equality')
        seen={(0,2)}; todo=deque(seen)
        while todo:
            state=todo.popleft()
            for i in range(2):
                for j in range(2+w):
                    y=list(state);y[i]^=1<<j;y=tuple(y)
                    if y in states and y not in seen:seen.add(y);todo.append(y)
        need((1,3) not in seen and len(seen)==1<<(2*w),'arbitrarily supplied tested neutral bits cannot change base projection')
        need(sum(a!=b for i in range(2) for a,b in zip(path(2+w,(0,2)[i]),path(2+w,(1,3)[i])))==4,
             'simultaneous two-slot repair costs exactly four owner clocks')
        rows.append(dict(workspace_bits=w,states=len(states),one_slot_component=len(seen)))
    # Base Gray cycle cannot be repaired by calling two kernel decorations
    # different normalized resources: two slots forced to q0 still collide.
    need(len(assignments([{0:0,1:0},{0:0,1:0}]))==0,'kernel multiplicity creates no resource capacity')
    return rows

def main():
    result=dict(status='FINITE-EXACT PASS; all-size theorem proved in ONE_BIT_WORKSPACE.md',
                exhaustive=exhaustive_graphs(),kernels=seeded_kernels(),native=native(),native_kernels=native_kernels(),
                shared_owner_hostile=shared_owner_hostile())
    pins=['guarded_coverage_20261023/PORTS.md','guarded_coverage_20261023/CUBES.md',
          'guarded_coverage_20261023/SYNTHESIS.md',
          'session45_injective_multisource_moves_20261023/incoming-integration/WORKSPACE_BOUNDARY.md',
          'session45_injective_multisource_moves_20261023/injective-completion/THEOREM.md']
    result['source_sha256']={p:hashlib.sha256((REPO/'05-knowledge/results'/p).read_text(encoding='utf-8').replace('\r\n','\n').encode('utf-8')).hexdigest() for p in pins}
    result['hash_convention']='UTF-8 text with CRLF normalized to LF'
    result['checks']=CHECKS
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true',help='regenerate the frozen receipt after a deliberate change')
    args=parser.parse_args()
    receipt=HERE/'receipt.json'
    if args.write:
        receipt.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    else:
        expected=json.loads(receipt.read_text(encoding='utf-8'))
        if result!=expected:
            raise RuntimeError('frozen receipt differs; inspect change before --write')
    print(f'PASS {CHECKS:,}; '+('wrote frozen receipt' if args.write else 'frozen receipt matches'))

if __name__=='__main__':main()
