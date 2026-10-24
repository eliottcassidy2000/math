#!/usr/bin/env python3
"""Independent actual-record audit. Standard library; no producer imports or asserts."""
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from math import comb, factorial
from pathlib import Path
import argparse, hashlib, json

HERE = Path(__file__).resolve().parent
CHECKS = Counter()

def require(ok, tag):
    CHECKS[tag] += 1
    if not ok:
        raise RuntimeError('Independent audit failed: '+tag)

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def validate_inputs():
    kernel = HERE/'inputs/kernel'
    m = json.loads((kernel/'MANIFEST.json').read_text())
    actual = {str(p.relative_to(kernel)) for p in kernel.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name != 'MANIFEST.json'}
    require(actual == set(m['files']), 'kernel manifest closure')
    for name, value in m['files'].items():
        require(hashlib.sha256((kernel/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest() == value, 'kernel input digest')
    completion = HERE/'inputs/completion'
    m = json.loads((completion/'MANIFEST.json').read_text())
    actual = {str(p.relative_to(completion)) for p in completion.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name != 'MANIFEST.json'}
    require(actual == {x['path'] for x in m['files']}, 'completion manifest closure')
    for row in m['files']:
        p = completion/row['path']
        require(p.stat().st_size == row['bytes'] and digest(p) == row['sha256'], 'completion input digest')
    return {name: digest(HERE/'inputs'/name/'MANIFEST.json') for name in ('kernel','completion')}

P0 = (3,0,6,7,1,5,4,2)
P1 = (3,1,6,7,0,5,4,2)
G = (4,0,2,3,1,5,6,7)

def source(r):
    # Read the 28-bit seed in upper-triangle lexicographic order.
    bits = format(29558785, '028b')
    seed = [[0]*8 for _ in range(8)]
    for (a,b), bit in zip(combinations(range(8),2), bits):
        seed[a][b] = int(bit); seed[b][a] = 1-int(bit)
    n = 8*r+3
    S = [[0]*n for _ in range(n)]
    for a,b in combinations(range(n),2):
        if a//8 == b//8 and b < 8*r:
            value = seed[a%8][b%8]
        elif a >= 8*r:
            value = int((b-a)%3 == 1)
        else:
            value = 1
        S[a][b] = value; S[b][a] = 1-value
    return tuple(map(tuple,S))

def record(r, key, phase):
    out = []
    for j in range(r):
        out.extend(8*j+v for v in (P1 if key & (1<<j) else P0))
    out.extend(8*r+(phase+j)%3 for j in range(3))
    return tuple(out)

def K(S, P):
    E = {frozenset((u,v)) for u,v in zip(P,P[1:])}
    return tuple(tuple(0 if u==v else S[u][v] if frozenset((u,v)) in E else S[v][u] for v in range(len(S))) for u in range(len(S)))

def w(S,P):
    # Full path-coordinate matrix, not the producer's deleted-diagonal word.
    return tuple(S[P[i]][P[j]] for i,j in combinations(range(len(P)),2))

def transpose(S):
    return tuple(tuple(S[j][i] for j in range(len(S))) for i in range(len(S)))

def tile(S,P):
    T = K(S,P)
    return frozenset((w(S,P), w(T,P), w(transpose(S),P[::-1]), w(transpose(T),P[::-1])))

def hamming(x,y):
    return sum(a!=b for a,b in zip(x,y))

def connected_sizes(states, budget, cost=hamming):
    rest = set(states); sizes=[]
    while rest:
        component = {rest.pop()}; frontier=list(component)
        while frontier:
            x=frontier.pop()
            new={y for y in rest if cost(x,y)<=budget}
            rest.difference_update(new); component.update(new); frontier.extend(new)
        sizes.append(len(component))
    return sorted(sizes)

def native():
    rows=[]
    seed=tuple(row[:8] for row in source(1)[:8])
    degree=[sum(row) for row in seed]
    automorphisms=[p for p in permutations(range(8)) if all(degree[i]==degree[p[i]] for i in range(8)) and all(seed[i][j]==seed[p[i]][p[j]] for i in range(8) for j in range(8))]
    require(automorphisms==[tuple(range(8))], 'seed matrix has trivial automorphism group')
    for graph in (seed,K(seed,P0)):
        for root in range(8):
            found={root};todo=[root]
            while todo:
                v=todo.pop()
                for u in range(8):
                    if graph[v][u] and u not in found:found.add(u);todo.append(u)
            require(len(found)==8, 'seed source and seed partner strong connectivity')
    for r in range(1,6):
        S=source(r); n=len(S); keys=range(1<<r)
        paths={(b,k):record(r,b,k) for b in keys for k in range(3)}
        targets={key:K(S,P) for key,P in paths.items()}
        tiles={key:tile(S,P) for key,P in paths.items()}
        for (b,k),P in paths.items():
            require(set(P)==set(range(n)) and all(S[u][v] for u,v in zip(P,P[1:])), 'native Hamiltonian record')
            require(len(tiles[b,k])==4, 'four genuinely distinct F/R orientations')
            require(tiles[b,k]==tiles[b,0], 'phase lifts are one normalized resource')
            require(K(transpose(S),P[::-1])==transpose(targets[b,k]), 'literal F/R partner square')
            require(all(P[::-1][n-1-t]==P[t] for t in range(n)), 'all converse clocks reflect')
            if r<=4:
                tails={P[:-3]+t for t in permutations(range(8*r,n)) if all(S[u][v] for u,v in zip(P[:-3]+t,(P[:-3]+t)[1:]))}
                require(tails=={paths[b,j] for j in range(3)}, 'entire prefix-pinned menu')
            for delta in (1,2):
                perm=tuple(u if u<8*r else 8*r+(u-8*r+delta)%3 for u in range(n))
                Q=paths[b,(k+delta)%3]
                require(tuple(perm[u] for u in P)==Q and hamming(P,Q)==3, 'literal kernel transport and work')
                require(all(S[u][v]==S[perm[u]][perm[v]] for u in range(n) for v in range(n)), 'matrix kernel automorphism')
                require(all(targets[b,k][u][v]==targets[b,(k+delta)%3][perm[u]][perm[v]] for u in range(n) for v in range(n)), 'matrix kernel target transport')
            for j in range(r):
                gg=G if b&(1<<j)==0 else tuple(G.index(i) for i in range(8))
                perm=tuple(8*j+gg[u%8] if 8*j<=u<8*j+8 else u for u in range(n))
                require(all(targets[b,k][u][v]==targets[b^(1<<j),k][perm[u]][perm[v]] for u in range(n) for v in range(n)), 'endpoint-marked target transport')
        require(len(set(tiles.values()))==1<<r, 'distinct keys are distinct normalized units')
        require(len(set().union(*tiles.values()))==4*(1<<r), 'normalized reservoirs disjoint')
        # A capacity-one unit rejects two same-key slots even at different phases.
        require(tiles[0,0]==tiles[0,1] and tiles[0,0].isdisjoint(tiles[1,0]), 'capacity-one duplicate key hostile')
        rows.append({'r':r,'order':n,'records':len(paths),'unit_tiles':1<<r,'F_R_words':4*(1<<r)})
    S=source(1); A=record(1,0,0); B=record(1,1,0)
    KA,KB=K(S,A),K(S,B)
    require(S[3][0]==1 and S[3][1]==1, 'common-source fork orientations')
    require(KA[3][0]==1 and KB[3][0]==0 and KA!=KB, 'same target class not same named partner')
    require((3,0) in set(zip(A,A[1:])) and (3,1) in set(zip(B,B[1:])), 'two full observers require conflicting successors')
    return {'matrices':rows,'no_common_global_K_control':{'source_order':11,'keys':[0,1],'conflicting_required_successors':[[3,0],[3,1]],'same_named_partner':False,'private_extensions_cannot_repair':True}}

# Modular Gaussian elimination is intentionally different from the producer's graph propagation.
def affine_rank(n, edges, pins):
    rows=[]
    for i,j,d in edges:
        row=[0]*(n+1); row[i]=(row[i]-1)%3;row[j]=(row[j]+1)%3;row[-1]=d%3;rows.append(row)
    for i,k in pins:
        row=[0]*(n+1);row[i]=1;row[-1]=k%3;rows.append(row)
    rank=0
    for col in range(n):
        pivot=next((r for r in range(rank,len(rows)) if rows[r][col]),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        inv=1 if rows[rank][col]==1 else 2
        rows[rank]=[(v*inv)%3 for v in rows[rank]]
        for r in range(len(rows)):
            if r!=rank and rows[r][col]:
                c=rows[r][col];rows[r]=[(a-c*b)%3 for a,b in zip(rows[r],rows[rank])]
        rank+=1
    return None if any(not any(row[:n]) and row[-1] for row in rows) else rank

def relation_components(n, edges):
    parts=[{i} for i in range(n)]
    for a,b,_ in edges:
        left=next(x for x in parts if a in x);right=next(x for x in parts if b in x)
        if left is not right:parts.remove(left);parts.remove(right);parts.append(left|right)
    return sorted(parts,key=min)

def literal(n,edges,pins):
    r=max(1,(n-1).bit_length()); q=8*r
    paths=[[record(r,i,k) for k in range(3)] for i in range(n)]
    return tuple(ks for ks in product(range(3),repeat=n) if all(paths[i][ks[i]][q+d%3]==paths[j][ks[j]][q] for i,j,d in edges) and all(paths[i][ks[i]][q]==q+k%3 for i,k in pins))

def catalogue(n):
    pairs=list(combinations(range(n),2))
    for digits in product(range(4),repeat=len(pairs)):
        yield [(i,j,d-1) for (i,j),d in zip(pairs,digits) if d]

def voltage_and_packets():
    rows=[]; unique={}
    for n in range(1,5):
        local={};graphs=0;pin_cases=0
        for edges in catalogue(n):
            graphs+=1
            plist=[()] if n==4 else [tuple((i,d-1) for i,d in enumerate(v) if d) for v in product(range(4),repeat=n)]
            for pins in plist:
                pin_cases+=1;states=literal(n,edges,pins);rank=affine_rank(n,edges,pins)
                require(len(states)==(0 if rank is None else 3**(n-rank)), 'literal acceptance versus modular row rank')
                parts=relation_components(n,edges);fixed={i for i,_ in pins}
                free=[part for part in parts if not part&fixed]
                if states:
                    require(len(states)==3**len(free), 'one C3 torsor per unpinned component')
                    local.setdefault(states,(edges,pins,parts,free))
        # All indexed comparisons, including same-slot clock loops, are literal tests.
        r=max(1,(n-1).bit_length());q=8*r
        paths=[[record(r,i,k) for k in range(3)] for i in range(n)]
        for i,j in product(range(n),repeat=2):
            for t,u in product(range(q+3),repeat=2):
                actual={(a,b) for a,b in product(range(3),repeat=2) if (i!=j or a==b) and paths[i][a][t]==paths[j][b][u]}
                expected={(a,b) for a,b in product(range(3),repeat=2) if (i!=j or a==b) and ((paths[i][0][t]==paths[j][0][u]) if t<q and u<q else False if (t<q)!=(u<q) else (b-a-t+u)%3==0)}
                require(actual==expected, 'arbitrary indexed equality including self occurrence')
        for states,(edges,pins,parts,free) in local.items():
            base=states[0];pred=Counter({0:1})
            for part in free:
                nxt=Counter(pred)
                for d,count in pred.items():nxt[d+len(part)]+=2*count
                pred=nxt
            require(Counter(hamming(base,x) for x in states)==pred, 'joint slots and triple-clock enumerator')
            for budget in range(n+1):
                heavy=sum(len(part)>budget for part in free);light=len(free)-heavy
                require(connected_sizes(states,budget)==[3**light]*(3**heavy), 'exact heavy-component budget graph')
            # Deterministic direct routes never release a resource or revisit a departed record.
            for first,last in product(states,repeat=2):
                now=list(first);seen=[{k} for k in now]
                for part in free:
                    if now[min(part)]==last[min(part)]:continue
                    before=tuple(now)
                    for i in part:
                        require(last[i] not in seen[i], 'direct route no full-record return')
                        now[i]=last[i];seen[i].add(now[i])
                    require(tuple(now) in states and hamming(before,now)==len(part), 'direct route accepted indivisible step')
                require(tuple(now)==last, 'direct route exact endpoint')
        unique[n]=local
        rows.append({'slots':n,'simple_graphs':graphs,'graph_pin_instances':pin_cases,'nonempty_distinct_fibres':len(local)})
    # Loops/parallel constraints and conflicting pins are outside the simple graph catalogue.
    for edges,pins in [([(0,0,1)],()), ([(0,0,0)],()), ([(0,1,0),(0,1,1)],()), ([(0,1,2),(1,0,1)],()), ([],((0,0),(0,1)))]:
        states=literal(2,edges,pins);rank=affine_rank(2,edges,pins)
        require(len(states)==(0 if rank is None else 3**(2-rank)), 'loop parallel and duplicate pin controls')
    return rows,unique

def overlay_minimum(n,A,B,da,db,weights):
    _,pina,partsa,_=da;_,pinb,partsb,_=db
    abase=A[0];bbase=B[0]
    # Root phases relative to the actual observed base solutions, not propagated potentials.
    posa={i:c for c,part in enumerate(partsa) for i in part}
    posb={i:len(partsa)+c for c,part in enumerate(partsb) for i in part}
    edges=[(posa[i],posb[i],(abase[i]-bbase[i])%3) for i in range(n)]
    pins={(posa[i],0) for i,_ in pina}|{(posb[i],0) for i,_ in pinb}
    best=sum(weights)
    for keep in range(1<<n):
        removed=sum(weights[i] for i in range(n) if not (keep>>i)&1)
        if removed>best:continue
        if affine_rank(len(partsa)+len(partsb),[edges[i] for i in range(n) if (keep>>i)&1],pins) is not None:best=removed
    return best

def boundary_and_helper(unique):
    stats=[]
    for n,cat in unique.items():
        # Include all unpinned fibres through four slots and all pinned fibres through two;
        # three-slot pinned controls retain deterministic independent representatives.
        entries=list(cat.items())
        if n==3:entries=[(s,d) for s,d in entries if not d[1] or d[1] in ( ((0,0),), ((2,2),), ((0,0),(2,2)) )]
        total=positive=0
        for (A,da),(B,db) in product(entries,repeat=2):
            weights=tuple(1+i%3 for i in range(n))
            direct=min(sum(weights[i] for i in range(n) if a[i]!=b[i]) for a in A for b in B)
            deletion=overlay_minimum(n,A,B,da,db,weights)
            require(direct==deletion, 'weighted native boundary equals feasible edge deletion')
            require((direct==0)==bool(set(A)&set(B)), 'positive weights detect intersection')
            total+=1;positive+=direct>0
        stats.append({'slots':n,'fibres':len(entries),'ordered_pairs':total,'positive_boundaries':positive})
    S=source(3); paths={(i,k):record(3,i,k) for i in range(5) for k in range(3)}
    states=[]
    for h in (0,1):
        for a,b in product(range(3),repeat=2):states.append((a,(a+h)%3,b,(b+h)%3,h))
    require(len(states)==18, 'eighteen actual helper assemblies')
    for state in states:
        require(len({tile(S,paths[i,k]) for i,k in enumerate(state)})==5, 'helper five capacity-one units')
    require(connected_sizes(states,2)==[9,9] and connected_sizes(states,3)==[18], 'helper exact sharp budget three')
    clock_cost=lambda a,b:sum(hamming(paths[i,a[i]],paths[i,b[i]]) for i in range(5))
    require(connected_sizes(states,8,clock_cost)==[9,9] and connected_sizes(states,9,clock_cost)==[18], 'helper exact sharp nine owner clocks')
    for h in (0,1):
        mode=[x for x in states if x[-1]==h]
        require(connected_sizes(mode,2)==[9], 'each helper mode budget two connected')
        for i in range(4):require({x[i] for x in mode}==set(range(3)), 'helper modes all unary phases survive')
    left=[x for x in states if x[-1]==0];right=[x for x in states if x[-1]==1]
    require(min(hamming(a,b) for a in left for b in right)==3, 'exact old boundary plus helper cost')
    # Full pairwise projections of the five-variable relation cannot define the relation.
    pairprojections={(i,j):{(x[i],x[j]) for x in states} for i,j in combinations(range(5),2)}
    binary_closure=[x for x in product(range(3),range(3),range(3),range(3),range(2)) if all((x[i],x[j]) in values for (i,j),values in pairprojections.items())]
    require(set(states)<set(binary_closure), 'higher-arity condition not binary pairwise closure')
    # Any sequence of these fixed-key states has constant occupied tile set; zero releases.
    held=[tuple(tile(S,paths[i,k]) for i,k in enumerate(x)) for x in states]
    require(all(x==held[0] for x in held), 'continuous normalized occupancy no reacquisition')
    return {'overlay':stats,'helper':{'accepted_states':18,'pairwise_relaxation_states':len(binary_closure),'slot_budget_components':{'2':[9,9],'3':[18]},'clock_threshold':9,'continuously_occupied_units':5}}


def component_sets(states,budget,cost):
    rest=set(states);parts=[]
    while rest:
        reached={rest.pop()};todo=list(reached)
        while todo:
            x=todo.pop();new={y for y in rest if cost(x,y)<=budget}
            rest-=new;reached|=new;todo.extend(new)
        parts.append(frozenset(reached))
    return set(parts)

def conditional_charge_quotients(unique):
    entries=list(unique[3].items());patterns=[(0,0,0),(0,1,2),(0,2,1),(0,1,0)]
    chargecodes=[(0,0,0),(1,0,0),(1,2,0),(1,1,1),(2,1,2),(3,0,1),(3,3,3),(1,2,3)]
    pairs=[(0,7),(7,0),(1,18),(18,50),(50,108),(108,1),(7,50),(50,7)]
    weights=(1,2,1);cases=0;nonempty=0
    for ai,bi in pairs:
        charts=[entries[ai],entries[bi]]
        for codes in chargecodes:
            charge=[patterns[c] for c in codes]
            for target in range(4):
                states=[(h,)+x for h,(rows,_) in enumerate(charts) for x in rows if (h+sum(charge[i][x[i]] for i in range(3)))%4==target]
                fullcost=lambda a,b: int(a[0]!=b[0])+sum(weights[i] for i in range(3) if a[i+1]!=b[i+1])
                for budget in range(1,5):
                    visible=[]
                    for h,(rows,data) in enumerate(charts):
                        base=rows[0];free=data[3];keep=[]
                        for part in free:
                            Q={sum(charge[i][(base[i]+t)%3] for i in part)%4 for t in range(3)}
                            if len(Q)>1 or sum(weights[i] for i in part)>budget:keep.append(min(part))
                        visible.append(keep)
                    fibres=defaultdict(list)
                    for x in states:fibres[(x[0],)+tuple(x[i+1] for i in visible[x[0]])].append(x)
                    for fibre in fibres.values():
                        require(len(component_sets(fibre,budget,fullcost))==1, 'conditional charged quotient fibres connected')
                    qcost=lambda a,b:min(fullcost(x,y) for x in fibres[a] for y in fibres[b])
                    qparts=component_sets(list(fibres),budget,qcost)
                    lifted={frozenset(x for label in part for x in fibres[label]) for part in qparts}
                    require(lifted==component_sets(states,budget,fullcost), 'varying-chart full charged connectivity versus retained quotient')
                    cases+=1;nonempty+=bool(states)
    return {'cases':cases,'nonempty_cases':nonempty,'old_slots':3,'helper_modes':2,'weights':list(weights),'charge_group':'Z/4','scope':'Abstract exact additive sidecars, not an asserted electrical realization.'}

# Tournament codes and actual partner-matrix restrictions for independent global completion checks.
def matrix(n,code):
    S=[[0]*n for _ in range(n)]
    for j,(a,b) in enumerate(combinations(range(n),2)):
        S[a][b]=(code>>j)&1;S[b][a]=1-S[a][b]
    return tuple(map(tuple,S))

def restriction(S,C):
    return tuple(S[u][v] for u,v in combinations(C,2))

def forest_blocks(n,paths):
    edges=set(e for P in paths for e in zip(P,P[1:]))
    indeg=Counter(v for u,v in edges);outdeg=Counter(u for u,v in edges)
    if any(d>1 for d in indeg.values()) or any(d>1 for d in outdeg.values()):return None
    parent=list(range(n))
    def find(x):
        while parent[x]!=x:x=parent[x]
        return x
    for u,v in edges:
        a,b=find(u),find(v)
        if a==b:return None
        parent[a]=b
    successor=dict(edges);out=[]
    for v in range(n):
        if indeg[v]:continue
        b=[v]
        while b[-1] in successor:b.append(successor[b[-1]])
        out.append(tuple(b))
    return out

def predicted(n,paths,e,pins=()):
    B=forest_blocks(n,paths)
    if B is None:return Counter()
    f=e-len(B)+1
    require(f>=0,'uncovered-pair lower bound')
    poly=Counter()
    for order in permutations(B):
        P=sum(order,())
        if any(P[c]!=v for v,c in pins):continue
        forward=sum(a[-1]<b[0] for a,b in zip(order,order[1:]))
        for j in range(f+1):poly[j+forward]+=comb(f,j)
    return poly

def observer_completion():
    rows=[]
    for n in range(1,5):
        verts=tuple(range(n));pairs=list(combinations(verts,2))
        hosts=[]
        for code in range(1<<len(pairs)):
            S=matrix(n,code)
            for P in permutations(verts):
                if all(S[u][v] for u,v in zip(P,P[1:])):hosts.append((S,P,K(S,P),code))
        if n<=3:
            cellmenu=[C for k in range(n+1) for C in combinations(verts,k)]
            families=[tuple(cellmenu[j] for j in range(len(cellmenu)) if (mask>>j)&1) for mask in range(1<<len(cellmenu))]
        else:
            cellmenu=[C for k in (2,3,4) for C in combinations(verts,k)]
            families=[family for k in range(4) for family in combinations(cellmenu,k)]
            families += [((),(0,),(1,),(2,),(3,))]
        cases=feasible=marked=pin_cases=0
        for cells in families:
            covered={p for C in cells for p in combinations(C,2)};cv=sorted(covered);uncovered=set(pairs)-covered;e=len(uncovered)
            actual=defaultdict(list)
            for S,P,T,code in hosts:
                signature=tuple((restriction(S,C),restriction(T,C)) for C in cells)
                weight=sum(S[u][v] for u,v in uncovered)
                actual[signature].append((weight,P,code))
            for bits in product((0,1),repeat=len(cv)):
                S=[[0]*n for _ in range(n)]
                for (u,v),b in zip(cv,bits):S[u][v]=b;S[v][u]=1-b
                menus=[[P for P in permutations(C) if all(S[u][v] for u,v in zip(P,P[1:]))] for C in cells]
                for paths in product(*menus):
                    signature=[]
                    for C,P in zip(cells,paths):
                        adjacent={frozenset((u,v)) for u,v in zip(P,P[1:])}
                        local_K=tuple(S[u][v] if frozenset((u,v)) in adjacent else S[v][u] for u,v in combinations(C,2))
                        signature.append((restriction(S,C),local_K))
                    matches=actual[tuple(signature)]; got=Counter(w for w,P,code in matches)
                    expected=predicted(n,paths,e)
                    require(got==expected, 'full source and actual K restriction completion polynomial')
                    B=forest_blocks(n,paths)
                    require(bool(matches)==(B is not None), 'directed linear forest exact existence')
                    if B is not None:
                        require(len(matches)==factorial(len(B))*2**(e-len(B)+1), 'labelled marked completion count')
                    reverse=predicted(n,[P[::-1] for P in paths],e)
                    require(reverse==Counter({e-j:v for j,v in expected.items()}),'converse conditioned polynomial reciprocity')
                    # All single pins for n<=3; regular deterministic selected n=4 inputs.
                    if n<=3 or cases%29==0:
                        for v,c in product(verts,repeat=2):
                            pins=((v,c),)
                            require(predicted(n,paths,e,pins)==Counter(w for w,P,code in matches if P[c]==v),'absolute global clock single pins')
                            pin_cases+=1
                        if n>=2:
                            for c,d in product(verts,repeat=2):
                                pins=((0,c),(n-1,d))
                                require(predicted(n,paths,e,pins)==Counter(w for w,P,code in matches if P[c]==0 and P[d]==n-1),'absolute global clock simultaneous pins')
                                pin_cases+=1
                    cases+=1;feasible+=bool(matches);marked+=len(matches)
        rows.append({'n':n,'cell_families':len(families),'local_instances':cases,'feasible_instances':feasible,'marked_completions':marked,'clock_pin_cases':pin_cases})
    # Fixed-host versus free-completion hostile: blocks 0->1 and 2->3 cannot connect.
    S=[[0]*4 for _ in range(4)]
    orientations=[(0,1),(2,3),(2,1),(0,3),(0,2),(1,3)]
    for u,v in orientations:S[u][v]=1;S[v][u]=0
    paths=((0,1),(2,3));blocks=forest_blocks(4,paths)
    admissible=[sum(O,()) for O in permutations(blocks) if all(S[a[-1]][b[0]] for a,b in zip(O,O[1:]))]
    require(blocks is not None and not admissible,'fixed-host transition graph separate obligation')
    # Two distinct block orders share some source hosts: 16 marked, exactly 12 unmarked.
    source_codes=[]
    for code in range(64):
        T=matrix(4,code)
        if T[0][1] and T[2][3]:
            orders=[P for P in ((0,1,2,3),(2,3,0,1)) if all(T[u][v] for u,v in zip(P,P[1:]))]
            source_codes.extend([code]*len(orders))
    require(len(source_codes)==16 and len(set(source_codes))==12,'marked versus unmarked hostile')
    # Subsequence does not imply a retained local partner.
    S=matrix(3,7);P=(0,1,2);local=((0,2),)
    require(K(S,P)[0][2]!=S[0][2] and forest_blocks(3,local) is not None,'subsequence alone loses actual K restriction')
    return {'exhaustive_actual_K_controls':rows,'fixed_host_counterexample':{'blocks':[[0,1],[2,3]],'free_marked_completions':16,'fixed_host_observers':0},'marked_unmarked_control':{'marked_pairs':16,'distinct_sources':12}}

def hamiltonian_count_dp(S):
    n=len(S);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,value in enumerate(dp[mask]):
            if not value:continue
            for u in range(n):
                if not(mask>>u)&1 and S[v][u]:dp[mask|1<<u][u]+=value
    return sum(dp[-1])

def eulerian():
    rows=[];A=[1]
    for n in range(1,7):
        if n>1:
            old=A;A=[(k+1)*(old[k] if k<len(old) else 0)+(n-k)*(old[k-1] if k>0 else 0) for k in range(n)]
        E=n*(n-1)//2;observed=Counter()
        for code in range(1<<E):observed[code.bit_count()]+=hamiltonian_count_dp(matrix(n,code))
        predicted_poly=Counter()
        for k,c in enumerate(A):
            for j in range(E-n+2):predicted_poly[j+k]+=c*comb(E-n+1,j)
        require(observed==predicted_poly,'Eulerian identity versus all-source subset dynamic programming')
        require(sum(observed.values())==factorial(n)*2**(E-n+1),'empty-family total marked observers')
        rows.append({'n':n,'labelled_sources':1<<E,'Eulerian_ascents':A[:],'polynomial':[observed[k] for k in range(E+1)],'marked_total':sum(observed.values())})
    return rows

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    pins=validate_inputs();nr=native();vc,unique=voltage_and_packets();bc=boundary_and_helper(unique);cq=conditional_charge_quotients(unique);oc=observer_completion();ec=eulerian()
    out={'status':'PASS','producer_manifest_sha256':pins,'checks':sum(CHECKS.values()),'checks_by_kind':dict(sorted(CHECKS.items())),'native':nr,'voltage':vc,'boundaries':bc,'conditional_charge_quotients':cq,'completion':oc,'Eulerian':ec,'scope':'Independent native matrices and full K restrictions; no producer imports, no Git, no external runtime dependencies. All-size claims proved in AUDIT.md; finite controls are not extrapolated.'}
    s=json.dumps(out,indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(s)
    else:print(s,end='')

if __name__=='__main__':main()
