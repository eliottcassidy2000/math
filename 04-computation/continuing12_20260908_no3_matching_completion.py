#!/usr/bin/env python3
"""Exact finite controls for finite-direction matching and incidence branching."""
from collections import Counter, defaultdict, deque
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations
from pathlib import Path
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
HERE=Path(__file__).resolve()
GATES=Counter()


def check(name,condition):
    GATES[name]+=1
    if not condition:raise RuntimeError(name)


def boards(n):
    pairs=tuple(combinations(range(n),2));degree=[0]*n;rows=[]
    def visit(r):
        if r==n:
            if all(d==2 for d in degree):yield tuple(rows)
            return
        for pair in pairs:
            if any(degree[c]==2 for c in pair):continue
            for c in pair:degree[c]+=1
            if all(2-d<=n-r-1 for d in degree):
                rows.append(pair);yield from visit(r+1);rows.pop()
            for c in pair:degree[c]-=1
    yield from visit(0)


def points(rows):return tuple((r,c) for r,pair in enumerate(rows) for c in pair)


def lines(pts,slopes):
    out=[]
    for slope in slopes:
        fam=defaultdict(int)
        for i,(r,c) in enumerate(pts):fam[c-slope*r]|=1<<i
        out.extend((slope,label,mask) for label,mask in sorted(fam.items()))
    return out


def delete_exact(masks):
    triples=set()
    for mask in masks:
        bits=[1<<i for i in range(mask.bit_length()) if mask>>i&1]
        triples.update(sum(t) for t in combinations(bits,3))
    @lru_cache(None)
    def hit(remaining):
        if not remaining:return 0
        first=remaining[0];best=None
        while first:
            bit=first&-first;first-=bit
            cand=bit|hit(tuple(t for t in remaining if not t&bit))
            if best is None or cand.bit_count()<best.bit_count():best=cand
        return best
    mask=hit(tuple(sorted(triples)))
    return mask.bit_count(),mask


def tau(pts,slopes):return delete_exact([m for _,_,m in lines(pts,slopes)])


def graph_data(masks,ncells):
    incidence=[tuple(v for v,m in enumerate(masks) if m>>i&1) for i in range(ncells)]
    adj=[set() for _ in masks];edge_cell={};private=defaultdict(list)
    for i,owners in enumerate(incidence):
        if len(owners)==2:
            a,b=owners;adj[a].add(b);adj[b].add(a)
            check("distinct_shared_cell", (a,b) not in edge_cell)
            edge_cell[a,b]=i
        elif len(owners)==1:private[owners[0]].append(i)
    return incidence,adj,edge_cell,private


def matching(adj):
    @lru_cache(None)
    def solve(mask):
        if not mask:return ()
        bit=mask&-mask;v=bit.bit_length()-1;rest=mask^bit
        best=solve(rest)
        for u in sorted(adj[v]):
            if rest>>u&1:
                cand=((min(u,v),max(u,v)),)+solve(rest^(1<<u))
                if len(cand)>len(best):best=cand
        return best
    return solve((1<<len(adj))-1)


def independence(adj):
    best=0
    for mask in range(1<<len(adj)):
        if mask.bit_count()<=best:continue
        if all(not (mask>>a&1 and mask>>b&1) for a in range(len(adj)) for b in adj[a] if a<b):best=mask.bit_count()
    return best


def components(adj):
    seen=set();out=[]
    for v in range(len(adj)):
        if v in seen:continue
        stack=[v];seen.add(v);comp=[]
        while stack:
            a=stack.pop();comp.append(a)
            for b in adj[a]-seen:seen.add(b);stack.append(b)
        out.append(comp)
    return out


def boolean(masks):
    best=0
    for chosen in range(1<<len(masks)):
        union=0
        for j,m in enumerate(masks):
            if chosen>>j&1:union|=m
        best=max(best,union.bit_count()-2*chosen.bit_count())
    return best


def double_cover(adj):
    right={}
    def augment(a,seen):
        for b in sorted(adj[a]):
            if b in seen:continue
            seen.add(b)
            if b not in right or augment(right[b],seen):right[b]=a;return True
        return False
    for a in range(len(adj)):augment(a,set())
    left={a:b for b,a in right.items()}
    reach_l=set(range(len(adj)))-set(left);reach_r=set();queue=deque(reach_l)
    while queue:
        a=queue.popleft()
        for b in adj[a]:
            if left.get(a)==b or b in reach_r:continue
            reach_r.add(b)
            if b in right and right[b] not in reach_l:reach_l.add(right[b]);queue.append(right[b])
    cover_l=set(range(len(adj)))-reach_l;cover_r=reach_r
    check("double_cover_min_cut", len(cover_l)+len(cover_r)==len(left))
    check("double_cover_vertex_cover", all(a in cover_l or b in cover_r for a in range(len(adj)) for b in adj[a]))
    return left,cover_l,cover_r


def certify_triples(masks,ncells):
    inc,adj,edge_cell,private=graph_data(masks,ncells)
    check("triple_regime_supplier", all(m.bit_count()==3 for m in masks) and all(len(a)<=2 for a in inc))
    v=len(masks);mat=matching(adj);nu=len(mat);alpha=independence(adj)
    exact,deleted=delete_exact(masks)
    selected={edge_cell[e] for e in mat};matched={x for e in mat for x in e}
    for a in range(v):
        if a not in matched:
            options=private[a] or [edge_cell[tuple(sorted((a,b)))] for b in sorted(adj[a])]
            selected.add(options[0])
    integer_mask=sum(1<<i for i in selected)
    check("matching_cover_witness", len(selected)==v-nu and all((m&integer_mask).bit_count()>=1 for m in masks))
    check("matching_vs_exact_deletion", exact==v-nu)
    check("boolean_vs_independence", boolean(masks)==alpha)
    double,cl,cr=double_cover(adj);M=len(double)
    f={e:Q(int(double.get(e[0])==e[1])+int(double.get(e[1])==e[0]),2) for e in edge_cell}
    x=[Q(0)]*ncells
    for e,z in f.items():x[edge_cell[e]]=z
    for a in range(v):
        deficit=1-sum(f[tuple(sorted((a,b)))] for b in adj[a])
        check("fractional_matching_capacity", deficit>=0)
        if deficit:
            options=private[a] or [edge_cell[tuple(sorted((a,b)))] for b in sorted(adj[a])]
            x[options[0]]+=deficit
    y=[1-Q(int(a in cl)+int(a in cr),2) for a in range(v)]
    check("fractional_deletion_primal", all(0<=z<=1 for z in x) and all(sum(x[i] for i in range(ncells) if m>>i&1)>=1 for m in masks))
    check("fractional_deletion_dual", all(z>=0 for z in y) and all(sum(y[a] for a in owners)<=1 for owners in inc))
    check("fractional_strong_certificate", sum(x)==sum(y)==v-Q(M,2))
    return dict(vertices=v,edges=sorted(edge_cell),matching=mat,nu=nu,alpha=alpha,
                double_cover_matching=sorted(double.items()),tau=exact,beta=alpha,
                tau_lp=str(v-Q(M,2)),integer_delete_mask=integer_mask,
                fractional_delete=[str(z) for z in x],fractional_line_dual=[str(z) for z in y])


def graph_bmatching(edges,capacities):
    @lru_cache(None)
    def take(i,caps):
        if i==len(edges):return 0,0
        best=take(i+1,caps);a,b,cell=edges[i]
        if caps[a] and caps[b]:
            nxt=list(caps);nxt[a]-=1;nxt[b]-=1
            count,mask=take(i+1,tuple(nxt));cand=count+1,mask|(1<<cell)
            if cand[0]>best[0]:best=cand
        return best
    return take(0,tuple(capacities))


def exceptional_branch(masks,ncells):
    inc=[tuple(a for a,m in enumerate(masks) if m>>i&1) for i in range(ncells)]
    H=[i for i,a in enumerate(inc) if len(a)>=3]
    C=[i for i,a in enumerate(inc) if 1<=len(a)<=2]
    best=-1;bestmask=0;valid=0
    for choice in range(1<<len(H)):
        caps=[2]*len(masks);retained=0
        for j,i in enumerate(H):
            if choice>>j&1:
                retained|=1<<i
                for a in inc[i]:caps[a]-=1
        if min(caps,default=0)<0:continue
        valid+=1;edges=[]
        for i in C:
            owners=inc[i]
            if len(owners)==1:
                leaf=len(caps);caps.append(1);edges.append((owners[0],leaf,i))
            else:edges.append((*owners,i))
        value,mask=graph_bmatching(edges,caps)
        value+=choice.bit_count()
        if value>best:best=value;bestmask=mask|retained
    exact,deleted=delete_exact(masks)
    check("exceptional_branch_exact", len(H)+len(C)-best==exact)
    check("exceptional_retention_witness", all((m&bestmask).bit_count()<=2 for m in masks))
    return dict(exceptional_cells=H,valid_branches=valid,total_branches=1<<len(H),tau=exact,retained_mask=bestmask)


complete={};small_board_bank=[]
for n in range(2,6):
    stats=Counter()
    for rows in boards(n):
        pts=points(rows);over=[m for _,_,m in lines(pts,(1,-1,2)) if m.bit_count()>2]
        stats["all"]+=1
        branch=exceptional_branch(over,len(pts))
        stats["with_exceptional_cells"]+=bool(branch["exceptional_cells"])
        if n<=4:small_board_bank.append(pts)
        inc=[sum(m>>i&1 for m in over) for i in range(len(pts))]
        if all(m.bit_count()==3 for m in over) and max(inc,default=0)<=2:
            stats["regime"]+=1;data=certify_triples(over,len(pts))
            stats["fractional_gap"]+=Q(data["tau_lp"])<data["tau"]
            adj=[set() for _ in over]
            for a,b in data["edges"]:adj[a].add(b);adj[b].add(a)
            for comp in components(adj):
                if all(len(adj[a])==2 for a in comp):
                    stats["five_cycle"]+=len(comp)==5
                    stats["triangle_component"]+=len(comp)==3
    check("complete_grid_cardinality", stats["all"]=={2:1,3:6,4:90,5:2040}[n])
    check("complete_grid_regime_count", stats["regime"]=={2:1,3:6,4:73,5:1240}[n])
    check("no_five_cycle_before_six", stats["five_cycle"]==0)
    complete[n]=dict(stats)

named_rows={
    "triangle":((0,3),(1,2),(0,4),(3,4),(1,2)),
    "five_cycle":((0,3),(1,2),(0,4),(1,4),(3,5),(2,5)),
    "concurrent":((2,4),(0,1),(2,3),(3,4),(0,1)),
    "triangle_tail":((0,1),(0,2),(3,4),(1,3),(2,4)),
    "larger_line":((0,1),(0,1),(2,3),(2,4),(3,4)),
}
witnesses={}
for name,rows in named_rows.items():
    pts=points(rows);all_lines=lines(pts,(1,-1,2));over=[m for _,_,m in all_lines if m.bit_count()>2]
    record=dict(rows=rows,overfull=[(s,k,m) for s,k,m in all_lines if m.bit_count()>2],branch=exceptional_branch(over,len(pts)))
    if name in ("triangle","five_cycle","triangle_tail"):
        record["matching"]=certify_triples(over,len(pts))
    witnesses[name]=record
    full_four=[m for _,_,m in lines(pts,(1,-1,2,-2)) if m.bit_count()>2]
    record["four_direction_branch"]=exceptional_branch(full_four,len(pts))
    small_board_bank.append(pts)
check("five_cycle_exact_values", tuple(witnesses["five_cycle"]["matching"][k] for k in ("vertices","nu","beta","tau_lp","tau"))==(5,2,2,"5/2",3))
five=points(named_rows["five_cycle"])
old_lines=lines(five,(1,-1));unsafe=0
for _,_,m in old_lines:
    if m.bit_count()>2:unsafe|=m
J=sum(max(0,(m&~unsafe).bit_count()-2) for _,_,m in lines(five,(2,)))
check("five_cycle_beats_prior_bounds", tau(five,(1,-1))[0]==2 and J==0)
check("concurrent_hyperedge_hostile", witnesses["concurrent"]["branch"]["exceptional_cells"]==[4] and witnesses["concurrent"]["branch"]["tau"]==1)
check("triangle_tail_hostile", tuple(witnesses["triangle_tail"]["matching"][k] for k in ("beta","tau_lp","tau"))==(2,"2",2))
check("larger_line_outside_edge_cover", any(m.bit_count()>3 for _,_,m in witnesses["larger_line"]["overfull"]))

# Literal capacity correspondence on every retained subset of the five controls.
for rows in named_rows.values():
    pts=points(rows);over=[m for _,_,m in lines(pts,(1,-1,2)) if m.bit_count()>2]
    inc=[tuple(a for a,m in enumerate(over) if m>>i&1) for i in range(len(pts))]
    for retained in range(1<<len(pts)):
        Hkeep=[i for i,a in enumerate(inc) if len(a)>=3 and retained>>i&1]
        capacity=[2-sum(a in inc[i] for i in Hkeep) for a in range(len(over))]
        remaining=[i for i,a in enumerate(inc) if 1<=len(a)<=2 and retained>>i&1]
        branch_feasible=min(capacity,default=0)>=0 and all(sum(a in inc[i] for i in remaining)<=capacity[a] for a in range(len(over)))
        literal=all((m&retained).bit_count()<=2 for m in over)
        check("every_retained_subset_branch_equivalence", branch_feasible==literal)

# Independent abstract incidence universe: all labelled subcubic graphs on 0..5
# vertices. These are algebraic incidence controls, not claimed geometric boards.
abstract_counts={}
for n in range(6):
    pairs=list(combinations(range(n),2));number=0
    for chosen in range(1<<len(pairs)):
        edges=[e for j,e in enumerate(pairs) if chosen>>j&1];deg=Counter(x for e in edges for x in e)
        if max(deg.values(),default=0)>3:continue
        masks=[0]*n;cell=0
        for a,b in edges:masks[a]|=1<<cell;masks[b]|=1<<cell;cell+=1
        for a in range(n):
            for _ in range(3-deg[a]):masks[a]|=1<<cell;cell+=1
        data=certify_triples(masks,cell)
        check("abstract_incidence_graph_reconstruction", data["edges"]==edges)
        number+=1
    abstract_counts[n]=number


def safe_original(pts,slopes):
    unsafe=0
    for _,_,m in lines(pts,slopes):
        if m.bit_count()>2:unsafe|=m
    return tuple(p for i,p in enumerate(pts) if not unsafe>>i&1)


directions=(1,-1,2,-2)
orders=[tuple((s,) for s in order) for order in permutations(directions)]
orders+=sorted({(tuple(sorted(p[:2])),tuple(sorted(p[2:]))) for p in permutations(directions)})
block_count=0;dual_count=0
for pts in small_board_bank:
    full=tau(pts,directions)[0]
    for blocks in orders:
        previous=();U=pts;charged=set();cost=0
        for block in blocks:
            amount,deleted=tau(U,block);cost+=amount
            chosen={p for i,p in enumerate(U) if deleted>>i&1}
            previous+=block;next_U=safe_original(pts,previous)
            check("block_charge_disjoint_stratum", not chosen&set(next_U) and not chosen&charged)
            charged.update(chosen);U=next_U
        check("finite_direction_block_inequality", cost<=full and cost==len(charged))
        block_count+=1
for rows in named_rows.values():
    pts=points(rows);old=[m for _,_,m in lines(pts,(1,-1)) if m.bit_count()>2];full=tau(pts,directions)[0]
    for chosen in range(1<<len(old)):
        union=0
        for i,m in enumerate(old):
            if chosen>>i&1:union|=m
        rest=tuple(p for i,p in enumerate(pts) if not union>>i&1)
        bound=union.bit_count()-2*chosen.bit_count()+tau(rest,(2,-2))[0]
        check("two_new_direction_residual_dual", bound<=full);dual_count+=1

# Explicit fixed-three-direction family; no assertion of saturated bounding grids.
translations=[]
for k in range(1,5):
    pts=tuple((r+100*j,c+300*j) for j in range(k) for r,c in five)
    over=[m for _,_,m in lines(pts,(1,-1,2)) if m.bit_count()>2]
    # Component additivity avoids a needless 2^(5k) independence enumeration.
    inc,adj,edge_cell,private=graph_data(over,len(pts));comps=components(adj)
    check("separated_five_cycle_components", len(comps)==k and all(len(c)==5 for c in comps) and all(len(a)==2 for a in adj))
    check("occupied_rows_columns_two", all(z==2 for z in Counter(r for r,c in pts).values()) and all(z==2 for z in Counter(c for r,c in pts).values()))
    check("separated_exact_deletions", delete_exact(over)[0]==3*k)
    translations.append(dict(copies=k,cells=len(pts),tau=3*k,tau_lp=str(Q(5*k,2)),beta=2*k))

certificate=dict(status="FINITE-EXACT controls; uniform analytical statements in paired report",
    complete_grids=complete,named_witnesses=witnesses,abstract_subcubic_graphs=abstract_counts,
    block_bank=dict(boards=len(small_board_bank),ordered_partitions=len(orders),checks=block_count,residual_dual_checks=dual_count),
    translated_families=translations,gates=dict(sorted(GATES.items())),gate_total=sum(GATES.values()),
    source_sha256=hashlib.sha256(HERE.read_bytes()).hexdigest())
outdir=HERE.parent
if outdir.name=="04-computation":outdir=outdir.parent/"05-knowledge"/"results"
out=outdir/(HERE.stem+"_certificate.json")
out.write_bytes((json.dumps(certificate,indent=2,sort_keys=True)+"\n").encode())
print("PASS continuing12 finite-direction matching completion")
print("Complete grid census:",json.dumps(complete,sort_keys=True))
print("Abstract labelled subcubic graphs:",json.dumps(abstract_counts,sort_keys=True))
print("Block bank:",json.dumps(certificate["block_bank"],sort_keys=True))
print("Five-cycle: tau2=2, J3=0, beta3=2, tauLP=5/2, tau3=3")
print("Concurrent triple: one exceptional cell, exact repair1; triangle-tail: no gap")
print("Always-active gates:",sum(GATES.values()))
print("Certificate:",out.name)
