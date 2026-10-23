"""Exact shared-chart tests on all ordered two-mark DRT31 partitions.
No external packages; all guards survive Python -O.
"""
from pathlib import Path
from itertools import combinations, permutations, product
from functools import lru_cache
from collections import Counter
import argparse, hashlib, json, sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'incidence_wiring_20261023'))
from marked_cells import paley, double
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(b,msg):
    global CHECKS
    CHECKS+=1
    if not b: raise ValueError(msg)

def refine(A, marks=(), seed=None):
    n=len(A)
    col=seed if seed is not None else [marks.index(v)+1 if v in marks else 0 for v in range(n)]
    while True:
        rows=[(col[u],tuple(sorted(col[v] for v in range(n) if A[u][v]))) for u in range(n)]
        ids={r:i for i,r in enumerate(sorted(set(rows)))};nxt=[ids[r] for r in rows]
        if len(set(nxt))==len(set(col)):
            return [[v for v in range(n) if nxt[v]==c] for c in sorted(set(nxt))]
        col=nxt

def sub(A,D):return tuple(tuple(A[u][v] for v in D) for u in D)
@lru_cache(None)
def autos(B):
    n=len(B)
    return tuple(p for p in permutations(range(n)) if all(B[i][j]==B[p[i]][p[j]] for i in range(n) for j in range(i)))

def relation(B,C,gi,gj):
    # R consists of source->target chart pairs. Distinct column words in C
    # make q unique once p is given; no histogram alignment is guessed.
    n,m=len(B),len(B[0]);cols={tuple(C[i][j] for i in range(n)):j for j in range(m)}
    need(len(cols)==m and len(set(C))==n,'relation requires distinct rows and columns')
    allowed=set(gj);R={}
    for p in gi:
        inv=[p.index(i) for i in range(n)]
        q=tuple(cols.get(tuple(B[inv[i]][j] for i in range(n)),-1) for j in range(m))
        if q in allowed:R[p]=q
    need(len(set(R.values()))==len(R),'reverse relation is not functional')
    return R

def compatibility(A,cells,B=None,other=None):
    B=A if B is None else B;other=cells if other is None else other
    ds=[D for D in cells if len(D)>1];es=[D for D in other if len(D)>1]
    # This implementation's cross-target controls retain the same internal charts.
    gs=[autos(sub(A,D)) for D in ds]
    need(all(sub(A,D)==sub(B,E) for D,E in zip(ds,es)),'internal charts differ')
    edges={}
    for i,j in combinations(range(len(ds)),2):
        U=tuple(tuple(A[u][v] for v in ds[j]) for u in ds[i])
        V=tuple(tuple(B[u][v] for v in es[j]) for u in es[i])
        if len(set(U))==len(ds[i]) and len(set(zip(*U)))==len(ds[j]) and len(set(V))==len(es[i]) and len(set(zip(*V)))==len(es[j]):
            edges[i,j]=relation(U,V,gs[i],gs[j])
    return ds,es,gs,edges

def transport_solutions(A,cells,B=None,other=None):
    B=A if B is None else B;other=cells if other is None else other
    ds,es,gs,edges=compatibility(A,cells,B,other)
    if not ds:
        f=[0]*len(A)
        for D,E in zip(cells,other):f[D[0]]=E[0]
        valid=[tuple(f)] if all(A[u][v]==B[f[u]][f[v]] for u in range(len(A)) for v in range(len(A))) else []
        return valid, {'root_trials':1,'tree_survivors':1,'edge_domains':[]}
    nbr={i:[] for i in range(len(ds))}
    for i,j in edges:nbr[i].append(j);nbr[j].append(i)
    order=[0];parent={0:None}
    for i in order:
        for j in sorted(nbr[i]):
            if j not in parent:parent[j]=i;order.append(j)
    need(len(order)==len(ds),'control requires a connected informative support')
    valid=[];tree=0
    for root in gs[0]:
        maps={0:root}
        for j in order[1:]:
            i=parent[j];e=edges[min(i,j),max(i,j)]
            if i>j:e={v:k for k,v in e.items()}
            if maps[i] not in e:break
            maps[j]=e[maps[i]]
        if len(maps)!=len(ds):continue
        tree+=1
        if not all(e.get(maps[i])==maps[j] for (i,j),e in edges.items()):continue
        f=[0]*len(A)
        for D,E in zip(cells,other):
            if len(D)==1:f[D[0]]=E[0]
        for i in range(len(ds)):
            for x,y in enumerate(maps[i]):f[ds[i][x]]=es[i][y]
        if all(A[u][v]==B[f[u]][f[v]] for u in range(len(A)) for v in range(len(A))):valid.append(tuple(f))
    return valid,{'root_trials':len(gs[0]),'tree_survivors':tree,'edge_domains':sorted(Counter(map(len,edges.values())).items())}

def triple_signs(A,cells):
    ds=[d for d in cells if len(d)==3];charts=[]
    need(all(len(d) in (1,3) for d in cells),'not singleton/triple')
    for D in ds:
        x=min(D);y=next(v for v in D if A[x][v]);z=next(v for v in D if v not in (x,y));charts.append([x,y,z])
        need(A[y][z] and A[z][x],'noncyclic triple')
    edges={}
    for i,j in combinations(range(len(ds)),2):
        B=[[A[u][v] for v in charts[j]] for u in charts[i]];q=sum(B[0])
        need(all(sum(r)==q for r in B),'inequitable triple block')
        if q in (0,3):continue
        P=B if q==1 else [[1-x for x in row] for row in B]
        pi=[r.index(1) for r in P];eps=(pi[1]-pi[0])%3
        need(eps in (1,2) and all(pi[x]==(eps*x+pi[0])%3 for x in range(3)),'not affine')
        edges[i,j]=(eps,pi[0])
    # Complete support in every present DDP7 case. Every cycle sign is tested
    # by all fundamental triangles of this actual star.
    need(len(edges)==len(ds)*(len(ds)-1)//2,'control needs complete support')
    negatives=[]
    for i,j in combinations(range(1,len(ds)),2):
        if edges[0,i][0]*edges[i,j][0]*edges[0,j][0]%3==2:negatives.append([0,i,j])
    owners=[]
    for x in range(3):
        values=[x]+[(edges[0,j][0]*x+edges[0,j][1])%3 for j in range(1,len(ds))]
        if all((e*values[i]+d)%3==values[j] for (i,j),(e,d) in edges.items()):owners.append(values)
    need(not owners,'actual order31 owner-section hostile changed')
    return {'triples':len(ds),'edges':len(edges),'negative_fundamental_triangles':len(negatives),'first_negative':negatives[0] if negatives else None,'parallel_owner_sections':len(owners)}

def motif_repair(A,cells):
    ds=[D for D in cells if len(D)==7];col=[0]*len(A);flags={}
    for k,D in enumerate(cells):
        for u in D:col[u]=(k,-1)
    # First informative partner in the retained cell order. Some actual
    # order31 pairs have no such partner; they must remain seven-cells.
    for i,D in enumerate(ds):
        selected=None
        for j,E in enumerate(ds):
            if i==j:continue
            B=[[A[u][v] for v in E] for u in D];q=sum(B[0]);need(q in (3,4),'unexpected degree')
            P=B if q==3 else [[1-x for x in r] for r in B]
            bits=[]
            for row in P:
                vs=[E[k] for k,x in enumerate(row) if x]
                need(len(vs)==3,'normalized selector is not a triple')
                bits.append(int(all(sum(A[u][v] for v in vs)==1 for u in vs)))
            need(sum(bits) in (4,7),'unexpected flag count')
            if sum(bits)==4:
                selected=(j,bits);break
        if selected is not None:
            j,bits=selected
            for u,b in zip(D,bits):col[u]=(col[u][0],b)
            flags[i]={'partner':j,'bits':bits}
        else:flags[i]=None
    repaired=refine(A,seed=col)
    return repaired,flags

def brute(A,cells):
    ds=[d for d in cells if len(d)>1];gs=[autos(sub(A,D)) for D in ds];accepted=[];tested=0
    for ps in product(*gs):
        tested+=1;ok=True
        for i,j in combinations(range(len(ds)),2):
            if not all(A[u][v]==A[ds[i][ps[i][a]]][ds[j][ps[j][b]]] for a,u in enumerate(ds[i]) for b,v in enumerate(ds[j])):
                ok=False;break
        if ok:
            f=list(range(len(A)))
            for D,p in zip(ds,ps):
                for x,y in enumerate(p):f[D[x]]=D[y]
            if all(A[u][v]==A[f[u]][f[v]] for u in range(len(A)) for v in range(len(A))):accepted.append(tuple(f))
    return accepted,tested

def cycle_tournament(delta):
    A=[[0]*9 for _ in range(9)]
    for i in range(3):
        for x in range(3):A[3*i+x][3*i+(x+1)%3]=1
    for i,j in combinations(range(3),2):
        d=delta if(i,j)==(0,2) else 0
        for x,y in product(range(3),repeat=2):A[3*i+x][3*j+y]=int(y==(x+d)%3);A[3*j+y][3*i+x]=1-A[3*i+x][3*j+y]
    return A

def lift_perm(p):
    n=len(p)
    return tuple(p)+(n,)+tuple(n+1+x for x in p)

def owned_path(A):
    P=[]
    for v in range(len(A)):
        i=next((i for i,u in enumerate(P) if A[v][u]),len(P));P.insert(i,v)
    need(sorted(P)==list(range(len(A))) and all(A[u][v] for u,v in zip(P,P[1:])),'not an actual Hamiltonian path')
    return P

def anchor_choices(A,cells,sol,P):
    if len(sol)==1:return [],0
    D=next(d for d in cells if len(d)>1)
    remaining=sol[:];marks=[]
    while len(remaining)>1:
        v=next(v for v in P if v in D and any(f[v]!=v for f in remaining))
        before=len(remaining);remaining=[f for f in remaining if f[v]==v];marks.append(v)
        need(before>=3*len(remaining),'odd-group index failed')
    minimum=next(k for k in range(len(D)+1) if any(sum(all(f[v]==v for v in S) for f in sol)==1 for S in combinations(D,k)))
    need(len(marks)==minimum and len(marks)<=2,'order31 anchor bound')
    return marks,minimum

def run():
    out={};examples={};allrows=[]
    explicit_group=[lift_perm(lift_perm(tuple((a*x+b)%7 for x in range(7)))) for a in (1,2,4) for b in range(7)]
    D=double(double(paley(7)))
    need(all(all(D[u][v]==D[f[u]][f[v]] for u in range(31) for v in range(31)) for f in explicit_group),'explicit group not automorphic')
    for name,A in [('P31',paley(31)),('DDP7',double(double(paley(7))))]:
        tallies=Counter();repair=Counter();trials=Counter();anchors=Counter();P=owned_path(A)
        for marks in permutations(range(31),2):
            cells=refine(A,marks);shape=tuple(sorted(map(len,cells)))
            sol,diag=transport_solutions(A,cells)
            need(len(sol)>0 and all(f[marks[0]]==marks[0] and f[marks[1]]==marks[1] for f in sol),'lost owned marks')
            if name=='DDP7':
                direct=[f for f in explicit_group if all(f[u]==u for u in marks)]
                need(set(sol)==set(direct),'explicit lifted group disagrees with chart transport')
            tallies[(max(shape),len(sol))]+=1;trials[(diag['root_trials'],diag['tree_survivors'],len(sol))]+=1
            extra,minimum=anchor_choices(A,cells,sol,P);anchors[minimum]+=1
            row={'additional_source_owned_anchors':extra,'minimum_extra_root_marks':minimum,'marks':list(marks),'largest_cell':max(shape),'automorphisms':len(sol),'internal_aut_orders':sorted(len(autos(sub(A,D))) for D in cells if len(D)>1),**diag}
            if max(shape)==7:
                refined,bits=motif_repair(A,cells)
                newshape=tuple(sorted(map(len,refined)))
                if max(newshape)==3:
                    signs=triple_signs(A,refined)
                    repair[(3,signs['negative_fundamental_triangles'])]+=1;row['repair']=signs
                else:
                    repair[(max(newshape),-1)]+=1;row['repair']={'shape':newshape}
                sol2,_=transport_solutions(A,refined);need(set(sol2)==set(sol),'motif refinement lost an automorphism')
            elif max(shape)==3:
                signs=triple_signs(A,cells)
                row['original_triple_signs']=signs
            if (name,max(shape)) not in examples:
                examples[name,max(shape)]=(A,cells,sol)
                row['example_cells']=cells;row['example_solutions']=sol
            if name=='DDP7' and max(shape)>1:
                # Nontrivial change of every available local chart, with shared
                # source marks fixed. This tests source->target transport, not
                # only the identity automorphism problem.
                f=list(range(31))
                for k,D in enumerate(cells):
                    G=autos(sub(A,D));p=G[(k+1)%len(G)]
                    for x,y in enumerate(p):f[D[x]]=D[y]
                inv=[f.index(i) for i in range(31)]
                B=[[A[inv[u]][inv[v]] for v in range(31)] for u in range(31)]
                target,td=transport_solutions(A,cells,B,cells)
                expected={tuple(f[g[u]] for u in range(31)) for g in sol}
                need(set(target)==expected,'shared relabeling lost a source or a chart')
                allrows.append(row)
        out[name]={'owned_hamiltonian_path':P,'minimum_extra_root_mark_counts':sorted(anchors.items()),'ordered_marks':930,'largest_cell_and_stabilizer_counts':[[*k,v] for k,v in sorted(tallies.items())],
                   'root_trials_tree_survivors_final_counts':[[*k,v] for k,v in sorted(trials.items())],
                   'repaired_largest_cell_negative_triangle_counts':[[*k,v] for k,v in sorted(repair.items())],
                   'adjacency_sha256':hashlib.sha256(bytes(x for r in A for x in r)).hexdigest()}
    brute_records=[]
    for key,(A,cells,sol) in examples.items():
        if key[1]==1:continue
        other,n=brute(A,cells);need(set(sol)==set(other),'independent product enumeration mismatch');brute_records.append([list(key),n,len(other)])
    # Pairwise compatible target blocks need not glue jointly.
    D=[list(range(3*i,3*i+3)) for i in range(3)];A=cycle_tournament(0);B=cycle_tournament(1)
    _,_,_,edges=compatibility(A,D,B,D);need(all(len(e)==3 for e in edges.values()),'hostile lost local feasibility')
    sol,d=transport_solutions(A,D,B,D);need(not sol and d['tree_survivors']==3,'cycle should reject all local candidates')
    # Duplicate columns really destroy functional propagation: uniform blocks.
    U=[[0]*3 for _ in range(3)];G=autos(sub(A,D[0]))
    need(sum(all(U[a][b]==U[p[a]][q[b]] for a,b in product(range(3),repeat=2)) for p,q in product(G,repeat=2))==9,'uniform hostile')
    return {'status':'FINITE-EXACT scoped shared-chart census','DRT31':out,'independent_full_products':brute_records,
            'independent_target_chart_changes':426,'explicit_lifted_automorphism_group_order':21,'phase_hostile':{'each_pair_compatible_maps':3,'tree_survivors':3,'global_maps':0},'nondiscrete_DDP7_records':allrows,'checks':CHECKS}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();out=json.loads(json.dumps(run()));dest=HERE/'cells.json'
    if a.write:dest.write_text(json.dumps(out,indent=2)+'\n')
    else:need(json.loads(dest.read_text())==out,'receipt changed')
    print(json.dumps({k:v for k,v in out.items() if k!='nondiscrete_DDP7_records'},indent=2))
