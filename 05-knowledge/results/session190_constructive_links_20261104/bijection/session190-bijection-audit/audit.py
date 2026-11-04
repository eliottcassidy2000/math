#!/usr/bin/env python3
"""Independent Session190 audit. No imports from audited producers or repository.
Use adjacency matrices and permutation transport; exact rational arithmetic.
"""
import itertools, json, math, hashlib
from fractions import Fraction as F
from pathlib import Path

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'session190-bijection'

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

def parity(k):
    return 1 if k%2==0 else -1

def fwht(a):
    a=list(a); d=1
    while d<len(a):
        for b in range(0,len(a),2*d):
            for k in range(b,b+d):
                x,y=a[k],a[k+d]; a[k]=x+y; a[k+d]=x-y
        d*=2
    return a

def mul(a,b):
    return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]

def rank_exact(a):
    a=[[F(x) for x in row] for row in a]; rank=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(rank,len(a)) if a[i][col]),None)
        if pivot is None:continue
        a[rank],a[pivot]=a[pivot],a[rank]
        divisor=a[rank][col];a[rank]=[x/divisor for x in a[rank]]
        for i in range(rank+1,len(a)):
            m=a[i][col]
            a[i]=[x-m*y for x,y in zip(a[i],a[rank])]
        rank+=1
        if rank==len(a):break
    return rank

def audit(n):
    pairs=list(itertools.combinations(range(n),2)); pos={p:i for i,p in enumerate(pairs)}
    N=len(pairs); size=2**N; full=size-1
    permutations=list(itertools.permutations(range(n)))
    # Precompute each relabelling from actual image endpoint pairs, independently
    # of the packet's bit-permutation/XOR implementation.
    transport=[]
    for p in permutations:
        transport.append([(pos[tuple(sorted((p[u],p[v])))], p[u]<p[v]) for u,v in pairs])
    def graph_image(g, table):
        return sum(1<<j for i,(j,up) in enumerate(table) if (g>>i)&1)
    def tournament_image(t,table):
        return sum(1<<j for i,(j,up) in enumerate(table) if bool((t>>i)&1)==up)
    def sense(t,u,v):
        return bool((t>>pos[tuple(sorted((u,v)))])&1)==(u<v)
    def character(t,g):
        return parity((g&~t).bit_count())
    tclass=[-1]*size; gclass=[-1]*size; ts=[]; gs=[]; even=[]; signed=[]; automorphs=[]; g_aut=[]
    for mask in range(size):
        if tclass[mask]<0:
            orb={tournament_image(mask,tr) for tr in transport}
            k=len(ts); ts.append(mask)
            for t in orb:tclass[t]=k
            automorphs.append([p for p,tr in zip(permutations,transport) if tournament_image(mask,tr)==mask])
        if gclass[mask]<0:
            signs={}; aut=0; is_even=True
            for tr in transport:
                g=graph_image(mask,tr)
                sg=parity(sum(1 for i,(_,up) in enumerate(tr) if mask>>i&1 and not up))
                signs.setdefault(g,set()).add(sg)
                if g==mask:
                    aut+=1
                    if sg<0:is_even=False
            k=len(gs); gs.append(mask); g_aut.append(aut)
            for g in signs:gclass[g]=k
            if is_even:
                check(all(len(s)==1 for s in signs.values()),'transporter inconsistency')
                even.append(mask); signed.append({g:next(iter(s)) for g,s in signs.items()})
    rt=[tclass.count(k) for k in range(len(ts))]; rg=[len(s) for s in signed]
    # Independent all-label Walsh values from coefficient vectors, rather than
    # evaluating each signed-copy polynomial only on canonical representatives.
    W=[[0]*len(even) for _ in ts]
    for j,signs in enumerate(signed):
        a=[0]*size
        for g,sg in signs.items():a[g]=sg*parity(g.bit_count())
        vals=fwht(a)
        for k,t in enumerate(ts):W[k][j]=vals[t]
        check(all(vals[t]==W[tclass[t]][j] for t in range(size)),('label invariance',n,j))
    count=len(ts)
    check(count==len(even),'count equality')
    B=[[F(rt[i]*W[i][j]*W[i][j],size*rg[j]) for j in range(count)] for i in range(count)]
    check(all(sum(row)==1 for row in B),'row normalization')
    check(all(sum(col)==1 for col in zip(*B)),'column normalization')
    for j in range(count):
        for k in range(count):
            check(sum(F(rt[i],size)*W[i][j]*W[i][k] for i in range(count))==(rg[j] if j==k else 0),'weighted orthogonality')
    # Change every chosen graph representative via the reversing permutation,
    # and make independently selected direction choices on its edges.
    rev=tuple(reversed(range(n))); revtr=transport[permutations.index(rev)]
    for j,g in enumerate(even):
        gp=graph_image(g,revtr)
        oriented=[(v,u) if (i%2==0) else (u,v) for i,(u,v) in enumerate(pairs) if gp>>i&1]
        alt={}
        for p in permutations:
            im=0; sg=1
            for u,v in oriented:
                a,b=p[u],p[v]; im|=1<<pos[tuple(sorted((a,b)))]; sg*=1 if a<b else -1
            check(im not in alt or alt[im]==sg,'alternative orientation transporter')
            alt[im]=sg
        for i,t in enumerate(ts):
            v=sum(sg*character(t,h) for h,sg in alt.items())
            check(v*v==W[i][j]**2,'orientation/representative gauge')
    A=[[0]*count for _ in ts]
    for i,t in enumerate(ts):
        for bit in range(N):A[i][tclass[t^(1<<bit)]]+=1
    # Matrix powers and the continuous-time kernel at exp(-2t)=1/3.
    powA=[[int(i==j) for j in range(count)] for i in range(count)]
    transition_checks=0
    for m in range(5):
        for i in range(count):
            for j in range(count):
                expect=F(rt[j],size)*sum(F(W[i][k]*W[j][k],rg[k])*(N-2*g.bit_count())**m for k,g in enumerate(even))
                check(powA[i][j]==expect,('full flip powers',n,m,i,j));transition_checks+=1
        powA=mul(powA,A)
    for i,t in enumerate(ts):
        direct=[0]*count
        for u in range(size):direct[tclass[u]]+=2**(N-(t^u).bit_count())
        for j in range(count):
            expect=F(rt[j],size)*sum(F(W[i][k]*W[j][k],rg[k])*F(1,3)**g.bit_count() for k,g in enumerate(even))
            check(F(direct[j],3**N)==expect,('independent flip distribution',n,i,j));transition_checks+=1
        moment=sum(B[i][k]*g.bit_count() for k,g in enumerate(even))
        check(2*moment==N-A[i][i],'exit rate factor two')
    # Directed paths are found from adjacency predicates, not precomputed bit masks.
    D=[[0]*count for _ in ts]; paths_by_rep=[]; stabilizers=0; swaps=0
    for i,t in enumerate(ts):
        paths=[p for p in permutations if all(sense(t,u,v) for u,v in zip(p,p[1:]))]
        paths_by_rep.append(paths)
        for p in paths:
            path_edges={tuple(sorted(uv)) for uv in zip(p,p[1:])}
            path_mask=sum(1<<pos[e] for e in path_edges)
            u=t^(full^path_mask); j=tclass[u]; D[i][j]+=1
            common={e for b,e in enumerate(pairs) if ((t>>b)&1)==((u>>b)&1)}
            check(common==path_edges,'common arc path')
            common_aut=[q for q,tr in zip(permutations,transport) if tournament_image(t,tr)==t and tournament_image(u,tr)==u]
            check(common_aut==[tuple(range(n))],'pair stabilizer');stabilizers+=1
            if n>=3:
                for tr in transport:
                    check(not (tournament_image(t,tr)==u and tournament_image(u,tr)==t),'endpoint exchange stabilizer')
                swaps+=1
        aut=len(automorphs[i])
        check(all(v%aut==0 for v in D[i]),'aut quotient integrality')
        if n>=3:check((D[i][i]//aut)%2==0,'diagonal quotient parity')
    check(all(F(D[i][j],len(automorphs[i]))==F(D[j][i],len(automorphs[j])) for i in range(count) for j in range(count)),'incidence quotient symmetry')
    DW=mul(D,W)
    C=[[sum(F(rt[i]*W[i][j]*DW[i][k],size*rg[j]) for i in range(count)) for k in range(count)] for j in range(count)]
    check(mul(W,C)==DW,'path invariant operator')
    # Entire labelled path-Walsh matrix independently via direct double sums.
    # n<=4 only; higher dimensions are deliberately not claimed.
    labelled_entries=0
    if 2<=n<=4:
        undirected=[p for p in permutations if p<p[::-1]]
        L=[]
        for t in range(size):
            neighbours=[]
            for p in permutations:
                if all(sense(t,u,v) for u,v in zip(p,p[1:])):
                    pm=sum(1<<pos[tuple(sorted((u,v)))] for u,v in zip(p,p[1:]))
                    neighbours.append(t^(full^pm))
            L.append(neighbours)
        for g in range(size):
            for h in range(size):
                lhs=F(sum(character(t,g)*sum(character(u,h) for u in us) for t,us in enumerate(L)),size)
                rhs=F(0); diff=g^h
                for p in undirected:
                    edge_indices=[pos[tuple(sorted((u,v)))] for u,v in zip(p,p[1:])]
                    pm=sum(1<<e for e in edge_indices)
                    if diff&~pm or diff.bit_count()%2:continue
                    sg=parity((h&~pm).bit_count())
                    for a,b in zip(p,p[1:]):
                        bit=pos[tuple(sorted((a,b)))]
                        if diff>>bit&1:sg*=1 if a<b else -1
                    rhs+=F(sg,2**(n-2))
                check(lhs==rhs,('labelled path Walsh',n,g,h));labelled_entries+=1
    converse=[tclass[full^t] for t in ts]
    check(all(B[i]==B[converse[i]] for i in range(count)),'converse B invariance')
    if n in [5,6]:
        check(all((B[i]==B[j])==(j in {i,converse[i]}) for i in range(count) for j in range(count)),'converse-only fibres')
    rank=rank_exact(B);distinct=len({tuple(r) for r in B})
    degrees=[[sum(row[k] for k,g in enumerate(even) if g.bit_count()==d) for d in range(N+1)] for row in B]
    return_distinct=len({tuple(row) for row in degrees})
    # Static producer files are used only as comparison targets after all
    # independent constructions and checks above are completed.
    expected=json.loads((SOURCE/'results.json').read_text())['orders'][n-1]
    for key,value in [('tournament_masks',ts),('even_graph_masks',even),('Walsh_matrix',W),('coupling',[[str(x) for x in r] for r in B]),('full_flip_adjacency',A)]:
        check(expected[key]==value,('packet mismatch',n,key))
    if n>=2:
        pe=json.loads((SOURCE/'path_results.json').read_text())[n-2]
        check(pe['path_complement_adjacency']==D,('path packet mismatch',n))
        check(pe['Walsh_basis_operator']==[[str(x) for x in row] for row in C],('path basis mismatch',n))
    check(expected['coupling_rank']==rank and expected['coupling_distinct_rows']==distinct and expected['return_spectral_distinct_rows']==return_distinct,'rank/fibre count agreement')
    # Independently verify every supplied Birkhoff term and reconstruct B.
    reconstructed=[[F(0)]*count for _ in ts]
    for term in expected['birkhoff_decomposition']:
        w=F(term['weight']); p=term['permutation']
        check(w>0 and sorted(p)==list(range(count)),'Birkhoff invalid term')
        for i,j in enumerate(p):reconstructed[i][j]+=w
    check(reconstructed==B,'Birkhoff reconstruction')
    out={'n':n,'tournament_types':count,'even_graph_types':len(even),'all_graph_types':len(gs),'even_order_group_graph_types':sum(a%2==0 for a in g_aut),'labelled_invariance_checks':size*count,'full_transition_entries_checked':transition_checks,'path_stabilizer_pairs_checked':stabilizers,'endpoint_swap_exclusions_checked':swaps,'full_labelled_path_Walsh_entries_checked':labelled_entries,'path_counts':[len(p) for p in paths_by_rep],'coupling_rank':rank,'distinct_B_rows':distinct,'distinct_return_rows':return_distinct,'all_checks':'PASS'}
    if n==5:
        i,j,k=[ts.index(t) for t in [2,11,4]]
        check(B[i]==B[j] and A[i][k]==3 and A[j][k]==0 and converse[i]==j,'off-diagonal loss witness')
        out['off_diagonal_loss']={'source_masks':[2,11],'same_B_row':True,'converses':True,'distinct_target_mask':4,'single_flip_counts':[3,0]}
    if n==3:out['path_Walsh_operator']=[[str(x) for x in row] for row in C]
    if n==4:
        check(ts==[0,2,4,5] and even==[0,3,7,30],'n4 masks')
        i,j=ts.index(0),ts.index(4)
        check(B[i]==B[j] and len(paths_by_rep[i])==1 and len(paths_by_rep[j])==5,'n4 loss witness')
        out['n4_B']=[[str(x) for x in row] for row in B];out['n4_A']=A;out['n4_D']=D
        out['n4_path_witnesses']={str(t):[list(p) for p in paths_by_rep[ts.index(t)]] for t in [0,4]}
        for t in [0,4]:
            check(tclass[full^t]==tclass[t],'n4 self converse')
            check(len(automorphs[ts.index(t)])==1,'n4 trivial aut')
        reach=[{v for v in range(n) if v!=u and sense(4,u,v)} for u in range(n)]
        for _ in range(n):reach=[r|set().union(*(reach[v] for v in r)) for r in reach]
        check(all(len(r)==n for r in reach),'n4 strongly connected mask 4')
    return out

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--max-n',type=int,default=6);parser.add_argument('--out',default=str(HERE/'audit_results.json'));args=parser.parse_args()
    rows=[]
    for n in range(1,args.max_n+1):
        r=audit(n); rows.append(r); print(json.dumps(r),flush=True)
    Path(args.out).write_text(json.dumps({'status':'PASS','independent_from_producer':True,'rows':rows},indent=2)+'\n')
