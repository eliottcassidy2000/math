#!/usr/bin/env python3
"""Second dyadic layer: half-Laplacian pairing on the blue harmonic kernel."""
import json,itertools as it
from pathlib import Path
import sympy as sp
from check_dyadic_network import lap,need,rank2,components,tree,val2
from check_dyadic_lift import dyadic_smith
HERE=Path(__file__).resolve().parent

def nullbasis2(M):
    m=len(M);n=len(M[0]) if m else 0
    A=[[int(v)&1 for v in row] for row in M];pivs=[];k=0
    for j in range(n):
        i=next((i for i in range(k,m) if A[i][j]),None)
        if i is None:continue
        A[k],A[i]=A[i],A[k]
        for i in range(m):
            if i!=k and A[i][j]:A[i]=[x^y for x,y in zip(A[i],A[k])]
        pivs.append(j);k+=1
    basis=[]
    for j in range(n):
        if j in pivs:continue
        v=[0]*n;v[j]=1
        for i,p in enumerate(pivs):v[p]=A[i][j]
        basis.append(v)
    return basis

def pairing(Lb,Lk):
    A=Lb+2*Lk;vecs=nullbasis2(Lb.tolist());V=sp.Matrix(len(vecs),A.rows,lambda i,j:vecs[i][j])
    raw=V*A*V.T;need(all(int(v)%2==0 for v in raw))
    beta=raw.applyfunc(lambda v:(int(v)//2)%2)
    B=beta.tolist();dim=len(vecs);r=rank2(B)
    need(beta==beta.T)
    # Separate blue half-energy and black incidence Gram contribution.
    rawblue=V*Lb*V.T;need(all(int(v)%2==0 for v in rawblue))
    Bblue=rawblue.applyfunc(lambda v:(int(v)//2)%2)
    Bblack=(V*Lk*V.T).applyfunc(lambda v:int(v)%2)
    need((Bblue+Bblack).applyfunc(lambda v:int(v)%2)==beta)
    vals=dyadic_smith(A.tolist());need(dim-r==sum(v>=2 for v in vals))
    return dict(kernel_dimension=dim,pairing_rank=r,pairing_radical_dimension=dim-r,blue_part_rank=rank2(Bblue.tolist()),black_part_rank=rank2(Bblack.tolist()),blue_part_alternating=not any(Bblue[i,i] for i in range(dim)),full_pairing_alternating=not any(beta[i,i] for i in range(dim)))

def rewire_witness(row):
    # Same exact actual n=5 node/color degrees and black graph, alter two blue edges.
    n=row['nodes'];B=row['blue_edges'];K=row['paired_black_edges'];actual=lap(n,B)+2*lap(n,K)
    actualvals=dyadic_smith(actual[:-1,:-1].tolist());actualrank=sum(v>0 for v in actualvals)
    for e,f in it.combinations(range(len(B)),2):
        a,b,wa=B[e];c,d,wb=B[f]
        if len({a,b,c,d})<4:continue
        for swapped in [[(a,c,1),(b,d,1)],[(a,d,1),(b,c,1)]]:
            new=[(x,y,w-(i==e)-(i==f)) for i,(x,y,w) in enumerate(B) if w-(i==e)-(i==f)>0]+swapped
            if len(components(n,new+K))!=1:continue
            A=lap(n,new)+2*lap(n,K);vals=dyadic_smith(A[:-1,:-1].tolist());rank=sum(v>0 for v in vals)
            if rank!=actualrank:
                need([sum(w*(int(i==v)+int(j==v)) for i,j,w in B) for v in range(n)]==[sum(w*(int(i==v)+int(j==v)) for i,j,w in new) for v in range(n)])
                return dict(removed=[B[e],B[f]],added=swapped,actual_two_rank=actualrank,rewired_two_rank=rank,actual_tree_count=tree(actual),rewired_tree_count=tree(A),rewired_blue_edges=new,black_graph_unchanged=True)
    return None


def rewire_black_witness(row):
    # Keep the whole actual blue graph, and preserve every paired-black degree.
    n=row['nodes'];B=row['blue_edges'];K=row['paired_black_edges']
    A0=lap(n,B)+2*lap(n,K);v0=dyadic_smith(A0[:-1,:-1].tolist())
    for e,f in it.combinations(range(len(K)),2):
        a,b,wa=K[e];c,d,wb=K[f]
        if len({a,b,c,d})<4:continue
        for swapped in [[(a,c,1),(b,d,1)],[(a,d,1),(b,c,1)]]:
            new=[(x,y,w-(i==e)-(i==f)) for i,(x,y,w) in enumerate(K) if w-(i==e)-(i==f)>0]+swapped
            if len(components(n,B+new))!=1:continue
            A=lap(n,B)+2*lap(n,new);vals=dyadic_smith(A[:-1,:-1].tolist())
            if sum(v>=2 for v in vals)!=sum(v>=2 for v in v0):
                need(sum(v>0 for v in vals)==sum(v>0 for v in v0))
                need([sum(w*(int(i==v)+int(j==v)) for i,j,w in K) for v in range(n)]==[sum(w*(int(i==v)+int(j==v)) for i,j,w in new) for v in range(n)])
                basis=nullbasis2(lap(n,B)[:-1,:-1].tolist());V=sp.Matrix(basis)
                Q0=(V*A0[:-1,:-1]*V.T).applyfunc(lambda z:int(z)//2%2)
                Q1=(V*A[:-1,:-1]*V.T).applyfunc(lambda z:int(z)//2%2)
                need(list(Q0.diagonal())==list(Q1.diagonal()))
                return dict(actual_half_pairing=[[int(v) for v in row] for row in Q0.tolist()],rewired_half_pairing=[[int(v) for v in row] for row in Q1.tolist()],removed_paired_black_edges=[K[e],K[f]],added_paired_black_edges=swapped,blue_graph_unchanged=True,all_paired_black_degrees_unchanged=True,actual_two_primary_exponents=[v for v in v0 if v],rewired_two_primary_exponents=[v for v in vals if v],actual_tree_count=tree(A0),rewired_tree_count=tree(A),rewired_paired_black_edges=new,rewired_half_form=pairing(lap(n,B)[:-1,:-1],lap(n,new)[:-1,:-1]))
    return None

if __name__=='__main__':
 rows=json.loads((HERE/'dyadic_network_receipt.json').read_text());out=[]
 for row in rows:
    r=row['nodes'];Lb=lap(r,row['blue_edges']);Lk=lap(r,row['paired_black_edges']);ans=pairing(Lb[:-1,:-1],Lk[:-1,:-1]);ans['n']=row['n'];out.append(ans);print(ans,flush=True)
 witness=rewire_witness(rows[2]);need(witness is not None);print('BLUE HOSTILE',witness,flush=True)
 black_witness=rewire_black_witness(rows[2]);need(black_witness is not None);print('BLACK HOSTILE',black_witness,flush=True)
 (HERE/'bockstein_receipt.json').write_text(json.dumps(dict(rows=out,degree_preserving_rewiring=witness,blue_fixed_black_rewiring=black_witness),indent=2)+'\n')
 print('ALL CHECKS PASSED')
