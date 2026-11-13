"""Exact small controls for occupancy powers, tail-marked cofactors and growth.

All labelled tournaments through order5; all loopless digraphs through order3.
No external packages, floating point, or assertions disabled by python -O.
"""
from functools import lru_cache
from itertools import permutations,product
from math import factorial
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def need(p,message):
    if not p:raise RuntimeError(message)

def plus(a,b,scale=1):
    out=dict(a)
    for k,v in b.items():
        out[k]=out.get(k,0)+scale*v
        if not out[k]:del out[k]
    return out

def times(a,b):
    out={}
    for x,u in a.items():
        for y,v in b.items():
            if not x&y:out[x|y]=out.get(x|y,0)+u*v
    return {k:v for k,v in out.items() if v}

def determinant(M):
    n=len(M);out={}
    for p in permutations(range(n)):
        term={0:(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))}
        for i,j in enumerate(p):term=times(term,M[i][j])
        out=plus(out,term)
    return out

def numeric_det(M):
    return sum((-1)**sum(p[i]>p[j] for i in range(len(M)) for j in range(i+1,len(M)))
               *prod(M[i][j] for i,j in enumerate(p)) for p in permutations(range(len(M))))

def prod(items):
    out=1
    for a in items:out*=a
    return out

@lru_cache(None)
def inventories(A):
    n=len(A);end=[[0]*n for _ in range(1<<n)];start=[[0]*n for _ in range(1<<n)]
    for v in range(n):end[1<<v][v]=start[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask>>v)&1:continue
            old=mask^(1<<v)
            if old:
                end[mask][v]=sum(end[old][u]*A[u][v] for u in range(n))
                start[mask][v]=sum(start[old][u]*A[v][u] for u in range(n))
    P={0:1};P.update((mask,sum(row)) for mask,row in enumerate(end) if sum(row))
    E=[{mask:row[v] for mask,row in enumerate(end) if row[v]} for v in range(n)]
    S=[{mask:row[v] for mask,row in enumerate(start) if row[v]} for v in range(n)]
    return P,E,S

def occupancy_power(A):
    n=len(A);layer=[{1<<v:1} for v in range(n)];P={0:1}
    for _ in range(n):
        for row in layer:P=plus(P,row)
        nxt=[{} for _ in range(n)]
        for u in range(n):
            for v in range(n):
                if A[u][v]:nxt[v]=plus(nxt[v],times(layer[u],{1<<v:A[u][v]}))
        layer=nxt
    need(not any(layer),'occupancy power terminates')
    return P

def brute_paths(A):
    n=len(A);H=[[0]*n for _ in range(n)]
    for p in permutations(range(n)):
        if all(A[u][v] for u,v in zip(p,p[1:])):H[p[0]][p[-1]]+=1
    return H

def out_trees(A,r):
    n=len(A);vertices=[v for v in range(n) if v!=r]
    options=[[u for u in range(n) if A[u][v]] for v in vertices]
    total=0;sf={};leaves={}
    for choice in product(*options):
        parent=dict(zip(vertices,choice));good=True
        for v in vertices:
            seen=set()
            while v!=r and v not in seen:seen.add(v);v=parent[v]
            if v!=r:good=False;break
        if not good:continue
        total+=1;deg=[0]*n
        for u in choice:deg[u]+=1
        ell=sum(d==0 for d in deg);leaves[ell]=leaves.get(ell,0)+1
        need(ell==1+sum(max(0,d-1) for d in deg),'tree leaf/excess law')
        if max(deg,default=0)<=1:
            mask=sum((d==1)<<u for u,d in enumerate(deg));sf[mask]=sf.get(mask,0)+1
    return total,sf,leaves

def cofactor(A,r,tail=True,values=None):
    n=len(A)
    if values is not None:
        B=[[A[u][v]*values[u if tail else v] for v in range(n)] for u in range(n)]
        L=[[(sum(B[u][j] for u in range(n)) if i==j else 0)-B[i][j] for j in range(n)] for i in range(n)]
        return numeric_det([[L[i][j] for j in range(n) if j!=r] for i in range(n) if i!=r])
    L=[[{} for _ in range(n)] for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if A[u][v]:
                z={1<<(u if tail else v):A[u][v]}
                L[v][v]=plus(L[v][v],z);L[u][v]=plus(L[u][v],z,-1)
    return determinant([[L[i][j] for j in range(n) if j!=r] for i in range(n) if i!=r])

def path_determinants(A):
    n=len(A)
    def entry(i,j,c):return plus({0:1} if i==j else {},{1<<j:c} if c else {})
    minus=determinant([[entry(i,j,-A[i][j]) for j in range(n)] for i in range(n)])
    plusdet=determinant([[entry(i,j,A[i][j]) for j in range(n)] for i in range(n)])
    complement=determinant([[entry(i,j,1-A[i][j]) for j in range(n)] for i in range(n)])
    return minus,plusdet,complement

def walk_total(A,k):
    v=[1]*len(A)
    for _ in range(k):v=[sum(v[u]*A[u][w] for u in range(len(A))) for w in range(len(A))]
    return sum(v)

def all_graphs(n,tournament):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)] if tournament else [(i,j) for i in range(n) for j in range(n) if i!=j]
    for bits in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            b=(bits>>k)&1;A[i][j]=b
            if tournament:A[j][i]=1-b
        yield tuple(map(tuple,A))

def test(A,tournament,totals):
    n=len(A);full=(1<<n)-1;P,E,S=inventories(A);ends=brute_paths(A)
    need(sum(map(sum,ends))==P.get(full,0),'permutations versus subset DP')
    need(occupancy_power(A)==P,'matrix occupancy powers')
    minus,plusdet,complement=path_determinants(A)
    need(times(P,minus)==complement,'general path determinant identity')
    if tournament:
        all_subsets={mask:1 for mask in range(1<<n)}
        need(times(all_subsets,plusdet)==complement,'tournament factorization')
        need(all(P.get(mask,0)%2==1 for mask in range(1<<n)),'all induced tournament counts odd')
        C=[[times(E[u],S[v]).get(full,0) for v in range(n)] for u in range(n)]
        lam=[E[u].get(full,0)-S[u].get(full,0)+sum(C[u]) for u in range(n)]
        need(all(C[u][u]==0 for u in range(n)),'two paths cannot share a port')
        need(all(a%2==0 for a in lam),'birth linear coefficient is even')
        need(all((C[u][v]+C[v][u])%2==0 for u in range(n) for v in range(u+1,n)),
             'birth mixed coefficient is even')
        need(all(lam[u]==sum((C[u][v]+C[v][u])*A[v][u] for v in range(n))
                 for u in range(n)), 'inherited exact incoming-cut field')
        need(sum(lam)==sum(C[u][v]+C[v][u] for u in range(n) for v in range(u+1,n)),
             'source and sink birth both preserve H')
        totals['birth_coefficient_checks']+=1
    for r in range(n):
        arbs,sf,leaves=out_trees(A,r)
        tail=cofactor(A,r,True);head=cofactor(A,r,False)
        expected={full^(1<<v):ends[r][v] for v in range(n) if ends[r][v]}
        need(tail==sf==expected,'tail-weighted MatrixTree endpoint projection')
        need(head==({full^(1<<r):arbs} if arbs else {}),'head marking retains all trees')
        need(cofactor(A,r,values=[1]*n)==arbs,'ordinary MatrixTree cofactor')
        need(leaves.get(1,0)==sum(ends[r]),'one-leaf coefficient')
        if n<=3:
            for t in range(n):
                rest=full^(1<<t);value=0;sub=rest
                while True:
                    value+=(-1)**(n-1-sub.bit_count())*cofactor(A,r,values=[(sub>>v)&1 for v in range(n)])
                    if not sub:break
                    sub=(sub-1)&rest
                need(value==ends[r][t],'Boolean Mobius determinant extraction')
                totals['mobius_endpoint_checks']+=1
        totals['rooted_tree_checks']+=1
    if n>1:
        old=tuple(tuple(row[:-1]) for row in A[:-1]);P0,E0,S0=inventories(old)
        left={0:1};right={0:1}
        for v in range(n-1):
            left=plus(left,E0[v],A[v][-1]);right=plus(right,S0[v],A[-1][v])
        expected=dict(P0)
        for mask,value in times(left,right).items():expected[mask|(1<<(n-1))]=value
        need(expected==P,'full new-vertex insertion kernel')
        totals['growth_checks']+=1
    if n<=3 or tournament and n<=4:
        for u in range(n):
            for v in range(n):
                if u==v or tournament and not A[u][v]:continue
                B=[list(row) for row in A]
                if tournament:
                    B[u][v]=0;B[v][u]=1
                    expected=plus(plus(P,times(E[u],S[v]),-1),times(E[v],S[u]))
                else:
                    delta=1-2*A[u][v];B[u][v]+=delta
                    expected=plus(P,times(E[u],S[v]),delta)
                need(inventories(tuple(map(tuple,B)))[0]==expected,'exact edge change kernel')
                totals['edge_change_checks']+=1

def main():
    totals=dict(tournaments=0,general_digraphs=0,rooted_tree_checks=0,growth_checks=0,
                edge_change_checks=0,mobius_endpoint_checks=0,birth_coefficient_checks=0)
    for n in range(1,6):
        for A in all_graphs(n,True):test(A,True,totals);totals['tournaments']+=1
    for n in range(1,4):
        for A in all_graphs(n,False):test(A,False,totals);totals['general_digraphs']+=1
    named={
        'TT4':tuple(tuple(int(i<j) for j in range(4)) for i in range(4)),
        'C3_then_sink':((0,1,0,1),(0,0,1,1),(1,0,0,1),(0,0,0,0)),
        'C3_plus_two_isolates':((0,1,0,0,0),(0,0,1,0,0),(1,0,0,0,0),(0,0,0,0,0),(0,0,0,0,0)),
        'out_star3':((0,1,1),(0,0,0),(0,0,0)),
        'directed_chain3':((0,1,0),(0,0,1),(0,0,0)),
        'bidirected_edge':((0,1),(1,0)),
        'empty2':((0,0),(0,0)),
        'TT3':((0,1,1),(0,0,1),(0,0,0)),
        'known_one_missing_pair_H7':tuple(tuple(map(int,r)) for r in
            ('00001','10000','01000','11100','01110')),
        'known_one_missing_pair_H21':tuple(tuple(map(int,r)) for r in
            ('000011','100000','110000','111000','011100','001110')),
    }
    examples={}
    for name,A in named.items():
        n=len(A);P,_,_=inventories(A)
        examples[name]=dict(A=A,H=P.get((1<<n)-1,0),walks_n_minus_1=walk_total(A,n-1),
            walks_n=walk_total(A,n),out_arborescences=sum(out_trees(A,r)[0] for r in range(n)),
            root0_tail_cofactor_squarefree=cofactor(A,0,True),root0_head_cofactor=cofactor(A,0,False))
    need(examples['TT4']['H']==1 and examples['TT4']['walks_n']==0,'nilpotent hostile')
    need(examples['C3_then_sink']['H']==3 and examples['C3_then_sink']['walks_n_minus_1']==6,'repeated-walk hostile')
    need(examples['out_star3']['H']==0 and examples['out_star3']['out_arborescences']==1,'branching hostile')
    need(examples['TT3']['H']==1 and examples['TT3']['out_arborescences']==2,'tree parity hostile')
    defects=[]
    for name,target,pair,completion_counts in [
        ('known_one_missing_pair_H7',7,(0,2),(9,9)),
        ('known_one_missing_pair_H21',21,(1,5),(41,23))]:
        A=named[name];n=len(A);full=(1<<n)-1;P,E,S=inventories(A)
        absent=[(u,v) for u in range(n) for v in range(u+1,n) if not A[u][v] and not A[v][u]]
        need(absent==[pair] and all(not (A[u][v] and A[v][u]) for u in range(n) for v in range(u+1,n)),
             'exactly one missing pair and no reciprocal pair')
        need(P[full]==target==sum(map(sum,brute_paths(A))),'known one-defect literal path count')
        repairs=[]
        for (u,v),expected in zip((pair,pair[::-1]),completion_counts):
            B=[list(row) for row in A];B[u][v]=1;B=tuple(map(tuple,B))
            h=inventories(B)[0][full];delta=times(E[u],S[v]).get(full,0)
            need(h==expected==sum(map(sum,brute_paths(B))) and h==target+delta,
                 'known one-defect completion and exact arc kernel')
            repairs.append(dict(added_arc=[u,v],H=h,added_paths=delta))
        defects.append(dict(name=name,H=target,missing_pair=pair,completions=repairs))
    result=dict(status='PASS',universe='All labelled tournaments n1..5 and every loopless directed graph n1..3, plus named controls; loops excluded, digons allowed in general class',
                counts=totals,examples=examples,known_one_defect_controls=defects,
                known_one_defect_source='05-knowledge/results/kps_s18_ternary_chains_digraph_weights_20261109/out/one_defect_witness.out (THM-5071 (4c), already known)',
                independent_paths=['permutation endpoint counts','subset dynamic program','square-zero matrix powers',
                    'square-zero determinant expansion','parent-function spanning-tree enumeration','numeric cofactors and Boolean Mobius extraction'])
    (HERE/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',counts=totals,examples={k:{a:v[a] for a in ('H','walks_n_minus_1','out_arborescences')} for k,v in examples.items()}),indent=2))

if __name__=='__main__':main()
