"""Visited-vertex transfer, ordinary walks, and oriented-graph hostile controls."""
from itertools import permutations
import numpy as np
from typed_glue import need,save,decode
from context_head import ordinary_bank
from inflation_probe import end_table


def hp(T):return int(end_table(np.array(T,dtype=np.int64))[-1].sum())


def brute(T):
    return sum(all(T[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(len(T))))


def squarefree(T):
    # Multilinear coefficient propagation: a repeated vertex kills the monomial.
    n=len(T);state={(1<<v,v):1 for v in range(n)}
    for _ in range(n-1):
        nxt={}
        for (S,v),c in state.items():
            for w in range(n):
                if T[v]>>w&1 and not(S>>w&1):
                    key=S|1<<w,w;nxt[key]=nxt.get(key,0)+c
        state=nxt
    return sum(state.values())


def run():
    controls=0;digraphs={};n=4;arcs=[(i,j) for i in range(n) for j in range(n) if i!=j]
    for mask in range(1<<len(arcs)):
        T=[0]*n
        for k,(a,b) in enumerate(arcs):
            if mask>>k&1:T[a]|=1<<b
        h=brute(T);need(h==hp(T)==squarefree(T),'three independent path counts');controls+=1
        if h in (0,2,7):digraphs.setdefault(h,dict(rows=T,H=h))
    oriented={};raw_walks={}
    for (n,code),T in sorted(ordinary_bank(6).items()):
        if n<3:continue
        A=np.array([[(T[i]>>j)&1 for j in range(n)] for i in range(n)],dtype=object)
        power=np.linalg.matrix_power(A,n);walks=int(power.sum());h=hp(T)
        if walks in (0,7,21) or walks%2==0:
            raw_walks.setdefault(str(walks),dict(rows=T,n=n,H=h,length_n_walks=walks))
        for a in range(n):
            for b in range(n):
                if not(T[a]>>b&1):continue
                U=list(T);U[a]^=1<<b;v=hp(U)
                if v in (7,21) and v not in oriented:
                    need(v==brute(U)==squarefree(U),'oriented forbidden-count witness')
                    oriented[v]=dict(rows=U,H=v,completion_rows=T,completion_H=h,deleted_arc=[a,b])
    need(set(oriented)=={7,21},'both tournament gaps occur after one arc deletion')
    save('walk_controls',dict(status='FINITE-EXACT',scope='All4096 loopless digraphs on4 labelled vertices, three independent Hamilton-path counts. One-arc deletions from every ordinary tournament through6; witnesses, not an arbitrary oriented-graph census.',digraph_controls=controls,general_digraph_examples=digraphs,oriented_gap_witnesses=oriented,ordinary_walk_examples=raw_walks))
    print('DIGRAPH CONTROLS',controls,'ORIENTED',oriented,'RAW WALK HEAD',list(raw_walks.items())[:8])


if __name__=='__main__':run()
