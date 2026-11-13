#!/usr/bin/env python3
"""Exact labelled cyclic-triangle/root-star fibres, small complete universe."""
from itertools import combinations, permutations
from pathlib import Path
import json

def require(condition,message):
    if not condition:raise RuntimeError(message)

def inventory(T):
    n=len(T);dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v,c in enumerate(dp[mask]):
            rem=T[v]&~mask
            while c and rem:
                b=rem&-rem;rem-=b
                dp[mask|b][b.bit_length()-1]+=c
    return [1]+[sum(row)for row in dp[1:]],dp

def state(T,r=0):
    H,E=inventory(T);mask=len(H)-1
    opposite=tuple(mask^(1<<i)^z for i,z in enumerate(T));_,S=inventory(opposite)
    return H[-1],S[-1][r],H[mask^(1<<r)],sum(E[s][r]*H[mask^s]for s in range(mask+1))

def triangle_mask(T,triples):
    out=0
    for i,(a,b,c) in enumerate(triples):
        if ((T[a]>>b&1)==(T[b]>>c&1)==(T[c]>>a&1)):out|=1<<i
    return out

def triangles(T):
    return {s for s in combinations(range(len(T)),3) if triangle_mask(T,[s])}

def birth(T,r):
    n=len(T);out=list(T);out[r]|=1<<n
    out.append(((1<<n)-1)^(1<<r))
    return tuple(out)

def flip(T,e):
    a,b=e;out=list(T);out[a]^=1<<b;out[b]^=1<<a;return tuple(out)

def blowup(Q,blocks):
    offsets=[];n=0
    for B in blocks:offsets.append(n);n+=len(B)
    out=[0]*n
    for i,B in enumerate(blocks):
        extra=sum(((1<<len(C))-1)<<offsets[j] for j,C in enumerate(blocks)if Q[i]>>j&1)
        for a,row in enumerate(B):out[offsets[i]+a]=(row<<offsets[i])|extra
    return tuple(out)

def programs(T,U,ports,depth):
    """All named-port choices to this depth; fresh ports remain available."""
    require(triangles(T)==triangles(U),'program cyclic-triangle mismatch')
    require(all(T[r]==U[r]for r in ports),'program available-port star mismatch')
    require(inventory(T)[0][-1]==inventory(U)[0][-1],'program Hamilton count mismatch')
    if depth==0:return 1
    return 1+sum(programs(birth(T,r),birth(U,r),ports+(len(T),),depth-1)for r in ports)

def main():
    report=[];witness=None;birth_controls=0;program_controls=0
    for n in range(2,6):
        pairs=list(combinations(range(n),2));triples=list(combinations(range(n),3));bank={};ambiguous=0;ready_bank={};hbank={}
        for bits in range(1<<len(pairs)):
            T=[0]*n
            for i,(a,b)in enumerate(pairs):
                if bits>>i&1:T[a]|=1<<b
                else:T[b]|=1<<a
            key=(triangle_mask(T,triples),T[0]);z=state(T)
            hbank[tuple(T)]=z[0]
            if T[0].bit_count()==n-2:
                if key[0]in ready_bank:require(ready_bank[key[0]]==z,'indegree-one response mismatch')
                else:ready_bank[key[0]]=z
            child=birth(T,0)
            expected=triangles(T)|{(0,v,n)for v in range(1,n)if T[v]&1}
            require(triangles(child)==expected,'birth hypergraph formula mismatch')
            birth_controls+=1
            if key in bank and bank[key][0]!=z:
                ambiguous+=1
                if witness is None:witness=dict(n=n,key=key,left=bank[key],right=(z,T))
            else:bank.setdefault(key,(z,T))
            if n<=4 and bank[key][1]!=T:
                program_controls+=programs(tuple(T),tuple(bank[key][1]),(0,),3)
        neutral_pair_controls=0
        for T,h in hbank.items():
            edges=[e for e in pairs if hbank[flip(T,e)]==h]
            for e,f in combinations(edges,2):
                if set(e)&set(f):continue
                require(hbank[flip(flip(T,e),f)]==h,'smaller disjoint-neutral-flip obstruction')
                neutral_pair_controls+=1
        report.append(dict(n=n,graphs=1<<len(pairs),fibres=len(bank),ambiguous=ambiguous,indegree_one_triangle_only_fibres=len(ready_bank),disjoint_neutral_pair_controls=neutral_pair_controls))
    C=(2,4,1);A=blowup((2,0),(C,(0,)));B=blowup((0,1),(C,(0,)))
    U=blowup(C,(A,A,A));V=blowup(C,(B,A,A));r=7
    require(triangles(U)==triangles(V)and U[r]==V[r],'symmetry-change source mismatch')
    require(state(U,r)==state(V,r),'symmetry-change marked response mismatch')
    sym=dict(left=U,right=V,root=r,state=state(U,r),aut_orders_inherited=[81,27])
    T=(22,8,2,5,14);V=(22,8,2,21,6)
    require(T[0]==V[0]and state(T)==(9,3,3,9)and state(V)==(9,4,5,10),'H9 source or state mismatch')
    require(triangles(T)!=triangles(V),'H9 hostile failed to leave triangle fibre')
    NR=[state(birth(birth(W,0),len(W)),len(W))[0]for W in (T,V)]
    require(NR==[27,29],'H9 continuation mismatch')
    hostile=dict(left=T,right=V,left_state=state(T),right_state=state(V),left_triangles=sorted(triangles(T)),right_triangles=sorted(triangles(V)),NR=NR)
    T=(6,16,50,7,41,11);e=(0,1);f=(2,4);graphs=[T,flip(T,e),flip(T,f),flip(flip(T,e),f)]
    counts=[inventory(G)[0][-1]for G in graphs]
    direct=[sum(all(G[w[i]]>>w[i+1]&1 for i in range(5))for w in permutations(range(6)))for G in graphs]
    require(counts==direct==[37,37,37,31],'neutral-trade counting mismatch')
    for G in graphs:
        for v in range(6):
            reached=1<<v
            for _ in range(6):
                for u in range(6):
                    if reached>>u&1:reached|=G[u]
            require(reached==63,'neutral-trade graph is not strongly connected')
    trade=dict(rows=T,edges=[e,f],four_graphs=graphs,H=counts,independent_permutation_H=direct,all_four_strong=True,mixed_difference=-6,minimum_order_in_complete_labelled_universe=6)
    out=dict(status='FINITE-EXACT',universe='all labelled tournaments orders2..5, root0, full cyclic-triangle vertex-set hypergraph and root outgoing star; all allowed programs length<=3 for nontrivial fibres through order4; every vertex-disjoint neutral edge pair through order5',report=report,first_witness=witness,birth_formula_controls=birth_controls,program_node_controls=program_controls,symmetry_change=sym,escape_hostile=hostile,disjoint_neutral_trade=trade)
    Path(__file__).with_name('hypergraph_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
