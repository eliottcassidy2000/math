#!/usr/bin/env python3
"""Independent finite controls for fixed-increment, sparse-backarc and run-entropy arguments.
Standard library; no producer imports; failure checks remain active under -O.
"""
from fractions import Fraction
from itertools import combinations,product,permutations
from math import comb,factorial
import json

def require(c,msg):
    if not c:raise RuntimeError(msg)

def graph(n,back):
    back=set(back);out=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            a,b=(j,i) if (i,j) in back else (i,j);out[a]|=1<<b
    return out

def hp_subsets(out):
    n=len(out);size=1<<n;d=[[0]*n for _ in range(size)];H=[0]*size;H[0]=1
    for i in range(n):d[1<<i][i]=1
    for s in range(1,size):
        H[s]=sum(d[s])
        for i in range(n):
            v=d[s][i]
            if not v:continue
            todo=out[i]&~s
            while todo:
                bit=todo&-todo;todo-=bit;j=bit.bit_length()-1;d[s|bit][j]+=v
    return H

def components(out):
    n=len(out);reach=[out[i]|(1<<i) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    unused=(1<<n)-1;C=[]
    while unused:
        bit=unused&-unused;i=bit.bit_length()-1;s=0
        for j in range(n):
            if reach[i]>>j&1 and reach[j]>>i&1:s|=1<<j
        C.append(s);unused&=~s
    return C

def increment_audit():
    cases=positive=0;pairs=set();largest=0;pair_controls=0
    for n in range(1,6):
        arcs=list(combinations(range(n),2));allmask=(1<<n)-1
        for code in range(1<<len(arcs)):
            out=graph(n,[p for i,p in enumerate(arcs) if code>>i&1]);H=hp_subsets(out);cc=components(out)
            core=0
            for component in cc:
                if component.bit_count()>1:core|=component
            for sub in range(allmask+1):
                require(H[sub&core]==H[sub],'pair kernel lost smaller H')
                require(H[core]==H[allmask],'pair kernel lost larger H')
                gain=H[allmask]-H[sub]
                if gain>0:
                    require(core.bit_count()<=gain+2,'whole induced-pair core bound')
                    require((sub&core)!=core,'positive pair kernel not proper')
                    require(2*H[sub]<=factorial(gain+1),'induced refined factorial cap')
                    pair_controls+=1
            for x in range(n):
                before=H[allmask^(1<<x)];after=H[allmask];delta=after-before;cases+=1
                require(delta>=0 and delta%2==0,'deletion monotonicity/parity')
                if not delta:continue
                positive+=1;K=delta//2
                C=next(c for c in cc if c>>x&1);a=1
                for c in cc:
                    if c!=C:a*=H[c]
                require(delta==a*(H[C]-H[C^(1<<x)]),'SCC product difference')
                require(K%a==0 and a<=K,'external multiplier divisibility')
                k=K//a;m=C.bit_count();largest=max(largest,m)
                require(m<=2*k+2,'odd-cycle core bound')
                require(after<=K*factorial(2*K+2),'global factorial cap')
                require(2*before<=factorial(2*K+1),'refined old-endpoint factorial cap')
                if K==1:pairs.add((before,after))
    require(pairs=={(1,3),(3,5)},'unit increment classification')
    return {'strict_induced_pair_kernel_controls':pair_controls,'all_vertex_deletions':cases,'positive_increment_cases':positive,'plus_two_pairs':sorted(pairs),'maximum_core_order_tested':largest,'universe':'all labelled tournaments through order5, each vertex'}

def arc_audit():
    graphs=spans=nontrivial=0
    for n in range(2,8):
        arcs=list(combinations(range(n),2));allmask=(1<<n)-1
        for b in range(min(3,len(arcs))+1):
            for back in combinations(arcs,b):
                out=graph(n,back);H=hp_subsets(out);graphs+=1
                for l,r in back:
                    removed=set()
                    for a,z in back:
                        if (a,z)==(l,r) or not(l<=a<z<=r):continue
                        candidates=[v for v in (a,z) if v not in (l,r)]
                        require(candidates,'duplicate chosen arc');removed.add(candidates[0])
                    vertices=[v for v in range(l,r+1) if v not in removed];m=len(vertices)
                    require(m>=r-l-b+2,'vertex deletion invoice')
                    for i,j in combinations(vertices,2):
                        require(bool(out[j]>>i&1)==((i,j)==(l,r)),'unexpected remaining backarc')
                    mask=sum(1<<v for v in vertices)
                    require(H[mask]<=H[allmask],'induced path monotonicity')
                    if m>=3:require(H[mask]==(1<<(m-2))+1,'single long arc formula')
                    else:require(H[mask]==1,'two-vertex boundary')
                    d=r-l
                    if d>=b+1:
                        nontrivial+=1;require(H[allmask]>=(1<<(d-b))+1,'span inequality')
                    spans+=1
                cc=components(out)
                for C in cc:
                    if C.bit_count()<2:continue
                    V=[v for v in range(n) if C>>v&1];bc=[(i,j) for i,j in back if i in V and j in V]
                    # Prefix cut coverage is asserted within the induced presentation of a strong SCC.
                    for cut in range(1,len(V)):
                        L=set(V[:cut]);R=set(V[cut:]);require(any(i in L and j in R for i,j in bc),'strong cut uncovered')
    return {'sparse_presentations':graphs,'chosen_backarc_controls':spans,'nontrivial_log_span_controls':nontrivial,'universe':'all presentations through order7 with at most3 backward arcs'}

def cycle_cover_audit():
    pairs={1:set(),2:set()};presentations=deletions=0
    def literal(out,vertices):
        return sum(all(out[a]>>b&1 for a,b in zip(w,w[1:])) for w in permutations(vertices))
    for n in range(3,7):
        forced={(i,i+1) for i in range(n-1)}|{(n-1,0)}
        free=[(i,j) for i,j in combinations(range(n),2) if (i,j) not in forced and (j,i) not in forced]
        for code in range(1<<len(free)):
            out=[0]*n
            for a,b in forced:out[a]|=1<<b
            for k,(i,j) in enumerate(free):
                a,b=(i,j) if code>>k&1 else (j,i);out[a]|=1<<b
            F=literal(out,range(n));presentations+=1
            for x in range(n):
                h=literal(out,[v for v in range(n) if v!=x]);d=F-h;deletions+=1
                require(d>0 and d%2==0,'strong cycle-cover deletion gain')
                if d//2 in pairs:pairs[d//2].add((h,F))
    require(pairs[1]=={(1,3),(3,5)},'literal cycle-cover +2 classification')
    require(pairs[2]=={(1,5),(5,9)},'literal cycle-cover +4 classification')
    return {'strong_fixed_Hamilton_cycle_presentations':presentations,'literal_vertex_deletions':deletions,'plus_two_pairs':sorted(pairs[1]),'plus_four_pairs':sorted(pairs[2]),'scope':'All strong classes orders3..6 covered via Hamilton-cycle relabelling; literal permutations only.'}

def changes(w):return sum(a!=b for a,b in zip(w,w[1:]))
def run_audit():
    subsequences=0
    for n in range(7):
        for w in product(range(4),repeat=n):
            base=changes(w)
            for mask in range(1<<n):
                sub=tuple(w[i] for i in range(n) if mask>>i&1)
                require(changes(sub)<=base,'subsequence run increase')
                for projection in (0,1,2,3):
                    require(changes(tuple(c&projection for c in sub))<=changes(sub),'apex projection run increase')
                subsequences+=1
    counts=0
    for q in (1,2,4):
        for N in range(6):
            for R in range(N+1):
                actual=sum(changes(w)<=R for n in range(N+1) for w in product(range(q),repeat=n))
                expected=1+q*sum(comb(N,r+1)*(q-1)**r for r in range(min(R,N-1)+1))
                require(actual==expected,'hockey-stick word count');counts+=1
    # exp(7/2)>30 follows from a finite lower Taylor sum, so (log30+1)/5<9/10.
    lower=sum(Fraction(7,2)**k/factorial(k) for k in range(11))
    require(lower>30,'explicit entropy constant')
    return {'subsequence_cases_with_all4_bit_projections':subsequences,'exact_word_count_cases':counts,'exp_3_5_lower_bound':str(lower),'constant_conclusion':'(ln(30)+1)/5 < 9/10','scope':'finite deletion/projection/count controls and exact constant bound; asymptotic theorem requires stated core-size proof'}

def main():
    print(json.dumps({'status':'PASS','fixed_increment':increment_audit(),'strong_literal_cycle_cover':cycle_cover_audit(),'fixed_backarc':arc_audit(),'run_entropy':run_audit()},indent=2,sort_keys=True))
if __name__=='__main__':main()
