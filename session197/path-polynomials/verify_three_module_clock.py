#!/usr/bin/env python3
from verify import *
from math import isqrt

Q=(30,64,10,66,174,335,21,367,95)
TAU=(0,7,8,5,6,3,4,1,2)

def substitute_blocks(blocks):
    starts=[];n=0
    for block in blocks:starts.append(n);n+=len(block)
    out=[0]*n
    for i,block in enumerate(blocks):
        for k,row in enumerate(block):
            u=starts[i]+k;out[u]=row<<starts[i]
            for j,other in enumerate(blocks):
                if Q[i]>>j&1:out[u]|=((1<<len(other))-1)<<starts[j]
    validate(out);return tuple(out),starts

def graph(T,U,V):
    blocks=[(0,)]*9;blocks[0]=T;blocks[2]=U;blocks[8]=V
    return substitute_blocks(blocks)

def strong(a):
    for root in range(len(a)):
        seen=1<<root
        while True:
            nxt=seen
            for v in range(len(a)):
                if seen>>v&1:nxt|=a[v]
            if nxt==seen:break
            seen=nxt
        if seen!=(1<<len(a))-1:return False
    return True

def run():
    require(strong(Q),'source strong')
    require(all((Q[u]>>v&1)==(Q[TAU[v]]>>TAU[u]&1) for u in range(9) for v in range(9)),'source anti')
    grid=[]
    for m in (1,2,3):
        row=[]
        for s,t in ((1,1),(1,2),(2,1),(2,2)):
            a,starts=graph(transitive(m),transitive(s),transitive(t));h=h_subsets(a)[-1]
            require(h==18*3**m*2**s*2**t+1,'exact grid')
            row.append(h)
            # Check the uniform cross-cover descriptions at one representative clone.
            for v,cover in ((0,{starts[4],starts[6]}),(2,{starts[6]}),(8,{starts[4]})):
                r=starts[v];edges=[(u,w) for u in range(len(a)) for w in range(len(a)) if (a[r]>>u&1) and (a[w]>>r&1) and (a[u]>>w&1)]
                require(all(u in cover or w in cover for u,w in edges),('cross cover',m,s,t,v,edges))
        grid.append(row)
    require(grid==[[217,433,433,865],[649,1297,1297,2593],[1945,3889,3889,7777]],'frozen grid')
    inputs=[(0,),transitive(2),transitive(3),(2,4,1)]
    ps={a:bin_polynomial(a) for a in inputs};checks=0
    for T,U,V in product(inputs,repeat=3):
        a,_=graph(T,U,V);h=h_subsets(a)[-1]
        expected=18*value(ps[T],3)*value(ps[U],2)*value(ps[V],2)+value(ps[T],1)*value(ps[U],1)*value(ps[V],1)
        require(h==expected,('arbitrary input identity',T,U,V,h,expected));checks+=1
    anti_controls=[]
    for aexp,bexp in ((3,3),(3,4),(5,3),(5,5),(7,3)):
        m=bexp-2;s=(aexp-1)//2;blocks=[(0,)]*9;blocks[0]=transitive(m);blocks[2]=blocks[8]=transitive(s)
        a,starts=substitute_blocks(blocks);tau=[]
        for v,block in enumerate(blocks):
            tau.extend(starts[TAU[v]]+len(block)-1-k for k in range(len(block)))
        require(all((a[u]>>v&1)==(a[tau[v]]>>tau[u]&1) for u in range(len(a)) for v in range(len(a))),'lifted anti')
        require(strong(a),'strong module graph');h=h_subsets(a)[-1]
        require(h==2**aexp*3**bexp+1,'symmetric count')
        anti_controls.append({'a':aexp,'b':bexp,'n':len(a),'H':h})
    witness,_=graph(transitive(5),transitive(2),(0,));h=h_subsets(witness)[-1]
    require(len(witness)==14 and h==34993==7*4999 and strong(witness),'7p physical witness')
    primes=[p for p in range(2,isqrt(4999)+1) if all(p%d for d in range(2,isqrt(p)+1))]
    remainders=[(p,4999%p) for p in primes];require(all(r for p,r in remainders),'prime certificate')
    require(6**7+1==7**2*29*197,'composite diagonal hostile')
    print(json.dumps({'physical_interpolation_grid':grid,'arbitrary_input_triples':checks,'self_converse_controls':anti_controls,'seven_prime_witness':{'a':4,'b':7,'n':14,'rows':witness,'H':h,'prime_quotient':4999,'trial_division_remainders':remainders},'composite_diagonal':[6**7+1,7,7,29,197]},sort_keys=True,indent=2))
if __name__=='__main__':run()
