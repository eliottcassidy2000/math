#!/usr/bin/env python3
"""Complete classified odd mixed-prime symmetry sectors at n=11,12.
Uses proved affine orbit-image normalizers; no subgroup database.
The copied full-isomorphism/sign engine is pinned in ENGINE_PROVENANCE.md.
"""
from itertools import permutations,product
from collections import Counter,defaultdict
from functools import lru_cache
import json,sys,time
import fixed_cube_engine as m

def need(x,msg):m.need(x,msg)

@lru_cache(None)
def affine(p):
    return [tuple((a*i+b)%p for i in range(p)) for a in range(1,p) for b in range(p)]

@lru_cache(None)
def normalizer_cached(n,H,gens):
    remaining=set(range(n));blocks=[];fixed=[]
    while remaining:
        u=min(remaining);B=tuple(sorted({h[u] for h in H}));remaining-=set(B)
        (fixed if len(B)==1 else blocks).append(B[0] if len(B)==1 else B)
    sizegroups=defaultdict(list)
    for i,B in enumerate(blocks):sizegroups[len(B)].append(i)
    ans=[]
    swaps=list(product(*(list(permutations(v)) for v in sizegroups.values())))
    for blockimages in swaps:
        target={i:j for Is,Js in zip(sizegroups.values(),blockimages) for i,j in zip(Is,Js)}
        for maps in product(*(affine(len(B)) for B in blocks)):
            for rootimages in permutations(fixed):
                q=list(range(n))
                for i,(B,a) in enumerate(zip(blocks,maps)):
                    C=blocks[target[i]]
                    for k,u in enumerate(B):q[u]=C[a[k]]
                for u,v in zip(fixed,rootimages):q[u]=v
                q=tuple(q);iq=m.inverse(q)
                if all(m.compose(m.compose(q,g),iq) in H for g in gens):ans.append(q)
    return tuple(ans)

def normalizer(n,H,gens):return list(normalizer_cached(n,H,tuple(gens)))
m.normalizer=normalizer

# Work on pair-orbit coordinates, retaining the entire signed action. This
# removes only actions identical on both the tournament cube and its signed dual.
def relative_records(n,H,gens,tournament):
    pairs,index,orbits,base=m.pair_data(n,H);N=normalizer(n,H,gens)
    coord={x:i for i,x in enumerate(orbits)};actions=set()
    for p in N:
        r=m.prep(p,pairs,index)
        image=tuple(coord[m.move(x,r)] for x in orbits)
        shifted=m.move(base,r,True)^base
        flip=sum(1<<i for i,x in enumerate(orbits) if shifted&x)
        need(shifted==sum(orbits[i] for i in range(len(orbits)) if flip>>i&1),'cube affine image')
        signbits=sum(1<<i for i,x in enumerate(orbits) if m.sign(x,r)==-1)
        actions.add((image,flip,signbits))
    # Keep multiplicities so relative stabilizer orders remain full N orders.
    mult=Counter()
    for p in N:
        r=m.prep(p,pairs,index);image=tuple(coord[m.move(x,r)] for x in orbits)
        shifted=m.move(base,r,True)^base
        flip=sum(1<<i for i,x in enumerate(orbits) if shifted&x)
        signbits=sum(1<<i for i,x in enumerate(orbits) if m.sign(x,r)==-1)
        mult[image,flip,signbits]+=1
    need(set(mult)==actions,'action compression')
    acts=[]
    for (im,flip,sgn),weight in mult.items():
        acts.append((tuple(1<<i for i in im),flip,sgn,weight))
    def move(s,r):
        ans=0;x=s
        while x:
            b=x&-x;ans|=r[0][b.bit_length()-1];x-=b
        return ans^(r[1] if tournament else 0)
    seen=set();records=[]
    for s in range(1<<len(orbits)):
        if s in seen:continue
        copies=set();stab=0;even=True
        for r in acts:
            y=move(s,r);copies.add(y)
            if y==s:
                stab+=r[3]
                if not tournament and (s&r[2]).bit_count()%2:even=False
        seen|=copies
        need(len(copies)*stab==len(N),'orbit stabilizer with compressed action')
        if tournament or even:
            x=(base if tournament else 0)^sum(o for i,o in enumerate(orbits) if s>>i&1)
            records.append({'mask':x,'relative_stabilizer_order':stab,'orbit_size':len(copies)})
    need(len(seen)==1<<len(orbits),'fixed cube covered')
    return records,pairs,orbits,N
m.relative_records=relative_records

def invariant(A):
    ca=[sum(row) for row in A];n=len(A)
    while True:
        palette=sorted(set(ca))
        sig=[(ca[i],tuple(sum(A[i][j] for j in range(n) if ca[j]==c) for c in palette),tuple(sum(A[j][i] for j in range(n) if ca[j]==c) for c in palette)) for i in range(n)]
        keys={s:i for i,s in enumerate(sorted(set(sig)))};nc=[keys[s] for s in sig]
        # The final signature, not merely the list of compressed labels.
        if len(set(ca))==len(set(nc)):return tuple(sorted(sig))
        ca=nc

def classify(records,n,pairs,tournament):
    reps=[];types=[];buckets=defaultdict(list)
    for r in records:
        A=m.matrix(r['mask'],n,pairs,tournament);key=invariant(A)
        j=next((i for i in buckets[key] if m.iso(A,reps[i])),None)
        if j is None:
            j=len(reps);reps.append(A);types.append({'multiplicity':0,'representative_mask':r['mask']});buckets[key].append(j)
        types[j]['multiplicity']+=1
    return reps,types
m.classify=classify

def cases():
    ans=[]
    for n in (11,12):
        t=m.cyc(n,[list(range(7))]);a=tuple([2*i%7 for i in range(7)]+list(range(7,n)));b=m.cyc(n,[[7,8,9]])
        c5=m.cyc(n,[list(range(5))]);c3=m.cyc(n,[[5,6,7]]);d3=m.cyc(n,[[8,9,10]])
        t11=m.cyc(n,[list(range(11))]);a11=tuple([3*i%11 for i in range(11)]+list(range(11,n)))
        ans.extend([
            (n,'F21_fixed',[t,a]),(n,'C15_5plus3',[c5,c3]),
            (n,'C21_7plus3',[t,b]),(n,'F21_diagonal_7plus3',[t,m.compose(a,b)]),
            (n,'F21_times_C3',[t,a,b]),(n,'F55_fixed',[t11,a11]),
            (n,'C15_5plus3plus3',[c5,m.compose(c3,d3)]),
            (n,'C5_times_C3_squared',[c5,c3,d3])])
        if n==12:
            q=m.cyc(n,[list(range(7,12))]);ans.extend([(n,'C35_7plus5',[t,q]),(n,'F21_times_C5',[t,a,q])])
    return ans

def main():
    selected=cases()
    if '--small' in sys.argv:selected=[x for x in selected if x[1] not in ('F21_fixed','C15_5plus3')]
    rows=[]
    for n,name,gens in selected:
        start=time.monotonic();r=m.census_case(n,name,gens);rows.append(r)
        print(n,name,'R',r['relative_tournament_orbits'],'T',r['tau_T'],'E',r['tau_E'],'F_T',r['tournament_fusion_defect'],'F_E',r['even_graph_fusion_defect'],'Z',r['sign_escape_relative_orbits'],'seconds',round(time.monotonic()-start,2),file=sys.stderr,flush=True)
    print(json.dumps({'status':'PASS','scope':'all classified odd non-prime-power subgroup conjugacy classes in S11 and S12' if '--small' not in sys.argv else 'selected small cubes from the complete n=11,12 classification','checks':m.CHECKS,'class_counts_by_degree':{str(n):sum(x[0]==n for x in cases()) for n in (11,12)},'rows':rows},indent=2,sort_keys=True))
if __name__=='__main__':main()
