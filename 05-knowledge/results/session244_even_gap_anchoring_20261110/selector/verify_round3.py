#!/usr/bin/env python3
"""Efficient endpoint source: full profiles, stabilizers, continuation arithmetic."""
import importlib.util
import itertools
import json
from pathlib import Path
p=Path(__file__).with_name('verify_round2.py');sp=importlib.util.spec_from_file_location('second',p)
r=importlib.util.module_from_spec(sp);sp.loader.exec_module(r)
v=r.v


def source(n):
    if n==1:return v.graph_word('Z'),0
    a=v.graph_word(v.word(n));root=len(a)
    a,_=r.graph_step(a,0,'O');a,_=r.graph_step(a,0,'Z')
    return a,root


def stabilizers(a,root):
    n=len(a);others=[i for i in range(n) if i!=root];count=0
    deg=[x.bit_count() for x in a]
    for perm in itertools.permutations(others):
        f=list(range(n))
        for i,j in zip(others,perm):f[i]=j
        if any(deg[i]!=deg[f[i]] for i in range(n)):continue
        if all(((a[i]>>j)&1)==((a[f[i]]>>f[j])&1) for i in range(n) for j in range(n)):count+=1
    return count


def main():
    examples=[]
    for n in range(1,129):
        a,root=source(n);g=v.g(n);eps=n%2
        expected=(g+2*n,n,2*g+2*n-eps,g+n-eps)
        v.require(v.profile(a,root)==expected,('full profile',n))
        v.require(all(v.reachable(a,q) for q in range(len(a))),('strong',n))
        v.require(2*n+1<=expected[0]<=3*n and 3*expected[2]<=4*expected[0],('cost',n))
        v.require(len(a)==(3 if n==1 else (n-1).bit_length()+3),('order',n))
        v.require((expected[0]==2*n+1)==(n&(n-1)==0),('lower equality',n))
        v.require((expected[0]==3*n)==(n==1 or n>=3 and (n-1)&(n-2)==0),('upper equality',n))
        if n<=16:v.require(stabilizers(a,root)==1,('stabilizer',n))
        if n>=2:
            z=len(a)-1
            v.require(a[root]==1<<z and a[z]==1 and root==z-1,('fixed-chain proof',n))
        if n<=32:
            b=a;s=expected
            for t in range(1,4):
                b,_=r.graph_step(b,root,'Z');s=v.step(s,'Z')
                v.require(v.profile(b,root)==s,('literal tail',n,t))
                v.require(s[0]==(2**t+1)*(g+n)-eps,('tail formula',n,t))
        if n in [1,2,3,7,21,37,101]:examples.append({'N':n,'profile':expected,'vertices':len(a)})
    tails=set()
    for n in range(1,2002):
        t=1
        while True:
            h=(2**t+1)*(v.g(n)+n)-n%2
            if h>2001:break
            tails.add(h);t+=1
    v.require(all(h not in tails for h in [13,31,127,173]),'tail hostiles')
    # The parent-requested 257 endpoint is separately recounted on its literal graph.
    a,q=source(257);prof=v.profile(a,q)
    v.require(prof==(v.g(257)+514,257,2*v.g(257)+513,v.g(257)+256),'257 profile')
    examples.append({'N':257,'profile':prof,'vertices':len(a)})
    print(json.dumps({'status':'PROVED efficient endpoint construction and fixed-root continuation law; exact replay',
      'literal_full_profiles':129,'literal_Z_tail_profiles':96,'root_stabilizer_enumerations':16,
      'examples':examples,'all_N_fixed_root_Z_tail_image_through_2001':len(tails),
      'certified_tail_hostiles':[13,31,127,173],
      'tail_misses_through_101':[h for h in range(1,102,2) if h not in tails]},indent=2,sort_keys=True))


if __name__=='__main__':main()
