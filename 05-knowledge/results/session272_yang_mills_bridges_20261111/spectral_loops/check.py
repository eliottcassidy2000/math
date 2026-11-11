#!/usr/bin/env python3
"""Exact finite controls, not a reconstruction/source verification."""
from fractions import Fraction as Q
from itertools import product
from math import exp

def lattice(torus):
    n=2
    verts=range(n if torus else n+1)
    edges=[]
    for x in verts:
        for y in verts:
            for d in range(2):
                if torus or (x if d==0 else y)<n:
                    edges.append((x,y,d))
    idx={e:i for i,e in enumerate(edges)}
    nxt=lambda v: (v+1)%n if torus else v+1
    faces=[]
    for x in range(n):
        for y in range(n):
            faces.append(sum(1<<idx[e] for e in [(x,y,0),(nxt(x),y,1),(x,nxt(y),0),(x,y,1)]))
    def sign(bits,mask): return -1 if (bits&mask).bit_count()%2 else 1
    nums=[0]*16; Z=0; toron_num=0
    toron=sum(1<<idx[(x,0,0)] for x in range(n))
    for bits in range(1<<len(edges)):
        Bs=[sign(bits,p) for p in faces]
        weight=2**sum(b==1 for b in Bs) # q=1/3; factors 2 (B=+),1 (B=-)
        Z+=weight; toron_num+=weight*sign(bits,toron)
        for subset in range(16):
            boundary=0; face_prod=1
            for j,p in enumerate(faces):
                if subset>>j&1: boundary^=p; face_prod*=Bs[j]
            assert sign(bits,boundary)==face_prod # independent edge/face path
            nums[subset]+=weight*sign(bits,boundary)
    q=Q(1,3)
    for subset,num in enumerate(nums):
        a=subset.bit_count()
        expected=(q**a+q**(4-a))/(1+q**4) if torus else q**a
        assert Q(num,Z)==expected,(torus,subset,Q(num,Z),expected)
    if torus: assert toron_num==0
    print(('torus' if torus else 'disk'),f'edges={len(edges)} configurations={1<<len(edges)} subsets=16 PASS')
    if torus: print('torus a=1:',Q(nums[1],Z),'a=2:',Q(nums[3],Z),'noncontractible:',Q(toron_num,Z))

lattice(False); lattice(True)
# Fixed-word spectral controls: A swaps with epsilon; B gives diagonal signs.
count=0
for n in range(1,9):
    e=Q(1,2**n)
    for length in range(9):
        for word in product('AB',repeat=length):
            v=[Q(1),Q(0)]
            for a in reversed(word): v=[e*v[1],e*v[0]] if a=='A' else [v[0],-v[1]]
            assert abs(v[1])<=e
            for u in [e/2,e,2*e,Q(1)]:
                mass=v[1]**2 if e<=u else 0
                assert mass<=u*u
            count+=1
print('uniform-phase fixed monomial checks:',count,'PASS; gaps 1/2 through 1/256')
# Parity transfer: independently propagate probabilities.
q=Q(1,3); p=(1+q)/2
v=[Q(1),Q(0)]
for a in range(21):
    assert v[0]-v[1]==q**a
    v=[p*v[0]+(1-p)*v[1],(1-p)*v[0]+p*v[1]]
print('two-channel transfer powers 0..20 PASS')
# Floating diagnostic only: no exact claim depends on these numbers.
for t in [1,2,10,100]:
    c=2*(1-(1+t)*exp(-t))/(t*t)
    print(f'gapless continuous control t={t}: C(t)={c:.12g}, t^2*C(t)={t*t*c:.12g}')
print('Exact finite checks complete; YM hypotheses, continuum limits and mass gap NOT tested.')
