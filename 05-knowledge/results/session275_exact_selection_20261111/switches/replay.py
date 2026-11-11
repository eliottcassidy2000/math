#!/usr/bin/env python3
from carrier import *
from math import gcd
A=(6,16,2,7,13,31);B=(4,37,0,7,15,29)
X=(46,8,2,4,15,30);Y=(46,20,0,6,13,30)
MA,MB=matrix(profile(A)),matrix(profile(B))
D=[[MA[i][j]-MB[i][j] for j in range(4)] for i in range(4)]
assert D==[[0,6,2,4],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
I=matrix(profile((0,)))
def matword(bits):
 M=I
 for b in bits:M=mul(M,MA if b else MB)
 return M

def weights(m):
 def g(r):return 4 if r==0 else (162*64**r-44*5**r)//59
 return [5**i*g(m-1-i) for i in range(m)]
def count(M):return M[0][0]+M[0][3]
def decode(m,target):
 base=count(matword([0]*m));q=target-base;bits=[]
 for w in weights(m):
  if q>=w:q-=w;bits.append(1)
  else:bits.append(0)
 return bits if q==0 else None
receipt={'profiles':{k:profile(t) for k,t in [('A',A),('B',B),('X',X),('Y',Y)]}}
# Literal pair compositions through 11 vertices, no matrix counts as expectations.
literal=[]
for bits in product((0,1),repeat=2):
 T=fuse(A if bits[0] else B,A if bits[1] else B)
 assert hp(T)[-1]==count(matword(bits))
 literal.append({'word':bits,'H':hp(T)[-1]})
 hp.cache_clear()
# Noncommutative spacers: D C D must vanish; literal contexts through 15 vertices.
Z=[[0]*4 for _ in range(4)]
spacers=[(0,),(2,0),(0,1),(2,4,1)]
for C in spacers:assert mul(mul(D,matrix(profile(C))),D)==Z
context=[]
for C in spacers:
 values={}
 for a,b in product((0,1),repeat=2):
  T=fuse(fuse(A if a else B,C),A if b else B)
  values[a,b]=hp(T)[-1];hp.cache_clear()
 assert values[1,1]-values[1,0]-values[0,1]+values[0,0]==0
 context.append({'spacer':C,'counts':{str(k):v for k,v in values.items()}})
levels=[]
for m in range(1,9):
 W=weights(m);base=count(matword([0]*m));seen=set()
 assert all(W[i]>sum(W[i+1:]) for i in range(m))
 for bits in product((0,1),repeat=m):
  h=count(matword(bits));assert h==base+sum(w*b for w,b in zip(W,bits))
  assert decode(m,h)==list(bits);seen.add(h)
 assert len(seen)==2**m
 g=0
 for w in W:g=gcd(g,w)
 assert g==4
 # Every tested off-spectrum target rejected by exact decoder, no completeness claim.
 for h in range(base,min(base+500,max(seen)+1)):
  assert (decode(m,h) is not None)==(h in seen)
 levels.append({'m':m,'order':1+5*m,'baseline':base,'weights':W,'distinct':len(seen),'gcd':g})
# Elementary h=1 boundary: literal root-deletion-transitive words, independently counted.
boundary=0
for n in range(1,8):
 for bits in product((0,1),repeat=n):
  T=(0,)
  for b in bits:T=fuse(T,(2,0) if b else (0,1))
  p=profile(T)
  assert p['h']==1 and p['Q']==2**n
  assert p['V']==1+sum((1-b)*2**i for i,b in enumerate(bits))
  boundary+=1;hp.cache_clear()
receipt.update(literal_pairs=literal,spacer_controls=context,levels=levels,boundary_controls=boundary,
 scope='All-order theorems in PROOF.md; finite controls of new consequences, not Chan-Pak coverage')
Path(__file__).with_name('RESULTS.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'literal_pairs':len(literal),'literal_spacer_cubes':len(context),'affine_levels':len(levels),'boundary_controls':boundary}))
