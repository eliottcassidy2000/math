"""Verify physical 33-offset dyadic window and the next three genuine signed shortcut steps."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/"moonshine_prime_blocks"))
from replay import wrap as compose,reads,dp
def count(T):return sum(dp(T)[-1])
import csv,sys
def tt(n):
 return tuple((1<<n)-(1<<(i+1)) for i in range(n))
def step(x):
 s=1 if (3*x+1)%4==0 else -1
 z=3*x+s; v=(abs(z)&-abs(z)).bit_length()-1
 assert v>=2
 return s,v,z//(1<<v)
def orbit(B,J):
 v=0;b=B;yield 0,b,v
 for j in range(1,J+1):
  sig,w,b=step(b);v+=w;yield j,b,v
def spine(D,I,w):
 low=tt(w); n=len(D); M=w+n; apex=3|(I<<w)
 T=tuple(r|(((1<<n)-1)<<w) for r in low)+tuple(x<<w for x in D)
 return tuple(x|(((apex>>i&1)^1)<<M)for i,x in enumerate(T))+(apex,)
def spine_data(D,I):
 E=dp(D);f=(1<<len(D))-1; H=[1]+[sum(e)for e in E[1:]]
 R=sum(H[f^a]*sum(x for i,x in enumerate(E[a])if not(I>>i&1)) for a in range(1,f+1))
 A=tuple(x|(((I>>i&1)^1)<<len(D)) for i,x in enumerate(D))+(I,)
 return 3*(H[-1]+R),count(A)-2*H[-1]
def replay(path,depth=3):
 P={}
 with open(path) as f:
  for r in csv.DictReader(f,delimiter='\t'):
   j,b,s,G,I=[int(r[k])for k in ('level','alpha','scale','G','I')]
   D=tuple(map(int,r['outgoing_masks'].split(',')))
   if r['kind']=='spine':
    assert spine_data(D,I)==(G,s*b)
   else:
    h,g,q=reads(D,I)
    assert G==g and h+q-g==s*b
   a,rem=divmod(G,s*3**j)
   assert not rem and a>0 and not (a&(a-1))
   k=a.bit_length()-1
   for w in (1,2):
    T=spine(D,I,w+2) if r['kind']=='spine' else compose(D,I,tt(w))
    assert count(T)==s*(3**j*(1<<(k+w))+b)
   P[j,b,s]=(D,I,k)
 R=range(-33,32,2)
 needed={(j,b,s)for B in R for j,b,v in orbit(B,depth)for s in (1,7)}
 assert needed==set(P)
 m0=max(P[j,b,s][-1]+v+1 for B in R for j,b,v in orbit(B,depth)for s in (1,7))
 for m in (m0,m0+1,32):
  for B in R:
   trail=[3**j*(1<<(m-v))+b for j,b,v in orbit(B,depth)]
   for x,y in zip(trail,trail[1:]):
    s,k,z=step(x)
    assert z==y and 0<y<x
 p=[3**j*(1<<(16-v))+b for j,b,v in orbit(-17,2)]
 from math import isqrt
 assert p==[65519,49139,18427] and all(all(q%d for d in range(2,isqrt(q)+1)) for q in p)
 print(f'{len(P)} verified template rows, {2*len(P)} full graph controls, common m >= {m0}')
if __name__=='__main__':replay(sys.argv[1],int(sys.argv[2])if len(sys.argv)>2 else 3)
