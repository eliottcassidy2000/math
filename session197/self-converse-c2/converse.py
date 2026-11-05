"""Physical source pairs and literal converse maps; not a scalar completeness guess."""
from functools import lru_cache
from core import inventories,require,adjacency_map
from bilinear import construct,converse,value

@lru_cache(None)
def strong_sources(n):
 if n==1:return (((0,),0),)
 out=[]
 for m in range(1,n):
  w=n-m-1
  if m==1 and w==0:continue
  for V,c in strong_sources(m):
   for bits in range(1<<max(0,w-1)):
    A=list(V)+[0]*(w+1);r=m
    for x in range(m):
     if x==c:A[r]|=1<<x
     else:A[x]|=1<<r
     for y in range(m+1,n):A[x]|=1<<y
    for i in range(w):
     v=m+1+i
     if bits>>i&1:A[r]|=1<<v
     else:A[v]|=1<<r
     for j in range(i+1,w):A[v]|=1<<(m+1+j)
    out.append((tuple(A),r))
 return tuple(out)

@lru_cache(None)
def sources(n):
 out=[]
 for m in range(1,n+1):
  for V,c in strong_sources(m):
   A=list(V)+[0]*(n-m)
   for i in range(n):
    for j in range(max(i+1,m),n):A[i]|=1<<j
   out.append((tuple(A),c))
 return tuple(out)

def swap_map(kind,n,m):
 k=1 if kind=='N'else 2
 central=[m]if k==1 else[m+1,m]
 return list(range(m+k,m+k+n))+central+list(range(m))

def cycle_words(A):
 n=len(A);out=[]
 def visit(word,mask):
  if len(word)==n:
   if A[word[-1]]&1:out.append(tuple(word))
   return
  rem=A[word[-1]]&(((1<<n)-1)^mask)
  while rem:
   bit=rem&-rem;rem-=bit;v=bit.bit_length()-1
   visit(word+[v],mask|bit)
 visit([0],1);return out

def relabel(A,word):
 return tuple(sum(((A[x]>>y)&1)<<j for j,y in enumerate(word))for x in word)

def canonical(A,cycles=None):
 C=cycle_words(A)if cycles is None else cycles
 return min(relabel(A,w[i:]+w[:i])for w in C for i in range(len(A)))

def anti_maps(A,cycles=None):
 C=cycle_words(A)if cycles is None else cycles;n=len(A);first=C[0];out=[]
 # Any anti-isomorphism takes this directed Hamilton cycle to a reverse cycle.
 for w in C:
  for i in range(n):
   p=[0]*n
   for j,v in enumerate(first):p[v]=w[(i-j)%n]
   if adjacency_map(converse(A),A,p):out.append(tuple(p))
 return out

def profile(A,c):
 q=inventories(A);require(q['s'][c]==1,'source start count')
 return q['H'],q['e'][c],q['R'][c],q['allH'][((1<<len(A))-1)^(1<<c)]

def diagonal_inverse(target):
 from math import isqrt
 from bilinear import source_bank,build_source
 if target==15:return dict(kind='R5',rows=[6,12,24,17,3])
 if target<1:return None
 for v,recipe in source_bank(isqrt((target-1)//3)).items():
  if v==(1,1,1,1):continue
  for kind in('N','M'):
   if value(kind,v,v)==target:
    V,c=build_source(recipe);T,_=construct(kind,V,c,V,c)
    return dict(kind=kind,source=list(v),recipe=list(recipe),rows=T)
 return None
