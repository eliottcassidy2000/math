#!/usr/bin/env python3
"""Exact finite tests of the all-m binary affine-incidence theorem."""
from itertools import combinations, product
from collections import Counter
from math import factorial
CHECKS=0
def check(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise RuntimeError(msg)
def gaussian(n,k,q=2):
 if k<0 or n<k:return 0
 a=b=1
 for i in range(k):a*=q**(n-i)-1;b*=q**(k-i)-1
 check(a%b==0,'Gaussian coefficient integral')
 return a//b
def rank2(vs):
 piv={}
 for v in vs:
  while v:
   i=v.bit_length()-1
   if i in piv:v^=piv[i]
   else:piv[i]=v;break
 return len(piv)
def adim(S):
 o=S[0]
 return rank2(x^o for x in S[1:])
def flats2(m):
 n=1<<m;spaces={tuple(sorted((a,b,a^b))) for a,b in combinations(range(1,n),2)}
 flats=set()
 for row in spaces:
  a,b=row[:2]
  for c,d in product(range(2),repeat=2):
   flats.add(frozenset(x for x in range(n) if (x&a).bit_count()%2==c and (x&b).bit_count()%2==d))
 return flats
def response_counts(flats,t):
 return Counter(S for B in flats for S in combinations(sorted(B),t))
rows=[]
for m in [3,4,5]:
 n=1<<m;flats=flats2(m)
 check(len(flats)==4*gaussian(m,2),'flat count')
 check(all(len(B)==1<<(m-2) for B in flats),'flat size')
 response={t:response_counts(flats,t) for t in range(1,5)}
 for t in range(1,5):
  for S in combinations(range(n),t):
   check(response[t][S]==gaussian(m-adim(S),2),('incidence formula',m,S))
 if m>=4:
  plane=gaussian(m-2,2);independent=gaussian(m-3,2)
  check(plane>independent,'four-port rank separation')
  for o in range(n):
   # Reconstruction works at every named origin, not only the zero chart.
   add=[[None]*n for _ in range(n)]
   for a in range(n):add[a][a]=o;add[a][o]=add[o][a]=a
   for a,b in combinations([x for x in range(n) if x!=o],2):
    cs=[c for c in range(n) if c not in (o,a,b) and response[4][tuple(sorted((o,a,b,c)))]==plane]
    check(cs==[o^a^b],'unique recovered fourth point')
    add[a][b]=add[b][a]=cs[0]
   for a,b in product(range(n),repeat=2):check(add[a][b]==a^b^o,'based affine addition')
  # Explicit hostile relabelling: nonaffine transposition0<->1.
  perm=[1,0]+list(range(2,n))
  transported={frozenset(perm[x] for x in B) for B in flats}
  check(transported!=flats,'distinct labelled affine structures')
  for t in range(1,4):check(response_counts(transported,t)==response[t],'all lower-port rows agree')
  check(response_counts(transported,4)!=response[4],'four-port row separates hostile')
  # Construct minimal Boolean quadratic word bank independently: all products
  # of affine linear forms have exactly the enumerated codimension-two supports.
  words=set()
  for a,b in combinations(range(1,n),2):
   for c,d in product(range(2),repeat=2):
    words.add(frozenset(x for x in range(n) if ((a&x).bit_count()%2^c) and ((b&x).bit_count()%2^d)))
  check(words==flats,'quadratic products give every flat')
  generators=[(1<<n)-1]
  generators += [sum(1<<x for x in range(n) if x>>i&1) for i in range(m)]
  generators += [sum(1<<x for x in range(n) if (x>>i&1) and (x>>j&1)) for i,j in combinations(range(m),2)]
  rm_words=[0]
  for v in generators:rm_words += [w^v for w in rm_words]
  check(len(set(rm_words))==1<<len(generators),'full RM dimension')
  minimum=min(w.bit_count() for w in rm_words if w)
  check(minimum==1<<(m-2),'RM minimum weight')
  min_supports={frozenset(i for i in range(n) if w>>i&1) for w in rm_words if w.bit_count()==minimum}
  check(min_supports==flats,'complete RM minimum-word support bank')
 else:
  check(flats=={frozenset(s) for s in combinations(range(n),2)},'m3 sharp hostile: all pairs')
  check(not response[3] and not response[4],'m3 no higher-port data')
 gl=1
 for i in range(m):gl*=n-(1<<i)
 structures=factorial(n)//(n*gl)
 check(factorial(n)%(n*gl)==0,'exact labelled geometry quotient')
 bits=(structures-1).bit_length()
 rows.append((m,len(flats),[gaussian(m-i,2) for i in range(4)],structures,bits))
# Separate q=3,m=3 control: triples distinguish collinearity.
points=list(product(range(3),repeat=3))
def add3(a,b):return tuple((x+y)%3 for x,y in zip(a,b))
def scale3(c,a):return tuple(c*x%3 for x in a)
lines=set()
for a in points:
 for v in points[1:]:lines.add(frozenset(add3(a,scale3(t,v)) for t in range(3)))
check(len(lines)==117,'AG3(F3) lines')
triples=Counter(S for L in lines for S in combinations(sorted(L),3))
collinear=0
for S in combinations(points,3):
 a,b,c=S;v=tuple((y-x)%3 for x,y in zip(a,b))
 inline=c in {add3(a,scale3(t,v)) for t in range(3)}
 check(triples[S]==int(inline),'q3 third-order line detector')
 collinear+=inline
check(collinear==117,'117 collinear triples')
print('CHECKS',CHECKS)
print('m, codim2flat count, point/pair/triple/independent4 responses, labelled structures, minimum fixed-length bits')
for row in rows:print(row)
print('For m4/m5 every named origin reconstructs exactly its translated XOR law.')
print('Nonaffine label swap preserves all <=3 responses and changes four-port rows.')
print('m3: codim2flats are all pairs; all30 labelled affine geometries remain invisible even to complete incidence.')
print('q3,m3: triples detect collinearity (117 positive of2925 triples); no claim to recover scalar labels.')
