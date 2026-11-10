#!/usr/bin/env python3
"""Exact lawful reset certificate. Standard library only; no assertions."""
import itertools,json,subprocess,random,hashlib,sys
from fractions import Fraction
from pathlib import Path
BASE='0e140250d07c'
M={'O':[[1,0,0,0],[0,0,0,0],[0,0,2,0],[0,0,0,1]],'Z':[[1,-1,1,0],[0,0,0,1],[0,0,2,-1],[0,0,0,1]],'J':[[1,-1,1,0],[1,-1,0,0],[1,0,1,0],[1,0,0,0]]}
I=[[int(i==j) for j in range(4)]for i in range(4)]
t=(1,0,2,1)
def need(x,m):
 if not x:raise RuntimeError(m)
def mul(a,b):return [[sum(x*y for x,y in zip(row,col))for col in zip(*b)]for row in a]
def apply(a,v):return tuple(sum(x*y for x,y in zip(row,v))for row in a)
def matrix(w):
 a=I
 for c in w:a=mul(M[c],a)
 return a
def state(w):return apply(matrix(w),t)
def rank(a):
 a=[[Fraction(x)for x in r]for r in a];k=0
 for j in range(len(a[0])):
  p=next((i for i in range(k,len(a))if a[i][j]),None)
  if p is None:continue
  a[k],a[p]=a[p],a[k];v=a[k][j];a[k]=[x/v for x in a[k]]
  for i in range(len(a)):
   if i!=k:
    v=a[i][j];a[i]=[x-v*y for x,y in zip(a[i],a[k])]
  k+=1
 return k
def graph(w):
 a=[2,0];root=0
 for c in w:
  n=len(a);a.append(0)
  for v in range(n):
   if c!='O' and v==root:a[n]|=1<<v
   else:a[v]|=1<<n
  if c=='J':root=n
 return a
def count(a):
 n=len(a);dp={}
 for i in range(n):dp[(1<<i,i)]=1
 for mask in range(1,1<<n):
  for j in range(n):
   z=dp.get((mask,j),0)
   if z:
    for k in range(n):
     if not(mask>>k&1) and a[j]>>k&1:
      key=(mask|1<<k,k);dp[key]=dp.get(key,0)+z
 return sum(dp.get(((1<<n)-1,j),0)for j in range(n))
def words(n):return map(''.join,itertools.product('OZJ',repeat=n))
r=matrix('OJOJ');need(r==[[3,0,6,0],[1,0,2,0],[3,0,6,0],[1,0,2,0]],'reset outer product')
need(matrix('OJ')==matrix('JZ'),'full rewrite')
for m in range(6):
 for k in range(1,7):
  v=state('O'*(k-1)+'J')
  target=[[z,0,z*2**(m+1),0]for z in v]
  for start in ('OJ','JZ'):need(matrix(start+'Z'*m+'O'*k+'J')==target,'general reset matrix')
a=graph('OJOJ')
for start in range(len(a)):
 seen=1<<start
 while True:
  more=seen
  for i in range(len(a)):
   if seen>>i&1:more|=a[i]
  if more==seen:break
  seen=more
 need(seen==(1<<len(a))-1,'literal strongness')
hist=[];resets=[]
for n in range(1,5):
 d={}
 for w in words(n):
  z=rank(matrix(w));d[z]=d.get(z,0)+1
  if z==1:resets.append(w)
 hist.append(d)
need(resets==['OJOJ','OJJZ','JZOJ','JZJZ'],'shortest resets')
# Stronger certificate than a finite state sample: multiply exact symbolic rows.
for n in range(6):
 for w in words(n):
  a=matrix(w);need(mul(a,r)==[[z,0,2*z,0]for z in apply(a,(3,1,3,1))],'all sampled suffix matrices')
literal=[]
for u in ('','O','Z','J','OZ','ZZ'):
 for w in ('','J','Z'):
  word=u+'OJOJ'+w;h=state(u)[0]+2*state(u)[2]
  got=count(graph(word));expected=h*state('J'+w)[0]
  need(got==expected,'independent literal reset count');literal.append([word,got])
# Hostile: scalar equality without matching response does not propagate.
need(state('O')[0]==state('')[0] and state('OJ')[0]!=state('J')[0],'hostile equal scalar')
paths=subprocess.check_output(['git','ls-tree','-r','--name-only',BASE],text=True).splitlines()
terms=['renewal','continuant','confluence','controllability','semigroup','invariant-quotient','automaton']
pool=sorted(p for p in paths if p.startswith('01-canon/theorems/') and any(t in p.lower()for t in terms))
indices=random.Random(24620261110).sample(range(len(pool)),min(5,len(pool)))
sampled=[{'index':i,'path':pool[i],'sha256':hashlib.sha256(subprocess.check_output(['git','show',BASE+':'+pool[i]])).hexdigest()}for i in indices]
if '--fail' in sys.argv:need(False,'deliberate hostile failure')
print(json.dumps({'base':BASE,'seed':24620261110,'pool_size':len(pool),'sample':sampled,'rank_histograms':hist,'shortest_rank_one_words':resets,'literal_reset_rows':graph('OJOJ'),'reset_matrix':r,'literal_checks':literal,'status':'PASS'},indent=2,sort_keys=True))
