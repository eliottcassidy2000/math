#!/usr/bin/env python3
"""All-prefix periodic-counter and bounded-macro certificates.
Standard library only. Exact integer/rational arithmetic; safe under python -O.
"""
import itertools,json,hashlib,math
from fractions import Fraction as Q
from pathlib import Path

def require(ok,msg):
 if not ok:raise ValueError(msg)
TRANS=[[[]for _ in range(36)]for _ in range(4)];ACCEPT=[]
for o,x,y,z in itertools.product(range(2),range(3),range(3),range(2)):
 s=18*o+6*x+2*y+z
 if x!=2 and (y==1 or y==0 and o==0):ACCEPT.append(s)
 for c in range(4):
  p=(c>>o)&1;q=(c>>(1-o))&1
  TRANS[c][s].append(18*o+6*(2 if p else 1)+2*y+z)
  if y or p:TRANS[c][s].append(18*o+6*x+2*(2 if q else 1)+z)
  if z or q:TRANS[c][s].append(18*o+6*x+2*y+1)
INIT=(1,)+(0,)*17+(1,)+(0,)*17
F=tuple(int(i in ACCEPT)for i in range(36))
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def step(v,c):
 r=[0]*36
 for i,a in enumerate(v):
  for j in TRANS[c][i]:r[j]+=a
 return tuple(r)
def state(w):
 v=INIT
 for c in w:v=step(v,c)
 return v
def value(w):return dot(F,state(w))
def row_step(r,c):return tuple(sum(r[j]for j in TRANS[c][i])for i in range(36))
def row_word(r,w):
 for c in reversed(w):r=row_step(r,c)
 return r

def span_reachable(prime=None):
 basis={};vectors=[]
 def add(v):
  r=[x%prime for x in v]if prime else list(map(Q,v))
  for j,b in sorted(basis.items()):
   if r[j]:
    a=r[j];r=[(x-a*y)%prime for x,y in zip(r,b)]if prime else [x-a*y for x,y in zip(r,b)]
  p=next((i for i,a in enumerate(r)if a),None)
  if p is None:return False
  a=pow(r[p],prime-2,prime)if prime else 1/r[p]
  basis[p]=tuple(x*a%prime for x in r)if prime else tuple(x*a for x in r)
  vectors.append(v);return True
 add(INIT);i=0
 while i<len(vectors):
  for c in range(4):
   v=step(vectors[i],c)
   if prime:v=tuple(x%prime for x in v)
   add(v)
  i+=1
 return vectors

def inverse(M):
 n=len(M);a=[[Q(x)for x in r]+[Q(i==j)for j in range(n)]for i,r in enumerate(M)]
 for j in range(n):
  i=next(i for i in range(j,n)if a[i][j]);a[i],a[j]=a[j],a[i]
  z=a[j][j];a[j]=[x/z for x in a[j]]
  for i in range(n):
   if i!=j and a[i][j]:
    z=a[i][j];a[i]=[x-z*y for x,y in zip(a[i],a[j])]
 return [r[n:]for r in a]
INTERP={}
def integer_weights(L,index):
 if L not in INTERP:
  a=2**L;b=3**L
  W=inverse([[1,k,k*k,a**k,k*a**k,b**k]for k in range(3,9)])
  INTERP[L]=W
 r=INTERP[L][index];den=math.lcm(*(x.denominator for x in r))
 return tuple(int(x*den)for x in r),den

def response_rows(w):
 r=F;rows=[]
 for k in range(9):
  if k>=3:rows.append(r)
  r=row_word(r,w)
 return rows

def coefficient_row(rows,L,index):
 weights,den=integer_weights(L,index)
 r=tuple(sum(a*row[i]for a,row in zip(weights,rows))for i in range(36))
 return r,den

def graph(w):
 n=len(w)+2;r=[0]*n;r[0]=2
 for j,c in enumerate(w,2):
  for a in (0,1):
   if c>>a&1:r[a]|=1<<j
   else:r[j]|=1<<a
  for i in range(2,j):r[i]|=1<<j
 return r

def dp(r):
 n=len(r);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for mask in range(1,1<<n):
  for j in range(n):
   if mask>>j&1:
    old=mask^(1<<j)
    if old:D[mask][j]=sum(D[old][i]for i in range(n)if r[i]>>j&1)
 return sum(D[-1])

def check_common_flag():
 # Bounded local basis E,B,D=G-B; final local basis E,O.
 ss=list(itertools.product(range(3),range(3),range(2)))
 def changed_column(o,t):
  out=[0]*36;choices=[]
  for i,x in enumerate(t):
   choices.append([(x,1)] if i==2 or x==0 else ([(2,1)] if x==1 else [(1,1),(2,-1)]))
  for choice in itertools.product(*choices):
   x,y,z=(a for a,b in choice);out[18*o+6*x+2*y+z]+=math.prod(b for a,b in choice)
  return tuple(out)
 checks=0
 for o in range(2):
  for t in ss:
   col=changed_column(o,t)
   for c in range(4):
    bits=((c>>o)&1,(c>>(1-o))&1);tri={}
    for i in range(3):
     x=t[i];incoming=1 if i==0 else bits[i-1]
     targets=[]
     if i==2:
      if x or incoming:targets=[(1,1)]
     elif x!=2:
      alpha=1 if x==1 else incoming
      if alpha:
       targets=[(1,alpha)]
       if not bits[i]:targets.append((2,alpha))
     for y,a in targets:
      dest=list(t);dest[i]=y;dest=tuple(dest);tri[dest]=tri.get(dest,0)+a
    expect=[0]*36
    for dest,a in tri.items():
     require(all(x<=y for x,y in zip(t,dest)),'tensor order')
     z=changed_column(o,dest);expect=[x+a*y for x,y in zip(expect,z)]
    require(tuple(expect)==step(col,c),'explicit common basis identity')
    require(tri.get(t,0)==sum(x==1 for x in t),'common diagonal')
    checks+=1
 return checks

def main():
 flag_checks=check_common_flag()
 R=span_reachable();require(len(R)==28,'reachable dimension')
 # Modulo p closure gives a finite, all-word modular invariant certificate.
 primes=(2,3,5,7,11,13,17,19,23,31,43,73,127)
 RP={p:span_reachable(p)for p in primes}
 records=[];stats={};digest=hashlib.sha256()
 def certify(w,kind):
  L=len(w);rows=response_rows(w);r2,d2=coefficient_row(rows,L,2);a=dot(r2,INIT)
  # Every reachable physical vector has normalizing coordinate 1.
  require(all(dot(r2,v)==a for v in R),'prefix-independent quadratic coefficient')
  if a:
   rec={'word':''.join(map(str,w)),'kind':kind,'counter':'quadratic','coefficient':str(Q(a,d2))}
  else:
   r,d=coefficient_row(rows,L,1);g=math.gcd(*r);require(g,'nonzero linear row')
   r=tuple(x//g for x in r);h=dot(r,INIT);found=None
   for p in primes:
    if h%p and all(dot(r,v)%p==h%p for v in RP[p]):found=(p,h%p);break
   require(found is not None,'nonzero modular linear counter '+str(w))
   rec={'word':''.join(map(str,w)),'kind':kind,'counter':'linear','scale':str(Q(g,d)),
        'prime':found[0],'residue':found[1]}
   # Primitive integer row plus all-word modular proof, not a scalar finite fit.
   digest.update((str(w)+str(r)+str(Q(g,d))+str(found)+'\n').encode())
  records.append(rec);return rec
 for L in range(2,7):
  batch=[]
  for w in itertools.product(range(4),repeat=L):
   if len(set(w))>1:batch.append(certify(w,'all_blocks'))
  stats[str(L)]={'mixed_blocks':len(batch),'quadratic':sum(r['counter']=='quadratic'for r in batch),
   'linear':sum(r['counter']=='linear'for r in batch)}
 # These composition periods eliminate two distinct constant-letter branches.
 for L in range(1,7):
  for c,d in itertools.permutations(range(4),2):certify((c,)*L+(d,)*L,'constant_branch_composition')
 # Independent graph DP controls for complete small-word universe and selected longer words.
 checks=0
 for L in range(5):
  for w in itertools.product(range(4),repeat=L):
   require(value(w)==dp(graph(w)),'subset-DP control');checks+=1
 for w in ((1,2,0,1,2,0),(0,1,0,3,2,3),(2,3,3,0,0),(1,3,2,0,0),(2,1,3,0,0)):
  require(value(w)==dp(graph(w)),'longer DP control');checks+=1
 # Hostile and positive controls: constant periods can genuinely be affine.
 require([value((2,)+(0,)*k)for k in range(8)]==[2**(k+1)+1 for k in range(8)],'affine constant control')
 require([value((2,3,3,0)+(0,)*k)for k in range(8)]==[42*3**k-17 for k in range(8)],'pure ternary control')
 require(value((1,2,0))==13,'H13 hostile')
 # Full-vector constant-letter annihilator on the invariant reachable space.
 for c in range(4):
  for v in R:
   v1=step(v,c);v2=step(v1,c);v3=step(v2,c);v4=step(v3,c)
   require(all(a-6*b+11*d-6*e==0 for a,b,d,e in zip(v4,v3,v2,v1)),'constant control recurrence')
 rec={'status':'PASS','scope':'all prefixes; mixed periods through length 6; distinct constant-branch composition periods through length 12',
      'reachable_dimension':len(R),'common_flag_column_checks':flag_checks,'modular_reachable_dimensions':{str(p):len(RP[p])for p in primes},
      'block_counts':stats,'composition_blocks':72,'independent_DP_checks':checks,
      'linear_rows_sha256':digest.hexdigest(),'certificates':records}
 raw=(json.dumps(rec,sort_keys=True,indent=2)+'\n').encode();Path(__file__).with_name('receipt.json').write_bytes(raw)
 print(json.dumps({'status':'PASS','all_mixed_blocks':sum(s['mixed_blocks']for s in stats.values()),
  'composition_blocks':72,'independent_DP_checks':checks,'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__':main()
