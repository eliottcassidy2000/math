#!/usr/bin/env python3
import itertools,json

def tournament(n,b):
 r=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if b>>k&1:r[i]|=1<<j
   else:r[j]|=1<<i
   k+=1
 return r

def hp(r):
 n=len(r);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,1<<n):
  for i in range(n):
   if not d[s][i]:continue
   a=r[i]&~s
   while a:
    bit=a&-a;a-=bit;j=bit.bit_length()-1
    d[s|bit][j]+=d[s][i]
 return sum(d[-1])

def expand(r,v,m):
 labels=[i for i in range(len(r)) if i!=v]+[v]*m
 out=[0]*len(labels)
 for a,i in enumerate(labels):
  for b,j in enumerate(labels):
   if (i==j and a<b) or (i!=j and r[i]>>j&1):out[a]|=1<<b
 return out

def single_scc(r,v):
 reach=[x|(1<<i) for i,x in enumerate(r)]
 for k in range(len(r)):
  for i in range(len(r)):
   if reach[i]>>k&1:reach[i]|=reach[k]
 return not any(i!=v and reach[v]>>i&1 and reach[i]>>v&1 for i in range(len(r)))

def kernel(r,v):
 n=len(r);c=[0]*(n+1)
 for p in itertools.permutations([i for i in range(n) if i!=v]):
  poly=[1]
  for k in range(n):
   e=int(k==0 or k==n-1 or r[p[k-1]]>>p[k]&1)
   f=int((k==0 or r[p[k-1]]>>v&1) and (k==n-1 or r[v]>>p[k]&1))
   q=[0]*(len(poly)+1)
   for j,x in enumerate(poly):q[j]+=e*x;q[j+1]+=f*x
   poly=q
  for j,x in enumerate(poly):c[j]+=x
 return c

def ordered_covers_tt(m,j):
 return sum((-1)**k*__import__('math').comb(j,k)*(j-k)**m for k in range(j+1))

if hp([40,1,3,6,15,30])!=19:raise RuntimeError('positive attainable-target control')
roots=0;nontrivial=0;composites=0;max_c2=0
for n in range(1,6):
 for b in range(1<<(n*(n-1)//2)):
  r=tournament(n,b);h=hp(r)
  for v in range(n):
   roots+=1;c=kernel(r,v);single=single_scc(r,v)
   if c[1]!=h or (c[2]==0 if n>=2 else True)!=single:raise RuntimeError(('kernel',n,b,v,c,single))
   nontrivial+=not single;max_c2=max(max_c2,c[2] if n>=2 else 0)
   for m in range(1,4):
    actual=hp(expand(r,v,m));composites+=1
    predicted=sum(c[j]*ordered_covers_tt(m,j) for j in range(1,len(c)))
    if actual!=predicted:raise RuntimeError(('formula',n,b,v,m))
    if single and actual!=h:raise RuntimeError('singleton')
    if not single and actual<h+(2**m-2):raise RuntimeError('growth')
    if m==2 and (actual==h)!=single:raise RuntimeError('two-test')
padding_checks=0
for n in range(1,5):
 for b in range(1<<(n*(n-1)//2)):
  r=tournament(n,b);h=hp(r)
  def q2(rows):
   n=len(rows);total=0
   for mask in range(1<<n):
    hs=[]
    for subset in (mask,((1<<n)-1)^mask):
     ids=[i for i in range(n) if subset>>i&1]
     rr=[sum((1<<j) for j,w in enumerate(ids) if rows[v]>>w&1) for v in ids]
     hs.append(hp(rr) if rr else 1)
    total+=hs[0]*hs[1]
   return total
  q=q2(r)
  for m in range(4):
   rows=[x|(((1<<m)-1)<<n) for x in r]+[sum(1<<j for j in range(n+i+1,n+m)) for i in range(m)]
   if hp(rows)!=h or q2(rows)!=(1<<m)*q:raise RuntimeError('arbitrary padding')
   padding_checks+=1
print(json.dumps(dict(status='PASS',universe='all labeled tournaments of orders 1..5; every marked vertex; TT blowups of sizes 1..3',rooted_contexts=roots,nontrivial_scc_roots=nontrivial,literal_composites=composites,max_c2=max_c2,arbitrary_input_padding_checks=padding_checks,independent_checks=['vertex-subset Hamilton-path DP','permutation-and-gap kernel','transitive-closure SCC criterion']),indent=2))
