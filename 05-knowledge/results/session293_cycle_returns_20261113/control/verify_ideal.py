#!/usr/bin/env python3
import itertools,math,json
checks=0
def req(q):
 global checks;checks+=1
 if not q:raise RuntimeError('ideal check failed')
def rows_of(n,mask):
 r=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if mask>>k&1:r[i]|=1<<j
   else:r[j]|=1<<i
   k+=1
 return r
def profile(r):
 n=len(r);N=1<<n;H=[0]*N;H[0]=1;end=[[0]*n for _ in range(N)]
 for i in range(n):end[1<<i][i]=1
 for s in range(1,N):
  for v in range(n):
   if end[s][v]:
    q=r[v]&~s
    while q:
     b=q&-q;q-=b;end[s|b][b.bit_length()-1]+=end[s][v]
  H[s]=sum(end[s])
 F=[[0]*N for _ in range(n+1)];F[0][0]=1
 for k in range(1,n+1):
  for s in range(1,N):
   t=s
   while t:F[k][s]+=H[t]*F[k-1][s^t];t=(t-1)&s
 return H,[F[k][-1] for k in range(1,n+1)]
def exterior_coeff(m):
 n=2*m;a=[0]*(m+2)
 if m==0:return [0,1]
 for perm in itertools.permutations(range(n)):
  legal=[];required=[];ok=True
  for g in range(n+1):
   allowed=(g==0 or perm[g-1]%2==1) and (g==n or perm[g]%2==0)
   if allowed:legal.append(g)
   if 0<g<n and perm[g-1]>perm[g]:
    if not allowed:ok=False;break
    required.append(g)
  if ok:
   for k in range(max(1,len(required)),len(legal)+1):a[k]+=math.comb(len(legal)-len(required),k-len(required))
 return a
coeff=[exterior_coeff(m) for m in range(5)]
for m,row in enumerate(coeff):
 for k,a in enumerate(row[1:],1):req(a%math.factorial(k-1)==0)
 req(row[-1]==math.factorial(m))
examples=[];stats={}
for n in range(1,6):
 distinct=set()
 for mask in range(1<<(n*(n-1)//2)):
  r=rows_of(n,mask);H,F=profile(r);D=[math.factorial(k)*v for k,v in enumerate(F)]
  g=math.gcd(*D);req(H[-1]%g==0 and g%2==1)
  probes=[sum(row[k]*F[k-1] for k in range(1,min(len(row),n+1))) for row in coeff[:n]]
  req(math.gcd(*probes)==g)
  for p in [2,3,5,7]:
   req((g%p==0)==all(F[k-1]%p==0 for k in range(1,min(p,n+1))))
  distinct.add((H[-1],g))
  if n==5 and mask in [8,76]:
   aut=sum(all(((r[i]>>j)&1)==((r[s[i]]>>s[j])&1) for i in range(n) for j in range(n)) for s in itertools.permutations(range(n)))
   examples.append({'mask':mask,'rows':r,'H':H[-1],'ordered_covers':F,'scaled_covers':D,'context_ideal':g,'aut':aut,'probes':probes})
 stats[n]=len(distinct)
req(examples[0]['context_ideal']==3 and examples[0]['aut']==1)
req(examples[1]['context_ideal']==15 and examples[1]['aut']==5)
print(json.dumps({'checks':checks,'probe_coefficients':coeff,'examples':examples,'distinct_H_ideal_pairs':stats},indent=2))
