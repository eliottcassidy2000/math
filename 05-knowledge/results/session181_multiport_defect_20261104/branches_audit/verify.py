from itertools import combinations,product,permutations
import json

def ck(x,msg):
 if not x:raise RuntimeError(msg)
def hp(r):
 n=len(r);d=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for m in range(1,1<<n):
  for v in range(n):
   if d[m][v]:
    for w in range(n):
     if not(m>>w&1) and r[v]>>w&1:d[m|1<<w][w]+=d[m][v]
 return sum(d[-1]) if n else 1

def sub(r,V):return [sum(1<<j for j,v in enumerate(V) if r[u]>>v&1) for u in V]
def strong(r):
 n=len(r)
 for v in range(n):
  frontier={v};seen=set()
  while frontier:
   u=frontier.pop();seen.add(u)
   frontier.update(w for w in range(n) if r[u]>>w&1 and w not in seen)
  if len(seen)!=n:return False
 return True

def blow(q,sizes):
 b=[i for i,s in enumerate(sizes) for _ in range(s)]
 return [sum(1<<v for v,j in enumerate(b) if (i==j and u<v) or (i!=j and q[i]>>j&1)) for u,i in enumerate(b)]

mono=0;blows=0;qs=0
for n in range(1,5):
 es=list(combinations(range(n),2))
 for bits in range(1<<len(es)):
  r=[0]*n
  for k,(u,v) in enumerate(es):
   if bits>>k&1:r[u]|=1<<v
   else:r[v]|=1<<u
  h=hp(r)
  for m in range(1,1<<n):
   ck(h>=hp(sub(r,[v for v in range(n) if m>>v&1])),'induced bound');mono+=1
  if n>=3 and strong(r):
   qs+=1
   for sizes in product((1,2),repeat=n):
    b=blow(r,sizes);hb=hp(b)
    ck(strong(b),'blowup strong')
    ck(hb>=max(2**s+1 for s in sizes),'module bound');blows+=1
for a in range(1,7):
 r=blow([2,4,1],(a,1,1))
 brute=sum(all(r[p[i]]>>p[i+1]&1 for i in range(len(r)-1)) for p in permutations(range(len(r))))
 ck(brute==2**a+1,'triangle permutation')
finite=0
for K in range(3,31):
 for L in range(1,21):
  exact_sum=sum(j*2**(j*(j-1)//2)*L**(j+2) for j in range(3,K+1))
  ck(exact_sum<=K*K*2**(K*(K-1)//2)*L**(K+2),'template summation bound');finite+=1
for n in range(1,13):
 vals={19}
 for j in range(n):vals={3*v+d for v in vals for d in (2,4)}
 ck(len(vals)==2**n,'branch injection');ck(max(vals)==3**n*21-2,'branch max')
ck(hp(blow([2,0],(7,7)))==1,'nonstrong hostile')
print(json.dumps({'status':'PASS','induced_controls':mono,'strong_quotients':qs,'blowup_controls':blows,'triangle_permutation_sizes':list(range(1,7)),'finite_summation_bounds':finite,'branch_depths':12,'method':'independent forward subset DP; explicit non-removable checks'}))
