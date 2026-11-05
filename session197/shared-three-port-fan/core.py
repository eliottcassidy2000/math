from itertools import permutations

def require(v,msg):
 if not v:raise RuntimeError(msg)

def dp(rows):
 n=len(rows);N=1<<n;E=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for v in range(n):E[1<<v][v]=1
 for mask in range(1,N):
  for v in range(n):
   if mask>>v&1:
    rem=mask^(1<<v)
    if rem:E[mask][v]=sum(E[rem][u] for u in range(n) if rem>>u&1 and rows[u]>>v&1)
  H[mask]=sum(E[mask])
 return H,E

def profile(rows,r):
 n=len(rows);F=(1<<n)-1;H,E=dp(rows);rev=[F^(1<<v)^rows[v] for v in range(n)];_,S=dp(rev)
 R=L=0
 for mask in range(1<<n):
  if mask>>r&1:R+=S[mask][r]*H[F^mask];L+=E[mask][r]*H[F^mask]
 return H[F],E[F][r],R,H[F^(1<<r)],S[F][r],L

def literal_profile(rows,r):
 n=len(rows);h=e=s=R=L=0
 for word in permutations(range(n)):
  ok=[bool(rows[a]>>b&1) for a,b in zip(word,word[1:])]
  if all(ok):h+=1;e+=word[-1]==r;s+=word[0]==r
  for k in range(n+1):
   if not(all(ok[:max(k-1,0)]) and all(ok[k:])):continue
   if k<n and word[k]==r:R+=1
   if k>0 and word[k-1]==r:L+=1
 exterior=[i for i in range(n) if i!=r]
 d=sum(all(rows[a]>>b&1 for a,b in zip(w,w[1:])) for w in permutations(exterior))
 return h,e,R,d,s,L

def fan(X,a,m):
 require(m>=1,'positive chain multiplicity')
 k=len(X);b=k;out=X[:]+[0]*(m+1)
 for u in range(k):
  for z in range(k+1,k+m+1):out[u]|=1<<z
  if u!=a:out[u]|=1<<b
 out[b]=1<<a
 for z in range(k+1,k+m+1):out[z]|=1<<b;out[z]|=sum(1<<w for w in range(z+1,k+m+1))
 return out,b

def left(rows,ports,S,drop):
 n=len(rows);exceptions={ports[i] for i in range(2) if S>>i&1};out=rows[:]+[0]
 for v in range(n):
  if v in exceptions:out[v]|=1<<n
  else:out[n]|=1<<v
 ps=ports+[n];del ps[drop];return out,ps

def plus2(rows):
 rows,ports=left(rows,[0,1],0,1)
 return left(rows,ports,1,2)

def plus4(rows):return left(rows,[1,4],3,2)

def relabel(rows,mapping):
 out=[0]*len(rows)
 for a,r in enumerate(rows):
  for b in range(len(rows)):
   if r>>b&1:out[mapping[a]]|=1<<mapping[b]
 return out

def strong(rows):
 F=(1<<len(rows))-1
 for r in range(len(rows)):
  seen=1<<r;old=-1
  while seen!=old:
   old=seen
   for v in range(len(rows)):
    if seen>>v&1:seen|=rows[v]
  if seen!=F:return False
 return True
