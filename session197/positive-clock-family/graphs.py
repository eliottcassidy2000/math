from itertools import permutations,combinations
SEED=[2,12,145,149,99,143,175,19]
ROOT=7
PORTS=[0,1,4]
def require(x,s):
 if not x:raise RuntimeError(s)
def hp(rows):
 n=len(rows);N=1<<n;D=[[0]*n for _ in range(N)]
 for v in range(n):D[1<<v][v]=1
 for mask in range(1,N):
  for v in range(n):
   if mask>>v&1:
    rem=mask^(1<<v)
    if rem:D[mask][v]=sum(D[rem][u] for u in range(n) if rem>>u&1 and rows[u]>>v&1)
 return sum(D[-1])
def clone(rows,r,m):
 old=[v for v in range(len(rows)) if v!=r];n=len(old);out=[0]*(n+m)
 for i,v in enumerate(old):
  for j,w in enumerate(old):
   if rows[v]>>w&1:out[i]|=1<<j
 for j in range(m):
  z=n+j
  for i,v in enumerate(old):
   if rows[r]>>v&1:out[z]|=1<<i
   else:out[i]|=1<<z
  for k in range(j+1,m):out[z]|=1<<(n+k)
 return out

def literal_runs(rows,r):
 old=[v for v in range(len(rows)) if v!=r];a=[0]*(len(rows)+1)
 for w in permutations(old):
  for mask in range(1,1<<len(rows)):
   word=[]
   for i in range(len(rows)):
    if mask>>i&1:word.append(r)
    if i<len(w):word.append(w[i])
   if all(rows[u]>>v&1 for u,v in zip(word,word[1:])):a[mask.bit_count()]+=1
 return a

def run_dp(rows,r):
 old=[v for v in range(len(rows)) if v!=r];n=len(old);N=1<<n;K=n+2
 D=[[[0]*K for _ in range(n+1)] for _ in range(N)]
 D[0][n][1]=1
 for i in range(n):D[1<<i][i][0]=1
 for mask in range(N):
  for k in range(K):
   for i,v in enumerate(old):
    x=D[mask][i][k]
    if not x:continue
    for j,w in enumerate(old):
     if not mask>>j&1 and rows[v]>>w&1:D[mask|1<<j][j][k]+=x
    if k+1<K and rows[v]>>r&1:D[mask][n][k+1]+=x
   x=D[mask][n][k]
   if x:
    for j,w in enumerate(old):
     if not mask>>j&1 and rows[r]>>w&1:D[mask|1<<j][j][k]+=x
 return [sum(D[-1][i][k] for i in range(n+1)) for k in range(K)]

def birth(rows,ports,side,mask,drop):
 n=len(rows);special={ports[i] for i in range(len(ports)) if mask>>i&1};out=rows[:]+[0]
 for v in range(n):
  new_wins=(v in special) if side=='R' else (v not in special)
  if new_wins:out[n]|=1<<v
  else:out[v]|=1<<n
 ps=ports+[n];del ps[drop];return out,ps

def strong(rows):
 for v in range(len(rows)):
  S=1<<v
  while True:
   T=S
   for u in range(len(rows)):
    if S>>u&1:T|=rows[u]
   if T==S:break
   S=T
  if S!=(1<<len(rows))-1:return False
 return True
