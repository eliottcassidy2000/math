import itertools,collections,json
# mask bit i means apex i -> chain vertex

def make(k,arcs):
 states=list(itertools.product(range(3),repeat=k+1)); ix={s:i for i,s in enumerate(states)}
 perms=list(itertools.permutations(range(k))); transitions=[]; ends=[]
 for p in perms:
  tr=[]
  for mask in range(1<<k):
   edges=[]
   for si,s in enumerate(states):
    dest=[]
    for g in range(k+1):
     if s[g]==0 and g>0 and not(mask>>p[g-1]&1):continue
     t=list(s);t[g]=2 if g==k or not(mask>>p[g]&1) else 1
     dest.append(ix[tuple(t)])
    edges.append(dest)
   tr.append(edges)
  transitions.append(tr)
  ends.append([i for i,s in enumerate(states) if all(x==2 or x==0 and (g==0 or g==k or arcs[p[g-1]][p[g]]) for g,x in enumerate(s))])
 init=tuple(tuple(int(i==0) for i in range(len(states))) for p in perms)
 def step(v,m):
  out=[]
  for pi,row in enumerate(v):
   nxt=[0]*len(states)
   for i,c in enumerate(row):
    if c:
     for j in transitions[pi][m][i]:nxt[j]+=c
   out.append(tuple(nxt))
  return tuple(out)
 def val(v):return sum(sum(row[i] for i in ends[p]) for p,row in enumerate(v))
 return init,step,val

def hp(word,k=2):
 n=len(word)+k
 if n==0:return 1
 dp=[{} for _ in range(1<<n)]
 def edge(i,j):
  if i<k and j<k:return i<j
  if i<k:return word[j-k]>>i&1
  if j<k:return not(word[i-k]>>j&1)
  return i<j
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,1<<n):
  for i,c in dp[s].items():
   for j in range(n):
    if not(s>>j&1) and edge(i,j):dp[s|1<<j][j]=dp[s|1<<j].get(j,0)+c
 return sum(dp[-1].values())
if __name__=='__main__':
 init,step,val=make(2,[[0,1],[0,0]])
 tested=0
 for n in range(5):
  for w in itertools.product(range(4),repeat=n):
   v=init
   for m in w:v=step(v,m)
   if val(v)!=hp(w):raise RuntimeError((w,val(v),hp(w)))
   tested+=1
 print('Independent subset-DP controls passed:',tested)
