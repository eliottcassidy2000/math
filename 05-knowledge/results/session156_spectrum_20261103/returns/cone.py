from itertools import combinations

def dpc(rows):
 n=len(rows); dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   for u in range(n):
    if s>>u&1 and not s>>v&1 and rows[u]>>v&1:dp[s|1<<v][v]+=dp[s][u]
 return sum(dp[-1]),dp[-1]
def add_a(rows,p):
 n=len(rows);return tuple(row|((1<<n) if i!=p else 0) for i,row in enumerate(rows))+(1<<p,)
for n in range(1,6):
 pairs=list(combinations(range(n),2))
 for code in range(1<<len(pairs)):
  rows=[0]*n
  for i,(u,v) in enumerate(pairs):
   if code>>i&1:rows[u]|=1<<v
   else:rows[v]|=1<<u
  h,end=dpc(rows)
  for p in range(n):
   a=add_a(rows,p); ha,_=dpc(a); aa=add_a(a,n); haa,_=dpc(aa)
   if (ha-h)%2 or (haa-h)%2: raise AssertionError("parity")
   x=(ha-h)//2; y=(haa-h)//2-x
   if y-x!=end[p]: raise AssertionError((n,code,p,h,x,y,end[p]))
   if not (0<=x<=y and h<=x+y and h+x>=y): raise AssertionError((n,code,p,h,x,y))
 print('n',n,'passed',1<<len(pairs),'tournaments')
