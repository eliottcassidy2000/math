def require(condition, message):
    if not condition:
        raise RuntimeError(message)

import itertools,json

def H(n,arcs):
 adj=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
 for x,y in arcs:adj[x]^=1<<y;adj[y]^=1<<x
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   bits=adj[i]&~mask
   while bits:
    bit=bits&-bits;bits-=bit;j=bit.bit_length()-1
    dp[mask|bit][j]+=dp[mask][i]
 return sum(dp[-1]) if n else 1
count=0
for n in range(1,6):
 pairs=list(itertools.combinations(range(n),2))
 for mask in range(1<<len(pairs)):
  arcs=[p for j,p in enumerate(pairs) if mask>>j&1];r=len(arcs);h=H(n,arcs)
  require(h%2==1, 'check failed: h%2==1')
  covered=sorted({v for x,y in arcs for v in range(x,y+1)})
  ren={v:j for j,v in enumerate(covered)}
  require(H(len(covered),[(ren[x],ren[y]) for x,y in arcs])==h, 'check failed: H(len(covered),[(ren[x],ren[y]) for x,y in arcs])==h')
  for x,y in arcs:
   s=y-x
   if s>=r+1:require(h>=1+2**(s-r), 'check failed: h>=1+2**(s-r)')
  count+=1
for t in range(1,11):require(H(t+2,[(0,t+1)])==1+2**t, 'check failed: H(t+2,[(0,t+1)])==1+2**t')
print(json.dumps({'all_labelled_orders_1_to_5':count,'kernel_deletion_and_span_bound':'PASS','single_arc_interior_1_to_10':'PASS'}))
