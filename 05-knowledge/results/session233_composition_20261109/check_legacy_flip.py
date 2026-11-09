"""Exact regression for the two repaired THM-250 prose errors."""
import itertools,json,sys
if '--fail' in sys.argv: raise RuntimeError('deliberate legacy-check sentinel')
def count(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for s in range(1,1<<n):
  for i,v in enumerate(d[s]):
   rem=a[i]&~s
   while rem:
    b=rem&-rem;rem-=b;d[s|b][b.bit_length()-1]+=v
 return sum(d[-1])
cases=0
for n in range(2,10):
 vals=set()
 for a,b in itertools.combinations(range(n),2):
  adj=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
  adj[a]^=1<<b;adj[b]|=1<<a
  h=count(adj);target=1 if b==a+1 else 1+(1<<(b-a-1))
  if h!=target:raise RuntimeError((n,a,b,h,target))
  if b==a+1:
   p=list(range(n));p[a],p[b]=p[b],p[a]
   if not all(adj[p[i]]>>p[j]&1 for i in range(n) for j in range(i+1,n)):raise RuntimeError('new order not transitive')
  vals.add(h);cases+=1
 if 2 in vals or len(vals)!=n-1:raise RuntimeError(('corollary',n,vals))
print(json.dumps({'status':'PASS','orders':[2,9],'single_flip_cases':cases,'corrected_corollary':True},indent=2))
