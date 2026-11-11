import itertools,json,time
from collections import defaultdict

def hcount(rows):
 n=len(rows); dp=[[0]*n for _ in range(1<<n)]
 for i in range(n): dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i,v in enumerate(dp[mask]):
   if v:
    rem=rows[i]&~mask
    while rem:
     b=rem&-rem; rem-=b; dp[mask|b][b.bit_length()-1]+=v
 return sum(dp[-1])
def compute(n):
 edges=[(i,j) for i in range(n) for j in range(i+2,n)]; m=len(edges)
 vals=[]; canon=[]
 for mask in range(1<<m):
  rows=[1<<(i+1) if i<n-1 else 0 for i in range(n)]
  for k,(i,j) in enumerate(edges):
   if mask>>k&1: rows[j]|=1<<i
   else: rows[i]|=1<<j
  vals.append(hcount(rows))
  groups=defaultdict(list)
  for i,r in enumerate(rows): groups[r.bit_count()].append(i)
  blocks=[groups[d] for d in sorted(groups)]
  best=None
  for parts in itertools.product(*(itertools.permutations(b) for b in blocks)):
   perm=sum(parts,())
   code=0
   for a in range(n):
    for b in range(a+1,n): code=(code<<1)|((rows[perm[a]]>>perm[b])&1)
   if best is None or code<best:best=code
  canon.append(best)
 cuts=[]
 for i in range(n): cuts.append(sum(1<<k for k,e in enumerate(edges) if i in e))
 allcuts={0}
 for d in cuts: allcuts |= {x^d for x in list(allcuts)}
 tiles=[1<<k for k in range(m)]
 twotiles={x^y for x in tiles for y in tiles}|set(tiles)
 result={}
 for name,moves in [('tile',[1<<k for k in range(m)]),('gauge',list(set(cuts)-{0})),('both',list(set(cuts+tiles)-{0})),('all_gauges_and_tiles',list((allcuts|set(tiles))-{0})),('two_tiles',list(twotiles-{0})),('full',list((allcuts|twotiles|{(1<<m)-1})-{0})),('reframed_tile',tiles),('reframed_two_tiles',list(twotiles-{0})),('reframed_gauge',list(allcuts-{0})),('reframed_full',list((allcuts|twotiles|{(1<<m)-1})-{0}))]:
  parent=list(range(1<<m))
  def root(x):
   while x!=parent[x]:parent[x]=parent[parent[x]];x=parent[x]
   return x
  if name.startswith('reframed_'):
   seen={}
   for x,c in enumerate(canon):
    if c in seen: parent[root(x)]=root(seen[c])
    else: seen[c]=x
  ec={0:0,2:0}
  for x,h in enumerate(vals):
   for d in moves:
    y=x^d
    if y>x and abs(vals[y]-h)<=2:
     ec[abs(vals[y]-h)]+=1
     a,b=root(x),root(y);parent[a]=b
  comp=defaultdict(list)
  for x,h in enumerate(vals):comp[root(x)].append(x)
  summaries=[]
  for xs in comp.values():
   v=sorted(set(vals[x] for x in xs)); assert v==list(range(min(v),max(v)+1,2))
   summaries.append({'size':len(xs),'min':v[0],'max':v[-1],'values':len(v),'representative':xs[0]})
  summaries.sort(key=lambda z:(-z['values'],-z['size'],z['min']))
  result[name]={'components':len(comp),'edges':ec,'largest_intervals':summaries[:12], 'singletons':sum(len(x)==1 for x in comp.values())}
 return {'n':n,'vertices':len(vals),'edges':edges,'isomorphism_classes':len(set(canon)),'spectrum':sorted(set(vals)),'maxH':max(vals),'graphs':result}
if __name__=='__main__':
 out=[compute(n) for n in range(3,8)]
 print(json.dumps(out,indent=2))
