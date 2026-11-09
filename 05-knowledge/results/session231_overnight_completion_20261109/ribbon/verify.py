from itertools import combinations
import json

def require(ok, message):
 if not ok: raise RuntimeError(message)

def h(rows):
 n=len(rows); d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for m in range(1,1<<n):
  for v in range(n):
   if d[m][v]:
    rem=rows[v]&~m
    while rem:
     b=rem&-rem;rem-=b;w=b.bit_length()-1
     d[m|b][w]+=d[m][v]
 return sum(d[-1]) if n else 1

def make(n,arcs):
 rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
 for x,y in arcs:rows[x]^=1<<y;rows[y]^=1<<x
 return rows

def induced(rows,keep):
 return [sum(((rows[v]>>w)&1)<<j for j,w in enumerate(keep)) for v in keep]

checks=0
for n in range(2,8):
 pairs=list(combinations(range(n),2))
 for r in range(0,min(3,len(pairs))+1):
  for arcs in combinations(pairs,r):
   rows=make(n,arcs);H=h(rows)
   keep=sorted({i for x,y in arcs for i in range(x,y+1)})
   require(h(induced(rows,keep))==H,"neutral compression")
   if r:
    L=(H-1).bit_length()-1 if H>1 else -1
    # H=1 can include backward arcs of span one; bound uses X>=3.
    L=max(1,L)
    require(len(keep)<=r*(L+r+1),"kernel size")
   for x,y in arcs:
    removed={next(v for v in edge if v not in (x,y)) for edge in arcs if edge!=(x,y)}
    local=[v for v in range(x,y+1) if v not in removed]
    single=h(induced(rows,local))
    # Span-one arc alone is transitive; formula applies with interior nonempty.
    if len(local)>=3:require(single==1+2**(len(local)-2),"single arc formula")
    require(single<=H,"deletion monotonicity")
    if y-x>=r and y-x-r>=1:require(H>=1+2**(y-x-r),"span bound")
   checks+=1
print(json.dumps({'tournaments_tested':checks,'orders':'2..7','back_arcs':'0..3','checks':'deletion monotonicity, neutral compression, single-arc induced control'},indent=2))
