from itertools import combinations

def endpoints(rows):
 n=len(rows); dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,1<<n):
  for i in range(n):
   z=dp[s][i]
   if not z:continue
   rem=rows[i]&~s
   while rem:
    b=rem&-rem;rem-=b;j=b.bit_length()-1;dp[s|b][j]+=z
 return dp

def h(rows):return sum(endpoints(rows)[-1])
def reverse(rows):
 n=len(rows);return tuple(((1<<n)-1)^r^(1<<i) for i,r in enumerate(rows))
def covers(rows):
 n=len(rows);en=endpoints(rows);st=endpoints(reverse(rows));full=(1<<n)-1
 q=[[0]*n for _ in range(n)]
 for s in range(1,full):
  for u,a in enumerate(en[s]):
   if not a:continue
   for v,b in enumerate(st[full^s]):
    if b:q[u][v]+=a*b
 return q

def ext(rows,D):
 n=len(rows);return tuple(r | ((1<<n) if not D>>i&1 else 0) for i,r in enumerate(rows))+(D,)
def induced(rows,vertices):
 return tuple(sum(1<<j for j,v in enumerate(vertices) if rows[u]>>v&1) for u in vertices)
def ordinal(*blocks):
 rows=[];N=sum(map(len,blocks));offset=0
 for b in blocks:
  above=((1<<N)-1)^((1<<(offset+len(b)))-1)
  rows.extend((r<<offset)|above for r in b);offset+=len(b)
 return tuple(rows)
def strong(rows):
 n=len(rows)
 for start in range(n):
  seen=1<<start
  while True:
   x=seen
   for i in range(n):
    if seen>>i&1:x|=rows[i]
   if x==seen:break
   seen=x
  if seen!=(1<<n)-1:return False
 return True
if __name__=='__main__':
 Q=(8,1,3,6)
 for d in [2,3,6,7]:
  a=ext(Q,d);c2=sum(map(sum,covers(a)));print(d,a,h(a),c2,2*h(a)+c2,c2%3)

def inflate(quotient, blocks):
    starts = []; total = 0
    for block in blocks:
        starts.append(total); total += len(block)
    rows = []
    for i, block in enumerate(blocks):
        external = 0
        for j, target in enumerate(blocks):
            if quotient[i] >> j & 1:
                external |= ((1 << len(target)) - 1) << starts[j]
        rows.extend((r << starts[i]) | external for r in block)
    return tuple(rows)
