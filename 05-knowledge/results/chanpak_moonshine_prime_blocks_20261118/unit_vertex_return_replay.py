import itertools,json
def matrices(n):
 es=list(itertools.combinations(range(n),2))
 for code in range(1<<len(es)):
  A=[0]*n
  for k,(i,j) in enumerate(es):
   a,b=(i,j) if code>>k&1 else (j,i);A[a]|=1<<b
  yield A

def H(A):
 n=len(A);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for mask in range(1<<n):
  for i in range(n):
   rest=A[i]&~mask
   while rest:
    bit=rest&-rest;rest-=bit;j=bit.bit_length()-1;d[mask|bit][j]+=d[mask][i]
 return sum(d[-1])

out={};tested=0
for n in range(3,7):
 counts={};examples={}
 for T in matrices(n):
  h=H(T)
  for v in range(n):
   ss=[i for i in range(n) if i!=v];B=[sum((T[i]>>j&1)<<k for k,j in enumerate(ss)) for i in ss]
   d=h-H(B);tested+=1
   if d==2:
    assert h in [3,5]
    counts[h]=counts.get(h,0)+1;examples.setdefault(h,{'T':T,'v':v})
 out[n]={'unit_deletions_by_child_H':counts,'examples':examples}
print(json.dumps({'tested':tested,'orders':out},indent=2))
