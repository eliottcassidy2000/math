import itertools,collections,json,math

def code(A):
 return sum(((A[i]>>j)&1)<<k for k,(i,j) in enumerate(itertools.combinations(range(len(A)),2)))
def relabel(A,p):
 return [sum(((A[p[i]]>>p[j])&1)<<j for j in range(len(A))) for i in range(len(A))]
def switch(A,s):
 n=len(A)
 return [sum((bool(A[i]>>j&1)^bool(s[i]^s[j]))<<j for j in range(n) if i!=j) for i in range(n)]
def norm(A):
 s=[0]
 for i in range(len(A)-1):s.append(s[-1]^(1-(A[i]>>(i+1)&1)))
 return switch(A,s)
def H(A):
 n=len(A);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1<<n):
  for i in range(n):
   for j in range(n):
    if not mask>>j&1 and A[i]>>j&1:dp[mask|1<<j][j]+=dp[mask][i]
 return sum(dp[-1])
def tilings(n):
 free=[(i,j) for i in range(n) for j in range(i+2,n)]
 for b in range(1<<len(free)):
  A=[0]*n
  for i in range(n-1):A[i]|=1<<(i+1)
  for k,(i,j) in enumerate(free):
   u,v=(i,j) if b>>k&1 else (j,i);A[u]|=1<<v
  yield A
out={}
for n in range(2,6):
 perms=list(itertools.permutations(range(n)));X=list(tilings(n));data=[]
 for A in X:
  rel=[relabel(A,p) for p in perms]
  canon=min(map(code,rel));scanon=min(code(norm(B)) for B in rel)
  aut=sum(B==A for B in rel);gs=sum(norm(B)==A for B in rel)
  data.append((canon,scanon,aut,gs,H(A)))
 sheets=collections.defaultdict(list)
 for idx,d in enumerate(data):sheets[d[1]].append(idx)
 sheetkeys=sorted(sheets);flow=[[0]*len(sheets) for _ in sheets]
 for i,d in enumerate(data):
  si=sheetkeys.index(d[1]);sj=sheetkeys.index(data[len(X)-1-i][1]);flow[si][sj]+=1
 blueflow=[[0]*len(sheets) for _ in sheets]
 for i,A in enumerate(X):
  R=relabel([sum((not(A[u]>>v&1))<<v for v in range(n) if u!=v) for u in range(n)],list(reversed(range(n))))
  if R==A:
   si=sheetkeys.index(data[i][1]);sj=sheetkeys.index(data[len(X)-1-i][1]);blueflow[si][sj]+=1
 rows=[]
 for sk in sheetkeys:
  inds=sheets[sk];classes=collections.Counter(data[i][0] for i in inds)
  rec=[]
  for c,w in classes.items():
   d=next(data[i] for i in inds if data[i][0]==c)
   assert w*d[2]==d[4]
   rec.append({'canonical':c,'H':d[4],'Aut':d[2],'tiling_weight':w})
  gs=data[inds[0]][3]
  A=X[inds[0]]; cuts=collections.Counter()
  for b in range(1<<(n-1)):
   B=switch(A,[(b>>i)&1 for i in range(n)])
   cuts[min(code(relabel(B,p)) for p in perms)]+=1
  for r in rec:assert cuts[r['canonical']]*r['Aut']==gs
  assert len(inds)*gs==math.factorial(n)
  if n%2:
   signs=[(n-1-A[i].bit_count())%2 for i in range(n)]
   C=switch(A,signs)
   assert len({(n-1-C[i].bit_count())%2 for i in range(n)})==1
   assert sum(relabel(C,p)==C for p in perms)==gs
   for b in range(1<<(n-1)):
    B=switch(A,[(b>>i)&1 for i in range(n)])
    bs=[(n-1-B[i].bit_count())%2 for i in range(n)]
    assert switch(B,bs)==C
  rows.append({'canonical_switch':sk,'G_switch':gs,'sheet_size':len(inds),'classes':rec})
 assert all(sum(flow[i])==len(sheets[k]) for i,k in enumerate(sheetkeys))
 assert all(flow[i][j]==flow[j][i] for i in range(len(sheets)) for j in range(len(sheets)))
 out[n]={'sheets':rows,'directed_complement_flow':flow,'directed_blue_flow':blueflow}
print(json.dumps(out,indent=2))
