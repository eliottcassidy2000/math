import contextlib,io,json
with contextlib.redirect_stdout(io.StringIO()):
 from independent import graph,cover,subset_count

def extend(G,m):
 n=len(G);N=n+m+1;R=[[False]*N for _ in range(N)]
 for i in range(n):
  for j in range(n):R[i][j]=G[i][j]
  for j in range(n,n+m):R[i][j]=True
  R[N-1][i]=True
 for i in range(n,n+m):
  for j in range(i+1,N):R[i][j]=True
 return R
results=[]
for n,As in [(5,({0,1},{0,2})),(7,({0,2},{0,1}))]:
 Gs=[graph(n,A) for A in As];cs=[cover(G) for G in Gs]
 for m in range(4):
  vals=[subset_count(extend(G,m)) for G in Gs]
  expected=[c[1] if m==0 else (2**m+1)*c[1]+2*(2**m-1)*c[2] for c in cs]
  assert vals==expected,(n,m,vals,expected)
  results.append({'core_order':n,'tail_length':m,'counts':vals,'difference':vals[1]-vals[0]})
print(json.dumps({'status':'PASS','consumer_checks':results},indent=2))
