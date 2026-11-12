import contextlib,io,json
with contextlib.redirect_stdout(io.StringIO()):
 from independent import subset_count

def mat(rows):return [[bool(r>>j&1) for j in range(len(rows))] for r in rows]
def append(G,source):
 n=len(G);R=[row+[False] for row in G]+[[False]*(n+1)]
 for i in range(n):
  # source seam root beats new vertex; new beats all other old vertices.
  new_to_old=(i!=0) if source else (i==0)
  R[n][i]=new_to_old;R[i][n]=not new_to_old
 return R
def fans(G,m,n):
 for _ in range(m):G=append(G,True)
 for _ in range(n):G=append(G,False)
 return G
Gs=[mat(x) for x in [(52,13,24,49,34,6),(52,13,24,33,42,6)]]
p=(0,2,4,1,3,5)
assert all(G[a][b] for G in Gs for a,b in zip(p,p[1:]))
checks=[]
for m in range(4):
 vals=[subset_count(fans(G,m,0)) for G in Gs]
 assert vals[1]-vals[0]==(2 if m==0 else 2**(m+1)-2)
 checks.append([m,0,*vals])
for m in (1,2):
 for n in (1,2):
  vals=[subset_count(fans(G,m,n)) for G in Gs]
  assert vals[1]-vals[0]==-16*2**(m+n)+18*2**m+12*2**n-16
  checks.append([m,n,*vals])
print(json.dumps({'status':'PASS','common_frame':p,'source_and_bifan_counts':checks},indent=2))
