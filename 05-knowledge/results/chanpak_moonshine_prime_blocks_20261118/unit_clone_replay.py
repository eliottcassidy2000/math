import itertools,math,json
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
def coeff(A,v):
 n=len(A);c=[0]*(n+1)
 for p in itertools.permutations(i for i in range(n) if i!=v):
  for bits in range(1,1<<n):
   w=[]
   for k in range(n):
    if bits>>k&1:w.append(v)
    if k<n-1:w.append(p[k])
   if all(A[a]>>b&1 for a,b in zip(w,w[1:])):c[bits.bit_count()]+=1
 return c

def inflate(A,v,m):
 labels=[i for i in range(len(A)) if i!=v]+[v]*m;n=len(labels)
 return [sum((i<j if a==b else bool(A[a]>>b&1))<<j for j,b in enumerate(labels) if i!=j) for i,a in enumerate(labels)]
def surj(m,r):return sum((-1)**(r-j)*math.comb(r,j)*j**m for j in range(r+1))

rows={};examples=[]
for n in range(1,6):
 tested=ready=0
 for A in matrices(n):
  h=H(A)
  for v in range(n):
   tested+=1;B=inflate(A,v,2);hb=H(B);delta=hb-h
   if delta==2:
    ready+=1;assert h==3 and hb==5
    assert all(H(inflate(B,w,2))-hb!=2 for w in range(n+1))
    outside=[i for i in range(n) if i!=v]
    tris=[(i,j) for i in outside for j in outside if A[v]>>i&1 and A[i]>>j&1 and A[j]>>v&1]
    assert len(tris)==1
    i,j=tris[0];C={v,i,j};I=[k for k in range(n) if k not in C and A[k]>>v&1];O=[k for k in range(n) if k not in C and A[v]>>k&1]
    assert all(A[k]>>z&1 for k in I for z in C|set(O))
    assert all(A[z]>>k&1 for k in O for z in C)
    def induced(S):return [sum(((A[a]>>b)&1)<<l for l,b in enumerate(S)) for a in S]
    assert H(induced(I))==1 if I else True
    assert H(induced(O))==1 if O else True
 rows[n]={'rooted_sources_checked':tested,'unit_clones':ready,'all_unit_children_H5_without_ready_root':True}
transport=0
for n in range(1,5):
 for A in matrices(n):
  for v in range(n):
   c=coeff(A,v)+[0,0];B=inflate(A,v,2)
   for w in [n-1,n]:
    cp=coeff(B,w)
    assert all(cp[j]==j*c[j]+(j+1)*c[j+1] for j in range(1,n+2))
    transport+=1
rows['marked_transports']=transport
print(json.dumps(rows,indent=2))
