import itertools,json
V=[(1,0),(1,1),(0,1)]
I=(1,0,0,1);T=(1,1,1,0);R=(1,1,0,1)
def mul(A,B,mod=2):return tuple(sum(A[2*i+k]*B[2*k+j] for k in range(2))%mod for i in range(2) for j in range(2))
def act(A,v):return tuple(sum(A[2*i+j]*v[j] for j in range(2))%2 for i in range(2))
def matrix(p):return V[p[0]]+V[p[2]] # columns, transpose below

def mat(p):
 c1=V[p[0]];c2=V[p[2]];return(c1[0],c2[0],c1[1],c2[1])
perms=list(itertools.permutations(range(3)));mats={p:mat(p) for p in perms};checks=0
for p in perms:
 A=mats[p]
 if any(act(A,V[x])!=V[p[x]] for x in range(3)):raise ValueError('encoding')
 for q in perms:
  pq=tuple(p[q[x]] for x in range(3))
  if mul(A,mats[q])!=mats[pq]:raise ValueError('composition')
  checks+=1
if mats[(1,2,0)]!=T or mats[(0,2,1)]!=R:raise ValueError('generators')

def rank(rows):
 piv={}
 for x in rows:
  while x:
   j=x.bit_length()-1
   if j in piv:x^=piv[j]
   else:piv[j]=x;break
 return len(piv)
for n,edges in [(3,[(0,1),(1,2),(0,2)]),(4,[(0,1),(1,2),(2,3),(0,3)])]:
 for ps in itertools.product(perms,repeat=len(edges)):
  rows=[]
  for (u,v),p in zip(edges,ps):
   M=mats[p]
   for a in range(2):
    row=1<<(2*v+a)
    for b in range(2):
     if M[2*a+b]:row^=1<<(2*u+b)
    rows.append(row)
  d=2*n-rank(rows)
  owners=sum(all(p[x[u]]==x[v] for (u,v),p in zip(edges,ps)) for x in itertools.product(range(3),repeat=n))
  if d>2 or owners!=2**d-1:raise ValueError('nonzero sections')
  checks+=1
allmat=list(itertools.product(range(2),repeat=4));inv=[A for A in allmat if (A[0]*A[3]-A[1]*A[2])%2];rot=[I,T,mul(T,T)]
subgroups=[]
for bits in range(1<<6):
 H=[mats[p] for i,p in enumerate(perms) if bits>>i&1]
 if I not in H or any(mul(a,b) not in H for a in H for b in H):continue
 fixed=sum(all(act(h,v)==v for h in H) for v in V)
 comm=[A for A in allmat if all(mul(A,h)==mul(h,A) for h in H)]
 subgroups.append(dict(order=len(H),owner_sections=fixed,commutant_dimension=len(comm).bit_length()-1,invertible_commutants=sum(A in inv for A in comm),allowed_rotation_commutants=sum(A in rot for A in comm)))
clock=[]
for s in range(1,13):
 mod=2**s;A=I
 for k in range(1,3*2**(s-1)+1):
  A=mul(A,T,mod)
  if A==I:break
 if k!=3*2**(s-1):raise ValueError('dyadic Fibonacci matrix clock')
 clock.append(dict(s=s,order=k))
print(json.dumps(dict(status='PASS',checks=checks,subgroups=subgroups,dyadic_clocks=clock,scope='exact 3-owner permutation carrier; zero vector excluded and internal chart restrictions retained'),indent=2))
