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
checks=0;roots=0;profiles={}
for n in range(1,5):
 for A in matrices(n):
  for v in range(n):
   c=coeff(A,v);R=(n+1)//2
   assert all(z==0 for z in c[R+1:]);assert c[1]==H(A)
   triangle=any(A[v]>>i&1 and A[i]>>j&1 and A[j]>>v&1 for i in range(n) for j in range(n))
   assert (len(c)>2 and c[2]>0)==bool(triangle)
   hs=[]
   for m in range(1,6):
    z=H(inflate(A,v,m));assert z==sum(c[r]*surj(m,r) for r in range(1,R+1));hs.append(z);checks+=1
   recovered=[0]*(R+1)
   for r in range(1,R+1):
    rem=hs[r-1]-sum(recovered[j]*surj(r,j) for j in range(1,r))
    assert rem%math.factorial(r)==0;recovered[r]=rem//math.factorial(r)
   assert recovered==c[:R+1]
   roots+=1;profiles.setdefault(tuple(c[1:R+1]),{'A':A,'v':v,'H_m_1_through_5':hs})
print(json.dumps({'rooted_sources':roots,'literal_inflation_checks':checks,'profiles':[{'c':p,**data} for p,data in profiles.items()]},indent=2))

for A, base in (([24,1,3,6,14],11),([18,8,3,5,14],13)):
 c=coeff(A,1)
 if c[1:4] != [base,6,0]: raise ValueError(c)
 for x in range(5):
  S=1<<x
  for _ in range(5):
   for i in range(5):
    if S>>i&1:S|=A[i]
  if S!=31:raise ValueError('not strong')
 for m in range(1,7):
  if H(inflate(A,1,m))!=6*(2**m)+base-12:raise ValueError((base,m))
print('PASS: literal strong reader pair 6*2^m +/- 1 through m=6')
