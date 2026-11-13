from math import gcd
BANK=(
 (5,5,(8,1,3,6)),
 (19,8,(40,1,3,6,15,30)),
 (217,12,(50,12,9,65,174,78,23,111)),
 (223,12,(248,209,155,242,224,134,36,64)),
 (599,13,(216,53,49,198,72,89,134,54)),
 (761,14,(446,84,88,50,160,6,185,302,126)),
 (773,14,(262,488,2,197,239,141,37,69,252)),
)
R=((24,9,35,4,14,27),(108,97,2,6,79,92,12))
def info(A,slots=False):
 n=len(A);N=(1<<n)-1
 D=[[0]*n for _ in range(N+1)];H=[1]+[0]*N
 for j in range(n):D[1<<j][j]=1
 for s in range(1,N+1):
  for j,d in enumerate(D[s]):
   if not d:continue
   H[s]+=d;a=A[j]&~s
   while a:
    b=a&-a;a^=b;D[s|b][b.bit_length()-1]+=d
 return (H[N],sum(H[s]*H[N^s]for s in range(N+1)))if slots else H[N]
def strong(A):
 for i in range(len(A)):
  S=1<<i
  for _ in A:
   for j,a in enumerate(A):
    if S>>j&1:S|=a
  if S!=(1<<len(A))-1:return False
 return True
def tt(m):
 return tuple(((1<<m)-1)^((1<<(i+1))-1) for i in range(m))
def join(A,B):
 return tuple(x|(((1<<len(B))-1)<<len(A))for x in A)+tuple(x<<len(A)for x in B)
def source_A(l,c):return -3+6*2**l+(2+4*c)*3**l
def decode(N,l,c):
 A=source_A(l,c)
 m=(N//A).bit_length()
 if not m:return None
 d=A*2**m-N
 if not (0<d<=2**(m-1)):return None
 z=d;cost=0;ex=[]
 for h,b,_ in BANK:
  e=0
  while z%h==0:z//=h;e+=1;cost+=b
  ex.append(e)
 if z!=1 or cost>m:return None
 return m,d,m-cost,tuple(ex)
def compile_target(N):
 out=[]
 for c in range(2):
  l=1
  while source_A(l,c)<N:
   d=decode(N,l,c)
   if d:out.append((c,l,d))
   l+=1
 return out
def material(c,l,recipe):
 m,d,pad,ex=recipe
 V=()
 for (_,_,A),e in zip(BANK,ex):
  for _ in range(e):V=join(V,A)
 V=join(V,tt(pad))
 B=(tt(l),)+((0,),)*(len(R[c])-2)+(V,)
 offsets=[];pos=0
 for A in B:offsets.append(pos);pos+=len(A)
 out=[]
 for i,A in enumerate(B):
  ext=sum(((1<<len(F))-1)<<offsets[j] for j,F in enumerate(B)if R[c][i]>>j&1)
  out.extend(ext|(x<<offsets[i]) for x in A)
 return tuple(out)
for i,(h,m,A) in enumerate(BANK):
 if not strong(A) or info(A,True)!=(h,2**m) or any(gcd(h,g)!=1 for g,_,_ in BANK[:i]):raise ValueError(A)
for N,c,w in ((3821,0,19),(61223,0,217),(6907,1,5),(29,0,1),(53,1,1)):
 r=decode(N,1,c)
 if not r or r[1]!=w or (c,1,r) not in compile_target(N):raise ValueError((N,r))
 A=material(c,1,r)
 if info(A)!=N or not strong(A):raise ValueError((N,A))
 print(N,'context',c,'recipe',r)
