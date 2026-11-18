from itertools import combinations
import json

def require(x):
 if not x:raise AssertionError

def subcounts(out):
 n=len(out); a=[0]*(1<<n); d=[[0]*n for _ in a]; a[0]=1
 for i in range(n):d[1<<i][i]=1
 for m in range(1,1<<n):
  for v in range(n):
   if m>>v&1:
    for u in range(n):
     if m>>u&1 and out[u]>>v&1:d[m][v]+=d[m^(1<<v)][u]
  a[m]=sum(d[m])
 return a

def profile(out):
 a=subcounts(out); full=len(a)-1
 return a[-1],sum(a[m]*a[full^m] for m in range(full+1))

def build(A,C):
 r,s=len(A),len(C); n=r+s+2; b=r; x=n-1; out=[0]*n
 for i in range(r):out[i]=A[i] | sum(1<<j for j in range(r,r+s+1))
 out[b]=sum(1<<j for j in range(r+1,n))
 for i in range(s):out[r+1+i]=C[i]<<(r+1)
 out[x]=sum(1<<j for j in list(range(r))+list(range(r+1,r+s+1)))
 for i,j in combinations(range(n),2):require(((out[i]>>j)&1)+((out[j]>>i)&1)==1)
 return out

def tt(n):return [(1<<i)-1 for i in range(n)]
P={}; total=0
for n in range(1,6):
 pairs=list(combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  out=[0]*n
  for k,(i,j) in enumerate(pairs):
   if bits>>k&1:out[i]|=1<<j
   else:out[j]|=1<<i
  hp=profile(out);total+=1
  require(n==1 or hp[1]%4==0)
  P.setdefault(hp,out)
checked=0
for (h,p),A in P.items():
 for (k,q),C in P.items():
  require((abs(k*p-h*q)==2)==(sorted([len(A),len(C)])==[1,2]))
  if len(A)+len(C)>8:continue
  left,right=build(A,C),build(C,A)
  require(profile(left)[0]==k*(p+h));require(profile(right)[0]==h*(q+k));checked+=1
for r in range(1,8):
 for s in range(1,9-r):require(profile(build(tt(r),tt(s)))[0]==2**r+1)
print(json.dumps({'labeled_tournaments':total,'profiles':len(P),'literal_pair_tests':checked,'count3_rows':build(tt(1),tt(2)),'count5_rows':build(tt(2),tt(1)),'status':'PASS'},indent=2))
