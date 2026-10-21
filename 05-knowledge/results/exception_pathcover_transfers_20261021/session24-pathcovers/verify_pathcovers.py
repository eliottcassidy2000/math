import itertools, math, json, argparse, hashlib
from collections import defaultdict

def need(ok,msg):
 if not ok: raise RuntimeError(msg)

def triangle(c):
 m=len(c)-1; L=m+1; size=L**3
 old=[0]*size;old[0]=1;out=[0]*(3*m+1)
 shifts=[(1,0,0,1),(0,1,0,1),(0,0,1,1),(1,1,0,1),(0,1,1,1),(1,0,1,1),(1,1,1,3)]
 for k in range(1,3*m+1):
  new=[0]*size;total=0
  for i in range(L):
   for j in range(L):
    for h in range(L):
     if i+j+h<k:continue
     idx=(i*L+j)*L+h
     v=new[idx-L*L-L-1] if i and j and h else 0
     for a,b,d,w in shifts:
      if i>=a and j>=b and h>=d:v+=w*old[idx-a*L*L-b*L-d]
     new[idx]=v
     if i and j and h:total+=v*c[i]*c[j]*c[h]
  out[k]=total;old=new
 need(out[-1]==math.factorial(3*m),'factorial endpoint')
 return out

def htriangle(c):
 a=c+[0]
 return 3*sum(a[r]**3+a[r+1]*a[r]**2+a[r+1]**2*a[r] for r in range(1,len(c)))

def fixed_triangle(c,d):
 m=len(c)-1; S=(m-1)//2;old={(0,0):1};out=[0]*(3*m+1)
 # multiply initial V=(1+Y+Y²)/(1-XY²)
 old={(s,2*s+t):1 for s in range(S+1) for t in range(3) if 2*s+t<=m}
 shifts=[(1,0,1),(0,1,2),(1,1,2),(0,2,1),(1,2,3)]
 for l in range((3*m-1)//2+1):
  out[2*l+1]=sum(v*d[2*s+1]*c[r] for (s,r),v in old.items())
  new={}
  for s in range(S+1):
   for r in range(m+1):
    v=new.get((s-1,r-2),0)
    for a,b,w in shifts:v+=w*old.get((s-a,r-b),0)
    if v:new[s,r]=v
  old=new
 return out

def adj_triangle(B):
 m=len(B);return [[int(((j//m-i//m)%3==1) if i//m!=j//m else B[i%m][j%m]) for j in range(3*m)] for i in range(3*m)]

def subset_c(A):
 n=len(A);N=1<<n;hp=[0]*N
 dp=[[0]*n for _ in range(N)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,N):
  hp[mask]=sum(dp[mask])
  for v in range(n):
   if dp[mask][v]:
    for w in range(n):
     if not(mask>>w&1) and A[v][w]:dp[mask|1<<w][w]+=dp[mask][v]
 # ordered concatenation convolution; all possible first nonempty paths
 prev=[0]*N;prev[0]=1;out=[0]
 for k in range(1,n+1):
  cur=[0]*N
  for mask in range(1,N):
   sub=mask
   while sub:
    cur[mask]+=hp[sub]*prev[mask^sub];sub=(sub-1)&mask
  out.append(cur[-1]);prev=cur
 return out

def brute_fixed(A,tau):
 n=len(A);out=[0]*(n+1)
 # symmetric vertex permutations only: choose order and direction of tau-pairs
 center=next(i for i in range(n) if tau[i]==i);pairs=[(i,tau[i]) for i in range(n) if i<tau[i]]
 for pairorder in itertools.permutations(pairs):
  for signs in itertools.product(range(2),repeat=len(pairs)):
   half=[p[b] for p,b in zip(pairorder,signs)];p=half+[center]+[tau[v] for v in half[::-1]]
   # cuts must be mirror symmetric; all bad adjacency positions forced to cut
   for mask in range(1<<len(pairs)):
    cuts={i for i in range(len(pairs)) if mask>>i&1};cuts |= {n-2-i for i in list(cuts)}
    if all(i in cuts or A[p[i]][p[i+1]] for i in range(n-1)):out[len(cuts)+1]+=1
 return out

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',default='receipt.json');args=ap.parse_args()
 C=[0,1];D=[0,1];A=[[0]];tau=[0];rows=[]
 for k in range(1,5):
  n=3**k;H=htriangle(C);F=fixed_triangle(C,D);a=3**((n-1)//2);e=3**((n-1-2*k)//4)
  need(H%a==0 and F[1]%e==0,'orbit integrality');h=H//a;b=F[1]//e;need((h-b)%2==0,'paired parity')
  row=dict(n=n,H=H,d1=F[1],aut=a,centralizer=e,h=h,b=b,paired_black_degree=(h-b)//2)
  if k<=4:
   CC=triangle(C);need(CC[1]==H,'Hamilton independent formula');row['c']=CC;row['d']=F
  if k<=2:
   m=len(A);A=adj_triangle(A);tau=[(0 if i//m==0 else 3-i//m)*m+tau[i%m] for i in range(3*m)]
   need(CC==subset_c(A),'independent subset cover');need(F==brute_fixed(A,tau),'independent symmetric permutations and cuts')
  rows.append(row);print(n,H,F[1],flush=True)
  if k<=4:C,D=CC,F
 # H243 only needs c81; fixed243 requires d81
 n=243;k=5;H=htriangle(C);F=fixed_triangle(C,D);a=3**121;e=3**58
 need(H%a==0 and F[1]%e==0,'243 orbit integrality');h=H//a;b=F[1]//e;need((h-b)%2==0,'243 parity')
 rows.append(dict(n=n,H=H,d1=F[1],aut=a,centralizer=e,h=h,b=b,paired_black_degree=(h-b)//2))
 result=dict(status='FINITE-EXACT',controls='full c subset DP and full d symmetric permutation/cut enumeration at n=3,9; Hamilton formula agrees with coefficient transfer n=3,9,27,81; factorial endpoints n=3,9,27,81',rows=rows,script_sha256=hashlib.sha256(open(__file__,'rb').read()).hexdigest())
 
 with open(args.output,'w') as f:json.dump(result,f,indent=2)
if __name__=='__main__':main()
