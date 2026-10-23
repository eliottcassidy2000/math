import math,json

def need(x,m):
 if not x:raise RuntimeError(m)
def fa(n):return 0 if n==0 else (3*n-math.isqrt(5*n*n)-1)//2
# Sign of a+b*alpha, computed by comparison of rational part with sqrt5 part.
def sg(a,b):
 if b==0:return (a>0)-(a<0)
 c=2*a+3*b
 if c==0:return -1 if b>0 else 1
 if c*b<0:return 1 if c>0 else -1
 d=c*c-5*b*b
 return (1 if c>0 else -1)*(1 if d>0 else -1 if d<0 else 0)
Nmax=10000
A=set();B=set()
for k in range(1,Nmax+1):
 x=3*k-fa(k)-1
 if x<=Nmax:A.add(x)
 y=2*((k+math.isqrt(5*k*k))//2)+k
 if y<=Nmax:B.add(y)
need(not(A&B),'Beatty color collision')
C=set(range(1,Nmax+1))-A-B;maxA=A|B
w='C'
while len(w)<Nmax:w=''.join({'C':'CA','A':'B','B':'CA'}[c] for c in w)
need(all(w[i-1]==('A' if i in A else 'B' if i in B else 'C') for i in range(1,Nmax+1)),'independent Beatty code')
checks=0;fail=[];support_miss=[]
for N in range(2,2001):
 count=sum(a in A and N-a in A for a in range(1,N));distinct=sum(a in A and N-a in A for a in range(1,(N+1)//2))
 theta=(-fa(N),N)
 if sg(theta[0]-1,theta[1]+2)<0:L=(0,0)
 elif sg(theta[0]-1,theta[1]+1)<0:L=(theta[0]-1,theta[1]+2)
 else:L=(1-theta[0],-theta[1])
 err=(count-(N-1)*L[0],-(N-1)*L[1]);bd=16*(N-1).bit_length()
 need(sg(err[0]-bd,err[1])<=0 and sg(err[0]+bd,err[1])>=0,'tent discrepancy')
 if L==(0,0):need(count==0,'outside support')
 maxdistinct=sum(a in maxA and N-a in maxA for a in range(1,(N+1)//2))
 if not maxdistinct:fail.append(N)
 need(N<7 or maxdistinct>0,'Amax tail')
 checks+=1
for k in range(2,101):
 N=k*k;d=sum(a in A and N-a in A for a in range(1,(N+1)//2))
 theta=(-fa(N),N);support=sg(theta[0]-1,theta[1]+2)>0
 if support and not d:support_miss.append(k)
 need(not d or support,'distinct outside support')
# Fixed-target choices: independently brute force all compatible colors N<=13.
K={'A':'A','B':'AB','C':'BC'}
import itertools,collections
polytests=0
for N in range(2,14):
 # Histogram of AA-edge counts; target edges form a matching.
 brute=collections.Counter()
 choices=[K[w[i-1]] for i in range(1,N)]
 for selected in itertools.product(*choices):
  r=sum(selected[a-1]=='A' and selected[N-a-1]=='A' for a in range(1,(N+1)//2))
  brute[r]+=1
 pred=collections.Counter({0:1})
 for a in range(1,(N+1)//2):
  b=N-a;total=len(choices[a-1])*len(choices[b-1]);aa=int('A'in choices[a-1] and 'A'in choices[b-1]);nxt=collections.Counter()
  for r,v in pred.items():nxt[r]+=v*(total-aa);nxt[r+1]+=v*aa
  pred=nxt
 if N%2==0:
  for r in pred:pred[r]*=len(choices[N//2-1])
 need(+pred==brute,'fixedtarget polynomial');polytests+=1
print(json.dumps({'status':'PASS','independent_Beatty_prefix':Nmax,'AA_pair_targets':checks,'Amax_failures_through2000':fail,'square_support_but_no_distinct_A_pair_k_through100':support_miss,'fixed_target_AA_polynomial_bruteforce_targets':polytests},indent=2))
