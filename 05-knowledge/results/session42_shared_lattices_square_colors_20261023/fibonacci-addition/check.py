from math import isqrt,sqrt
from collections import Counter
import json
checks=0
def require(v,msg):
 global checks
 checks+=1
 if not v:raise RuntimeError(msg)
def floor_alpha(n):
 return 0 if n==0 else (3*n-isqrt(5*n*n)-1)//2
# Exact a+b*alpha, alpha=(3-sqrt(5))/2.
def add(x,y):return(x[0]+y[0],x[1]+y[1])
def sub(x,y):return(x[0]-y[0],x[1]-y[1])
def scale(n,x):return(n*x[0],n*x[1])
def sign(x):
 a,b=x;c=2*a+3*b
 if b==0:return(c>0)-(c<0)
 if c>=0 and b<0:return 1
 if c<=0 and b>0:return -1
 z=c*c-5*b*b
 return ((c>0)-(c<0))*((z>0)-(z<0))
def low(x,y):return x if sign(sub(x,y))<=0 else y
def high(x,y):return x if sign(sub(x,y))>=0 else y
I={'B':((0,0),(1,-2)),'C':((1,-2),(1,-1)),'A':((1,-1),(1,0))}
def phase(n):return(-floor_alpha(n),n)
def color(n):
 if floor_alpha(n+1)>floor_alpha(n):return 'A'
 if floor_alpha(n+2)>floor_alpha(n+1):return 'C'
 return 'B'
def from_phase(n):
 t=phase(n)
 for c,(l,u) in I.items():
  if sign(sub(t,l))>=0 and sign(sub(t,u))<0:return c
 raise RuntimeError('no phase color')
def conv(c,d,N):
 theta=phase(N);l,u=I[c];p,q=I[d];ans=(0,0)
 for shift in (-1,0,1):
  lo=high(l,add(sub(theta,q),(shift,0)))
  hi=low(u,add(sub(theta,p),(shift,0)))
  if sign(sub(hi,lo))>0:ans=add(ans,sub(hi,lo))
 return ans
def tent(N):
 t=phase(N)
 return high((0,0),sub(low((1,0),add(t,(0,1))),high((1,-1),t)))
maxn=250**2
colors=['?']+[color(n) for n in range(1,maxn+1)]
for n in range(1,maxn+1):require(colors[n]==from_phase(n),'rotation code')
# Independent substitution prefix.
w='C'
while len(w)<maxn:w=''.join({'C':'CA','A':'B','B':'CA'}[x] for x in w)
require(w[:maxn]==''.join(colors[1:]),'substitution')
for N in range(2,1501):
 counts=Counter((colors[a],colors[N-a])for a in range(1,N))
 for c in I:
  for d in I:
   v=conv(c,d,N);err=sub((counts[c,d],0),scale(N-1,v));bound=16*(N-1).bit_length()
   require(sign(sub(err,(bound,0)))<=0 and sign(add(err,(bound,0)))>=0,'uniform bound')
 require(conv('A','A',N)==tent(N),'tent')
 require(sum(counts.values())==N-1,'partition')
square_stats=[];averages=Counter();positive=0;possible=0;diag_total=0;max_error=(0.0,None)
for k in range(2,251):
 N=k*k;counts=Counter((colors[a],colors[N-a])for a in range(1,N))
 diag=int(N%2==0 and colors[N//2]=='A');distinct=(counts['A','A']-diag)//2
 require(counts['A','A']>=diag and (counts['A','A']-diag)%2==0,'diagonal')
 L=tent(N);support=sign(L)>0
 require(support or counts['A','A']==0,'outside support')
 positive+=distinct>0;possible+=support;diag_total+=diag
 for c in I:
  for d in I:
   v=conv(c,d,N);err=sub((counts[c,d],0),scale(N-1,v));bound=16*(N-1).bit_length()
   require(sign(sub(err,(bound,0)))<=0 and sign(add(err,(bound,0)))>=0,'square uniform bound')
   averages[c+d]+=counts[c,d]/N
 if k in (10,25,50,100,250):square_stats.append(dict(k=k,positive_distinct_A_targets=positive,phase_support_targets=possible,denominator=k-1))
require(colors[2]==colors[7]=='A' and isqrt(2+7)**2==9,'forced edge')
require(all(not(colors[a]==colors[b]=='A' and isqrt(a+b)**2==a+b) for b in range(1,7) for a in range(1,b)),'minimal forced A edge')
# Branch extremes and exact one-target polynomial.
for N in range(2,1501):
 H=N-1;mx=sum(colors[a]!='C' and colors[N-a]!='C' for a in range(1,N))
 cc=sum(colors[a]==colors[N-a]=='C' for a in range(1,N))
 require(mx==H-2*floor_alpha(H+2)+cc,'max branch identity')
 distinct=sum(colors[a]!='C' and colors[N-a]!='C' for a in range(1,(N+1)//2))
 if N>=7:require(distinct>0,'max branch tail')
 if N==6:require(distinct==0,'last failure6')
for N in range(2,31):
 H=N-1;B=[a for a in range(1,N) if colors[a]=='B'];idx={v:i for i,v in enumerate(B)}
 pairs=[(a,N-a)for a in range(1,(N+1)//2)]
 raa=sum(colors[a]==colors[b]=='A' for a,b in pairs)
 rab=sum(sorted((colors[a],colors[b]))==['A','B'] for a,b in pairs)
 rbb=sum(colors[a]==colors[b]=='B' for a,b in pairs)
 u=len(B)-rab-2*rbb;require(u>=0,'free choices')
 polynomial=[0]*raa+[2**u]
 for factor,count in (([1,1],rab),([3,1],rbb)):
  for _ in range(count):
   new=[0]*(len(polynomial)+1)
   for i,v in enumerate(polynomial):
    for j,w in enumerate(factor):new[i+j]+=v*w
   polynomial=new
 actual=Counter()
 for mask in range(1<<len(B)):
  A=lambda v:colors[v]=='A' or(colors[v]=='B' and mask>>idx[v]&1)
  actual[sum(A(a) and A(b) for a,b in pairs)]+=1
 require([actual[i]for i in range(len(polynomial))]==polynomial,'branch polynomial')
 require(sum(polynomial)==1<<len(B),'branch mass')
a=(3-sqrt(5))/2
print(json.dumps({'status':'PASS','checks':checks,'exact_rotation_prefix':maxn,'all_target_pair_tests_through':1500,'square_targets_tested':249,'square_statistics':square_stats,'finite_mean_R_cd_over_square':{k:v/249 for k,v in sorted(averages.items())},'limiting_AA_mean':a*a,'limiting_distinct_A_square_target_fraction':2*a,'first_forced_A_edge':[2,7],'scope':'finite checks supplement the analytic all-size proofs; no prime weights'},indent=2))
