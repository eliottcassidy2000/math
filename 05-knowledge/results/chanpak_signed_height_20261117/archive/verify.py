from pathlib import Path
from itertools import product
from math import ceil
import sys,json
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'pubsub_branch'))
from replay import decode,h_of,table

def coefficients(t):return 12*t-3,15*t-4

def profile(m):return 2**m+1,4*2**m+2*(3**m-1)

def output(t,m,s):
 a,b=coefficients(t);H,P=profile(m)
 return a*H+b*P*2**s

def inverse(N):
 ans=[];m=1
 while True:
  H,P=profile(m)
  if 9*H+11*P>N:break
  s=0;Q=P
  while 9*H+11*Q<=N:
   t,r=divmod(N+3*H+4*Q,12*H+15*Q)
   if r==0 and t>=1:
    word=decode(t)
    if word is not None:ans.append((t,m,s))
   s+=1;Q*=2
  m+=1
 return ans

def require(x):
 if not x:raise ArithmeticError('control failed')

nchecks=0
for length in range(3):
 for letters in product('ABC',repeat=length):
  t=h_of(''.join(letters));a,b=coefficients(t)
  require(5*a-4*b==1)
  for m in range(1,5):
   for s in range(4):
    H,P=profile(m);P*=2**s;N=a*H+b*P;J=4*H+5*P
    require((5*N-b*J,a*J-4*N)==(H,P))
    require((t,m,s) in inverse(N));nchecks+=1
for word in ['ABCABCABCABCABCABC','C'*18,'A'*35]:
 t=h_of(word);N=output(t,7,9)
 require((t,7,9) in inverse(N));nchecks+=1
for m in range(1,7):
 n=m+2;rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
 rows[0]^=1<<(n-1);rows[-1]=1
 dp=table(tuple(rows));hs=[1]+[sum(r) for r in dp[1:]];full=len(hs)-1
 require((hs[-1],sum(hs[x]*hs[full^x] for x in range(full+1)))==profile(m))
checks=0
for t in [1,61,4397,71857]:
 a,b=coefficients(t);threshold=25*b*(a+4*b)
 for N in range(threshold|1,(threshold|1)+200,2):
  lower=(24*N+a+4*b-1)//(a+4*b);upper=(5*N-1)//b
  J=lower+(12-lower)%24
  require(J<=upper)
  H=5*N-b*J;P=a*J-4*N
  require(H>0 and H%2 and P>=4*H and P%8==0 and (P-H)%3==0);checks+=1
require(output(1,1,0)==159 and output(1,2,0)==397)
print(json.dumps({'exact_target_decodes':nchecks,'necessary_gate_interval_checks':checks,'literal_extreme_sources':6,'sample_changed_source_outputs':[159,397]},indent=2))
