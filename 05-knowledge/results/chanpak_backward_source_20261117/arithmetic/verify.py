from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'pubsub_branch'))
from replay import h_of,table

def candidate(N,t):
 a=12*t-3;b=15*t-4
 low=(24*N+48+a+4*b-1)//(a+4*b)
 high=(5*N-1)//b
 J=high-(high-12)%24
 if J<low:return None
 return J,5*N-b*J,a*J-4*N

def require(b):
 if not b:raise ArithmeticError('control failed')

count=0
for word in ['', 'A','B','C','ABC','CAB','A'*80,'BC'*20]:
 t=h_of(word);a=12*t-3;b=15*t-4
 N=a*187+b*3736
 got=candidate(N,t)
 if N<24*b*(a+4*b):require(got==(19428,187,3736))
 require(a*197+b*3728==N+2)
 for shift in range(-100,101,2):
  target=N+shift
  if target<24*b*(a+4*b):
   z=candidate(target,t)
   if z is not None:
    J,H,P=z
    require(J%24==12 and H>0 and H%2 and P>=4*H+48 and P%8==0 and (P-H)%3==0)
    require(a*H+b*P==target)
   count+=1
# Direct profile controls on all order-four and order-five tournaments,
# then 2-adically pad pairs to orders beyond any finite raw census.
profiles={}
for n in [4,5]:
 pset=set()
 for bits in range(1<<(n*(n-1)//2)):
  rows=[0]*n;i=0
  for u in range(n):
   for v in range(u+1,n):
    if bits>>i&1:rows[u]|=1<<v
    else:rows[v]|=1<<u
    i+=1
  d=table(rows);h=[1]+[sum(r) for r in d[1:]];M=len(h)-1
  H=h[-1];P=sum(h[s]*h[M^s] for s in range(M+1))
  if H>1:pset.add((H,P))
 profiles[n]=pset
vchecks=0
for n,ps in profiles.items():
 for h,p in ps:
  for k,q in ps:
   if p==q:continue
   for s in [0,1,2,10,1000]:
    diff=abs(q-p)*2**s;v=(diff&-diff).bit_length()-1
    require(n+s<=max(h,k)+v-2);vchecks+=1
print(json.dumps({'tall_candidate_controls':count,'padded_fine_pair_bound_controls':vchecks,'native_source_J':19428},indent=2))
