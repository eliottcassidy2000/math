"""Small exact controls for the graph/norm interface; no external dependencies."""
import json,math
from pathlib import Path
from itertools import combinations,permutations

def endpoints(T):
 n=len(T);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for mask in range(1,1<<n):
  for i,c in enumerate(d[mask]):
   q=T[i]&~mask
   while q:
    bit=q&-q;q-=bit;d[mask|bit][bit.bit_length()-1]+=c
 return d

def profile(T):
 d=endpoints(T);h=[sum(x) for x in d];h[0]=1;full=len(h)-1
 return h[full],sum(h[m]*h[full^m] for m in range(len(h)))
def cyclic(A,B):
 a=len(A);b=len(B);z=a+b;out=[]
 for x in A:out.append(x|(((1<<b)-1)<<a))
 for x in B:out.append((x<<a)|(1<<z))
 out.append((1<<a)-1);return out

def tt(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
def graph(n,mask):
 out=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  if mask>>k&1:out[i]|=1<<j
  else:out[j]|=1<<i
 return out

def need(q,msg):
 if not q:raise RuntimeError(msg)
small=[graph(n,m) for n in range(1,4) for m in range(1<<math.comb(n,2))]
pairs=0
for A in small:
 y,q=profile(A)
 need(q>=(len(A)+1)*y,'cut bound')
 for B in small:
  z,r=profile(B);N=profile(cyclic(A,B))[0]
  need(N==(q-y)*(r-z)+2*y*z,'cyclic formula');pairs+=1
hits=[];classes={};classified=0
for n in range(1,6):
 for mask in range(1<<math.comb(n,2)):
  T=graph(n,mask);y,q=profile(T);classified+=1
  need(q>=(n+1)*y,'cut bound exhaustive')
  if q in (5*y-3,5*y+3):
   item=(n,y,q);classes[item]=classes.get(item,0)+1
for t in range(1,13):
 x=2**t-1;N=x*x+2;v=(N-1)//2;m=math.isqrt(v)
 hits.append(dict(order=t,x=x,y=1,N=N,distinguished_m=m if m*m==v else None))
C3=[2,4,1];y,q=profile(C3);need((y,q)==(3,12),'C3 profile');need(profile(cyclic(C3,C3))[0]==99,'99')
need(profile(cyclic(tt(3),tt(3)))[0]==51,'51')
# Independent permutation count for the seven-vertex 51 and99 witnesses.
for T,target in [(cyclic(C3,C3),99),(cyclic(tt(3),tt(3)),51)]:
 need(sum(all(T[p[i]]>>p[i+1]&1 for i in range(len(p)-1)) for p in permutations(range(len(T))))==target,'permutation')
# Scalar factor3 switches and norm equality, positive odd y.
for y in range(1,501,2):
 for sign in (-1,1):
  x=4*y+3*sign;m=3*y+2*sign
  need(x*x+2*y*y==1+2*m*m,'factor3 norm switch')
result=dict(status='PASS',all_pair_checks=pairs,factor3_universe='all labelled tournaments on1..5 vertices',factor3_universe_count=classified,finite_factor3_profiles=[dict(n=n,H=y,Q2=q,labelled_count=c) for (n,y,q),c in sorted(classes.items())],transitive_diagonal=hits,C3=dict(H=3,Q2=12,N=99),factor3_odd_y_checks=500)
Path(__file__).with_name('norm_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
