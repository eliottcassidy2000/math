import json,math,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from decode import decode

def check(p,m):
 if not p:raise RuntimeError(m)
def g(n):
 z=1
 while n:
  u=n//2
  if n%2==0:z+=u+u%2
  n=u
 return z
def f(n):
 bits=[i for i in range(n.bit_length()) if n>>i&1]
 return sum((1<<(i-j-1))+(i==j+1) for i in bits for j in bits if i>j)
def rev(n,k):return int(f'{n:0{k}b}'[::-1],2) if k else 0
def value(w):return int(w,2) if w else 0
def poly(vals):
 y2,y4,y8=vals;s1=(y4-y2)//2;s2=(y8-y4)//4
 check((y4-y2)%2==0 and (y8-y4)%4==0 and (s2-s1)%6==0,'integrality')
 A=(s2-s1)//6;B=s1-6*A;D=y2-4*A-2*B
 return A,B,D
def hp(w):
 n=len(w)+1;out=[0]*n
 for i in range(1,n):
  for j in range(i+1,n):out[i]|=1<<j
  if w[i-1]=='1':out[0]|=1<<i
  else:out[i]|=1
 dp=[{} for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i,num in dp[mask].items():
   avail=out[i]&~mask
   while avail:
    low=avail&-avail;j=low.bit_length()-1;avail-=low
    dest=dp[mask|low];dest[j]=dest.get(j,0)+num
 return sum(dp[-1].values())

halfchecks=graphs=0
for l in range(0,10):
 B=1<<l
 for q in range(B):
  x=B-1-q;left=f'{x:0{l}b}' if l else '';p=q%2
  for beta in [None,0,1]:
   w=left+('' if beta is None else str(beta))+left[::-1]
   d=1-q*q-2*f(q)-p if beta is None else 1-2*(q*q+f(q)+(1-beta)*(q+p))
   for s in [1,2,5]:
    source='1'*s+w+'1'*s;N=(1<<(2*s+len(w)))+d
    check(g(value(source+'0'))==N,'half-word')
    check(N>0 and N%2==1,'positive odd')
    halfchecks+=1
   if l<=3:
    source='1'+w+'10';check(hp(source)==(1<<(2+len(w)))+d,'literal graph');graphs+=1
# The half-word candidate is exact, not just sufficient within a fixed q-range.
for odd in [False,True]:
 image={}
 for q in range(100):
  if odd:
   for b in [0,1]:image.setdefault(q*q+f(q)+(1-b)*(q+q%2),set()).add((q,b))
  else:image.setdefault(q*q+2*f(q)+q%2,set()).add((q,None))
 for T in range(5001):
  q=math.isqrt(T);got=set()
  if odd:
   for b in [0,1]:
    if q*q+f(q)+(1-b)*(q+q%2)==T:got.add((q,b))
  elif q*q+2*f(q)+q%2==T:got.add((q,None))
  check(got==image.get(T,set()),'sqrt inverse')
words=['']+[f'{r:0{k}b}' for k in range(1,4) for r in range(1<<k)]
templates=rowcases=0
for U in ['']+[bin(x)[2:] for x in range(1,8)]:
 P=value(U)
 for W in words:
  M=1<<len(W);r=value(W);R=rev(r,len(W));C=M-1-R
  for S in words:
   MS=1<<len(S);CS=MS-1-rev(value(S),len(S))
   for beta in [0,1]:
    for gamma in [0,1]:
     vs=[value(U+str(beta)*s+W+str(gamma)*s+S) for s in [1,2,3]]
     va,vb,vd=poly(vs);ga,gb,gd=poly([g(v) for v in vs]);templates+=1
     v4=value(U+str(beta)*4+W+str(gamma)*4+S)
     check(v4==va*256+vb*16+vd and g(v4)==ga*256+gb*16+gd,'polynomial test')
     for a in range(1,16,2):
      for c in range(1,17):
       K=c*MS+a*CS;actual=(a*ga+c*va,a*gb+c*vb)
       if beta==gamma==1:expected=(K*M*(P+1),a*C*(P+1)-K*(M-1-r))
       elif beta==gamma==0:expected=((a+K)*M*P,(a+K)*r-a*P*R)
       elif beta==1:expected=((a+K)*M*(P+1),-a*(R+1)*(P+1)-(a+K)*(M-r))
       else:expected=(K*M*P,a*P*(1+C)+K*(r+1))
       check(actual==expected,'coefficient law');A,L=actual;rowcases+=1
       if A>0 and A&(A-1)==0 and L==0:
        check(beta==gamma,'mixed exclusion')
        if a>1:check(W==str(beta)*len(W),'bank rigidity')
# All transitive bank rows are redundant in scalar support.
redundancy=0
for j in range(1,7):
 for v in range(1000):
  source=(1<<(j+1))*v+(1<<(j+1))-2
  check(g(source)==g(v)+(1<<j)*(v+1),'transitive support');redundancy+=1
# Complete finite image from the all-order exponent bound, independently enumerated.
image=set()
for m in range(2,16):
 for q in range(1<<max(0,(m-2)//2)):
  if m%2==0:
   if m>=2*q.bit_length()+2:image.add((1<<m)+1-q*q-2*f(q)-q%2)
  elif m>=2*q.bit_length()+3:
   for b in [0,1]:image.add((1<<m)+1-2*(q*q+f(q)+(1-b)*(q+q%2)))
for N in range(1,20001):
 cert=decode(N);check((cert is not None)==(N in image),'scalar decoder')
 if cert:check(g(value(cert['chain_word']))==N,'decoded witness')
receipt={'scalar_decoder_targets':20000,'status':'PASS','half_word_identity_checks':halfchecks,'literal_graph_DP_checks':graphs,'sqrt_membership_targets':10002,'two_clock_templates':templates,'coefficient_row_checks':rowcases,'transitive_support_checks':redundancy,'scope':'all-order bank rigidity proof; palindrome rays inherited; no completeness claim'}
print(json.dumps(receipt,indent=2,sort_keys=True))
