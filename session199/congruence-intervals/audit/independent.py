"""Independent integer audit: imports no primary functions or research packages."""
from pathlib import Path
from hashlib import sha256
from math import gcd,isqrt
import json,sys,subprocess,random
HERE=Path(__file__).resolve().parent
BASE=HERE.parent
# A copied release can place this audit/ beside PROOF.md and verify.py.
if not (BASE/'PROOF.md').exists():BASE=Path('session199-cross-audit/congruence-intervals').resolve()
PROOF='aa62ce77eaeb039a702d8af784e41d6480d1867efa7d2919c9003dcbbef55bff'
SCRIPT='243918e4a0978e5a8ec1670cc0a748f5ad9558d479d2efb25a862f7ec7a2bf46'
checks=0

def need(c,m):
 global checks
 checks+=1
 if not c:raise RuntimeError(m)
def g(U):
 # Direct four-coordinate O/Z physical compiler from leading binary bit.
 D,E,A,d=1,0,2,1
 for bit in bin(U-1)[3:]:
  if bit=='1':D,E,A,d=D,0,2*A,d
  else:D,E,A,d=D-E+A,d,2*A-d,d
 need(A==U and E==U%2 and d==1,'source dictionary')
 return D

def step(n,F,b,M):
 return (2*n+b)%M,(F+(n+(n&1)if b==0 else 0))%M

def password(n,F,w,M):
 for b in w:n,F=step(n,F,b,M)
 return n,F

def bits(n,k):return [(n>>j)&1 for j in range(k-1,-1,-1)]

def parameters(M,compressed):
 k=0;q=M
 while q%2==0:k+=1;q//=2
 b=(M-1).bit_length();a=2
 while (2**a-1)%q:a+=1
 T=(M//2-1)*(a+k)
 if compressed and q==1:T=0 if k==1 else (5*k*k-3*k-6)//2
 return k,q,b,a,T

def construct(P,F,u,z,M,compressed=True):
 k,q,b,a,T=parameters(M,compressed)
 reset=bits((-P*2**b)%M,b)
 n0,F0=password(P%M,F%M,reset,M)
 need(n0==0,'prefix reset')
 end=bits((u-1)%M,b)
 _,delta_end=password(0,0,end,M)
 required=(z-F0-delta_end)%M
 loops=[]
 if compressed and q==1:
  j=(-required//2)%(M//2)
  for i in range(k-1):
   if j>>i&1:
    loops.extend([1]*(i+1 if i else 2)+[0]*k)
    if i:loops.extend([1,1]+[0]*k)
 else:
  W=[1]*a+[0]*k
  nz,delta=password(0,0,W,M)
  if nz or gcd(delta,M)!=2:raise RuntimeError('loop increment')
  # Different solver: explicitly enumerate the short cyclic group.
  j=next(j for j in range(M//2)if (j*delta-required)%M==0)
  loops=W*j
 if len(loops)>T:raise RuntimeError('budget')
 word=reset+[0]*(T-len(loops))+loops+end
 if password(P%M,F%M,word,M)!=((u-1)%M,z%M):raise RuntimeError('residue finish')
 return word

def main():
 need(sha256((BASE/'PROOF.md').read_bytes()).hexdigest()==PROOF,'proof hash')
 need(sha256((BASE/'verify.py').read_bytes()).hexdigest()==SCRIPT,'script hash')
 state_tests=0;integer_tests=0;approximations=0
 for M in range(2,22,2):
  for P in range(M):
   for F in range(1,M,2):
    for u in range(M):
     for z in range(1,M,2):
      w=construct(P,F,u,z,M);state_tests+=1
      need(password(P,F,w,M)==((u-1)%M,z%M),'all modular endpoints')
 rng=random.Random(92173)
 for k in range(1,21):
  for _ in range(25):
   P=rng.getrandbits(100)+1;M=1<<k;u=rng.randrange(M);z=2*rng.randrange(M//2)+1
   F=g(P+1);w=construct(P,F,u,z,M)
   n=P
   for bit in w:n=2*n+bit
   U=n+1
   need(U%M==u and g(U)%M==z,'actual integer suffix')
   need((n>>len(w))==P,'exact prefix cylinder');integer_tests+=1
 for m in range(1,35):
  M=2*m;_,_,b,_,T=parameters(M,True);B=1<<(2*b+T)
  for t in [9,10,15,16,17,25,63,64,65,100,1000,10**6,10**20]:
   N=t*B+1;x=isqrt(N//B);Q=1<<(x.bit_length()-1)
   P=N//(B*(Q+1))-1
   need(P>=1,'positive prefix')
   w=construct(P,g(P+1),1,(N-Q)%M,M);n=P
   for bit in w:n=2*n+bit
   U=n+1;D=g(U);H=D+Q*U
   need(H<=N and (N-H)%M==0,'lower and congruence')
   need((N-H)**2<9*B*N,'strict error')
   need((Q+1)*U<=N<(Q+1)*(U+2*B),'prefix inequalities')
   vertices=(U-1).bit_length()+Q.bit_length()+1
   need(vertices<=N.bit_length()+2,'order bound');approximations+=1
 primary=[]
 for optimized in [False,True]:
  p=subprocess.run([sys.executable]+(['-O']if optimized else [])+[str(BASE/'verify.py')],cwd=BASE,text=True,capture_output=True,check=True)
  primary.append(p.stdout)
 need(primary[0]==primary[1],'original normal optimized byte match')
 report=dict(status='PASS',proof_sha256=PROOF,verifier_sha256=SCRIPT,independent_modular_endpoints=state_tests,independent_large_integer_prefixes=integer_tests,independent_approximations=approximations,checks=checks,primary_output_sha256=sha256(primary[0].encode()).hexdigest())
 print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':main()
