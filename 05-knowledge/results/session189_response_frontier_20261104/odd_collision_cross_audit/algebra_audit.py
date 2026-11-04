"""Independent weighted-path subset polynomials and ordinary Stirling transform."""
from pathlib import Path
import math,json,hashlib
D=Path(__file__).resolve().parent;B=64;MASK=(1<<B)-1

def need(c,s):
 if not c:raise RuntimeError(s)
def tournament(n,code):
 R=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if code>>k&1:R[i]|=1<<j
   else:R[j]|=1<<i
   k+=1
 return R

def paths(R):
 n=len(R);size=1<<n;G=[[0]*(n*n) for _ in range(size)];C=[[0]*size for _ in R];A=[[0]*size for _ in R];H=[0]*size
 need(math.factorial(n)*(1<<n)*6**n<1<<B,'packed coefficient no-carry bound')
 for S in range(1,size):
  for v in range(n):
   if not(S>>v&1):continue
   M=S^(1<<v)
   if not M:G[S][v*n+v]=1;continue
   for u in range(n):
    if not(M>>u&1):continue
    value=0
    for w in range(n):
     if M>>w&1:
      p=G[M][u*n+w];value+=p<<B
      if R[w]>>v&1:value+=p
    G[S][u*n+v]=value
  for u in range(n):
   for v in range(n):
    p=G[S][u*n+v];A[u][S]+=p;C[v][S]+=p;H[S]+=p
 return G,C,A,H

def full2(a,b,full):
 ans=0;S=full
 while S:
  ans+=a[S]*b[full^S];S=(S-1)&full
 return ans

def full3(a,b,c,full):
 ans=0;S=full
 while S:
  rest=full^S;U=rest
  if a[S]:
   while U:
    ans+=a[S]*b[U]*c[rest^U];U=(U-1)&rest
  S=(S-1)&full
 return ans

def decode(z,n):return [(z>>(B*j))&MASK for j in range(n)]
def kernel(R,p,q,r,cache):
 n=len(R);full=(1<<n)-1;G,C,A,H=cache
 def mij(i,j):return [g[i*n+j] for g in G]
 terms=[(C[p],mij(q,r),1),(C[q],mij(p,r),1),(mij(r,p),A[q],-1),(mij(r,q),A[p],-1)]
 X=[0]*(n-1);Y=[0]*(n-2)
 for a,b,sgn in terms:
  for j,x in enumerate(decode(full2(a,b,full),n-1)):X[j]+=sgn*x
  for j,y in enumerate(decode(full3(H,a,b,full),n-2)):Y[j]+=sgn*y
 for a,b in [(C[p],A[q]),(C[q],A[p])]:
  for c,sgn in [(C[r],1),(A[r],-1)]:
   for j,y in enumerate(decode(full3(a,b,c,full),n-2)):Y[j]+=sgn*y
 need(Y==[(j+1)*X[j+1] for j in range(n-2)],'formal two/three-channel derivative for arbitrary interior')
 return X,Y

def ordinary(m,r):
 L=m-2;u=r-1;p=[0]*(L+1)
 for j in range(u+1):p[L-u+j]+=math.comb(u,j)
 for j in range(L-u+1):p[u+j]-=math.comb(L-u,j)
 return p

def transform(p):
 n=len(p);S=[[0]*n for _ in range(n)];S[0][0]=1
 for d in range(1,n):
  for j in range(1,d+1):S[d][j]=S[d-1][j-1]+j*S[d-1][j]
 return [math.factorial(j)*sum(p[d]*S[d][j] for d in range(n)) for j in range(n)]
def weights(m):
 w={m-1:1}
 for r in range(2,(m-1)//2+1):w[r]=(-1)**r*math.comb(m-2,r-1) if m%2 else (-1)**r*(math.comb(m-3,r-1)-math.comb(m-3,r-2))
 return w

def main():
 controls=[]
 for n in [3,4]:
  for code in range(1<<(n*(n-1)//2)):
   R=tournament(n,code);cache=paths(R)
   for p in range(n):
    for q in range(p+1,n):
     for r in range(n):
      if r in [p,q]:continue
      X,Y=kernel(R,p,q,r,cache);controls.append(dict(rows=R,marks=[p,q,r],X=X,derivative_channel=Y))
 for n in range(5,10):
  R=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)];cache=paths(R)
  for r in range(2,n):
   X,Y=kernel(R,0,1,r,cache);expected=[2*x for x in transform(ordinary(n,r))];need(X==expected,'transitive kernel via independent Stirling transform');controls.append(dict(rows=R,marks=[0,1,r],X=X,derivative_channel=Y,ordinary_polynomial=ordinary(n,r)))
 combos=[]
 for m in range(3,102):
  out=[0]*(m-1)
  for r,w in weights(m).items():
   for d,c in enumerate(ordinary(m,r)):out[d]+=w*c
  want=[1]+[0]*(m-2) if m%2 else [1,2]+[0]*(m-3);need(out==want,'exact ordinary binomial cancellation')
  kernel_combined=[2*x for x in transform(out)];derivative=[(i+1)*kernel_combined[i+1] for i in range(len(kernel_combined)-1)]
  need(kernel_combined==([2]+[0]*(m-2) if m%2 else [2,4]+[0]*(m-3)),'finite transform gives odd/even kernel');combos.append(dict(m=m,weights=weights(m),ordinary_combination=out,complete_kernel_combination=kernel_combined,complete_derivative=derivative))
 outputs={'kernel_controls.json':controls,'binomial_controls.json':combos}
 for name,obj in outputs.items():
  p=D/name;txt=json.dumps(obj,separators=(',',':'))+'\n'
  if p.exists():need(p.read_text()==txt,'independent algebra data stable')
  else:p.write_text(txt)
 return dict(status='PASS',arbitrary_internal_kernels=792,additional_transitive_kernels=len(controls)-792,symbolic_binomial_orders=[3,101],all_even_combinations_retain_derivative_four=True,source_sha256=hashlib.sha256((D/'algebra_audit.py').read_bytes()).hexdigest(),output_pins={name:hashlib.sha256((D/name).read_bytes()).hexdigest() for name in outputs})
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
