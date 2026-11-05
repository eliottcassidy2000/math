"""Complete fixed-coefficient full-gadget inverse; only Python standard library."""
from functools import lru_cache
from pathlib import Path
import hashlib
import json
@lru_cache(None)
def c(L):
 if L==1:return 1,1,1
 if L%2==0:
  g,Q,B=c(L//2);return g,Q,2*B
 ell=(L+1)//2;g,Q,B=c(ell);return g+ell-ell%2,Q+B,2*B

def rev(x,k):
 y=0
 for _ in range(k):y=2*y+1-(x&1);x>>=1
 return y

def value(v,L):
 A,B,C,D,E=v;g,Q,Z=c(L);return A*g+B*L+C*Q+D*Z-E*(L%2)

def inverse(v,T):
 A,B,C,D,E=v
 if not(A>=1 and B>=A and 0<=C<=4*A and D>=0 and 0<=E<=A):raise ValueError('cone')
 if T<1:return [],0
 out=[1]if value(v,1)==T else[];ev=1
 if T<A+2*B+C+2*D:return out,ev
 K=T//B
 for k in range(1,(K-1).bit_length()+1):
  lo=(1<<(k-1))+1;hi=min(1<<k,K)
  if lo>hi:continue
  if k<=4:
   for L in range(lo,hi+1):
    ev+=1
    if value(v,L)==T:out.append(L)
   continue
  q=(k+4)//2;p=k-q;Z=1<<k;B0=1<<q;H=1<<p
  table=[0,1]+[0]*(B0-1)
  for n in range(2,B0+1):
   ell=(n+1)//2
   table[n]=table[n//2]if n%2==0 else table[ell]+ell-ell%2
  reversals=[rev(x,p)for x in range(H)]
  for b in range(B0):
   xlo=max(0,((lo-1-b)+B0-1)//B0);xhi=(hi-1-b)//B0
   if xlo>xhi:continue
   qb=rev(b,q);par=qb%2
   constant=A*(table[b+1]-1)+B*(b+1)+C*(1+H*qb)+D*Z-E*((b+1)%2)
   slope=A*qb+B*B0
   def f(x):return A*table[x+1]+slope*x+A*par*(x%2)+C*reversals[x]+constant
   while xlo<xhi:
    mid=(xlo+xhi)//2;ev+=1
    if f(mid)<T:xlo=mid+1
    else:xhi=mid
   ev+=1
   if f(xlo)==T:out.append(B0*xlo+b+1)
 return sorted(out),ev


def gadget(L):
 m=(L-1).bit_length();T=[0]*(m+1)
 for i in range(m):
  x=i+1
  T[x if (L-1)>>i&1 else 0]|=1<<(0 if (L-1)>>i&1 else x)
  for j in range(i+1,m):T[x]|=1<<(j+1)
 return T

def physical(T):
 n=len(T);N=1<<n;starts=[[0]*n for _ in range(N)];count=[1]+[0]*(N-1)
 for S in range(1,N):
  for v in range(n):
   if S>>v&1:
    U=S^(1<<v)
    starts[S][v]=1 if not U else sum(starts[U][w]for w in range(n)if U>>w&1 and T[v]>>w&1)
  count[S]=sum(starts[S])
 return count[-1],sum(starts[S][0]*count[(N-1)^S]for S in range(N)if S&1)

def main():
 P=Path(__file__).resolve().parent
 dependencies={'ROOTED.md':'ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e','SQRT_FIXED_PREFIX_INVERSE.md':'084b7c372524f41e99a4087179e0fd30baf35a0de58a0166e3bdcd84c5e407bf'}
 for name,digest in dependencies.items():
  if hashlib.sha256((P/name).read_bytes()).hexdigest()!=digest:raise RuntimeError(('dependency',name))
 identities=split_controls=block_controls=monotone=cases=physical_controls=0
 for L in range(1,8193):
  k=(L-1).bit_length();g,Q,Z=c(L)
  if Q!=1+rev(L-1,k)or Z!=1<<k:raise RuntimeError(('minimal reversal',L))
  identities+=1
 for k in range(5,13):
  q=(k+4)//2;p=k-q;B0=1<<q;H=1<<p
  for b in range(B0):
   u=rev(b,q)
   for x in range(H//2,H):
    L=B0*x+b+1;g,Q,Z=c(L)
    if Q!=1+H*u+rev(x,p):raise RuntimeError(('Q split',L))
    split_controls+=1
    if g!=c(x+1)[0]+u*x+(u%2)*(x%2)+c(b+1)[0]-1:raise RuntimeError(('g block',L))
    block_controls+=1
   for A in (1,3):
    for B in (A,A+1):
     for C in (0,4*A):
      v=A,B,C,1,A
      vals=[value(v,B0*x+b+1)for x in range(H//2,H)]
      for old,new in zip(vals,vals[1:]):
       if not old<new:raise RuntimeError(('strict branch',v,k,b))
       monotone+=1
 for A in range(1,5):
  for B in (A,A+1):
   for C in (0,A,4*A):
    for D in (0,1,3):
     for E in sorted({0,A//2,A}):
      v=A,B,C,D,E
      for T in range(251):
       got,_=inverse(v,T);want=[L for L in range(1,T//B+1)if value(v,L)==T]
       if got!=want:raise RuntimeError(('inverse',v,T,got,want))
       cases+=1
 for L in range(1,257):
  got=physical(gadget(L));g,Q,Z=c(L)
  if got!=(g,Q):raise RuntimeError(('physical gadget',L,got,(g,Q)))
  physical_controls+=1
 large=[]
 for v,L in [((9,9,30,30,3),1234567),((39,54,156,216,6),9876543),((100,300,400,1234,100),100000007)]:
  T=value(v,L);got,ev=inverse(v,T)
  if got!=[L]:raise RuntimeError(('large exact inverse',v,T,got))
  large.append(dict(coefficients=v,target=T,solutions=got,evaluations=ev))
 invalid=0
 for v in ((0,1,0,0,0),(2,1,0,0,0),(1,1,-1,0,0),(1,1,5,0,0),(1,1,0,-1,0),(1,1,0,0,-1),(1,1,0,0,2)):
  try:inverse(v,1)
  except ValueError:invalid+=1
  else:raise RuntimeError(('accepted invalid',v))
 # Huge coefficients are rejected before any target-sized table is allocated.
 if inverse((10**100,10**100,0,10**200,0),1)!=([],1):raise RuntimeError('large early reject')
 result=dict(status='PASS',dependencies=dependencies,Q_minimal_width_controls=identities,Q_split_controls=split_controls,g_block_controls=block_controls,strict_branch_controls=monotone,exhaustive_inverse_controls=cases,physical_gadget_controls=physical_controls,invalid_inputs=invalid,large=large)
 text=json.dumps(result,indent=2,sort_keys=True)+'\n';(P/'full_gadget_inverse_certificate.json').write_text(text);print(text,end='')

if __name__=='__main__':main()
