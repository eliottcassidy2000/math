"""Exact controls for the separately review-pending TT3 degree identities."""
from pathlib import Path
from itertools import permutations
from math import comb,factorial
import hashlib,json,sys
ROOT=Path(__file__).resolve().parent

def need(v,s):
 if not v:raise RuntimeError(s)
def rows(n,c):
 r=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   u,v=(i,j)if c>>k&1 else(j,i);r[u]|=1<<v;k+=1
 return r
def matrix(R):
 n=len(R);N=1<<n;G=[[[0]*N for _ in R]for _ in R]
 for s in range(n):
  d=[[0]*n for _ in range(N)];d[1<<s][s]=1
  for M in range(1,N):
   for v,z in enumerate(d[M]):
    if not z:continue
    G[s][v][M]=z;A=R[v]&~M
    while A:
     b=A&-A;A-=b;d[M|b][b.bit_length()-1]+=z
 return G
def plus(*xs):return [sum(v)for v in zip(*xs)]
def mul(A,B):
 N=len(A);out=[0]*N
 for M in range(N):
  S=M
  while True:
   out[M]+=A[S]*B[M^S]
   if not S:break
   S=(S-1)&M
 return out
def full(A,B):return sum(A[s]*B[len(A)-1-s]for s in range(len(A)))
def H(R):
 n=len(R);N=1<<n;d=[[0]*n for _ in range(N)]
 for i in range(n):d[1<<i][i]=1
 for M in range(1,N):
  for i,z in enumerate(d[M]):
   if z:
    q=R[i]&~M
    while q:
     b=q&-q;q-=b;d[M|b][b.bit_length()-1]+=z
 return sum(d[-1])
def attach(T,a,c):
 n=len(T);R=T[:]+[1<<(n+1),0]
 for i in range(n):
  if a>>i&1:R[i]|=1<<n
  else:R[n]|=1<<i
  if c>>i&1:R[i]|=1<<(n+1)
  else:R[n+1]|=1<<i
 return R
def blow(Q,blocks):
 starts=[];n=0
 for B in blocks:starts.append(n);n+=len(B)
 T=[0]*n
 for i,B in enumerate(blocks):
  for u,row in enumerate(B):
   T[starts[i]+u]=row<<starts[i]
   for j,C in enumerate(blocks):
    if Q[i]>>j&1:T[starts[i]+u]|=((1<<len(C))-1)<<starts[j]
 return T,starts
def setup(E,ao,mo):
 Q=[2,4,1]+[0]*len(E)
 for i,row in enumerate(E):
  Q[i+3]|=row<<3
  if ao>>i&1:Q[0]|=1<<(i+3)
  else:Q[i+3]|=1
  if mo>>i&1:Q[1]|=1<<(i+3);Q[2]|=1<<(i+3)
  else:Q[i+3]|=6
 return Q
def tt(n):return [sum(1<<j for j in range(i+1,n))for i in range(n)]
def tq(n,d):
 a=[0]*(d+1);a[0]=1
 for _ in range(n):a=[0]+[k*(a[k]+a[k-1])for k in range(1,d+1)]
 return a

def kernel_inputs(R,out):
 n=len(R);N=1<<n;G=matrix(R);E=[plus(*(G[u][v]for u in range(n)))for v in range(n)];S=[plus(*G[v])for v in range(n)];z=[0]*N
 inn=[v for v in range(n)if not out>>v&1];oo=[v for v in range(n)if out>>v&1]
 al=plus(z,*(E[v]for v in inn));be=plus(z,*(S[v]for v in oo));al[0]=be[0]=1
 ta=plus(z,*(G[u][v]for u in oo for v in inn));ell=[plus(z,*(G[i][v]for v in inn))for i in range(n)];rho=[plus(z,*(G[u][i]for u in oo))for i in range(n)]
 return G,E,S,al,be,ta,ell,rho

def quartic(T,p,q,r,i):
 G=matrix(T);n=len(T);C=[plus(*(G[u][v]for u in range(n)))for v in range(n)];R=[plus(*G[v])for v in range(n)]
 return sum(full(mul(C[u],G[v][w]),R[z])for u,v in [(p,q),(q,p)]for w,z in [(r,i),(i,r)])+sum(full(mul(C[w],G[z][u]),R[v])for u,v in [(p,q),(q,p)]for w,z in [(r,i),(i,r)])

def check_tomography():
 records=[]
 for n in range(1,4):
  for co in range(1<<(n*(n-1)//2)):
   R=rows(n,co)
   for out in range(1<<n):
    G,E,S,al,be,ta,ell,rho=kernel_inputs(R,out);ab=mul(al,be)
    for ac in range(8):
     A=rows(3,ac);T=[row|sum(1<<(3+i)for i in range(n)if out>>i&1)for row in A]+[(r<<3)|(0 if out>>i&1 else 7)for i,r in enumerate(R)]
     for r in range(3):
      p,q=[v for v in range(3)if v!=r];d=sum(A[v]>>r&1 for v in range(3));fac=[2*z for z in ta];fac[0]+=d;fac2=[2*z for z in ta];fac2[0]+=2-d
      for i in range(n):
       expected=full(mul(al,S[i]),fac)+full(mul(be,E[i]),fac2)+2*full(ab,plus(ell[i],rho[i]));actual=quartic(T,p,q,r,3+i);need(actual==expected,'three-plus-one identity');records.append([n,co,out,ac,r,i,actual])
 return records

def skew(E,Aout,Mout,bx,by):
 n=len(E);N=1<<n;nx=len(bx)-1;ny=len(by)-1;zero=(0,0,0);one={zero:1}
 def varx():return {(0,1,0):1}
 def vary():return {(0,0,1):1}
 def add(*ps):
  out={}
  for p in ps:
   for e,c in p.items():out[e]=out.get(e,0)+c
  return {e:c for e,c in out.items() if c}
 def sc(c,p):return {e:c*d for e,d in p.items() if c*d}
 def mul(*ps):
  out=one
  for p in ps:
   z={}
   for (M,i,j),a in out.items():
    for (S,u,v),b in p.items():
     if M&S or i+u>nx or j+v>ny:continue
     key=(M|S,i+u,j+v);z[key]=z.get(key,0)+a*b
   out={e:c for e,c in z.items() if c}
  return out
 G=matrix(E)
 def paths(starts,ends):return {(mask,0,0):sum(G[u][v][mask] for u in starts for v in ends) for mask in range(1,N) if sum(G[u][v][mask] for u in starts for v in ends)}
 ia=[v for v in range(n) if not Aout>>v&1];oa=[v for v in range(n) if Aout>>v&1];im=[v for v in range(n) if not Mout>>v&1];om=[v for v in range(n) if Mout>>v&1]
 al=add(one,paths(range(n),ia));be=add(one,paths(oa,range(n)));a=add(one,paths(range(n),im));b=add(one,paths(om,range(n)));sg=paths(om,im);mu=paths(oa,im);ka=paths(om,ia);ta=paths(oa,ia);x=varx();y=vary();h=add(x,y,mul(x,y));U=add(mul(mu,h),x,mul(x,y));V=add(mul(ka,h),y,mul(x,y));sh=mul(sg,h);powers=[one]
 for r in range(1,n+1):powers.append(mul(powers[-1],sh))
 ds=[None]+[add(*(sc(comb(r+k-1,k-1),p) for r,p in enumerate(powers))) for k in range(1,4)]
 Ax=add(mul(al,ds[1]),mul(a,V,ds[2])); Bx=add(mul(be,ds[1]),mul(b,U,ds[2]))
 kx=add(mul(b,x,add(one,y),Ax),sc(-1,mul(a,x,Bx)))
 ky=add(mul(b,y,Ax),sc(-1,mul(a,y,add(one,x),Bx)))
 value=lambda p:sum(c*bx[i]*by[j] for (M,i,j),c in p.items() if M==N-1)
 return [[i,j,c] for (M,i,j),c in sorted(kx.items()) if M==N-1 and i and j],[[i,j,c] for (M,i,j),c in sorted(ky.items()) if M==N-1 and i and j]

def check_binomial():
 profiles=[tq(n,n)for n in range(1,21)]+[[0,3,6,6]]+[[0]+[comb(n,k)for k in range(n+1)]for n in range(1,16)]
 rec=[]
 for qi,q in enumerate(profiles):
  need(all(q[k]**2>=q[k-1]*q[k+1]for k in range(1,len(q)-1)),'profile logconcavity')
  val=lambda i:q[i]if 0<=i<len(q)else 0
  for c in range(6):
   for m in range(16):
    lhs=sum(comb(m,k)*(2*val(k+c+1)*val(m-k+c+1)-val(k+c+2)*val(m-k+c)-val(k+c)*val(m-k+c+2))for k in range(m+1));N=m+2
    rhs=sum(comb(N,j)*(N-(2*j-N)**2)*val(c+j)*val(c+N-j)for j in range(N+1));need(lhs*N*(N-1)==rhs and lhs>=0,'binomial identity/sign');rec.append([qi,c,m,lhs])
 return rec

def check_skew():
 rec=[]
 for n in range(4):
  dim=n+3;profiles=[tq(k,dim)for k in range(1,9)]+[[0,3,6,6]+[0]*(dim-3)]
  for co in range(1<<(n*(n-1)//2)):
   E=rows(n,co)
   for ao in range(1<<n):
    for mo in range(1<<n):
     px,py=skew(E,ao,mo,[1]*(dim+1),[1]*(dim+1));M=blow(setup(E,ao,mo),[[0],[0],[0]]+[[0]]*n)[0]
     # The M placeholder keeps only E and one vertex with common M incidences.
     EM=[row<<1|(0 if mo>>i&1 else 1)for i,row in enumerate(E)];EM=[sum(1<<(i+1)for i in range(n)if mo>>i&1)]+EM;hm=H(EM)
     for qi,q in enumerate(profiles):
      x=sum(z*q[i]*q[j]for i,j,z in px);y=sum(z*q[i]*q[j]for i,j,z in py);lower=2*q[1]**2*hm;need(x-y>=lower>0,'strict skew');rec.append([n,co,ao,mo,qi,x,y,lower])
 return rec

def check_physical():
 T,_=blow([2,4,1],[tt(3),[0],[0]]);U,_=blow([4,1,2],[tt(3),[0],[0]]);n=len(T);N=1<<n;F=[];G=[]
 for a in range(N):
  for c in range(N):F.append(H(attach(T,a,c)));G.append(H(attach(U,a,c)))
 co=F[:];cg=G[:]
 for bit in range(2*n):
  for M in range(1<<(2*n)):
   if M>>bit&1:co[M]-=co[M^(1<<bit)];cg[M]-=cg[M^(1<<bit)]
 masks=[x for x in range(N)if x.bit_count()<=2];bank=[[a,c,co[(a<<n)|c],cg[(a<<n)|c]]for a in masks for c in masks]
 need(all(not z for m,z in enumerate(co)if (m>>n).bit_count()>2 or(m&(N-1)).bit_count()>2),'full bidegree')
 for r in range(3):
  p,q=[i for i in range(3)if i!=r]
  for i in [3,4]:need(co[((1<<p|1<<q)<<n)|(1<<r|1<<i)]==quartic(T,p,q,r,i),'literal quartic')
 return {'T':T,'U':U,'source_values':F,'target_values':G,'complete_named_bank':bank}

def check_boundaries():
 records=[]
 for T,ex,ey in [([2,380,377,369,353,321,1,127,193],386,404),([2,444,441,433,417,385,63,65,193],806,806)]:
  n=len(T);G=matrix(T);E=[plus(*(G[u][v]for u in range(n)))for v in range(n)];S=[plus(*G[v])for v in range(n)];ks=[full(E[0],S[i])-full(E[i],S[0])for i in range(n)];kx=ks[1];ky=sum(ks[2:6]);need((kx,ky)==(ex,ey),'unequal profile boundary');records.append({'rows':T,'virtual_z':0,'X':[1],'Y':[2,3,4,5],'skew_vector':ks,'K_X':kx,'K_Y':ky})
 return records

def main():
 t=check_tomography();b=check_binomial();s=check_skew();p=check_physical()
 result={'status':'PASS','proof_sha256':hashlib.sha256((ROOT/'PROOF.md').read_bytes()).hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'tomography_controls':len(t),'binomial_controls':len(b),'equal_profile_controls':len(s),'tomography':t,'binomial':b,'equal_profile':s,'physical':p,'unequal_profile_boundaries':check_boundaries(),'scope':'Sections1--3 exact controls; Section4 integration remains review-pending.'}
 print(json.dumps(result,sort_keys=True,separators=(',',':')))
if __name__=='__main__':main()
