from functools import lru_cache
from pathlib import Path
import hashlib,json
P=Path(__file__).resolve().parent

@lru_cache(None)
def gadget(L):
 if L==1:return 1,1,1
 if L%2==0:
  g,Q,B=gadget(L//2);return g,Q,2*B
 l=(L+1)//2;g,Q,B=gadget(l);return g+l-l%2,Q+B,2*B

def step(t,L):
 h,e,R=t;g,Q,B=gadget(L)
 return g*h-(L%2)*e+L*R,h-(L==1)*e,Q*h+B*R
@lru_cache(None)
def inverse(h,e,R):
 if h==1:return () if (h,e,R)==(1,1,1) else (1,)if(h,e,R)==(1,0,2)else None
 if h==3:return (2,)if(e,R)==(1,3)else None
 if h<3 or h%2==0 or e<1 or not h<=R<=2*(h-e):return None
 if 3*e>h or 3*e==h and R>h:
  old=(R-h+e,R-h,h-e)
  if not 1<=old[0]<h or 3*old[0]>2*h:return None
  w=inverse(*old)
  return None if w is None else w+(1,)
 if e%2==0:return None
 for L in range(2,R//e):
  g,Q,B=gadget(L);rem=R-Q*e
  if rem%B:continue
  R0=rem//B
  if R0<e:continue
  if L%2:
   e0=g*e+L*R0-h
   w=inverse(e,e0,R0)
   if w is not None:return w+(L,)
  elif h==g*e+L*R0:
   for e0 in range(e+1):
    w=inverse(e,e0,R0)
    if w is not None:return w+(L,)
 return None

def forward(N):
 Q=[(1,1,1)];seen=set(Q)
 for t in Q:
  for L in range(1,N+1):
   w=step(t,L)
   if w[0]>N:continue
   if w not in seen:seen.add(w);Q.append(w)
 return seen

def build(recipe):
 a=[0];root=0
 for L in recipe:
  n=len(a);m=(L-1).bit_length();a +=[0]*(m+1)
  for i in range(n):a[i]|=((1<<(m+1))-1)<<n
  a[root]^=1<<n;a[n]|=1<<root
  for i in range(m):
   x=n+1+i
   if (L-1)>>i&1:a[x]|=1<<n
   else:a[n]|=1<<x
   for j in range(i+1,m):a[x]|=1<<(n+1+j)
  root=n
 return a,root

def profile(a,r):
 n=len(a);N=1<<n;st=[[0]*n for _ in range(N)];en=[[0]*n for _ in range(N)];h=[1]+[0]*(N-1)
 for S in range(1,N):
  for v in range(n):
   if not S>>v&1:continue
   U=S^(1<<v)
   if U==0:st[S][v]=en[S][v]=1
   else:
    for z in range(n):
     if U>>z&1:
      if a[v]>>z&1:st[S][v]+=st[U][z]
      else:en[S][v]+=en[U][z]
  h[S]=sum(st[S])
 F=N-1;return (h[F],en[F][r],sum(st[S][r]*h[F^S]for S in range(1,N)if S>>r&1)),st[F][r],h[F^(1<<r)]

def deletion(t):
 h,e,R=t
 return 1 if h==1 else R-h+e if 3*e>h or 3*e==h and R>h else e

def oracle(D,E,U):
 """Map each possible deletion count to (core recipe, sink count)."""
 if D<1 or E<0 or E>D or U<1:return {}
 if E:
  if D==1:return {1:((),0)}if(E,U)==(1,1)else{}
  w=inverse(D,E,U)
  return {}if w is None else{deletion((D,E,U)):(w,0)}
 if D==1:
  return {1:((),U.bit_length()-1)}if U>=2 and U&(U-1)==0 else{}
 if D%2==0:return{}
 k=0;V=U
 while V>=2*D:
  if V%2:return{}
  V//=2;k+=1
 if k<1 or V<D:return{}
 out={}
 for endpoint in range(1,D+1):
  w=inverse(D,endpoint,V)
  if w is not None:
   out.setdefault(deletion((D,endpoint,V)),(w,k))
 return out

def sinks(a,k):
 a=list(a)
 for _ in range(k):
  n=len(a)
  a=[row|(1<<n)for row in a]+[0]
 return a

def adjacent(a,b):
 """New last vertex is the first clone of a new sole-out-neighbor vertex."""
 a=list(a);n=len(a)
 for i in range(n):
  if i!=b:a[i]|=(1<<n)|(1<<(n+1))
 a.extend([1<<b,(1<<b)|(1<<n)])
 return a,n+1

def main():
 dependencies={'BACKWARD_PROFILE_INVERSE.md':'8a07f1bb59f7ee3af9c172a1cfec8551c87125747eeb95dddef7335df33df1fb','verify_backward_profiles.py':'3f11b79a5470736e01ce1596a5d7b5fe9b75be9a78489cc285a9cfc2e22cf4b4','ROOTED.md':'ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e'}
 for name,digest in dependencies.items():
  if hashlib.sha256((P/name).read_bytes()).hexdigest()!=digest:raise RuntimeError(('dependency',name))
 N=61;bank=forward(N);expected={};core_count=0
 for t in sorted(bank):
  if t==(1,0,2):continue
  a,r=build(inverse(*t));got,start,d=profile(a,r)
  if got!=t or start!=1:raise RuntimeError(('core',t))
  D,E,U=t
  if oracle(*t).keys()!={d}:raise RuntimeError(('strong oracle',t,d))
  core_count+=1
  for k in range(1,4):
   expected.setdefault((D,0,(1<<k)*U),set()).add(d)
 weak_tests=0
 for D in range(1,N+1):
  for U in range(1,8*D+1):
   got=oracle(D,0,U)
   if set(got)!=expected.get((D,0,U),set()):raise RuntimeError(('weak membership',D,U,got,expected.get((D,0,U))))
   weak_tests+=1
 literal_graphs=literal_roots=0;witnesses={}
 for n in range(1,6):
  pairs=[(i,j)for i in range(n)for j in range(i+1,n)]
  for mask in range(1<<len(pairs)):
   a=[0]*n
   for q,(i,j)in enumerate(pairs):
    if mask>>q&1:a[i]|=1<<j
    else:a[j]|=1<<i
   literal_graphs+=1
   for r in range(n):
    t,start,d=profile(a,r)
    if start!=1:continue
    got=oracle(*t)
    if d not in got:raise RuntimeError(('literal source',a,r,t,d,got))
    literal_roots+=1
    for deletion_count,recipe in got.items():witnesses[t,deletion_count]=recipe
 for(t,d),(recipe,k)in sorted(witnesses.items()):
  a,r=build(recipe);a=sinks(a,k);p,s,q=profile(a,r)
  if(p,s,q)!=(t,1,d):raise RuntimeError(('oracle witness',t,d,recipe,k,p,s,q))
 collisions=[]
 for a,r,d,expected_output in [([2,28,25,17,1],0,1,(61,17,89)),([18,28,25,1,8],3,3,(57,15,87))]:
  core=profile(a,r);weak=sinks(a,1);p,s,q=profile(weak,r)
  if(p,s,q)!=((9,0,18),1,d):raise RuntimeError('weak collision')
  output,z=adjacent(weak,r);out,start,outd=profile(output,z)
  if out!=expected_output or start!=2:raise RuntimeError(('adjacent collision',out,start))
  collisions.append(dict(core_rows=a,root=r,core_profile=core,weak_rows=weak,deletion=d,output_rows=output,output_root=z,output_profile=out,output_start=start))
 if set(oracle(9,0,18))!={1,3}:raise RuntimeError('complete collision deletion set')
 result=dict(status='PASS',core_H_bound=N,core_profiles=core_count,weak_U_bound='8D',weak_integer_controls=weak_tests,labelled_graphs=literal_graphs,literal_start_one_roots=literal_roots,reconstructed_oracle_witnesses=len(witnesses),collision_deletions=sorted(oracle(9,0,18)),collisions=collisions,dependencies=dependencies)
 text=json.dumps(result,indent=2,sort_keys=True)+'\n';(P/'weak_start_one_certificate.json').write_text(text);print(text,end='')

if __name__=='__main__':main()
