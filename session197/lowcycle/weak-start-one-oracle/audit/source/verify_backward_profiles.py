from functools import lru_cache
from collections import deque
import json,hashlib
from pathlib import Path
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
N=61;bank=forward(N);checked=0
for h in range(1,N+1,2):
 for e in range(h+1):
  for R in range(1,2*h+1):
   w=inverse(h,e,R);checked+=1
   if (w is not None)!=((h,e,R)in bank):raise RuntimeError(('membership',h,e,R,w))
   if w is not None:
    t=(1,1,1)
    for L in w:t=step(t,L)
    if t!=(h,e,R):raise RuntimeError(('recipe',t,(h,e,R),w))
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
physical=0
for t in sorted(bank):
 a,r=build(inverse(*t));p,s0,d0=profile(a,r)
 h,e,R=t;predicted=1 if h==1 else R-h+e if 3*e>h or 3*e==h and R>h else e
 if d0!=predicted:raise RuntimeError(('deletion readout',t,d0,predicted))
 if p!=t or s0!=1:raise RuntimeError(('physical recipe',t,p,s0))
 physical+=1
bad=[8,1,3,6];p,s0,bad_d=profile(bad,2)
if p!=(5,1,9) or s0!=2 or bad_d!=3 or step(p,1)!=(13,4,14):raise RuntimeError('wrong-root hostile')
if inverse(5,1,7)!=(3,) or inverse(3,2,4)is not None:raise RuntimeError('formal inverse hostile')
expected={'ROOTED.md':'ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e','ARITHMETIC_INVERSE.md':'4ac1727eb27c85f8a39f5598eb9f0e1b020bc4020cbca0a94e6288cb0676c3b7','C1_COMPLETE.md':'d38a53b1c475df3578d58e83d9e9f97785d9455084c0a094aaed752577d8b289'}
for name,sha in expected.items():
 if hashlib.sha256((P/name).read_bytes()).hexdigest()!=sha:raise RuntimeError('dependency '+name)
large=[]
for original in ([4,17,1,1,4500],[96,3,4061],[21,2,1,2,4869],[3,1,18,1,9628]):
 t=(1,1,1)
 for L in original:t=step(t,L)
 recipe=inverse(*t)
 if recipe is None:raise RuntimeError('large inverse missing')
 replay=(1,1,1)
 for L in recipe:replay=step(replay,L)
 if replay!=t:raise RuntimeError('large inverse replay')
 large.append({'profile':t,'original_recipe':original,'inverse_recipe':recipe})
out={'large_full_profile_controls':large,'status':'PASS','bound':N,'forward_profiles':len(bank),'tested_triples':checked,'recursive_states':inverse.cache_info().currsize,'physical_recipe_checks':physical,'physical_deletion_readouts':physical,'strong_H_support':sorted({h for h,e,R in bank if h>1}),'wrong_root_hostile':{'rows':bad,'root':2,'profile':p,'start':s0,'actual_deletion_count':bad_d,'singleton_output':step(p,1)},'dependencies':expected}
text=json.dumps(out,sort_keys=True,indent=2)+'\n';(P/'backward_profile_certificate.json').write_text(text);print(text,end='')
