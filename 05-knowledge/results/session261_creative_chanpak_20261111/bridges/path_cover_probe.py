"""Exact small hostile controls for new bridge statements, not a spectrum search."""
import itertools,json
from functools import lru_cache
from pathlib import Path

def rows(n,b):
 r=[0]*n
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  if b>>k&1:r[i]|=1<<j
  else:r[j]|=1<<i
 return tuple(r)
def sub(r,vs):
 return tuple(sum(((r[u]>>v)&1)<<j for j,v in enumerate(vs)) for u in vs)
@lru_cache(None)
def H(r):
 n=len(r)
 if not n:return 1
 dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   if not dp[s][v]:continue
   a=r[v]&~s
   while a:
    q=a&-a; a-=q; w=q.bit_length()-1
    dp[s|q][w]+=dp[s][v]
 return sum(dp[-1])
def clone(r,z,m):
 n=len(r); vs=[v for v in range(n) if v!=z]; out=[0]*(n-1+m)
 for i,u in enumerate(vs):
  for j,v in enumerate(vs):out[i]|=((r[u]>>v)&1)<<j
  for j in range(n-1,n-1+m):
   if r[u]>>z&1:out[i]|=1<<j
   else:out[j]|=1<<i
 for i in range(n-1,n-1+m):
  for j in range(i+1,n-1+m):out[i]|=1<<j
 return tuple(out)
def graph(r,z):
 vs=[v for v in range(len(r)) if v!=z]; A=[v for v in vs if r[z]>>v&1]; B=[v for v in vs if not r[z]>>v&1]
 es=[(a,b) for a in A for b in B if r[a]>>b&1]
 d=0
 for s in range(1<<len(vs)):
  S={v for i,v in enumerate(vs) if s>>i&1}
  if all(v in S or any((a==v and b in S) or (b==v and a in S) for a,b in es) for v in vs):d+=1
 return A,B,es,d

def main():
 tests=0; collisions={}; witness=None; biclique=0; sharp=0
 for n in range(1,6):
  for bits in range(1<<(n*(n-1)//2)):
   r=rows(n,bits); h=H(r)
   for z in range(n):
    U=sub(r,[v for v in range(n) if v!=z]); h0=H(U); A,B,es,d=graph(r,z)
    vals=[H(clone(r,z,m)) for m in (1,2,3)]
    for m,hm in enumerate(vals,1):
     if hm<h0+2**(m-1)*(h-h0):raise RuntimeError(('clone bound',r,z,m))
    key=(n,h,h0,d)
    if key in collisions and collisions[key]['clone2']!=vals[1] and witness is None:
     witness={'key_n_H_Hdelete_d':key,'first':collisions[key],'second':{'rows':r,'root':z,'clone2':vals[1],'clone3':vals[2]}}
    collisions.setdefault(key,{'rows':r,'root':z,'clone2':vals[1],'clone3':vals[2]})
    if A and B and len(es)==len(A)*len(B):
     biclique+=1
     if h-h0<d-1:raise RuntimeError('biclique Q')
     if h-h0==d-1:sharp+=1
    tests+=1
 return {'scope':'All labelled tournaments orders 1..5, every root; actual clone sizes 1..3','rooted_cases':tests,'complete_biclique_cases':biclique,'sharp_biclique_cases':sharp,'scalar_continuation_collision':witness}
if __name__=='__main__':
 out=main(); p=Path(__file__).with_name('path_cover_receipt.json');p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

# Finite matrix certificate proves the all-k formula, independently of fitting.
def alternating_certificate():
 from fractions import Fraction as F
 states=list(itertools.product(range(3),repeat=3)); idx={s:i for i,s in enumerate(states)}
 def step(s,c):
  out=[]
  for i in range(3):
   if i>0 and s[i]==0 and c==2:continue
   t=list(s);t[i]=c;out.append(idx[tuple(t)])
  return out
 M=[[0]*27 for _ in states]
 for s,j in idx.items():
  for a in step(s,1):
   for b in step(states[a],2):M[b][j]+=1
 f=[int(s[0]!=1 and s[1]==2) for s in states]
 def rowmul(v):return [sum(v[i]*M[i][j] for i in range(27)) for j in range(27)]
 q=f[:]
 for lam in (1,1,4,4,9):
  mm=rowmul(q);q=[mm[i]-lam*q[i] for i in range(27)]
 if any(q):raise RuntimeError('all-order row annihilator failed')
 v=[int(s==(0,0,0)) for s in states]; seq=[]
 for k in range(13):
  a2=sum(f[i]*v[i] for i in range(27))
  formula=F(9**(k+1),100)+F((60*k+8)*4**k,675)-F(12*k+11,108)
  if a2!=formula:raise RuntimeError(('closed formula',k,a2,formula))
  h1=F(4**(k+1)+6*k+5,9); h2=h1+2*a2
  if k<=3:
   # Outside vertices 0..2k-1 in their transitive order, A even, B odd.
   nr=2*k+1; rr=[0]*nr
   for i in range(2*k):
    for j in range(i+1,2*k):rr[i]|=1<<j
    if i%2==0:rr[-1]|=1<<i
    else:rr[i]|=1<<(nr-1)
   if H(tuple(rr))!=h1 or H(clone(tuple(rr),nr-1,2))!=h2:raise RuntimeError('literal family')
  seq.append({'k':k,'a2':a2,'H_two_apex':int(h2)})
  v=[sum(M[i][j]*v[j] for j in range(27)) for i in range(27)]
 return {'identity':'f(M-I)^2(M-4I)^2(M-9I)=0, exact 27-coordinate integer check','family':seq}

if __name__=='__main__':
 p=Path(__file__).with_name('path_cover_receipt.json'); out=json.loads(p.read_text());out['alternating_certificate']=alternating_certificate();p.write_text(json.dumps(out,indent=2)+'\n');print('All-order matrix certificate and literal family controls passed.')

def component_controls():
 import random,math
 rng=random.Random(2611705); checks=[]
 for sizes in [[(1,1)]*k for k in range(1,5)]+[[(1,2),(2,1)],[(2,2),(1,1)],[(1,1),(1,2),(1,1)]]:
  av=[];bv=[];n=0
  for p,q in sizes:
   av.append(list(range(n,n+p)));n+=p
   bv.append(list(range(n,n+q)));n+=q
  A=sum(av,[]);B=sum(bv,[]);which={v:i for i,bb in enumerate(av) for v in bb};which.update({v:i for i,bb in enumerate(bv) for v in bb})
  E=[((2**p-1)*(2**q-1)+1)//2 for p,q in sizes]
  bound=sum(2*math.factorial(s.bit_count())*math.prod(E[i] for i in range(len(E)) if s>>i&1) for s in range(1,1<<len(E)))
  vals=[]
  for rep in range(8):
   r=[0]*(n+1)
   for u,v in itertools.combinations(range(n),2):
    if (u in A)==(v in A):forward=rng.randrange(2)
    else:
     a,b=(u,v) if u in A else (v,u);winner=a if which[a]==which[b] else b;forward=(winner==u)
    if forward:r[u]|=1<<v
    else:r[v]|=1<<u
   for a in A:r[n]|=1<<a
   for b in B:r[b]|=1<<n
   delta=H(tuple(r))-H(sub(r,list(range(n))))
   if delta<bound:raise RuntimeError(('component bound',sizes,r,delta,bound))
   vals.append(delta)
  checks.append({'component_sizes':sizes,'factorial_lower_bound':bound,'sample_deltas':vals})
 return {'seed':2611705,'trials_per_shape':8,'cases':checks}

if __name__=='__main__':
 p=Path(__file__).with_name('path_cover_receipt.json'); out=json.loads(p.read_text());out['component_controls']=component_controls();p.write_text(json.dumps(out,indent=2)+'\n');print('Component-block controls passed.')
