"""Exact checks of a restricted all-order physical residual replacement."""
from itertools import permutations
from pathlib import Path
import hashlib,json
from core import inventories,cycles,adjacency_map
HERE=Path(__file__).resolve().parent
GEN={
'O':((1,0,1,0),(1,2,1,2),(0,0,0,0),(0,0,1,1)),
'Z':((2,1,1,1),(0,1,0,0),(0,0,1,1),(0,0,0,0)),
'J':((0,1,0,1),(0,0,0,0),(1,0,1,1),(0,0,1,0))}

def req(ok,msg):
 if not ok:raise RuntimeError(msg)
def mv(M,x):return tuple(sum(a*b for a,b in zip(row,x))for row in M)
def old(x):
 a,f,e,p=x;return a+2*e+p,e,a+f+2*e+2*p,e+p
def new(x):
 D,E,U,d=x;return D-E-d,U-D+E-d,E,d-E
def p2(n):return n>0 and n&(n-1)==0
def append(T,c,l):
 n=len(T);A=list(T)+[0]
 for v in range(n):
  if l in('Z','J')and v==c:A[n]|=1<<v
  else:A[v]|=1<<n
 return A,n if l=='J'else c
def word(T,c,w):
 for l in w:T,c=append(T,c,l)
 return T,c
def profile(T,c):
 q=inventories(T);return q['H'],q['e'][c],q['R'][c],q['allH'][((1<<len(T))-1)^(1<<c)]
def tournaments(n):
 pairs=[(i,j)for i in range(n)for j in range(i+1,n)]
 for mask in range(1<<len(pairs)):
  T=[0]*n
  for k,(a,b)in enumerate(pairs):
   if mask>>k&1:T[a]|=1<<b
   else:T[b]|=1<<a
  yield T
def induced(T,V):return[sum(((T[u]>>v)&1)<<i for i,v in enumerate(V))for u in V]
def strip_tail(T,c,k):
 kept=list(range(len(T)));tail=[]
 for _ in range(k):
  mask=sum(1<<v for v in kept);sinks=[v for v in kept if not T[v]&mask]
  req(len(sinks)==1 and sinks[0]!=c,'canonical sink tail')
  z=sinks[0];kept.remove(z);tail.append(z)
 W=induced(T,kept);r=kept.index(c);D,E,U,d=profile(W,r)
 req(inventories(W)['s'][r]==1 and D==U,'balanced inverse core')
 X,b=word(W,r,'O'*k)
 req(adjacency_map(X,T,kept+list(reversed(tail)))and b==r,'all-edge tail reconstruction')
 return W,r

def main():
 cap=2000;base=(0,0,0,1);bank={base};q=[base]
 for x in q:
  for M in GEN.values():
   y=mv(M,x)
   if old(y)[2]<=cap and y not in bank:bank.add(y);q.append(y)
 req(len(bank)==32548,'inherited complete physical bank size')
 dyadic=set();core_images=set()
 for x in bank:
  D,E,U,d=old(x)
  if E:req(U<=2*(D-E)<2*D,'strong core strict ratio bound')
  if U%D==0 and p2(U//D)and U>D:
   k=(U//D).bit_length()-1
   req(E==0,'dyadic source weak')
   core=(D,d,D,d)
   req(core==(1,1,1,1)or new(core)in bank,'physical balanced core')
   dyadic.add((D,E,U,d))
 cores=[(1,1,1,1)]+[old(x)for x in bank if old(x)[0]==old(x)[2]]
 for D,E,U,d in cores:
  k=1
  while 2**k*D<=cap:
   core_images.add((D,0,2**k*D,d));k+=1
 req(dyadic==core_images,'entire capped dyadic profile iff')
 graphs=ratio_graphs=0
 for n in range(1,6):
  for T in tournaments(n):
   q=inventories(T)
   for c in range(n):
    if q['s'][c]!=1:continue
    D,E,U,d=profile(T,c)
    if U%D==0 and U>D and p2(U//D):
     strip_tail(T,c,(U//D).bit_length()-1);graphs+=1
    if n>1 and D==U and E and D%E==0 and p2(D//E-1)and D//E-1>=4:
     k=(D//E-1).bit_length()-2
     # Invert W=J(O(V)) canonically, then remove V's k sinks.
     cs=[v for v in range(n)if T[c]>>v&1];req(len(cs)==1,'balanced root sole neighbor');a=cs[0]
     kept=[v for v in range(n)if v!=c];mask=sum(1<<v for v in kept)
     zs=[v for v in kept if not T[v]&mask];req(len(zs)==1,'balanced inverse sink')
     z=zs[0];kept.remove(z);V=induced(T,kept);r=kept.index(a)
     C,b=strip_tail(V,r,k)
     X,b=word(C,b,'O'*(k+1)+'J')
     req(profile(X,b)==(D,E,U,d),'ratio-fiber response reconstruction')
     ratio_graphs+=1
 literal=[];cores=[([0],0)]+[word([2,0],0,w+'OJ')for w in('','O','Z','J')]
 for T,c in cores:
  h,E,U,d=profile(T,c);req(h==U,'balanced literal core')
  for k in range(1,4):
   A,a=word(T,c,'O'*(k+1)+'JO'+'Z'*k+'J')
   B,b=word(T,c,'O'*(2*k+3)+'J');qa,qb=inventories(A),inventories(B)
   target=(2**(2*k+3)+1)*h
   req(len(A)==len(B)==len(T)+2*k+4,'same order')
   req(qa['H']==qb['H']==target,'same H')
   req(cycles(A)==cycles(B)==1,'both C1')
   req(not any(qa['s'][v]==1 and qa['R'][v]==target for v in range(len(A))),'no old balanced root')
   req(qb['s'][b]==1 and qb['R'][b]==target,'new balanced root')
   literal.append({'core_profile':[h,E,U,d],'k':k,'order':len(A),'H':target})
 for k in range(1,41):
  wa='O'*(k+1)+'JO'+'Z'*k+'J';wb='O'*(2*k+3)+'J'
  req(len(wa)==len(wb)==2*k+4,'symbolic same order')
  for seed,expected in [((1,0,0,0),2**(2*k+3)+1),((-2,0,1,0),0)]:
   for w in(wa,wb):
    x=seed
    for l in w:x=mv(GEN[l],x)
    req(old(x)[0]==expected,'exact balanced-input Hamilton row')
  # The replacement's cover and Hamilton rows agree on the balanced plane.
  for seed in[(1,0,0,0),(-2,0,1,0)]:
   x=seed
   for l in wb:x=mv(GEN[l],x)
   req(old(x)[0]==old(x)[2],'symbolic balanced replacement')
 result={'status':'PASS','source_U_cap':cap,'nontrivial_source_profiles':len(bank),
  'dyadic_source_profiles':len(dyadic),'canonical_dyadic_graph_controls_through_order5':graphs,
  'canonical_ratio_graph_controls_through_order5':ratio_graphs,'literal_graph_pairs':literal,
  'symbolic_k_range':[1,40],
  'files_sha256':{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest()for n in('PROOF.md','verify.py','core.py')}}
 print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
