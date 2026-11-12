"""Independent literal graph controls for§6, no fusion-matrix evaluator."""
from check import H,hp
from charge_probe import count,f
from collections import Counter
from itertools import permutations
import json

def context(n,a):
 return [a<<1]+[((((1<<n)-1)^((1<<(i+1))-1))<<1)|(1-(a>>i&1)) for i in range(n)]
def fuse(A,B):
 n=len(A);out=[0]*(n+len(B)-1)
 for T,ids in [(A,list(range(n))),(B,[0]+list(range(n,len(out))))]:
  for i in range(len(T)):
   for j in range(len(T)):
    if T[i]>>j&1:out[ids[i]]|=1<<ids[j]
 for i in range(1,n):
  for j in range(n,len(out)):out[i]|=1<<j
 return out
T41=[56,13,25,48,34,6];T43=[52,13,24,49,34,6];T45=[52,13,24,33,42,6]
def V(t):
 n=len(t);N=1<<n
 hpall=[0]*N;hpall[0]=1;end=[[0]*n for _ in range(N)]
 for i in range(n):end[1<<i][i]=1
 for mask in range(1,N):
  for i in range(n):
   v=end[mask][i]
   if not v:continue
   rem=t[i]&~mask
   while rem:b=rem&-rem;rem^=b;end[mask|b][b.bit_length()-1]+=v
  hpall[mask]=sum(end[mask])
 return sum(end[mask][0]*hpall[(N-1)^mask] for mask in range(1,N))

def main():
 checks=Counter();trans=[]
 for n in range(6):
  for a in range(1<<n):
   C=context(n,a)
   left=H(fuse(C,T45))-H(fuse(C,T43));right=H(fuse(T45,C))-H(fuse(T43,C))
   pair=sum(f(j-i) for i in range(n) for j in range(i+1,n) if a>>i&1 and a>>j&1)
   assert left==2+4*pair
   assert (left==2)==(a.bit_count()<=1)
   assert (right==2)==(a==(1<<n)-1)
   checks['all_context_left_right_DP']+=1
   if a.bit_count()<=1:
    t=n-a.bit_length() if a else None
    shift=0 if t is None else 48 if t==0 else 57*(1<<t)
    counts=[H(fuse(C,T)) for T in [T41,T43,T45]]
    assert counts==[41+shift,43+shift,45+shift]
    checks['whole_three_value_translation_DP']+=1
    if a==1 or not a:trans.append(dict(n=n,cut=a,counts=counts))
 # All unit-deletion sandwiches up to two vertices on either side.
 for nl in range(3):
  for al in range(1<<nl):
   L=context(nl,al);q=1<<nl;A=sum((1<<i) for i in range(nl) if al>>i&1)
   for nr in range(3):
    for ar in range(1<<nr):
     R=context(nr,ar);B=sum((1<<i) for i in range(nr) if not(ar>>i&1))
     X=fuse(fuse(L,T43),R);Y=fuse(fuse(L,T45),R)
     dv=V(Y)-V(X);assert dv==2*A-16*q*B
     assert (dv==0)==(A==0 and B==0)
     if not dv:assert H(X)==43 and H(Y)==45
     checks['sandwich_receiver_formula_direct_DP']+=1
 frames=[(0,4,1,2,3,5),(0,2,4,1,3,5)]
 for lead in range(3):
  for t in range(4):
   n=lead+1+t;a=1<<lead;C=context(n,a)
   for (A,B),p in zip([(T41,T43),(T43,T45)],frames):
    AA=fuse(C,A);BB=fuse(C,B);old=[0 if i==0 else n+i for i in p]
    if t:frame=list(range(1,n+1))+old
    else:frame=list(range(1,n))+[0,n]+old[1:]
    assert sorted(frame)==list(range(n+6))
    assert all(AA[u]>>v&1 and BB[u]>>v&1 for u,v in zip(frame,frame[1:]))
    checks['transported_common_frames']+=1
 R159=[210,252,241,181,96,129,168,16]
 G159=fuse(context(2,1),T45);G157=fuse(context(2,1),T43)
 assert H(G159)==H(R159)==159 and V(G159)==215 and V(R159)==V(G157)==213
 for A,B,frame,expected in [(G159,R159,[1,2,5,7,4,6,3,0],2),(G157,R159,[0,1,2,4,6,3,5,7],3)]:
  diff=[(i,j) for i in range(8) for j in range(i+1,8) if (A[i]^B[i])>>j&1]
  assert len(diff)==expected
  assert all(A[u]>>v&1 and B[u]>>v&1 for u,v in zip(frame,frame[1:]))
  checks['finite_repair_frame_controls']+=1
 print(json.dumps(dict(checks=checks,translations=trans,status='PASS; all-order classification proved in PROOF.md, finite controls independently count literal graphs'),indent=2))
if __name__=='__main__':main()
