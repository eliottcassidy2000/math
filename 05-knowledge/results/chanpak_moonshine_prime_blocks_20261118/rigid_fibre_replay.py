from unit_gain_replay import K,L,R,state,forward,bits
from math import gcd

def profile(p,word):
 h,e,t,d=p
 for b in word:
  if b:h,e,t,d=h+t-e,d,2*t-d,d
  else:h,e,t,d=h,0,2*t,d
 return h,e,t,d
assert state(tuple(a&15 for a in L[:4]),3)==(3,0,12,3)
for j in range(1,13):
 M=16**j;B=4*M;u=(7*M-2)//5;V=(299*M+36-10*j)//10
 wl=B+4*B*V;wr=u+4*B*V
 rows=[];graphs=[]
 for A,p,w in ((L,(9,1,23,5),wl),(R,(9,1,25,5),wr)):
  word=bits(w);assert word[-2:]==[1,1]
  h,e,t,d=profile(p,word[:-1]);N=h+t-e
  rows.append((N,h-e,h+t,h))
  G=forward(A,3,word);graphs.append(G);n=len(G)
  assert [i for i,a in enumerate(G) if a.bit_count()==1]==[n-1]
  assert [i for i,a in enumerate(G) if a.bit_count()==n-2]==[4]
 a,b=rows;N=a[0]
 assert N==b[0] and a[2]==b[2]==N+5 and N%5==4
 assert gcd(N,N+5)==1
 assert a[3]-b[3]==32*M and gcd(a[3],b[3])==1
 assert 2*a[3]<N and 2*b[3]<N
 D4=3+12*wl-3*K(wl)
 assert 2*D4>N
 if j<=2:
  def isos(A,B,limit=2):
   n=len(A)
   if sorted(a.bit_count() for a in A)!=sorted(a.bit_count() for a in B):return 0
   match={};used=set();answer=0
   def rec():
    nonlocal answer
    if answer>=limit:return
    if len(match)==n:answer+=1;return
    best=None
    for i in range(n):
     if i in match:continue
     cs=[k for k in range(n) if k not in used and A[i].bit_count()==B[k].bit_count() and all(((A[i]>>u)&1)==((B[k]>>v)&1) for u,v in match.items())]
     if not cs:return
     if best is None or len(cs)<len(best[1]):best=i,cs
    i,cs=best
    for k in cs:
     match[i]=k;used.add(k);rec();used.remove(k);del match[i]
   rec();return answer
  X,Y=graphs
  assert isos(X,X)==isos(Y,Y)==1
  Yrev=tuple(((1<<len(Y))-1)^a^(1<<i) for i,a in enumerate(Y))
  assert isos(X,Y)==isos(X,Yrev)==0
  Xrev=tuple(((1<<len(X))-1)^a^(1<<i) for i,a in enumerate(X))
  assert isos(X,Xrev)==0
  print('literal structure',j,len(graphs[0]),N,'rigid/nonisomorphic/nonconverse PASS')
print('all-order fibre signature controls: PASS')
