from collections import Counter
from math import gcd

def cc(A,B,d,N,with_edges=False):
    S=[i for i in range(N+1)]
    def root(i):
        while i!=S[i]:S[i]=S[S[i]];i=S[i]
        return i
    E=[];cycles=0
    for x in range(1,N+1):
        for tx in (d*A,d*B,d*(A+B)):
            y=tx-x
            if y<=0 or y>=x:continue
            p,q=root(x),root(y)
            if p==q:cycles+=1
            else:S[p]=q
            E.append((x,y,tx))
    out=Counter(root(i) for i in range(1,N+1))
    p=Counter(i for x,y,t in E for i in (x,y))
    return sorted(out.values()),cycles,p,E

def exact(A,B,d,N):
    def e(T):return max(0,min(N,(T-1)//2)-max(1,T-N)+1)
    beta=max(0,(d-1)//2-max(1,d*B-N)+1)
    return N-sum(e(T) for T in (d*A,d*B,d*(A+B)))+beta,beta

trials=0
for B in range(2,14):
 for A in range(1,B):
  if gcd(A,B)!=1 or (A+B)%2!=1:continue
  for d in range(1,8):
   for N in range(1,d*(A+B)+3):
    p,r,_,_=cc(A,B,d,N)
    assert (len(p),r)==exact(A,B,d,N)
    trials+=1
   p,r,deg,ed=cc(A,B,d,d*B-1)
   assert p==sorted([B-1]+[(2*B)]*((d-1)//2)+([B] if d%2==0 else []))
   assert not max(deg.values(),default=0)>2
   for r in range(1,(d-1)//2+1):
    f=[(x,y,t) for x,y,t in ed if x%d in (r,d-r)]
    assert Counter(t for x,y,t in f)=={d*A:A,d*B:B,d*(A+B):B-A}
    bits=[int((r-j*d*A)%(d*B)<d*A) for j in range(B)]
    assert sum(bits)==A
    for q in range(B):
     for ell in range(1,B+1):assert abs(B*sum((bits*2)[q:q+ell])-ell*A)<B
A,B,d=64,225,9
for N in (2016,2017,2018,2019,2020,2021,2022,2023,2024,2600,2700):
  p,beta,g,e=cc(A,B,d,N)
  assert (len(p),beta)==exact(A,B,d,N)
  print(N,'components',len(p),'cycle rank',beta,'component range',(p[0],p[-1]))
  if N==2016:assert p==[218]*4+[224]+[230]*4
  if 2016<=N<=2020:assert len(p)==2025-N and beta==0
  if N==2024:
   assert p==[224,450,450,450,450]
   assert {x for x,k in g.items() if k==1}=={288,576}
for x in (6,12,18,24,60,600):
 U,V=(x+7)*(x+3),(x+5)*(x+1); R,S=(x+7)*(x+5),(x+3)*(x+1)
 assert U*V==R*S and gcd(U,V)==gcd(R,S)==1
 a,b,c=(U*U-V*V)//2,U*V,(U*U+V*V)//2
 z,h,w=(R*R-S*S)//2,R*S,(R*R+S*S)//2
 assert a*a+b*b==c*c and z*z+h*h==w*w and a<z<b==h
print('PASS',trials,'all-height tiny graph controls plus 2024, balanced words, and collision families')
