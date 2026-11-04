"""Exact integer four-variable expansion, without a computer-algebra dependency."""
from math import comb
C={(6,6):-120,(6,5):-240,(6,4):-120,(5,6):-120,(5,5):-336,(5,4):-288,(5,3):-96,(4,4):396,(4,3):480,(4,2):84,(3,4):396,(3,3):1136,(3,2):720,(3,1):80,(2,3):384,(2,2):701,(2,1):270,(2,0):13,(1,2):113,(1,1):146,(1,0):28,(0,1):8,(0,0):8}
def need(x,msg):
 if not x:raise RuntimeError(msg)
def q(n,k):return sum((-1)**(k-i)*comb(k,i)*i**n for i in range(k+1))
def value(a,b):return sum(c*q(a,s)*q(b,t) for (s,t),c in C.items())
def expand():
 def var(i):a=[0]*4;a[i]=1;return {tuple(a):1}
 def add(*ps):
  z={}
  for p in ps:
   for e,c in p.items():z[e]=z.get(e,0)+c
  return {e:c for e,c in z.items() if c}
 def mul(*ps):
  z={(0,0,0,0):1}
  for p in ps:
   out={}
   for e,c in z.items():
    for f,d in p.items():
     g=tuple(a+b for a,b in zip(e,f))
     if g[0]<=3 and g[1]<=3:out[g]=out.get(g,0)+c*d
   z={e:c for e,c in out.items() if c}
  return z
 def sc(k,p):return {e:k*c for e,c in p.items() if k*c}
 one={(0,0,0,0):1};x,y,P,R=map(var,range(4));h=add(x,y,mul(x,y));al=mul(add(one,P),add(one,R));be=one;a=add(one,R,mul(P,R));b=add(one,P,mul(P,R));sg=mul(P,R);ka=add(P,mul(P,R));D=add(one,sc(-1,mul(sg,h)));V=y;W=add(mul(ka,h),x)
 Nxx=mul(x,add(one,sc(-1,mul(sg,y))));Nyy=mul(y,add(one,sc(-1,mul(sg,x))));Nxy=mul(x,y,add(one,sg));Nyx=mul(x,y,sg);Cx=x;Cy=mul(y,add(one,x));Rx=mul(x,add(one,y));Ry=y
 AA=add(mul(al,D),mul(a,add(Cx,mul(ka,h))));BB=add(mul(be,D),mul(b,Ry));TT=add(Nyx,mul(ka,Ry));Ex=mul(a,Cx);Ey=mul(a,Cy);Sx=mul(b,Rx);Sy=mul(b,Ry);rx=Nyx;ry=Nyy;lx=add(Nxx,mul(ka,Rx));ly=add(Nyx,mul(ka,Ry));direct=add(sc(2,mul(AA,BB,add(Nxy,sc(-1,Nyx)))),sc(-1,mul(add(D,sc(2,TT)),add(mul(Ex,Sy),sc(-1,mul(Ey,Sx))))),mul(AA,add(mul(rx,Sy),sc(-1,mul(ry,Sx)))),mul(BB,add(mul(Ex,ly),sc(-1,mul(Ey,lx)))))
 L=add(sc(2,mul(al,be)),sc(-1,mul(al,b)),sc(-1,mul(a,be)));J=add(mul(al,b,V),mul(a,be,W),mul(a,b,add(sc(-1,mul(ka,h)),sc(-1,mul(x,y)))));rational=mul(x,y,add(mul(L,D,D),mul(J,D),sc(2,mul(a,b,V,W))));need(direct==rational,'complete relevant four-term numerator identity')
 power=one;inv={}
 for r in range(5):inv=add(inv,sc(comb(r+2,2),power));power=mul(power,sg,h)
 polynomial=mul(direct,inv);weights=[0,1,2,2];out={}
 for (s,t,i,j),c in polynomial.items():out[i,j]=out.get((i,j),0)+c*weights[s]*weights[t]
 out={e:c for e,c in out.items() if c};need(out==C,'complete physical C3 profile extraction')
 p=sum(c*q(6,s)*q(6,t) for (s,t),c in C.items() if c>0);n=sum(c*q(6,s)*q(6,t) for (s,t),c in C.items() if c<0 and s>=4 and t>=4);need(p==2080753736 and n==-2560896000 and value(6,6)==-573454264,'exact domination constants')
 return {'terms':[[s,t,c] for (s,t),c in sorted(C.items())],'positive_at_6_6':p,'retained_negative_at_6_6':n,'strict_bound_numerator':p+n,'exact_at_6_6':value(6,6)}
