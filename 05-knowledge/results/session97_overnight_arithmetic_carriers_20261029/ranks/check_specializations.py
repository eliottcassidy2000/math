"""Exact specialization audit of the public formulas; does not prove generic rank."""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import json
from check_characters import inv,require
P=Path(__file__).resolve().parent

def family(i,t):
 if i==273:
  L=-10340482*t*t+48716518*t-4436170;p=1345905*t*t+6071604*t-43997;q=169*(1757*t*t+169808*t+10560)
  D=11340016558200*t**4+100399938810750*t**3+94651330477309*t*t+1291548172478*t-5344262336985
  E=5408*(6409666485*t**4+59220491056*t**3-72530579204*t*t+37829024608*t+5489594880)
 elif i==302:
  L=446667*t*t+471466*t+239031;p=318552*t*t+368554*t-72570;q=733413*t*t-45082*t-14960
  D=5*(7174492962*t**4-7114589515*t**3-22069002960*t*t+3909144679*t-205134150)
  E=882769396002*t**4+811447034567*t**3-1174040743*t*t-32493137198*t-2386325360
  # Put minus-sign presentation into common positive-sign convention.
  L=-L
 elif i==4239:
  L=466548*t*t-828066*t+32416;p=-478380*t*t-1076788*t+103008;q=82612*t*t-180205*t+81622
  D=-26666192028*t**4-91335649060*t**3+11566653808*t*t+239872553168*t-53374213248
  E=-47423161625*t**4+118525002574*t**3-104259187849*t*t+76235905664*t-24372268764
 else:raise ValueError(i)
 B=p*E+q*D+p*q*(L-p-q)
 return [L,D+E,B,D*E,Q(0)]

def fourth(q):
 n=isqrt(isqrt(q.numerator));d=isqrt(isqrt(q.denominator));r=Q(n,d);require(r**4==q);return r
out=[]
for i,t in [(273,Q(12473,32946)),(302,Q(164518,924945)),(4239,Q(-2868768,362209))]:
 a=family(i,t);b=list(map(Q,json.loads((P/f'curve{i}.json').read_text())['ainvs']));_,_,_,c4,c6,D=inv(a);_,_,_,d4,d6,E=inv(b)
 u=fourth(c4/d4);require(c6==u**6*d6);require(D==u**12*E)
 # x=u² X+r, y=u³ Y+s u² X+w sends minimal model to supplied family.
 r=(u*u*(b[0]**2+4*b[1])-(a[0]**2+4*a[1]))/12
 s=(u*b[0]-a[0])/2;w=(u**3*b[2]-a[0]*r-a[2])/2
 pts=json.loads((P/f'curve{i}.json').read_text())['points']
 for X,Y in pts:
  X,Y=Q(X),Q(Y);x=u*u*X+r;y=u**3*Y+s*u*u*X+w
  require(y*y+a[0]*x*y+a[2]*y==x**3+a[1]*x*x+a[3]*x+a[4])
 out.append({'curve':i,'T':str(t),'u':str(u),'r':str(r),'s':str(s),'w':str(w),'transported_points':len(pts)})
print(json.dumps(out,indent=2))
