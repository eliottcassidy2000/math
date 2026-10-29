"""Independent coefficient arrays, discriminants of cubics, exhaustive x,y enumeration."""
import json
import sympy as s
T,X=s.symbols('T X');p=233
# Descending coefficients of L,p,q,D,E, transcribed from ICARM family equations.
DATA={
273:[[-10340482,48716518,-4436170],[1345905,6071604,-43997],[169*1757,169*169808,169*10560],[11340016558200,100399938810750,94651330477309,1291548172478,-5344262336985],[5408*z for z in [6409666485,59220491056,-72530579204,37829024608,5489594880]]],
302:[[-446667,-471466,-239031],[318552,368554,-72570],[733413,-45082,-14960],[5*z for z in [7174492962,-7114589515,-22069002960,3909144679,-205134150]],[882769396002,811447034567,-1174040743,-32493137198,-2386325360]],
4239:[[466548,-828066,32416],[-478380,-1076788,103008],[82612,-180205,81622],[-26666192028,-91335649060,11566653808,239872553168,-53374213248],[-47423161625,118525002574,-104259187849,76235905664,-24372268764]]}
def require(v,msg):
 if not v:raise ValueError(msg)
out={}
for id,rows in DATA.items():
 L,a,b,D,E=[s.Poly.from_list(r,T).as_expr() for r in rows];B=s.expand(a*E+b*D+a*b*(L-a-b))
 co=[s.Poly(z,T)for z in [L,D+E,B,D*E,0]]
 # Complete the square: Y=2y+Lx+B, Y²=F. disc(F)=16 Delta.
 F=4*X*(X+D)*(X+E)+(L*X+B)**2
 delta=s.Poly(s.discriminant(F,X)/16,T)
 require(delta.degree()==24,'degree over Q')
 require(s.gcd(delta,delta.diff()).degree()==0,'multiple singular fiber over Q')
 dm=s.Poly(delta.as_expr(),T,modulus=p)
 require(dm.degree()==24 and s.gcd(dm,dm.diff()).degree()==0,'not etale discriminant')
 # Explicit c4 check, redundant but guards against cusp and bad model.
 c4=s.Poly((L*L+4*(D+E))**2-24*(L*B+2*D*E),T,modulus=p)
 require(s.gcd(dm,c4).degree()==0,'cuspidal fiber')
 vals=[]
 for t in range(p+1):
  v=[int(c.eval(t))%p for c in co] if t<p else [int(c.nth(w))%p for c,w in zip(co,[2,4,6,8,12])]
  a1,a2,a3,a4,a6=v;count=1
  for x in range(p):
   rhs=(x*x*x+a2*x*x+a4*x+a6)%p;by=(a1*x+a3)%p
   count+=sum((y*y+by*y-rhs)%p==0 for y in range(p))
  vals.append(count)
 expected={273:58509,302:58404,4239:58374}[id]
 require(sum(vals)==expected,'surface count mismatch')
 out[id]={'count':sum(vals),'infinity_count':vals[-1],'finite_rational_bad_fibers':[t for t in range(p)if dm.eval(t)%p==0],'delta_infinity':int(dm.LC())%p,'fibers':vals,'H2_trace':sum(vals)-1-p*p}
print(json.dumps({'prime':p,'checks':'Q degree24 squarefree; mod233 degree24 squarefree and coprime c4; all x,y points; infinity included','surfaces':out},indent=2))
