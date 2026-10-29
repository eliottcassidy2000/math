"""Exact minimum trace-clock carriers for the units in Ramanujan's examples."""
from pathlib import Path
from math import gcd
import json
HERE=Path(__file__).resolve().parent
CHECKS=0

def need(test,message):
 global CHECKS
 CHECKS+=1
 if not test:raise RuntimeError(message)

def mul(x,y,a):
 r,s=x;u,v=y
 return r*u+s*v,r*v+s*u+a*s*v

def power(a,e):
 x=(1,0)
 for _ in range(e):x=mul(x,(0,1),a)
 return x

def tr(x,a):return 2*x[0]+a*x[1]
def norm(x,a):return x[0]**2+a*x[0]*x[1]-x[1]**2

def record(a,es,m):
 bs=[power(a,e) for e in es]
 B=gcd(*(b[1] for b in bs));delta=gcd(2,a);disc=a*a+4
 invariants=(delta,disc*B//delta)
 expected=gcd(m,invariants[0])*gcd(m,invariants[1])
 image=set();kernel=[]
 for r in range(m):
  for s in range(m):
   x=(r,s);y=(tr(x,a)%m,)+tuple(tr(mul(b,x,a),a)%m for b in bs)
   image.add(y)
   if not any(y):kernel.append(x)
   for b in bs:
    y0=tr(x,a);y1=tr(mul(b,x,a),a);y2=tr(mul(b,mul(b,x,a),a),a)
    need(y2==tr(b,a)*y1-norm(b,a)*y0,'Cayley-Hamilton recurrence')
    need(y1*y1-tr(b,a)*y0*y1+norm(b,a)*y0*y0==-(a*a+4)*b[1]**2*norm(x,a),'quadratic conserved form')
 need(len(kernel)==expected,'exact Smith kernel count')
 need(len(image)*expected==m*m,'minimal observation cells')
 return {'a':a,'powers':es,'m':m,'smith':invariants,'fiber':expected,'cells':len(image),'kernel':kernel}

def main():
 rows=[]
 for a in (1,2,5,12):
  for es in ((1,),(2,),(6,),(12,),(4,6),(6,12),(5,7)):
   for m in range(1,31):record(a,es,m)
 for n,a,e in ((22,2,12),(37,12,6),(58,5,12)):
  A,B=power(a,e);delta=gcd(2,a)
  rows.append({'Ramanujan_n':n,'a':a,'power':e,'unit_power':[A,B],
               'discriminant':a*a+4,'smith':[delta,(a*a+4)*B//delta]})
 controls=[record(2,(12,),8),record(12,(6,),37),record(5,(12,),29),
           record(5,(12,),13),record(5,(1,),13),record(1,(2,3),6),
           record(1,(4,6),6)]
 need(tuple(v%13 for v in power(5,12))==(1,0),'58 twelfth-power clock is identity modulo13')
 need(tuple(v%11 for v in power(2,12))==(10,0),'22 twelfth-power clock is minus identity modulo11')
 for a,e,factors in ((2,12,2**5*3**2*5*7*11),(12,6,2**4*3**2*5*7**2*29*37),
                     (5,12,2**4*3**4*5*7*11**2*13*29)):
  need((a*a+4)*power(a,e)[1]==factors,'trace determinant prime factorization')
 out={'status':'FINITE-EXACT','universe':{'a':[1,2,5,12],'m':[1,30],
       'power_menus':[[1],[2],[6],[12],[4,6],[6,12],[5,7]]},
       'Ramanujan_clocks':rows,'controls':controls,'checks':CHECKS}
 (HERE/'unit_clocks.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({'checks':CHECKS,'Ramanujan_clocks':rows,
                    'controls':[{k:v for k,v in r.items() if k!='kernel'} for r in controls]},indent=2))
if __name__=='__main__':main()
