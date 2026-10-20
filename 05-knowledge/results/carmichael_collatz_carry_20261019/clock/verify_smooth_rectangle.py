import math,cmath,json
rows=[]
def req(ok,msg):
 if not ok:raise RuntimeError(msg)
for D in [37,101,1009,10007,1000003,1000000007,1000000000001,1000000000000000003]:
 if math.gcd(D,6)!=1:continue
 L=0
 while 6**(2*(L+1))<=D:L+=1
 theta=math.pi/(L+1);b={i:math.sqrt(2/(L+1))*math.sin(i*theta) for i in range(1,L+1)}
 f={2**i*3**j:b[i]*b[j] for i in b for j in b}
 req(len(f)==L*L and max(f)<D,'frequency uniqueness')
 U={2*a%D:v for a,v in f.items()}
 V={2*(a//3)%D:v*cmath.exp(-2j*math.pi*(a//3)/D) for a,v in f.items()}
 V0={2*(a//3)%D:v for a,v in f.items()}
 req(all(a%3==0 for a in f),'integer quotient')
 def diff(x,y):return sum(abs(x.get(a,0)-y.get(a,0))**2 for a in set(x)|set(y))
 e=diff(U,f);o=diff(V,f);raw=diff(V0,f);delta=2*math.pi*6**L/(3*D)
 req(abs(e-2*(1-math.cos(theta)))<1e-12,'E energy')
 req(abs(raw-2*math.sin(theta)**2)<1e-12,'unphased O energy')
 req(diff(V,V0)<=delta*delta+1e-14,'phase error')
 ray=(e+o)/8;upper=(2*(1-math.cos(theta))+(math.sqrt(2)*math.sin(theta)+delta)**2)/8
 coarse=5*math.pi**2/(8*(L+1)**2)+math.pi**2*6**(2*L)/(9*D**2)
 req(ray<=upper+1e-12 and upper<=coarse+1e-12,'Rayleigh bounds')
 rows.append({'D':D,'L':L,'rayleigh_quotient':ray,'exact_upper_formula':upper,'coarse_upper_formula':coarse})
print(json.dumps({'status':'PASS','scope':'Integer frequency distinctness exact; trigonometric energy checks numerical at tolerance1e-12, general theorem proved separately','cases':rows},indent=2))
