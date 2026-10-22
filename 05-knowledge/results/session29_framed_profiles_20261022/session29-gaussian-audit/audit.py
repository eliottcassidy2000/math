from fractions import Fraction as Q
from math import comb, factorial
import importlib.util, hashlib, json
from pathlib import Path
P=Path(__file__).parent
spec=importlib.util.spec_from_file_location('frozen',P/'producer_frozen.py'); m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def endpoints(z): return Q(z.l,m.S),Q(z.h,m.S)
def contains(z,a,b=None):
 l,h=endpoints(z); assert l<=a<=(a if b is None else b)<=h,(l,a,b,h)
# Exact rational arithmetic controls, including signed reciprocal intervals.
vals=[Q(i,j) for i in range(-5,6) for j in (1,3,7)]
cases=0
for a in vals:
 for b in vals:
  if a>b:continue
  z=m.I(a,b);contains(z,a,b)
  if not a<=0<=b:contains(z.inv(),1/b,1/a)
  for c,d in [(Q(-7,3),Q(-1,7)),(Q(-2,7),Q(5,3)),(Q(1,7),Q(9,7))]:
   w=m.I(c,d);contains(z+w,a+c,b+d)
   terms=[a*c,a*d,b*c,b*d];contains(z*w,min(terms),max(terms));cases+=1
# Independent long alternating exp(-x) sums at rational arguments.
for x in (Q(0),Q(1),Q(6),Q(7),Q(14),Q(16)):
 s=sum(((-x)**j/Q(factorial(j)) for j in range(181)),Q(0))
 nxt=(-x)**181/Q(factorial(181)); contains(m.expminus(m.I(x)),min(s,s+nxt),max(s,s+nxt))
# Independent Delta from Euler product, not E4/E6 identity.
N=m.N
def conv(a,b):
 out=[0]*(N+1)
 for i,ai in enumerate(a):
  for j,bj in enumerate(b[:N+1-i]):out[i+j]+=ai*bj
 return out
prod=[1]+[0]*N
for n in range(1,N+1):
 fac=[0]*(N+1)
 for j in range(min(24,N//n)+1):fac[j*n]=(-1)**j*comb(24,j)
 prod=conv(prod,fac)
delta=[0]+prod[:-1];assert delta==m.delta
# Rebuild F by repeated multiplication independently using product Delta.
pows4=[[1]+[0]*N]
for j in range(11):pows4.append(conv(pows4[-1],m.e4))
powsd=[[1]+[0]*N]
for j in range(3):powsd.append(conv(powsd[-1],delta))
f=[0]*(N+1)
for j,c in enumerate((1,-2640,1813680,-244992000)):
 term=conv(powsd[j],pows4[11-3*j]);f=[a+c*b for a,b in zip(f,term)]
assert f==m.f
# Derivatives independently via Leibniz rule, not recurrence.
for order,abc in enumerate(m.DER):
 for k in range(N-3):
  source=(m.C[k],m.B[k],m.A[k]); derived=[]
  for degree in range(3):
   v=sum(comb(order,j)*(-k)**(order-j)*factorial(degree+j)//factorial(degree)*source[degree+j] for j in range(min(order,2-degree)+1))
   derived.append(v)
  assert tuple(reversed(derived))==tuple(a[k] for a in abc)
# Explicit derivative-majorant constants. n>=85; each quotient decreases in n.
n=Q(85);a=n**3+5*n*n;b=45*(n*n+5*n);c=45*n
p0=256*a+16*b+c;p1=32*a+b;p2=2*a
constants=[]
for d in range(4):
 bound=n**d*p0+(d*n**(d-1)*p1 if d else 0)+(comb(d,2)*n**(d-2)*p2 if d>=2 else 0)
 constants.append(str(bound/n**6));assert bound<20000*n**6
assert p0<2000*n**3 and n*p0+p1<3000*n**4
assert Q(1,5)*(Q(86,85))**6<Q(1,4)
assert 2*Q(1,factorial(61))<Q(1,2**250)
# Numerator tails: x>=6 so V+45U/x <= n²+7.5n.
assert n*n+Q(15,2)*n<20000*n**6
out={'status':'PASS','source_sha256':hashlib.sha256((P/'producer_frozen.py').read_bytes()).hexdigest(),'rational_interval_control_cases':cases,'exp_alternating_controls':[0,1,6,7,14,16],'independent_delta_and_theta_degree':N,'independent_stationary_derivatives_orders':[0,1,2,3],'tail_n85_derivative_constants_over_n6':constants,'common_tail_bound':str(m.errHigher)}
(P/'audit.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
