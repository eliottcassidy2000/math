"""Standard-library replay: exact prime collision and quadratic recurrence bank."""
from pathlib import Path
import json,math,sys
ROOT=Path(__file__).resolve().parents[2]
for name in ('joint_placement','chanpak_joint_placement_20261117'):
 if (ROOT/name).is_dir():sys.path.insert(0,str(ROOT/name/'arithmetic'))
from collisions import family,H
CERT=json.loads(Path(__file__).with_name('PRIME_CERTIFICATE.json').read_text())
def prime(p):
 if p<=100:
  return p>=2 and all(p%d for d in range(2,math.isqrt(p)+1))
 c=CERT[str(p)];fs={int(q):e for q,e in c['factors'].items()};a=c['base']
 if math.prod(q**e for q,e in fs.items())!=p-1:return False
 return all(prime(q) for q in fs) and pow(a,p-1,p)==1 and all(math.gcd(pow(a,(p-1)//q,p)-1,p)==1 for q in fs)
A=1267674777822284490394719748096
B=135388433329826992581696
C=3614891648909077
def polynomial(j):return A*4**j+B*2**j+C
def decode(N):
 d=B*B+4*A*(N-C)
 if d<0:return None
 r=math.isqrt(d)
 if r*r!=d or (r-B)%(2*A):return None
 x=(r-B)//(2*A)
 if x<1 or x&(x-1):return None
 return x.bit_length()-1
P=42284279505304223759
u,w,v,q=family(4,24)
if not prime(P) or H(u,w)!=P or H(v,q)!=P:raise ArithmeticError('prime collision')
if H(v,q+2)-H(u,w+2)!=8212:raise ArithmeticError('response')
for j in range(150):
 n=H(*family(4,2**(22+j))[:2])
 if n!=polynomial(j) or decode(n)!=j or decode(n+2) is not None:raise ArithmeticError(j)
D=B*B-4*A*C
if D<=0 or math.isqrt(D)**2==D:raise ArithmeticError('irreducibility')
g=math.gcd(polynomial(0),polynomial(1),polynomial(2))
if g!=1:raise ArithmeticError('fixed divisor')
print(json.dumps({'prime_collision':P,'addresses':[u,w,v,q], 'continuations':[H(u,w+2),H(v,q+2)],'prime_certificate_nodes':len(CERT),'exact_dyadic_inverse_controls':150,'fixed_divisor':g,'discriminant_floor_sqrt':math.isqrt(D)},indent=2))
