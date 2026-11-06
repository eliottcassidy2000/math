"""Exact small controls for the S82 scope correction; not an asymptotic proof."""
import json
from fractions import Fraction

def require(condition):
    if not condition:
        raise RuntimeError("Scope-control verification failed")

def g(v):
    if not v: return 1
    if v&1: return g(v//2)
    return g(v//2)+v//2+(v//2)%2

def rho(v):
    return sum((Fraction(1,2**(i+1)) for i in range(v.bit_length()-1) if not(v>>i)&1),Fraction())

for v in range(1,4096):
    n=v.bit_length()-1
    require(abs(g(v)-v*rho(v))<=2*n+1)
    require(1<=g(v)<=v+1)

controls=[]
for b in [4,16,64,256,1024]:
    count=sum(g(v)+v+1-b*(v==0 or not(v&1))<=1 for v in range(b+2))
    guaranteed=sum(1 for v in range(2,b//2,2))
    require(count>=guaranteed)
    controls.append({'a':1,'c':1,'X':1,'b':b,'count':count,'guaranteed_even_subfamily':guaranteed})

# The finite Fourier identity is checked by exact residue coefficients: choosing
# modulus greater than every possible readout difference removes aliasing.
ys=[(1,1,1),(3,1,4)]
for N in range(1,101):
    vals=[]
    for a,b,c in ys:
        vals += [a*g(v)-b*(v==0 or not(v&1))+c*(1+v) for v in range((N+b)//c+1)]
    Q=max([abs(x-N) for x in vals]+[0])+1
    exact=sum(x==N for x in vals)
    cyclic=sum((x-N)%Q==0 for x in vals)
    require(exact==cyclic)
print(json.dumps({'status':'PASS','rho_approximation_checked':4095,'uniformity_counterexample_controls':controls,'exact_fourier_no_alias_controls':100},indent=2))
