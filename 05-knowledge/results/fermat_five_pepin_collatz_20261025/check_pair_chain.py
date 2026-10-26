# THM-4581's pair chain (k,e): u_n = 3^k v_n + e_n with v = x = 2^(2L+1)+1, u = T(y) = 3x-1, start (k,e)=(1,-1).
# Verify the chain against the direct orbits and report the state at Terras times 2L+1, 2L+2 for L mod 16 classes.
from fractions import Fraction as Fr
import json
def T(x): return x>>1 if x%2==0 else (3*x+1)>>1
def step(k,e,beta):
    # sigma = parity of 3^max(0,-k) e
    f=e*Fr(3)**max(0,-k); assert f.denominator==1
    sigma=int(f.numerator)%2
    if sigma==0 and beta==0: return (k, e/2)
    if sigma==0 and beta==1: return (k, (3*e+1-Fr(3)**k)/2)
    if sigma==1 and beta==0: return (k+1, (3*e+1)/2)
    return (k-1, (e-Fr(3)**(k-1))/2)
out={}
for L in range(1,34):
    x=2**(2*L+1)+1; y=2**(2*L+2)+1
    v=x; u=T(y); k,e=1,Fr(-1)
    assert u==3**k*v+e
    states={}
    absorbed=None
    for n in range(0, 6*L+300):
        if (k,e)==(0,0): absorbed=n; break
        if n in (2*L+1, 2*L+2, 2*L+3): states[n]=(k,str(e))
        beta=v%2
        k,e=step(k,e,beta)
        v=T(v); u=T(u)
        assert u==Fr(3)**k*v+e, (L,n)
        if v<=2: break
    out[L]={'L_mod_4':L%4,'state_at_2L+1':states.get(2*L+1),'state_at_2L+2':states.get(2*L+2),'state_at_2L+3':states.get(2*L+3),'absorbed_at_terras':absorbed}
print(json.dumps(out,indent=0))
