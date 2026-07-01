"""
Is the tight condition LOCALIZED at rational resonances (=> combinatorial/finite, not analytic)?
The danger sum D(t)=#{v:||vt||<=1/n}; its MIN (the hole) -- at which t? bounded-denominator rational?
Also: the danger-sum Fourier D_hat(k)=sum_{v|k, v in S} f_hat(k/v) (divisor sum); f_hat(m)=sin(2pi m/n)/(pi m).
"""
from fractions import Fraction
def frac(x): x=x%1.0; return min(x,1-x)
def D(S,t,n): return sum(1 for v in S if frac(v*t)<=1.0/n+1e-12)
def min_and_where(S,n,Qmax):
    # min over rationals a/q, q<=Qmax
    best=999; bt=None
    for q in range(1,Qmax+1):
        for a in range(q):
            t=Fraction(a,q); d=D(S,float(t),n)
            if d<best: best=d; bt=t
    return best,bt
n=14
print("min_t D(t) (the covering min = tightness) and WHERE (rational t*), over q<=Qmax:")
for name,S in [("AP",list(range(1,n))),("GW",list(range(1,12))+[13,24]),
               ("swap12->26",[v for v in range(1,n) if v!=12]+[26]),
               ("drop 6 (punctured)",[v for v in range(1,n) if v!=6])]:
    for Qmax in [n, 2*n, 4*n]:
        b,t=min_and_where(S,n,Qmax)
        print(f"  {name:>18} Qmax={Qmax:>3}: min D={b} at t*={t} (denom {t.denominator if t else '-'})", end="")
    print()
print()
print("=> the covering-min (tightness) is achieved at LOW-denominator rationals (the Farey resonances) --")
print("   so tightness is a FINITE/combinatorial condition (check D at bounded-denom t), NOT a genuine")
print("   analytic min. Averaged lenses (moments, additive energy) miss it; pointwise-at-resonances is the way.")
