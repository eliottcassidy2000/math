"""
monad-explorer-2026-06-07-S5   (HYP-2298 Q2: characterize the 6+6k unit count)
============================================================================
CLAIM (to verify, exact integer arithmetic):  for the rank-4 bridge lattice
    L_t = Z[zeta6] (+) Z[zeta6]*w_t ,   w_t = ((2t-1)+i sqrt(4t-1))/(2t), |w_t|=1,
the number of UNIT VECTORS (lattice points z with |z|=1) is EXACTLY

        #units(L_t) = 12 + r_E(t),

where r_E(t) = #{ alpha in Z[zeta6] : N(alpha)=t } is the Eisenstein representation
count of t  (r_E(t) = 6*(d_{1,3}(t) - d_{2,3}(t)), the mod-3 divisor excess).

MECHANISM (proved by hand, checked here):
  * 6 vectors  : the triangular rosette  zeta6^j           (beta=0)        [always]
  * 6 vectors  : the w_t rosette         zeta6^j * w_t     (alpha=0)        [always]
  * r_E(t)     : z = alpha*(1 - w_t),  N(alpha)=t          (BOTH nonzero)
    because  1 - w_t = (1 - i sqrt(4t-1))/(2t)  has  |1-w_t|^2 = (1+(4t-1))/(4t^2) = 1/t,
    so |z|^2 = N(alpha)/t = 1  iff N(alpha)=t.
  No other unit vectors exist (alpha,beta forced collinear with rational ratio
  lambda; only lambda=-1 keeps beta in Z[zeta6] at the required norm).

The count depends on t ONLY through its splitting in Q(sqrt-3) -- the FIRST CM
field (triangular side).  The glued second direction sqrt(4t-1) does NOT affect
the count.  That is the sharpening of the seed's 'splitting in Q(sqrt-3,sqrt(4t-1))'.
"""
from itertools import product
from math import isqrt

def is_square(n):
    if n<0: return False
    r=isqrt(n); return r*r==n

def rE(t):
    """#{Eisenstein integers of norm t} = 6*(#div=1 mod3 - #div=2 mod3)."""
    if t<=0: return 0
    exc=0
    for dvd in range(1,t+1):
        if t%dvd==0:
            if dvd%3==1: exc+=1
            elif dvd%3==2: exc-=1
    return 6*exc

def count_units_exact(t, R):
    """EXACT count of lattice points of L_t on the unit circle.
       basis of values: {1, sqrt3, sqrt(d), sqrt(3d)}, d=4t-1. scale coords by S=4t."""
    d=4*t-1; r1,r2,r12=3,d,3*d
    MT={(0,0):(1,0),(0,1):(1,1),(0,2):(1,2),(0,3):(1,3),
        (1,1):(r1,0),(1,2):(1,3),(1,3):(r1,2),(2,2):(r2,0),(2,3):(r2,1),(3,3):(r12,0)}
    def mul(x,y):
        r=[0,0,0,0]
        for i in range(4):
            if x[i]==0: continue
            for j in range(4):
                if y[j]==0: continue
                a,b=(i,j) if i<=j else (j,i); c,idx=MT[(a,b)]; r[idx]+=x[i]*y[j]*c
        return r
    S=4*t; SS=S*S; m=2*t-1
    # S*(re,im) integer tuples over basis {1,sqrt3,sqrtd,sqrt3d}:
    #  1        -> re=(S,0,0,0)      im=0
    #  zeta6    -> S*zeta6=2t + 2t i sqrt3 /2 = 2t*( (1+i sqrt3)/2 )*? : S=4t, S*zeta6=4t*(1/2,..)=2t + 2t*sqrt3 i? 
    #     zeta6=(1+i sqrt3)/2 ; S*zeta6 = 2t + 2t i sqrt3  -> re=(2t,0,0,0) im=(0,2t,0,0)
    #  w_t=(m+i sqrtd)/(2t); S*w_t=2*(m+i sqrtd)=2m + 2 i sqrtd -> re=(2m,0,0,0) im=(0,0,2,0)
    #  zeta6*w_t : S*(zeta6 w_t)=4t*zeta6*w_t = 2*zeta6*(2t w_t)=... compute zeta6*(m+i sqrtd)= (1+i sqrt3)/2*(m+i sqrtd)
    #     = [ (m - sqrt(3d)) + i (sqrtd + m sqrt3) ] /2 ; times (S/(2t))=2 :
    #     = (m - sqrt3d) + i(sqrtd + m sqrt3)  -> re=(m,0,0,-1) im=(0,m,1,0)
    RE=[(S,0,0,0),(2*t,0,0,0),(2*m,0,0,0),(m,0,0,-1)]
    IM=[(0,0,0,0),(0,2*t,0,0),(0,0,2,0),(0,m,1,0)]
    TARGET=[SS,0,0,0]
    cnt=0; both=0
    for v in product(range(-R,R+1),repeat=4):
        if v==(0,0,0,0): continue
        re=[0,0,0,0]; im=[0,0,0,0]
        for k,(rr,ii) in zip(v,zip(RE,IM)):
            if k==0: continue
            for q in range(4):
                re[q]+=k*rr[q]; im[q]+=k*ii[q]
        nrm=[a+b for a,b in zip(mul(re,re),mul(im,im))]
        if nrm==TARGET:
            cnt+=1
            a,b,c,e=v
            if not(a==0 and b==0) and not(c==0 and e==0): both+=1
    return cnt, both

print("="*82)
print("VERIFY  #units(L_t) = 12 + r_E(t)   (exact integer arithmetic, box R)")
print("="*82)
print(f"  {'t':>3} {'d=4t-1':>7} {'deg?':>5} {'#units':>7} {'12+rE':>7} {'rE(t)':>6} {'both':>5} {'match':>6}")
allok=True
for t in range(2,32):
    d=4*t-1
    degen = (d%3==0 and is_square(d//3))   # w_t falls into Q(sqrt-3)
    if degen:
        print(f"  {t:>3} {d:>7} {'YES':>5} {'rank2':>7} {'--':>7} {'--':>6} {'--':>5}  (w_t in Q(sqrt-3))")
        continue
    R = 5 if t<=14 else 7
    cnt,both = count_units_exact(t,R)
    pred = 12 + rE(t)
    ok = (cnt==pred) and (both==rE(t))
    allok = allok and ok
    print(f"  {t:>3} {d:>7} {'no':>5} {cnt:>7} {pred:>7} {rE(t):>6} {both:>5} {('OK' if ok else 'FAIL!'):>6}")
print()
print(f"ALL MATCH: {allok}")
print()
print("k in (#units = 6 + 6k):  k = 1 + r_E(t)/6 = 1 + (#zeta6-orbits of norm-t Eisenstein ints).")
print(" r_E(t)=0  (t has a prime =2 mod3 to ODD power; t NOT Loeschian) -> 12 units (k=1)")
print(" r_E(t)=6  (t Loeschian, single orbit: 3,4,9,12,...)             -> 18 units (k=2)")
print(" r_E(t)=12 (t has a SPLIT prime =1 mod3, e.g. 13,21,31)          -> 24 units (k=3)")
print(" => the unit-vector count is a TRIANGULAR-SIDE (Q(sqrt-3)) invariant of t,")
print("    independent of the second CM direction sqrt(4t-1).")
