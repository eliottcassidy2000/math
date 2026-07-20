# Independent verification of kp-S128c120 THM-1605: on {-1,0,1}, M=1, constant a,b,c,
#   m * E_r[psi_m] = s^m * He_m(b/s),  s = sqrt(-2ac),
# with psi_m = [t^m] log v,  v = 1 + b t v + (r a c) t^2 v^2,  E_r[r^k] = k!.
# He_m = probabilists' Hermite (He_0=1, He_1=x, He_{m+1}=x He_m - m He_{m-1}).
from fractions import Fraction as Fr
import cmath
def psi_Er(a,b,c,M):
    # v as power series in t up to t^M, coeffs are polynomials in r (list by r-power), rationals
    # represent a series coeff as dict {r_power: Fraction}
    def padd(p,q):
        r=dict(p)
        for k,v in q.items(): r[k]=r.get(k,Fr(0))+v
        return {k:v for k,v in r.items() if v!=0}
    def pscale(p,s):
        return {k:v*s for k,v in p.items()} if s!=0 else {}
    W_r={1:Fr(a*c)}  # r a c  (a,c ints here) => coeff of r^1
    B={0:Fr(b)}
    # v = sum_{n>=0} v_n t^n, v_0=1; recursion from v = 1 + b t v + (rac) t^2 v^2
    v=[{0:Fr(1)}]  # v_0
    for n in range(1,M+1):
        # [t^n] of (b t v): b*v_{n-1}; [t^n] of (rac t^2 v^2): rac * sum_{i+j=n-2} v_i v_j
        term=pscale(v[n-1],Fr(b))
        conv={}
        for i in range(0,n-1):
            j=n-2-i
            if 0<=j<len(v):
                # v_i * v_j (poly mult in r), then * rac (r^1)
                for k1,c1 in v[i].items():
                    for k2,c2 in v[j].items():
                        kk=k1+k2+1; conv[kk]=conv.get(kk,Fr(0))+c1*c2*Fr(a*c)
        vn=padd(term,conv)
        v.append(vn)
    # log v = sum_{m>=1} psi_m t^m ; compute via L = log(1+ (v-1)) series, or L' v = v'
    # Use: psi_m from log series: log v = sum, with v_0=1. Newton: L_m = v_m - (1/m) sum_{k=1}^{m-1} k L_k v_{m-k}
    L=[{}]  # L_0 unused
    for m in range(1,M+1):
        Lm=dict(v[m])
        for k in range(1,m):
            # subtract (k/m) L_k v_{m-k}
            for k1,c1 in L[k].items():
                for k2,c2 in v[m-k].items():
                    kk=k1+k2; Lm[kk]=Lm.get(kk,Fr(0))-Fr(k,m)*c1*c2
        Lm={k:val for k,val in Lm.items() if val!=0}
        L.append(Lm)
    # E_r[psi_m] = sum_k [r^k] psi_m * k!
    from math import factorial
    Er=[None]
    for m in range(1,M+1):
        Er.append(sum(val*factorial(k) for k,val in L[m].items()))
    return Er
def He(m,x):
    if m==0: return 1
    h0,h1=1,x
    for k in range(1,m): h0,h1=h1,x*h1-k*h0
    return h1
for (a,b,c,tag) in [(1,1,1,"klein row (ac>0, s imaginary)"),(1,1,-1,"SIGN-MIXED a=1,b=1,c=-1"),
                    (1,3,1,"b=3"),(1,2,-1,"sign-mixed b=2")]:
    Er=psi_Er(a,b,c,8)
    s=cmath.sqrt(-2*a*c)
    print(f"\n[a={a},b={b},c={c}] {tag}: s=sqrt(-2ac)={s:.4f}")
    ok=True; anyzero=False
    for m in range(1,9):
        lhs=complex(m*Er[m])
        rhs=(s**m)*He(m,b/s)
        match=abs(lhs-rhs)<1e-9*(1+abs(rhs))
        ok=ok and match
        if abs(complex(Er[m]))<1e-12: anyzero=True
    print(f"  identity m*E_r[psi_m]==s^m He_m(b/s) for m=1..8: {ok}")
    print(f"  any E_r[psi_m]==0 (would allow nullcone)? {anyzero}  => one-sided forced: {not anyzero}")
