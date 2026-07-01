"""The NEXT STEP in the iota-odd certificate. (1) fundamental unit of Q(sqrt21) + its NORM (genus theory of
21=3*7); (2) the NARROW class group Z/2 = the residual odd degree; (3) Hasse-Davenport lift; (4) the compositum
where sqrt21 lives (both the tight-AP 7 and the covering 3)."""
import numpy as np, math
def pell_fund(d):  # fundamental unit of Q(sqrt d), handling (a+b sqrt d)/2 for d=1 mod4
    a0=math.isqrt(d); m,q,a=0,1,a0; h1,h=1,a0; k1,k=0,1
    for _ in range(1000):
        # test half-integer units (a+b sqrt d)/2 with a=b mod2 and (a^2-d b^2)=+-4
        for (A,B) in [(h,k)]:
            for nrm in (1,-1,4,-4):
                pass
        m=q*a-m; q=(d-m*m)//q; a=(a0+m)//q
        h1,h=h,a*h+h1; k1,k=k,a*k+k1
        for scale in [1]:
            val=h*h-d*k*k
            if val in (1,-1): return h,k,val,1
            if d%4==1 and (h*h-d*k*k) in (4,-4): return h,k,(h*h-d*k*k)//4,2
    return None
for d in [21,3,7,61]:
    r=pell_fund(d)
    if r:
        x,y,nrm,den=r
        eps = (x + y*math.sqrt(d))/den
        print(f"  Q(sqrt{d}): fundamental unit ~ ({x}+{y}sqrt{d})/{den} = {eps:.4f}, NORM = {nrm}  ({'+1 => NARROW class# = 2*wide (residual Z/2!)' if nrm==1 else '-1 => narrow=wide (no residual)'})")
# genus theory: narrow 2-rank = (#ramified primes in disc) - 1
def genus_2rank(d):
    # disc of Q(sqrt d)
    D = d if d%4==1 else 4*d
    fac=[]; n=abs(D); p=2
    while p*p<=n:
        if n%p==0:
            fac.append(p)
            while n%p==0: n//=p
        p+=1
    if n>1: fac.append(n)
    return len(fac)-1, D, fac
for d in [21,3,7]:
    r,D,fac=genus_2rank(d)
    print(f"  Q(sqrt{d}): disc={D}, ramified primes {fac}, narrow class group 2-RANK = {r}  ({'Z/2 residual' if r>=1 else 'trivial'})")
# Hasse-Davenport lift: Gauss over F_{p^2} = -(Gauss over F_p)^2 = -(i sqrt p)^2 = +p
print(f"\n  Hasse-Davenport: g_2(chi_p)=-g(chi_p)^2 = -(i sqrt p)^2 = +p (the iota-odd i sqrt p lifts to the EVEN prime p).")
# compositum: sqrt21 needs both 7 (tight AP N=14) and 3 (covering Phi6=183)
print(f"  COMPOSITUM: sqrt-7 in Q(zeta_14) [tight AP], sqrt-3 in Q(zeta_183) [covering]; sqrt21 in Q(zeta_lcm(14,183))=Q(zeta_{np.lcm(14,183)}).")
print(f"  => sqrt21 surfaces only when the tight-AP heptagon (7) and the covering Eisenstein (3) are BOTH present.")
