"""The '2' in p = delta/(2-delta) comes from |Delta^j P(0)| <= 2^j, i.e. from
the binary alphabet.  Replace it by q and re-solve the stationarity system:
    p = delta/(q-delta);  gamma = -H'(delta)/[H(p) + (1+H'(delta))(1-p)];
    alpha = H(delta)/[H(p)+1-p]  must equal (q-delta)/[(1+gamma)/gamma - p].
Does the root produce metallic ratios  x^2 = n x + 1 ?"""
from mpmath import mp, mpf, log, sqrt, findroot
mp.dps=30; L2=log(2)
H=lambda p: mpf(0) if p<=0 or p>=1 else (-p*log(p)-(1-p)*log(1-p))/L2
Hp=lambda d: log((1-d)/d)/L2
def resid(delta,q):
    d=mpf(delta); p=d/(q-d); hp=Hp(d)
    g=-hp/(H(p)+(1+hp)*(1-p))
    return H(d)/(H(p)+1-p) - (q-d)/((1+g)/g - p)
def solve(q, guess):
    d=findroot(lambda x: resid(x,mpf(q)), mpf(guess))
    p=d/(q-d); hp=Hp(d); g=-hp/(H(p)+(1+hp)*(1-p))
    return d,p,g
print(" q    delta*            p*                gamma*            1+gamma*")
res={}
for q,guess in ((2,"0.618"),(3,"0.5"),(4,"0.45"),(5,"0.4")):
    try:
        d,p,g=solve(q,guess); res[q]=(d,p,g)
        print(f" {q}  {mp.nstr(d,15):18s} {mp.nstr(p,15):18s} {mp.nstr(g,15):18s} {mp.nstr(1+g,15)}")
    except Exception as e:
        print(f" {q}  no root ({type(e).__name__})")
print()
print("metallic ratios  x = (n + sqrt(n^2+4))/2  and 1/x:")
for n in (1,2,3,4):
    x=(n+sqrt(n*n+4))/2
    print(f"  n={n}: x={mp.nstr(x,15)}  1/x={mp.nstr(1/x,15)}  sqrt(n^2+4)={mp.nstr(sqrt(n*n+4),15)}")
print()
for q,(d,p,g) in res.items():
    print(f"q={q}: delta*={mp.nstr(d,12)}   1/p*={mp.nstr(1/p,12)}   check delta*=1/phi_n?")
    for n in (1,2,3,4):
        x=(n+sqrt(n*n+4))/2
        if abs(d-1/x)<mpf("1e-8"): print(f"     MATCH delta* = 1/metallic({n})")
        if abs(1/p-sqrt(n*n+4))<mpf("1e-8"): print(f"     MATCH 1/p* = sqrt({n}^2+4)")
