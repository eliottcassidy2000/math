"""The closure equation reduces to delta^2 = 1-delta.  Check the algebra."""
from mpmath import mp, mpf, log, sqrt
mp.dps=40; L2=log(2)
H  = lambda d: (-d*log(d)-(1-d)*log(1-d))/L2
Hp = lambda d: log((1-d)/d)/L2

# Claim: with p = delta/(2-delta) and gamma = -H'(d)/D,  D = H(p)+(1+H')(1-p),
# the closure alpha1 = alpha2 is EQUIVALENT to  H(d) = -H'(d)*(2-d).
def alpha1(d):
    p=d/(2-d); return H(d)/(H(p)+1-p)
def alpha2(d):
    p=d/(2-d); hp=Hp(d); D=H(p)+(1+hp)*(1-p); g=-hp/D
    return (2-d)/((1+g)/g - p)
def reduced(d):
    return H(d) + Hp(d)*(2-d)          # = 0 iff closure holds
print("d      alpha1-alpha2        H(d)+H'(d)(2-d)")
for d in ("0.55","0.618033988749895","0.70"):
    d=mpf(d); print(f"{d}  {mp.nstr(alpha1(d)-alpha2(d),8):>18s}  {mp.nstr(reduced(d),8)}")
print()
# and H(d) = -H'(d)(2-d)  <=>  2 log2 d = log2(1-d)  <=>  d^2 = 1-d
phi=(1+sqrt(5))/2; d=1/phi
print("at d = 1/phi :  d^2 - (1-d) =", mp.nstr(d*d-(1-d),8))
print("               H(d)+H'(d)(2-d) =", mp.nstr(reduced(d),8))
print("               p = d/(2-d) - 1/sqrt5 =", mp.nstr(d/(2-d)-1/sqrt(5),8))
print("               -H'(d) - log2(phi) =", mp.nstr(-Hp(d)-log(phi)/L2,8))
# denominator collapse D = (1/2)log2 5
p=d/(2-d); D=H(p)+(1+Hp(d))*(1-p)
print("               D - (1/2)log2 5 =", mp.nstr(D-log(5)/(2*L2),8))
g=-Hp(d)/D
print("               gamma - log_5(phi^2) =", mp.nstr(g-2*log(phi)/log(5),8))
print()
# general alphabet q:  H(d) = -H'(d)(q-d)  <=>  d^q = (1-d)^(q-1)
print("general q:  closure  <=>  d^q = (1-d)^(q-1)")
from mpmath import findroot
for q in (2,3,4,5):
    r=findroot(lambda x: q*log(x)-(q-1)*log(1-x), mpf("0.6"))
    print(f"  q={q}: delta* = {mp.nstr(r,15)}   (earlier numeric: "
          f"{'0.618034 / 1/phi' if q==2 else '0.5701?/0.549700/0.538597'[ (q-3)*9 : (q-3)*9+8 ] if q>2 else ''})")
