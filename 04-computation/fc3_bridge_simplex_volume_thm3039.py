"""THE FC(3) BRIDGE, with the normalisation done carefully.

THM-3018: for HOMOGENEOUS f of degree d in n variables,
   L(f^m) = (dm+n-1)! * int_{Delta_{n-1}} g^m dA,   g = f|_Delta,
and the referee script for n=3 uses Delta_2 = {u,v>=0, u+v<=1} with dA = du dv.
So Area(Delta_2) = 1/2 (NOT 1 -- the n=2 case had Area(Delta_1) = length[0,1] = 1,
which is why the FC(2) forced level was exactly 1).

A counterexample g (nonzero, all moments m>=1 vanishing) therefore gives
   int_{Delta_2} e^{g} dA = sum_m (1/m!) int g^m dA = Area + 0 = 1/2.

  ==>  FC(3)-homog  <==  [ int_{Delta_2} e^{Q} dA != 1/2  for every nonconstant Q in Qbar[u,v] ].

MODEL with all moments m>=1 zero: for (u,v) uniform on Delta_2, the marginal density of u
is 2(1-u), so T = 1-(1-u)^2 is UNIFORM on [0,1].  Hence h(u,v) = exp(2 pi i (2u - u^2))
satisfies int_Delta h^m dA = Area * int_0^1 e^{2 pi i m t} dt = 0 for every m >= 1.
"""
import mpmath as mp
mp.mp.dps = 30

def tri(f):
    """int over Delta_2 = {u,v>=0, u+v<=1} of f(u,v) du dv"""
    return mp.quad(lambda u: mp.quad(lambda v: f(u, v), [0, 1-u]), [0, 1])

print("Area(Delta_2) =", mp.nstr(tri(lambda u, v: mp.mpf(1)), 20), " (expect 0.5)")
print()
h = lambda u, v: mp.e**(2j*mp.pi*(2*u - u**2))
print("model h(u,v) = exp(2 pi i (2u - u^2));  moments int_Delta h^m dA:")
for m in range(0, 6):
    val = tri(lambda u, v, m=m: h(u, v)**m)
    print(f"   m={m}: {mp.nstr(val, 14)}")
print()
I = tri(lambda u, v: mp.e**(h(u, v)))
print(f"int_Delta e^h dA = {mp.nstr(I, 20)}    Area = 0.5    diff = {mp.nstr(abs(I - mp.mpf(1)/2), 5)}")
print()
print("=> FORCED LEVEL for FC(3) is Area(Delta_2) = 1/2, NOT 1.")
print("   FC(3)-homog follows from:  int_{Delta_2} e^Q dA != 1/2  (nonconstant Q in Qbar[u,v]).")
print()
print("CONTROL -- THM-3018 sec 4's S_3 eigenvector g = u + w v + w^2 (1-u-v), w = e^{2 pi i/3}:")
w = mp.e**(2j*mp.pi/3)
g = lambda u, v: u + w*v + w**2*(1-u-v)
for m in (1, 2, 3, 4, 6):
    val = tri(lambda u, v, m=m: g(u, v)**m)
    print(f"   m={m}: {mp.nstr(val, 12)}   {'(vanishes: 3 does not divide m)' if m % 3 else '(SURVIVES -- so g is NOT a counterexample)'}")
print("   -> confirms THM-3018 sec 4: the S_3 mechanism kills 3-not-dividing-m only; the")
print("      surviving m=3k moments are nonzero, so this g is a CONTROL, not a witness.")
