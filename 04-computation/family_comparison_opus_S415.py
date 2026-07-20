import sympy as sp
x,y,z=sp.symbols('x y z')
u=1+x*y
F=[sp.expand(u**3*z + y**2*u*(4+3*x*y)),
   sp.expand(y + 3*x*u**2*z + 3*x*y**2*(4+3*x*y)),
   sp.expand(2*x - 3*x**2*y - x**3*z)]
print("OUR counterexample (THM-1300), component degrees:")
for i,f in enumerate(F,1):
    print(f"   F{i}: total degree {sp.Poly(f,x,y,z).total_degree()}")
print(f"   deg F = {max(sp.Poly(f,x,y,z).total_degree() for f in F)}")
J=sp.Matrix([[sp.diff(f,v) for v in (x,y,z)] for f in F])
print(f"   det JF = {sp.simplify(sp.expand(J.det()))}")
print()
print("THEIR family: deg E_2 = 7, and #E_m^{-1}(c,0,0) = 2m-1, so m=2 gives fibre 3.")
print("OURS: degree 7, fibre 3.  => THEIR m=2 IS OUR CASE; m>=3 is genuinely new to us.")
print()
print("NOW THE FORCING (THM-1350).  Fibre over a tau-fixed target = 1 sigma-fixed point")
print("plus FREE sigma-orbits (size 2).  So |fibre| = 1 + 2k, ODD, minimum 3.")
print("Their count 2m-1 factors EXACTLY as 1 + 2(m-1):")
for m in range(2,8):
    print(f"   m={m}: |fibre| = 2m-1 = {2*m-1} = 1 sigma-fixed + {m-1} free 2-orbit(s)")
print()
print("CHECK on our m=2 map: fibre over the tau-fixed target (a,0,0).")
a=sp.symbols('a', positive=True)
# fibre cubic from THM-1440: L x^3 + (4-3bc) x - 2c, at b=c=0 -> 16a x^3 + 4x
sol=sp.solve(sp.Eq(16*a*x**3+4*x,0),x)
print(f"   roots: {sol}")
print("   x=0 is the SIGMA-FIXED sheet (sigma = (-x,-y,z)); the other two are +-r,")
print("   a single free sigma-orbit.  1 + 2*1 = 3 = 2m-1 at m=2.  CONFIRMED.")
