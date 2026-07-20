import sympy as sp
x,y,z=sp.symbols('x y z')
def Fmap(x,y,z):
    u=1+x*y
    return (u**3*z + y**2*u*(4+3*x*y),
            y + 3*x*u**2*z + 3*x*y**2*(4+3*x*y),
            2*x - 3*x**2*y - x**3*z)
F=sp.Matrix(Fmap(x,y,z)); J=F.jacobian([x,y,z])
print("THE GAP I LEFT IN THM-1350: det JF const bounds only the PRODUCT of blocks.")
print("sigma acts as (-1,-1,+1), so on Fix(sigma)={x=y=0} the tangent space splits:")
print("  +1 eigenspace = z-direction (tangent to Fix);  -1 eigenspace = (x,y)-plane.")
print()
Jf=J.subs({x:0,y:0})
print("J at (0,0,z):"); sp.pprint(sp.simplify(Jf))
print()
# tangential block: dF1/dz restricted (F1 is the Fix(tau) coordinate)
tang=sp.simplify(sp.diff(Fmap(0,0,z)[0], z))
# normal block: 2x2 in (x,y) of the anti-invariant coords (F2,F3)
norm=sp.Matrix([[sp.diff(F[1],x),sp.diff(F[1],y)],[sp.diff(F[2],x),sp.diff(F[2],y)]]).subs({x:0,y:0})
nd=sp.simplify(norm.det())
print(f"  tangential determinant (d/dz of F|Fix) = {tang}")
print(f"  normal 2x2 block det at (0,0,z)        = {nd}")
print(f"  product = {sp.simplify(tang*nd)}   det JF = {sp.expand(J.det())}")
print()
print("THE REPAIR (the argument I should have written):")
print("  On Fix, det JF = det(tangential) * det(normal), both POLYNOMIALS in z.")
print("  Their product is the nonzero constant -2.  In the domain C[z], a product")
print("  of polynomials equal to a nonzero constant forces EACH factor to be a")
print("  nonzero constant (degrees add).  Hence the tangential block is a nonzero")
print("  constant, i.e. F|Fix IS Keller in dimension dim Fix(sigma).  The claim")
print("  stands; the proof needed this step.")
print()
print(f"  check both are constants here: tangential {tang} (deg 0? {sp.degree(sp.Poly(tang,z)) if tang.free_symbols else 0}),"
      f" normal {nd} (deg 0? {sp.degree(sp.Poly(nd,z)) if nd.free_symbols else 0})")
