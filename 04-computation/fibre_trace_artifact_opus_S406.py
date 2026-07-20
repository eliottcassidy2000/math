import sympy as sp
x,y,z = sp.symbols('x y z')
u = 1+x*y
F = [sp.expand(u**3*z + y**2*u*(4+3*x*y)),
     sp.expand(y + 3*x*u**2*z + 3*x*y**2*(4+3*x*y)),
     sp.expand(2*x - 3*x**2*y - x**3*z)]
def sum_of_fibre_x(maps, target):
    """eliminate to a univariate poly in x; -c2/L = sum of the fibre's x-coordinates"""
    I = [sp.expand(m - t) for m, t in zip(maps, target)]
    Gb = sp.groebner(I, z, y, x, order='lex')
    uni = [g for g in Gb.exprs if g.free_symbols <= {x}]
    if not uni: return None, None
    P = sp.Poly(uni[0], x)
    co = P.all_coeffs()
    return P.degree(), (sp.Rational(-co[1], co[0]) if len(co) > 1 else 0)
print("DECISIVE NUMERIC TEST -- is sum(fibre x-coords) = 0 an artifact?")
print("Compute the eliminated univariate polynomial in x at NUMERIC targets.\n")
for tgt in [(sp.Integer(1), sp.Integer(0), sp.Integer(0)),
            (sp.Integer(2), sp.Integer(1), sp.Integer(1)),
            (sp.Rational(3,2), sp.Integer(-2), sp.Integer(5))]:
    d, s = sum_of_fibre_x(F, tgt)
    print(f"  ORIGINAL F, target {tuple(map(str,tgt))}: deg={d}, sum of fibre x = {s}")
print()
lam = sp.Integer(1)
phi = {x: x + lam*y*z}
G = [sp.expand(f.subs(phi, simultaneous=True)) for f in F]
print("  composed G = F o phi, phi = (x + y*z, y, z)  [Keller, sigma/tau-equivariant]")
for tgt in [(sp.Integer(1), sp.Integer(0), sp.Integer(0)),
            (sp.Integer(2), sp.Integer(1), sp.Integer(1)),
            (sp.Rational(3,2), sp.Integer(-2), sp.Integer(5))]:
    d, s = sum_of_fibre_x(G, tgt)
    print(f"  COMPOSED G, target {tuple(map(str,tgt))}: deg={d}, sum of fibre x = {s}")
print()
print("READING: targets (1,0,0) are tau-FIXED (b=c=0) -- sum MUST be 0 there by THM-1445-B.")
print("  Targets with b,c != 0 are the real test: if the sum is nonzero for G but zero")
print("  for F, then 'sum = 0 everywhere' is a NORMALISATION ARTIFACT of the owner's F,")
print("  and only tau-ODDNESS (vanishing on the fixed locus) is structural.")
