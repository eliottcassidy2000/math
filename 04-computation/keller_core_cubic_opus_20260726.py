# The generic x-eliminant "core cubic" of the sporadic Keller map F, symbolically.
# Goal: its leading coefficient L(a,b,c) = the x-escape (non-properness) locus,
# and its discriminant D(a,b,c) = the (would-be) collision locus, exactly.
# Settles whether the non-properness set S_F has a codim-1 component or is only
# the three exact drop curves found in keller_plane_anatomy_opus_20260726.py.
# opus 2026-07-26.

import sympy as sp

x, y, z, a, b, c = sp.symbols('x y z a b c')
u = 1 + x*y
F1 = u**3*z + y**2*u*(4 + 3*x*y)
F2 = y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
A = [u**3, 3*x*u**2, -x**3]
Bv = [sp.expand(F1 - A[0]*z), sp.expand(F2 - A[1]*z), sp.expand(F3 - A[2]*z)]

G1 = sp.expand((a - Bv[0])*A[1] - (b - Bv[1])*A[0])
H1 = sp.expand(sp.cancel(G1 / u**2))
H2 = sp.expand((a - Bv[0])*A[2] - (c - Bv[2])*A[0])

print("deg_y H1 =", sp.Poly(H1, y).degree(), " deg_y H2 =", sp.Poly(H2, y).degree())
print("computing Res_y(H1,H2) over QQ[a,b,c,x] ...")
R = sp.resultant(sp.Poly(H1, y, domain=sp.QQ[a, b, c, x]),
                 sp.Poly(H2, y, domain=sp.QQ[a, b, c, x]))
R = sp.Poly(R.as_expr(), x)
print("total eliminant degree in x:", R.degree())
fl = sp.factor_list(R.as_expr())
print("factor structure:")
core = None
for f, m in fl[1]:
    fp = sp.Poly(f, x)
    dep = f.has(a) or f.has(b) or f.has(c)
    print(f"  deg_x={fp.degree()} mult={m} has-params={dep}"
          + (f"  factor={f}" if fp.degree() <= 1 or not dep else ""))
    if dep and fp.degree() >= 2:
        core = fp
if core is not None:
    print("\nCORE eliminant degree:", core.degree())
    LC = sp.factor(core.LC())
    print("leading coeff L(a,b,c) =", LC)
    print("L on {a=0}:", sp.factor(LC.subs(a, 0)))
    print("L on {c=0}:", sp.factor(LC.subs(c, 0)))
    D = sp.factor(sp.discriminant(core.as_expr(), x))
    print("disc D(a,b,c) =", D)
    # x^0 coefficient (product of roots up to sign/L): where a root is 0
    print("constant coeff =", sp.factor(core.as_expr().subs(x, 0)))
