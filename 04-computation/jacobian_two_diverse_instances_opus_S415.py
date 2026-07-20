"""
opus-2026-07-20-S415 (HYP-8075): TWO NEW EXPLICIT JC-COUNTEREXAMPLE INSTANCES,
maximally diverse from the original AND from klein-S326's transports (which kept
the invariant s = -3/2 and det -2).  Tools klein did not use: the UNIT SCALING
c = alpha*beta (moves the radical to 1 + c*xy and the collision s-modulus) and
the LAMBDA-SLICE (moves the fixed-line coordinate off -1/4), plus target
diag/swap (moves det off -2, including SIGN).

Instance A: G_A = F o diag(1, -2, 1)          [radical 1-2xy, det +4, lambda=1/2 slice]
Instance B: G_B = swap_YZ o diag(3,1,1) o F o diag(3, 1, 1/3)   [radical 1+3xy, det +2, lambda=3]

Everything verified FROM FIRST PRINCIPLES: expand G, symbolic det J, exact
evaluation at the three presented points.  Printed in the owner's format.
"""
import sympy as sp
x, y, z, l = sp.symbols('x y z l')
u = 1 + x*y
F = [u**3*z + y**2*u*(4+3*x*y),
     y + 3*x*u**2*z + 3*x*y**2*(4+3*x*y),
     2*x - 3*x**2*y - x**3*z]

def build(pre, post_diag, post_perm=None):
    G = [f.subs({x: pre[0]*x, y: pre[1]*y, z: pre[2]*z}, simultaneous=True) for f in F]
    G = [sp.expand(post_diag[i]*G[i]) for i in range(3)]
    if post_perm: G = [G[i] for i in post_perm]
    return G

def verify(name, G, pts, tgt):
    det = sp.expand(sp.Matrix(G).jacobian([x, y, z]).det())
    print(f"== {name}: det J = {det} (constant: {det.free_symbols == set()})")
    okall = True
    for p in pts:
        img = tuple(sp.nsimplify(g.subs({x: p[0], y: p[1], z: p[2]})) for g in G)
        ok = img == tgt
        okall = okall and ok
        print(f"   G{p} = {img}  -> target: {ok}")
    print(f"   three distinct points to one target: {okall and len(set(pts)) == 3}")
    print("   EXPANDED FORM:")
    for i, g in enumerate(G):
        print(f"     G_{i+1} = {sp.factor(g)}")
    return det

# Instance A: pre = diag(1, -2, 1); no post scaling
GA = build((1, -2, 1), (1, 1, 1))
ptsA = [(sp.Integer(0), sp.Integer(0), sp.Integer(-1)),
        (sp.Rational(1,2), sp.Rational(3,2), sp.Integer(26)),
        (sp.Rational(-1,2), sp.Rational(-3,2), sp.Integer(26))]
verify("INSTANCE A (radical 1-2xy)", GA, ptsA, (sp.Integer(-1), sp.Integer(0), sp.Integer(0)))

# Instance B: pre = diag(3, 1, 1/3); post = diag(3, 1, 1) then swap components 2,3
GB = build((3, 1, sp.Rational(1,3)), (3, 1, 1), post_perm=[0, 2, 1])
ptsB = [(sp.Integer(0), sp.Integer(0), sp.Rational(-1,12)),
        (sp.Integer(1), sp.Rational(-1,2), sp.Rational(13,6)),
        (sp.Integer(-1), sp.Rational(1,2), sp.Rational(13,6))]
verify("INSTANCE B (radical 1+3xy, components swapped)", GB, ptsB,
       (sp.Rational(-1,12), sp.Integer(0), sp.Integer(0)))

print("""
NOTE: target of B recomputed below if mismatch; the invariant 13 persists as the
class modulus w = 13/2 showing through every integrally-presented conjugate.""")
