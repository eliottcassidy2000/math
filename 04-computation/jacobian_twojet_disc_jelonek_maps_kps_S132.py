import sympy as sp

x, y, z, u, v, w, T = sp.symbols('x y z u v w T')

def nonsquare_over_C(expr):
    """rational function: square in C(u,v,w) iff num*den is a square in C[u,v,w];
    parity of QQ(i)-irreducible multiplicities decides (irreducibles stay squarefree
    and pairwise coprime over C). Returns True if certified NON-square."""
    num, den = sp.fraction(sp.cancel(expr))
    _, fl = sp.factor_list(sp.expand(num * den), gaussian=True)
    return any(m % 2 == 1 for f, m in fl)

print('--- MAP III: D4 vs C4 vs V4 (biquadratic z^4 + p z^2 + r, p=-v/u, r=w/u) ---')
p0 = -v / u
r0 = w / u
print('r = w/u nonsquare over C(u,v,w):', nonsquare_over_C(r0))
print('r(p^2-4r) =', sp.cancel(r0 * (p0**2 - 4 * r0)),
      ' nonsquare over C(u,v,w):', nonsquare_over_C(r0 * (p0**2 - 4 * r0)))
print('=> quartic irreducible (linear in w) + r nonsquare + r(p^2-4r) nonsquare => Gal = D4')

print()
print('--- MAP II: resolvent cubic irreducibility over C (linear-in-w gcd certificate) ---')
Rcl = sp.expand(T**3*u**4 + 2*T**2*u**3*v - 4*T*u**2*v**2 + 4*T*u**2*w - T*u*v - 8*u*v**3 + 8*u*v*w - 2*v**2 + w)
B = sp.expand(sp.diff(Rcl, w))            # coeff of w (Rcl is linear in w)
A = sp.expand(Rcl - w * B)
print('Rcl = A + w*B;  B =', B)
g = sp.gcd(A, B)
print('gcd(A, B) =', g, ' => resolvent irreducible over C(u,v,w):', g.is_constant() if hasattr(g,'is_constant') else g)

print()
print('--- MAP II: branch locus certificate ---')
F2 = (x, x*z**2 + y, y*z + y**2)
dJ = sp.Matrix(F2).jacobian([x, y, z]).det()
print('det J =', sp.expand(dJ))
ycrit = sp.solve(sp.Eq(dJ, 0), y)[0]
print('critical surface: y =', ycrit, ' (covers all of {detJ=0}: 1-4xz=0 forces detJ=-z/2 != 0)')
Fu, Fv, Fw = [sp.cancel(c.subs(y, ycrit)) for c in F2]
Delta = 256*u**3*v**2*w**2 - 256*u**3*w**3 + 128*u**2*v**3*w - 96*u**2*v*w**2 + 16*u*v**4 + 24*u*v**2*w - 27*u*w**2 + 4*v**3
val = sp.cancel(Delta.subs({u: Fu, v: Fv, w: Fw}))
print('Delta(F(critical point)) simplifies to:', sp.simplify(val))
print('u(F) on critical surface = x (not identically 0), w(F) =', sp.factor(Fw),
      '(not identically 0) => u, w from elimination are extraneous; branch = V(Delta)')
print('Delta irreducible over QQ:', len(sp.factor_list(Delta)[1]) == 1)
print('critical set is the irreducible rational graph y=2xz^2/(1-4xz) => image closure irreducible')

print()
print('--- MAP I: branch locus over C ---')
F1 = (x**2 - z**2 + y, x*z, y)
for s in (sp.I, -sp.I):
    Fi = [sp.expand(c.subs(x, s * z)) for c in F1]  # x^2+z^2=(x-iz)(x+iz)
    # image: u = y-2z^2? compute:
    print(' x =', s, '* z  ->  F =', Fi, '  eliminates to u - w', '-' if s == sp.I else '+', '2i v = 0 check:',
          sp.simplify((Fi[0] - Fi[2]) + (2 * s * sp.I) * Fi[1]) )
print('=> branch locus = {(u-w)^2 + 4v^2 = 0}  (the pair of conjugate planes u-w=±2iv)')
print('fiber over (5,0,1) (on {v=0}, off branch):')
sols = sp.solve([F1[0]-5, F1[1]-0, F1[2]-1], [x, y, z], dict=True)
print(sols)
dJ1 = sp.Matrix(F1).jacobian([x,y,z]).det()
print('detJ at those points:', [sp.simplify(dJ1.subs(s)) for s in sols])
print('=> 4 distinct unramified points though disc_z(5,0,1)=0: v^2 factor = z-primitive-element defect, NOT branch')

print()
print('--- swallowtail-normalized master identities (L = u = Jelonek polynomial) ---')
# Map III
q3m_disc = sp.cancel(16*u*w*(4*u*w - v**2)**2 / u**6)
print('Map III: u^5 * disc(monic quartic) =', sp.factor(sp.cancel(u**5 * q3m_disc)))
# Map II
q2_disc = sp.expand(u**3 * Delta)
q2m_disc = sp.cancel(q2_disc / (u**2)**6)
print('Map II : u^9 * disc(monic quartic) =', sp.expand(sp.cancel(u**9 * q2m_disc)) == Delta, '(equals Delta)')
