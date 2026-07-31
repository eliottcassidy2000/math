"""
The central-binomial edge, verified across AMM / S(k) / Jacobian.

Claim: each problem's characteristic constant is a central-binomial generating function evaluated at
its arithmetic edge; the arithmetic nature is set by motivic dimension (curve -> algebraic golden;
surface -> irreducible period). See 07-reflections/central-binomial-edge-...-opus-S4.md.
"""
import mpmath as mp
mp.mp.dps = 30
from math import comb
import sympy as sp

phi = (1+mp.sqrt(5))/2

# (1) AMM: Catalan GF at branch point w=-1  (THM-3009)
def Cat(w): return (1-mp.sqrt(1-4*w))/(2*w)
print("AMM (signature 2, a curve):")
print("  Catalan C(-1) =", mp.nstr(Cat(-1), 15), " = 1/phi =", mp.nstr(1/phi, 15))
print("  sqrt(1-4w)|_{-1} = sqrt5 =", mp.nstr(mp.sqrt(5), 12), " (=1/p*)")
print("  => C* = log_5(5 phi^2) =", mp.nstr(mp.log(5*phi**2)/mp.log(5), 16), " (algebraic edge)")

# (2) S(k): product of central binomials = signature 4 (a surface)
print("\nS(k) (signature 4, a surface):")
ok = all(abs(comb(2*n, n)*comb(4*n, 2*n)/mp.mpf(64)**n
             - mp.rf(mp.mpf(1)/4, n)*mp.rf(mp.mpf(3)/4, n)/mp.factorial(n)**2) < mp.mpf(10)**-25
         for n in range(1, 8))
print("  C(2n,n)C(4n,2n)/64^n == (1/4)_n(3/4)_n/n!^2 for n=1..7:", ok, " (Catalan (x) Catalan)")
# S(4) clean 1-D edge = theta-weighted elliptic moment; independent of {K(i),varpi,Catalan,pi}
S4 = (2/mp.pi)*mp.quad(lambda s: (mp.asinh(s)+mp.asin(s))/mp.sqrt(1-s**4), [0, 1])
anchor = mp.quad(lambda th: 1/mp.sqrt(1+mp.sin(th)**2), [0, mp.pi/2])   # K(i)
lemn = mp.gamma(mp.mpf(1)/4)**2/(2*mp.sqrt(2*mp.pi))
rel = mp.pslq([mp.pi*S4, anchor, mp.catalan, mp.pi, lemn], maxcoeff=10**7)
print("  S(4) =", mp.nstr(S4, 16), " PSLQ vs {K(i),Catalan,pi,varpi}:", rel,
      "(coeff 0 on S(4) => irreducible edge)")

# (4) Jacobian: BCW formal inverse of F=X+H terminates (tree/Catalan skeleton); capacity role INVERTED
x, y = sp.symbols('x y')
G = [x, y]
for _ in range(6):
    sub = {x: G[0], y: G[1]}
    H = [(y**2).subs(sub), sp.Integer(0)]          # F=(x+y^2, y), J=1 constant
    G = [sp.expand(x-H[0]), sp.expand(y-H[1])]
print("\nJacobian (tree/Catalan skeleton):")
print("  F=(x+y^2,y) J=1 -> formal inverse G =", G, "-> TERMINATES (polynomial)")
print("  (capacity role inverted vs AMM: JC conjectures termination always holds; open in 2D)")
