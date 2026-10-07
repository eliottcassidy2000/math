"""C2: THM-1300's map at the prime 2 (det JF = -2): the special fibre and 2-adic merges.

(a) Symbolic: F mod 2, the surface y = xz collapsed onto the a-axis.
(b) Image of F mod 2^k on (Z/2^k)^3: #image / 8^k -> Haar measure of F(Z_2^3).
    Area formula: int_{Z_2^3} |det JF| = 1/2 = int N(w) dw, N(w) = #(F^{-1}(w) cap Z_2^3).
    So F is injective a.e. on Z_2^3 iff mu(F(Z_2^3)) = 1/2.
(c) Dominance of F mod 2: exhaustive image size over F_{2^k}^3 for small k.
"""
import sys
import numpy as np
import sympy as sp

x, y, z = sp.symbols('x y z')
u = 1 + x*y
F = [sp.expand(u**3*z + y**2*u*(4 + 3*x*y)),
     sp.expand(y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y)),
     sp.expand(2*x - 3*x**2*y - x**3*z)]

print("(a) F mod 2:")
Fb = [sp.Poly(f, x, y, z, modulus=2) for f in F]
for i, f in enumerate(Fb):
    print("   F%d mod 2 =" % (i+1), f.as_expr())
sub = {y: x*z}
red = [sp.Poly(sp.expand(f.subs(sub)), x, y, z, modulus=2).as_expr() for f in F]
print("   on the surface y = xz:", red)

# (b) image of F mod 2^k
def F_mod(X, Y, Z, m):
    U = (1 + X*Y) % m
    f1 = (U*U % m * U % m * Z + Y*Y % m * U % m * ((4 + 3*X*Y) % m)) % m
    f2 = (Y + 3*X % m * (U*U % m) % m * Z + 3*X % m * (Y*Y % m) % m * ((4 + 3*X*Y) % m)) % m
    f3 = (2*X - 3*X*X % m * Y - X*X % m * X % m * Z) % m
    return f1, f2, f3

print("\n(b) #F((Z/2^k)^3) / 8^k  (-> mu(F(Z_2^3)); injective a.e. iff limit = 1/2)")
for k in range(1, 8):
    m = 2**k
    r = np.arange(m, dtype=np.int64)
    X, Y, Z = np.meshgrid(r, r, r, indexing='ij')
    X, Y, Z = X.ravel(), Y.ravel(), Z.ravel()
    f1, f2, f3 = F_mod(X, Y, Z, m)
    code = (f1 * m + f2) * m + f3
    uniq, counts = np.unique(code, return_counts=True)
    hist = np.bincount(counts)
    print("   k=%d  image=%d  ratio=%.6f  fibre-size histogram (size:count) %s"
          % (k, len(uniq), len(uniq)/m**3,
             {s: int(c) for s, c in enumerate(hist) if c and s <= 64}))
    sys.stdout.flush()
