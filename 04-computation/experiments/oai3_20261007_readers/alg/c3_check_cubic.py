import sympy as sp
x, y, z, a, b, c = sp.symbols('x y z a b c')
u = 1 + x*y
F1 = sp.expand(u**3*z + y**2*u*(4 + 3*x*y))
F2 = sp.expand(y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y))
F3 = sp.expand(2*x - 3*x**2*y - x**3*z)
phi = -2*y**3 + 3*b*y**2 - 18*a*y + (18*a*b - b**3 - 27*a**2*c)
s = {a: F1, b: F2, c: F3}
print("phi(y) after substituting F:", sp.expand(phi.subs(s)))
print("x*(y^2 - b y + 3a) - (b - y):", sp.expand((x*(y**2 - b*y + 3*a) - (b - y)).subs(s)))
print("x^3 z - (2x - 3x^2 y - c):", sp.expand((x**3*z - (2*x - 3*x**2*y - c)).subs(s)))
