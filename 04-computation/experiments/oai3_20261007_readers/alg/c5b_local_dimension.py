import itertools, sympy as sp
x, y, z = sp.symbols('x y z')
u = 1 + x*y
F = [sp.expand(u**3*z + y**2*u*(4 + 3*x*y)), sp.expand(y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y)), sp.expand(2*x - 3*x**2*y - x**3*z)]
# (1) mod-2 inseparability identity: F2*y^2 + F2^3 + F1^2*F3 == 0 in F_2[x,y,z]  (reduction of the fibre cubic)
idt = sp.Poly(sp.expand(F[1]*y**2 + F[1]**3 + F[0]**2*F[2]), x, y, z, modulus=2)
print("F2*y^2 + F2^3 + F1^2*F3 mod 2 is zero:", idt.is_zero)
# reconstruction mod 2: x*(y^2 + F2*y + F1) = F2 + y ; x^3*z = x^2*y + F3
print("x(y^2+F2 y+F1) - (F2+y) mod 2 zero:", sp.Poly(sp.expand(x*(y**2 + F[1]*y + F[0]) - (F[1] + y)), x, y, z, modulus=2).is_zero)
print("x^3 z - (x^2 y + F3) mod 2 zero:", sp.Poly(sp.expand(x**3*z - (x**2*y + F[2])), x, y, z, modulus=2).is_zero)
# (2) local dimension of the mod-2 fibre at each point of F_2^3: growth of dim F_2[x,y,z]/(I + m^K)
def std_count(G, B):
    lms = [sp.Poly(g, x, y, z, modulus=2).monoms(order='grevlex')[0] for g in G.exprs]
    return sum(1 for e in itertools.product(range(B), repeat=3) if not any(all(e[i] >= l[i] for i in range(3)) for l in lms))
def Fint(P):
    X, Y, Z = P; U = 1 + X*Y
    return ((U**3*Z + Y**2*U*(4 + 3*X*Y)) % 2, (Y + 3*X*U**2*Z + 3*X*Y**2*(4 + 3*X*Y)) % 2, (2*X - 3*X**2*Y - X**3*Z) % 2)
for Pb in itertools.product(range(2), repeat=3):
    wb = Fint(Pb)
    sub = {x: x + Pb[0], y: y + Pb[1], z: z + Pb[2]}
    I = [sp.expand(f.subs(sub) - c) for f, c in zip(F, wb)]
    vals = []
    for K in (4, 8, 12, 16):
        mK = [x**a*y**b*z**c for a in range(K+1) for b in range(K+1) for c in range(K+1) if a+b+c == K]
        G = sp.groebner(I + mK, x, y, z, order='grevlex', modulus=2)
        vals.append(std_count(G, K + 2))
    print("Pb =", Pb, " dim F_2[x,y,z]/(I + m^K), K = 4,8,12,16:", vals,
          "-> isolated (length %d)" % vals[-1] if vals[-1] == vals[-2] else "-> grows: Pb lies on a curve of the mod-2 fibre")
