# A nonzero LND of A = P/(H): in the Laurent chart A[1/p] = C[p^{+-1},x,y,z] take D = p^5 d/dy.
from sympy import symbols, expand, div, Poly, reduced
p, s, u, F, J = symbols('p s u F J')
x = s**2 + u**3 + p**2*F
H = expand(x**2*F - (1 + 2*s*x)*J - p**2*J**2 - p*u)
Du = p**2*x
Ds = p**5 + 3*p**2*x**2*u**2
DF = -2*s*p**3 - 6*s*x**2*u**2 - 3*x*u**2
DJ = -x*p**3 - 3*x**3*u**2
def D(e):
    return expand(e.diff(s)*Ds + e.diff(u)*Du + e.diff(F)*DF + e.diff(J)*DJ)
print("D(x) = 0:", D(x) == 0)
DH = D(H)
q, r = reduced(DH, [H], s, u, F, J, p)
print("D(H) mod H == 0:", r == 0, "; quotient:", q[0])
# local nilpotence on A: iterate on generators modulo H (it is LND on A[1/p] = C[p^{+-1},x,y,z] as p^5 d/dy)
for name, g in [('s', s), ('u', u), ('F', F), ('J', J)]:
    e, k = g, 0
    while True:
        e = D(e); k += 1
        rr = reduced(e, [H], s, u, F, J, p)[1] if e != 0 else 0
        if rr == 0 or k > 12:
            break
    print(f"D^{k}({name}) == 0 mod H:", rr == 0)
