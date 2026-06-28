"""
S75d: the three A..G recursion modes are CYCLOTOMIC-FACTORED characteristic polynomials (x-1)^depth * Phi_d.
EVEN A+B-C=(x-1)^2; FULL A+B+C-D-E-F+G=(x-1)^3; ODD A+B+D-C-E-F+G=(x-1)^3*Phi_2.
Extension: the full family = the cyclotomic factorization of x^n-1; the apex = real Phi_7 (de Moivre cubic).
moment-order depth = the (x-1)-multiplicity; character = the Phi_d factor.
"""
import sympy as sp
x=sp.symbols('x')

modes={
 'EVEN  A+B-C         ': [1,-2,1],         # f(n)=2f(n-1)-f(n-2)
 'FULL  A+B+C-D-E-F+G ': [1,-3,3,-1],      # f(n)=3f(n-1)-3f(n-2)+f(n-3)
 'ODD   A+B+D-C-E-F+G ': [1,-2,0,2,-1],    # f(n)=2f(n-1)-2f(n-3)+f(n-4)
}
print("THE THREE MODES as characteristic polynomials (cyclotomic factorization):")
for name,coef in modes.items():
    p=sum(c*x**(len(coef)-1-i) for i,c in enumerate(coef))
    print(f"  {name}: {sp.factor(p)}")

print("\nCYCLOTOMIC factorization of x^n-1 (the FULL family of relevant recursions):")
for n in (7,14):
    print(f"  x^{n}-1 = {sp.factor(x**n-1)}")
print(f"  real Phi_7 = de Moivre cubic = {sp.minimal_polynomial(2*sp.cos(2*sp.pi/7),x)}")
print(f"  Phi_1={sp.cyclotomic_poly(1,x)}  Phi_2={sp.cyclotomic_poly(2,x)}  Phi_3={sp.cyclotomic_poly(3,x)}")

print("\nSYNTHESIS: mode = (x-1)^depth * Phi_d.")
print("  (x-1)-multiplicity = MOMENT-ORDER DEPTH (kps (p+1)/2: n=14->4); Phi_d = the CHARACTER (chi_2,chi_3,chi_7).")
print("  cap char poly = (x-1)^depth * Phi_2 * Phi_7 * Phi_14; the apex Phi_7 (de Moivre) is the irreducible hardness.")
print("  EVEN/ODD/FULL = first cyclotomic factors; the apex chi_7 (real Phi_7) is the extension that carries the difficulty.")
