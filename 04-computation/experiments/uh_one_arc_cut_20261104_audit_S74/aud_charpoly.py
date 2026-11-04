# Characteristic polynomials of the order-5 strong tournament classes (out-rows from aud_routes classes 8),
# and checks of the owner's claimed polynomials / spectra and the PG(2,3) collinearity.
import sympy as sp
lam = sp.symbols('lam')
classes = {
    'Q5 (scores 1,2,2,2,3; H=15; hc=3)': [16, 5, 9, 3, 14],
    'UH a (H=9)': [8, 1, 3, 22, 7],
    'UH b (H=9)': [4, 1, 10, 19, 7],
    'UH c (H=11)': [2, 20, 9, 3, 13],
    'H=13, hc=2': [2, 12, 17, 5, 11],
    'RT5 (regular; H=15; hc=2)': [6, 24, 10, 17, 5],
}
def mat(rows):
    n = len(rows)
    return sp.Matrix(n, n, lambda i, j: (rows[i] >> j) & 1)
for name, rows in classes.items():
    A = mat(rows)
    assert all(A[i, j] + A[j, i] == (0 if i == j else 1) for i in range(5) for j in range(5))
    p = sp.expand(A.charpoly(lam).as_expr())
    print(name, '| charpoly:', p, '| factor:', sp.factor(p))
    print('   eigenvalues:', [sp.N(r, 6) for r in sp.Poly(p, lam).nroots()])
# Circ(Z5,{1,2}) directly
C = sp.Matrix(5, 5, lambda i, j: 1 if (j - i) % 5 in (1, 2) else 0)
print('Circ(Z5,{1,2}) charpoly:', sp.factor(C.charpoly(lam).as_expr()))
# owner's claims
for q in [lam**5 - 2*lam**2 - 2*lam - 1, lam**2*(lam**3 - lam - 1)]:
    print('owner polynomial', sp.expand(q), 'coefficient of lam^3:', sp.Poly(q, lam).coeff_monomial(lam**3),
          '(tournaments: tr A = tr A^2 = 0 force coefficients of lam^4 and lam^3 to vanish; lam^2 coefficient = -(# cyclic triples))')
# spectrum {2, i, -i, phi, -1/phi}: sum of squares
phi = (1 + sp.sqrt(5)) / 2
print('sum of squares of {2, i, -i, phi, -1/phi} =', sp.simplify(4 + (sp.I)**2 + (-sp.I)**2 + phi**2 + (1/phi)**2), '(must be tr A^2 = 0 for a tournament)')
# PG(2,3)
M = sp.Matrix([[0, 1, 0], [1, 1, 1], [1, 2, 1]])
print('det[(0:1:0),(1:1:1),(1:2:1)] =', M.det(), '-> collinear over every field')
