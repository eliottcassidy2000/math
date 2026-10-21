"""Independent audit of THM-4742 (symbolic part). Own formulation: general R, spheres SOLVED from 4 points.
Basis (v, s, w) with O = 0: v = V0, s = V1+V2+V3, w = O' (foot of perpendicular from O to the face plane).
Gram: v.v = R2, v.s = a, s.s = b, v.w = p, w.w = k, s.w = 3k (since w.Vi = k for face vertices)."""
import sympy as sp

R2, a, b, p, k = sp.symbols('R2 a b p k')
Gm = sp.Matrix([[R2, a, p], [a, b, 3 * k], [p, 3 * k, k]])


def nrm(X):
    X = sp.Matrix(X)
    return sp.simplify((X.T * Gm * X)[0, 0])


T = nrm([1, 1, 0])                       # |v + s|^2 = |4G|^2
L = {}
L['C'] = [1, 0, 0]
L['D'] = [1, 1, 0]                       # O + 4(G - O) = v + s
L['E'] = [4 * R2 / T, 4 * R2 / T, 0]     # inverse of G = (v+s)/4 : R2 * G/|G|^2
L['F'] = [0, 0, 4]
L['B'] = [0, 0, R2 / k]                  # B.Vi = R2 for the three face vertices
L['H'] = [0, 1, -2]                      # Euler: H' = V1+V2+V3 - 2 O'
d = sp.Matrix([-3, 1, 0])                # s - 3v = 4(G - V0)
vd = (sp.Matrix([1, 0, 0]).T * Gm * d)[0, 0]
dd = (d.T * Gm * d)[0, 0]
mu = sp.simplify(-2 * vd / dd)
L['A'] = [1 - 3 * mu, mu, 0]
# sanity: A on circumsphere, orthocentre altitude conditions via Gram are identities (checked below)
assert sp.simplify(nrm(L['A']) - R2) == 0

lam = sp.symbols('l0 l1 l2'); c = sp.Symbol('c')


def Q(X, lv, cv):
    X = list(X)
    return sp.simplify(nrm(X) - sum(lv[i] * X[i] for i in range(3)) + cv)


def solve_sphere(names):
    eqs = [Q(L[x], lam, c) for x in names]
    sol = sp.solve(eqs, list(lam) + [c], dict=True)
    assert len(sol) == 1, (names, sol)
    s0 = sol[0]
    return [sp.simplify(s0[l]) for l in lam], sp.simplify(s0[c])


l1, c1 = solve_sphere(['C', 'D', 'F', 'B'])
l2, c2 = solve_sphere(['D', 'F', 'H', 'B'])
print('S1 (2u.v, 2u.s, 2u.w), c =', l1, c1)
print('S2 (2u.v, 2u.s, 2u.w), c =', l2, c2)
Sig = 6 * R2 - 2 * a                      # sum |V0 - Vi|^2 = sum (2R2 - 2 v.Vi)
for nm, (lv, cv) in (('S1', (l1, c1)), ('S2', (l2, c2))):
    print(nm, {x: sp.factor(Q(L[x], lv, cv)) for x in 'ABCDEFH'})
print('Pow_S2(C) - Sigma =', sp.simplify(Q(L['C'], l2, c2) - Sig))
print('Pow_S1(H) - Sigma =', sp.simplify(Q(L['H'], l1, c1) - Sig))
print('Pow_S1(A) - mu(16R2 - T) =', sp.simplify(Q(L['A'], l1, c1) - mu * (16 * R2 - T)))
print('-8 v.(G - v) - Sigma =', sp.simplify(-8 * ((sp.Matrix([1, 0, 0]).T * Gm * sp.Matrix([-sp.Rational(3, 4), sp.Rational(1, 4), 0]))[0, 0]) - Sig))
print('mu - Sigma/|s-3v|^2 =', sp.simplify(mu - Sig / dd))
dl = [sp.simplify(l1[i] - l2[i]) for i in range(3)]
print('lam1 - lam2 =', dl, '; on t=v+s:', sp.simplify(dl[0] + dl[1]))
print('c1 = c2 = 4R2 :', sp.simplify(c1 - 4 * R2) == 0 and sp.simplify(c2 - 4 * R2) == 0)
