#!/usr/bin/env python3
# audit J1 (THM-4910): symbolic re-derivation in the Gram basis (v, s, w), all n: S2 DEFINED by the four corner powers (not by THM-4742's coefficients), then all eight labels, the defect, the sign facts and the quadruple determinants.
import sympy as sp

n, R2, vs, vw, ss, k = sp.symbols('n R2 vs vw ss k')
# Gram matrix of v = V0-O, s = sum_{i>=1}(Vi-O), w = O'-O ;  s.w = n k because (Vi - O').w = 0
Gm = sp.Matrix([[R2, vs, vw], [vs, ss, n*k], [vw, n*k, k]])
def dot(a, b): return sp.expand((sp.Matrix(a).T*Gm*sp.Matrix(b))[0])
a, b, c, kap = sp.symbols('a b c kap')        # a = 2u.v, b = 2u.s, c = 2u.w, kap = constant term
def f(x):                                      # power function relative to O, x given by basis coefficients
    return dot(x, x) - (a*x[0] + b*x[1] + c*x[2]) + kap
O = [0, 0, 0]; G = [1/(n+1), 1/(n+1), 0]; Gp = [0, 1/n, 0]; Op = [0, 0, 1]
eqs = [sp.Eq(f(O), (n+1)*R2), sp.Eq(f(G), n*(R2 - dot(G, G))),
       sp.Eq(f(Gp), (n-1)*(R2 - dot(Gp, Gp))), sp.Eq(f(Op), n*(R2 - dot(Op, Op)))]
sol = sp.solve(eqs, [a, b, c, kap], dict=True)
assert len(sol) == 1
sol = {key: sp.simplify(val) for key, val in sol[0].items()}
print("S2 from corner powers: 2u.v =", sol[a], "; 2u.s =", sol[b], "; 2u.w =", sol[c], "; kappa =", sol[kap])
thm4742 = {a: 2*vs + (2 - n)*R2, b: ss + 2*n*R2, c: (n + 1)*k + R2, kap: (n + 1)*R2}
print("equals THM-4742's S2 coefficients:", all(sp.simplify(sol[key] - thm4742[key]) == 0 for key in thm4742))
def F2(x): return sp.simplify(sp.together(f(x).subs(sol)))

T = dot([1, 1, 0], [1, 1, 0])                          # |t|^2, t = v + s = (n+1)(G - O)
gam2 = dot([0, 1/n, -1], [0, 1/n, -1])                 # |G' - O'|^2
rp2 = R2 - k                                           # facet circumradius^2 (|OO'|^2 = k)
mu = -2*dot([1, 0, 0], [-n, 1, 0]) / dot([-n, 1, 0], [-n, 1, 0])
labels = {
    "A":  [1 - n*mu, mu, 0],                           # v + mu (s - n v): second point of line V0G on S
    "B":  [0, 0, R2/k],                                # pole of the facet hyperplane
    "C*": [-(n - 1), 1, 0],                            # (n+1)G - n V0 = s - (n-1) v
    "D":  [1, 1, 0],                                   # O + (n+1)(G - O) = t
    "E":  [(n + 1)*R2/T, (n + 1)*R2/T, 0],             # inverse of G in S
    "F":  [0, 0, n + 1],                               # O + (n+1)(O' - O)
    "H'": [0, 1, -(n - 1)],                            # O' + n(G' - O') = s - (n-1) w
    "X":  [0, (rp2/gam2)/n, 1 - rp2/gam2],             # O' + r'^2 (G'-O')/|G'-O'|^2
}
# sanity: A is on the circumsphere and on the line V0 G
print("|A|^2 - R^2 =", sp.simplify(dot(labels["A"], labels["A"]) - R2))
for name, x in labels.items():
    print(f"Pow_S2({name}) =", F2(x))
Sigma = 2*n*R2 - 2*vs                                  # sum_i |V0 - Vi|^2 (|Vi - O| = R)
print("Pow_S2(V0) - Sigma =", F2([1, 0, 0]) - Sigma, "(Sigma = 2nR^2 - 2 v.s)")
# V0A * V0C*: A - V0 = mu (s - n v), C* - V0 = s - n v  -> signed product mu |s - n v|^2
print("V0A*V0C* (signed) - Sigma =", sp.simplify(mu*dot([-n, 1, 0], [-n, 1, 0]) - Sigma))
# corner inversion radii: powers at the corners are c_P (R^2 - |OP|^2) > 0 when the corner is inside S
for name, P, cp in (("O", O, n + 1), ("G", G, n), ("G'", Gp, n - 1), ("O'", Op, n)):
    print(f"Pow_S2({name}) - c(R^2-|OP|^2) =", sp.simplify(F2(P) - cp*(R2 - dot(P, P))))
# coplanarity of the 'top' (homothety) and 'bottom' (inversion/chord) quadruples, in basis coordinates
def det4(P1, P2, P3, P4):
    M = sp.Matrix([[P2[i] - P1[i] for i in range(3)], [P3[i] - P1[i] for i in range(3)], [P4[i] - P1[i] for i in range(3)]])
    return sp.factor(sp.simplify(M.det()))
print("det(top {D,F,C*,H'}) in basis coords =", det4(labels["D"], labels["F"], labels["C*"], labels["H'"]))
dbot = det4(labels["E"], labels["B"], labels["A"], labels["X"])
print("det(bottom {E,B,A,X}) identically zero?", sp.simplify(dbot) == 0)
# numeric instance of the bottom determinant at a random admissible Gram point (n = 3)
import random
random.seed(1)
inst = {n: 3, R2: 1, vs: sp.Rational(-1, 3), vw: sp.Rational(1, 7), ss: sp.Rational(5, 2), k: sp.Rational(1, 5)}
print("bottom det at a sample point (n=3):", sp.nsimplify(dbot.subs(inst)))
