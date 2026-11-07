#!/usr/bin/env python3
# audit J1 (THM-4910): exact rational checks in actual coordinates, n = 3..7: S2 from the corner powers, eight labels, Vieta positions and signs on each side, defect, distinctness; coincidence examples inside the regime.
import random
from fractions import Fraction as Fr
import sympy as sp

def sub(a, b): return [x - y for x, y in zip(a, b)]
def add(a, b): return [x + y for x, y in zip(a, b)]
def mul(c, a): return [c * x for x in a]
def dt(a, b): return sum(x * y for x, y in zip(a, b))
def solve(M, rhs):
    m = len(M); A = [list(r) + [b] for r, b in zip(M, rhs)]
    for col in range(m):
        piv = next((r for r in range(col, m) if A[r][col] != 0), None)
        if piv is None: raise ZeroDivisionError("singular")
        A[col], A[piv] = A[piv], A[col]
        for r in range(m):
            if r != col and A[r][col] != 0:
                fct = A[r][col] / A[col][col]; A[r] = [x - fct * y for x, y in zip(A[r], A[col])]
    return [A[i][m] / A[i][i] for i in range(m)]
def det3(Gm):
    return (Gm[0][0]*(Gm[1][1]*Gm[2][2]-Gm[1][2]*Gm[2][1]) - Gm[0][1]*(Gm[1][0]*Gm[2][2]-Gm[1][2]*Gm[2][0])
            + Gm[0][2]*(Gm[1][0]*Gm[2][1]-Gm[1][1]*Gm[2][0]))
def circumcentre(P):
    # point of aff(P) equidistant from all points of P
    p0 = P[0]; D = [sub(p, p0) for p in P[1:]]
    lam = solve([[dt(x, y) for y in D] for x in D], [dt(x, x) / 2 for x in D])
    return add(p0, [sum(l * d[j] for l, d in zip(lam, D)) for j in range(len(p0))])

def analyse(V, zero=0):
    N = len(V) - 1
    O = circumcentre(V); R2 = dt(sub(V[0], O), sub(V[0], O))
    G = mul(Fr(1, N + 1), [sum(c) for c in zip(*V)])
    Fc = V[1:]; Op = circumcentre(Fc); Gp = mul(Fr(1, N), [sum(c) for c in zip(*Fc)])
    rp2 = dt(sub(Fc[0], Op), sub(Fc[0], Op))
    v, w = sub(V[0], O), sub(Op, O); s = sub(mul(N, Gp), mul(N, O))
    basis = [v, s, w]; Gb = [[dt(x, y) for y in basis] for x in basis]
    regime = det3(Gb) != 0
    if not regime: return None
    # labels
    dV = sub(G, V[0]); mu = -2 * dt(sub(V[0], O), dV) / dt(dV, dV); A = add(V[0], mul(mu, dV))
    B = add(O, mul(R2 / dt(w, w), w))
    Cs = add(mul(N + 1, G), mul(-N, V[0]))
    D = add(O, mul(N + 1, sub(G, O))); g = sub(G, O); E = add(O, mul(R2 / dt(g, g), g))
    F = add(O, mul(N + 1, w)); H = add(Op, mul(N, sub(Gp, Op)))
    gg = sub(Gp, Op); X = add(Op, mul(rp2 / dt(gg, gg), gg))
    L = {"A": A, "B": B, "C*": Cs, "D": D, "E": E, "F": F, "H'": H, "X": X}
    # S2 from corner powers: f(x) = |x-O|^2 - 2 sum_i al_i (b_i.(x-O)) + kap
    corners = [(O, N + 1), (G, N), (Gp, N - 1), (Op, N)]
    rows, rhs = [], []
    for P, cP in corners:
        y = sub(P, O)
        rows.append([-2 * dt(bb, y) for bb in basis] + [Fr(1)])
        rhs.append(cP * (R2 - dt(y, y)) - dt(y, y))
    sol = solve(rows, rhs); al, kap = sol[:3], sol[3]
    def pw(P):
        y = sub(P, O); return dt(y, y) - 2 * sum(a_ * dt(bb, y) for a_, bb in zip(al, basis)) + kap
    def in_flat(P):
        y = sub(P, O); cf = solve(Gb, [dt(bb, y) for bb in basis])
        return add(O, [sum(ci * bb[j] for ci, bb in zip(cf, basis)) for j in range(len(O))]) == P
    res = {}
    res["flat"] = all(in_flat(P) for P in L.values())
    res["on_S2"] = all(pw(P) == 0 for P in L.values())
    res["real_radius"] = (dt(cen := [sum(a_ * bb[j] for a_, bb in zip(al, basis)) for j in range(len(O))], cen) - kap) > 0
    # side lines with their parameters (unit = |direction| not needed: use parameter tau, line P0 + tau*d)
    def tau(P0, d, P):
        # P = P0 + tau d (assert collinear)
        t = dt(sub(P, P0), d) / dt(d, d); assert add(P0, mul(t, d)) == P, "not on line"; return t
    sides = {"OG": (O, sub(G, O), ("D", "E"), (0, 1)), "OO'": (O, sub(Op, O), ("F", "B"), (0, 1)),
             "O'G'": (Op, sub(Gp, Op), ("H'", "X"), (0, 1)), "V0-median": (V[0], sub(Gp, V[0]), ("C*", "A"), (Fr(N, N + 1), 1))}
    ok_sides = True; ok_beyond = True
    for nm, (P0, d, (p, q), (t_first, t_second)) in sides.items():
        tp, tq = tau(P0, d, L[p]), tau(P0, d, L[q])
        # quadratic of the power along the line: |d|^2 t^2 + beta t + gamma; check roots are exactly tp, tq
        dd = dt(d, d); gam = pw(P0); beta = pw(add(P0, d)) - dd - gam
        ok_sides &= (dd * tp * tq == gam) and (-dd * (tp + tq) == beta)
        # both points beyond both corners of the side (corners at parameters t_first, t_second, t_first < t_second)
        ok_beyond &= tp > t_second and tq > t_second
    res["sides_vieta"] = ok_sides; res["beyond_far_corner"] = ok_beyond
    # corner powers all positive (all four corners outside S2)
    res["corners_outside"] = all(pw(P) > 0 for P, _ in corners)
    Sigma = sum(dt(sub(V[0], Vi), sub(V[0], Vi)) for Vi in V[1:])
    d_med = sub(Gp, V[0])
    prod_signed = tau(V[0], d_med, A) * tau(V[0], d_med, Cs) * dt(d_med, d_med)
    res["defect"] = (pw(V[0]) == Sigma == prod_signed) and Sigma > 0
    res["distinct"] = len({tuple(P) for P in L.values()})
    res["mod2"] = True
    return res, L

random.seed(20261107 + 1)
summary = {}
for N, trials in ((3, 200), (4, 60), (5, 30), (6, 15), (7, 8)):
    done = 0; agg = {}
    while done < trials:
        V = [[Fr(random.randint(-7, 7), random.randint(1, 3)) for _ in range(N)] for _ in range(N + 1)]
        try:
            out = analyse(V)
        except ZeroDivisionError:
            continue
        if out is None: continue
        res, _ = out; done += 1
        for key, val in res.items():
            if key == "distinct": agg.setdefault("distinct8", 0); agg["distinct8"] += (val == 8)
            else: agg.setdefault(key, 0); agg[key] += bool(val)
    print(f"n={N}: {trials} random exact rational simplices in the regime:", agg)

# ---------- coincidences inside the regime ----------
print("\nCoincidences inside the regime (labels on one side may coincide; the side is then tangent to S2):")
# (i) session23b's exceptional example (sqrt 3): D = E and B = F
r3 = sp.sqrt(3)
Vx = [sp.Matrix([0, 1, 0]), sp.Matrix([r3/2, 0, sp.Rational(1, 2)]), sp.Matrix([0, r3/2, sp.Rational(1, 2)]), sp.Matrix([0, -r3/2, sp.Rational(1, 2)])]
O = sp.Matrix([0, 0, 0]); R2 = (Vx[0] - O).dot(Vx[0] - O)
assert all(sp.simplify((P - O).dot(P - O) - R2) == 0 for P in Vx)
G = sum(Vx, sp.zeros(3, 1)) / 4; Op = sp.Matrix([0, 0, sp.Rational(1, 2)]); Gp = sum(Vx[1:], sp.zeros(3, 1)) / 3
assert all(sp.simplify((P - Op).dot(P - Op) - (Vx[1] - Op).dot(Vx[1] - Op)) == 0 for P in Vx[1:])
v, s, w = Vx[0] - O, 3 * (Gp - O), Op - O
print("  session23b example: det[v s w] =", sp.simplify(sp.Matrix.hstack(v, s, w).det()), "(nonzero => in THM-4910's regime)")
D = O + 4 * (G - O); E = O + R2 * (G - O) / (G - O).dot(G - O); F = O + 4 * w; B = O + R2 * w / w.dot(w)
print("  D == E:", sp.simplify(D - E) == sp.zeros(3, 1), "; B == F:", sp.simplify(B - F) == sp.zeros(3, 1))

# (ii) rational example with C* = A: right-triangle facet, V0 on the sphere |V0 - G'|^2 = (r'^2 - gamma^2)/(n-1)
V = [[Fr(2), Fr(2, 3), Fr(4, 3)], [Fr(0), Fr(0), Fr(0)], [Fr(4), Fr(0), Fr(0)], [Fr(0), Fr(2), Fr(0)]]
res, L = analyse(V)
print("  C*=A example V0=(2,2/3,4/3), facet (0,0,0),(4,0,0),(0,2,0):", res, " C*==A:", L["C*"] == L["A"])

# (iii) H' = X: facet on a circle with 3|O'G'|^2 = r'^2 (n = 3); search rational points on the unit circle
def circ(t): t = Fr(t); return (Fr(1) - t*t) / (1 + t*t), 2*t / (1 + t*t)
cands = sorted({Fr(p, q) for p in range(-12, 13) for q in range(1, 13)})
found = None
for i, t1 in enumerate(cands):
    for j in range(i + 1, len(cands)):
        for l in range(j + 1, len(cands)):
            P = [circ(t1), circ(cands[j]), circ(cands[l])]
            sx, sy = sum(p[0] for p in P), sum(p[1] for p in P)
            if sx * sx + sy * sy == 3:
                found = P; break
        if found: break
    if found: break
print("  H'=X facet search (unit circle, |sum|^2 = 3):", found)
if found:
    for V0 in ([Fr(1, 3), Fr(1, 5), Fr(1)], [Fr(2, 7), Fr(-1, 3), Fr(3, 2)]):
        V = [V0] + [[p[0], p[1], Fr(0)] for p in found]
        out = analyse(V)
        if out:
            res, L = out
            print("  H'=X example V0 =", V0, res, " H'==X:", L["H'"] == L["X"])
            break
print("  (no rational example with an axis-aligned facet can exist: 3 r'^2 = |sum of facet vertices - 3 O'|^2 is never a sum of two rational squares)")
# exact algebraic example of H' = X: facet on the unit circle at angles 0, +-alpha with (1 + 2 cos alpha)^2 = 3
ca = (r3 - 1) / 2; sa = sp.sqrt(1 - ca**2)
W = [sp.Matrix([sp.Rational(1, 3), sp.Rational(1, 5), 1]), sp.Matrix([1, 0, 0]), sp.Matrix([ca, sa, 0]), sp.Matrix([ca, -sa, 0])]
Op = sp.Matrix([0, 0, 0]); rp2 = 1                                  # facet circumcentre and circumradius^2
z0 = sp.nsimplify((W[0].dot(W[0]) - rp2) / (2 * W[0][2]))          # circumcentre O = (0, 0, z0)
O = sp.Matrix([0, 0, z0]); R2 = sp.simplify((W[0] - O).dot(W[0] - O))
assert all(sp.simplify((P - O).dot(P - O) - R2) == 0 for P in W)
G = sum(W, sp.zeros(3, 1)) / 4; Gp = sum(W[1:], sp.zeros(3, 1)) / 3
v, s, w = W[0] - O, 3 * (Gp - O), Op - O
H = Op + 3 * (Gp - Op); X = Op + rp2 * (Gp - Op) / (Gp - Op).dot(Gp - Op)
print("  algebraic H'=X example: det[v s w] =", sp.nsimplify(sp.simplify(sp.Matrix.hstack(v, s, w).det())),
      "; H' == X:", sp.simplify(H - X) == sp.zeros(3, 1), "; G != O:", sp.simplify(G - O) != sp.zeros(3, 1))
