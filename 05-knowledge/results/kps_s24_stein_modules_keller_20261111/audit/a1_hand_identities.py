"""a1_hand_identities.py  (independent audit of kind-pasteur-2026-11-11-S24; does not import the author's scripts)

Symbolic checks (sympy, exact) of the hand identities used in sections 2, 4 and in the draft Theorem 5.2:

 A. complex Hermite generating function:  exp(sV + tU - st) = sum H_(n,m)(U,V) s^m t^n/(m! n!)
    with H_(n,m)(U,V) = sum_j (-1)^j j! C(n,j) C(m,j) U^(n-j) V^(m-j)          (n <-> U, m <-> V)
 B. the four operators  s + d/dt, t + d/ds, d/ds, d/dt on C[s,t] satisfy the relations of U, V, d_u, d_v
 C. the Leibniz step of Theorem 4.1 on three etale maps (an automorphism, the toy on C^* x C, the slice of THM-1300):
       d_u^m d_v^n (V h) = V d_u^m d_v^n h + n d_u^m d_v^(n-1) h,   and the same with (U, s) for (V, t)
    together with the defining relations of the dual fields (d_u u = 1, d_u v = 0, d_v v = 1, d_v u = 0, [d_u, d_v] = 0)
 D. the frame identity  (d - dW ^)(alpha du + beta dv) = (T_u beta - T_v alpha) du ^ dv  (Proposition 2.1)
 E. the relations behind the family module of the draft Theorem 5.2 and the identification G = B:
       [f - d_v, g - d_u] = 0,  [d_u, f - d_v] = 1,  [d_v, g - d_u] = 1,  [d_u, g - d_u] = 0,  [d_v, f - d_v] = 0
 F. Proposition 2.3: the functional U^i V^j -> i! [i = j] kills T_u, T_v of every monomial (i, j <= 8)
"""
import sympy as sp
from math import comb, factorial

s, t, U, V = sp.symbols("s t U V")
x, y = sp.symbols("x y")
ok = True


def report(name, cond):
    global ok
    ok = ok and bool(cond)
    print("%-78s %s" % (name, "OK" if cond else "FAIL"))


# ---------------------------------------------------------------- A
def H(n, m):
    return sum((-1) ** j * factorial(j) * comb(n, j) * comb(m, j) * U ** (n - j) * V ** (m - j) for j in range(min(n, m) + 1))


N = 7
gen = sp.series(sp.exp(s * V + t * U - s * t), s, 0, N).removeO()
bad = 0
for m in range(N):
    cm = sp.expand(gen.coeff(s, m))                 # coefficient of s^m: a function of t, U, V
    ser = sp.series(cm, t, 0, N).removeO()
    for n in range(N):
        c = sp.expand(ser.coeff(t, n) * factorial(m) * factorial(n))
        if sp.expand(c - H(n, m)) != 0:
            bad += 1
report("A. generating function of H_(n,m), n, m <= 6 (49 cases)", bad == 0)
report("A'. H_(0,m) = V^m and H_(n,0) = U^n", all(H(0, m) == V ** m and H(m, 0) == U ** m for m in range(6)))

# ---------------------------------------------------------------- B
phi = sp.Function("phi")(s, t)
Uh = lambda F: s * F + sp.diff(F, t)
Vh = lambda F: t * F + sp.diff(F, s)
Ds = lambda F: sp.diff(F, s)
Dt = lambda F: sp.diff(F, t)
comm = lambda A, B, F: sp.simplify(A(B(F)) - B(A(F)))
report("B. [d/ds, s + d/dt] = 1", comm(Ds, Uh, phi) - phi == 0)
report("B. [d/dt, t + d/ds] = 1", comm(Dt, Vh, phi) - phi == 0)
report("B. [d/ds, t + d/ds] = 0", comm(Ds, Vh, phi) == 0)
report("B. [d/dt, s + d/dt] = 0", comm(Dt, Uh, phi) == 0)
report("B. [s + d/dt, t + d/ds] = 0", comm(Uh, Vh, phi) == 0)
report("B. [d/ds, d/dt] = 0", comm(Ds, Dt, phi) == 0)
# (s + d/dt)^i (t + d/ds)^j 1 = s^i t^j + lower terms
low = True
for i in range(5):
    for j in range(5):
        F = sp.Integer(1)
        for _ in range(j):
            F = Vh(F)
        for _ in range(i):
            F = Uh(F)
        F = sp.Poly(sp.expand(F), s, t)
        top = [mon for mon in F.monoms() if sum(mon) == i + j]
        low = low and top == [(i, j)] and F.coeff_monomial((i, j)) == 1
report("B. (s + d/dt)^i (t + d/ds)^j 1 = s^i t^j + lower degree (i, j <= 4)", low)


# ---------------------------------------------------------------- C, D, E on three etale maps
def dual_fields(u, v):
    J = sp.simplify(sp.diff(u, x) * sp.diff(v, y) - sp.diff(u, y) * sp.diff(v, x))
    du = lambda h: sp.simplify((sp.diff(h, x) * sp.diff(v, y) - sp.diff(h, y) * sp.diff(v, x)) / J)
    dv = lambda h: sp.simplify((sp.diff(u, x) * sp.diff(h, y) - sp.diff(u, y) * sp.diff(h, x)) / J)
    return J, du, dv


def power(op, h, k):
    for _ in range(k):
        h = op(h)
    return h


w = 1 / (1 + x * y)
maps = [
    ("automorphism (x + y^2, y + (x + y^2)^2)", x + y ** 2, y + (x + y ** 2) ** 2, [x, x * y + y ** 3]),
    ("toy (x^2, y/x) on C^* x C", x ** 2, y / x, [1 / x ** 2, x * y, y ** 2 / x]),
    ("slice of THM-1300 (c = 1)", y + 3 * x * w, x * w + x * w ** 2 - x ** 3 * w ** 3, [w, x, y * w ** 2]),
]
a0, b0 = sp.Integer(2), sp.Integer(3)
for name, u, v, tests in maps:
    J, du, dv = dual_fields(u, v)
    print("   map: %s ;  Jacobian = %s" % (name, sp.factor(J)))
    report("C. d_u u = 1, d_u v = 0, d_v v = 1, d_v u = 0", du(u) == 1 and du(v) == 0 and dv(v) == 1 and dv(u) == 0)
    UU, VV = u - a0, v - b0
    good = True
    for h in tests:
        good = good and sp.simplify(du(dv(h)) - dv(du(h))) == 0
        for m in range(0, 3):
            for n in range(0, 3):
                lhs = power(du, power(dv, VV * h, n), m)
                rhs = VV * power(du, power(dv, h, n), m) + (n * power(du, power(dv, h, n - 1), m) if n else 0)
                good = good and sp.simplify(lhs - rhs) == 0
                lhs = power(du, power(dv, UU * h, n), m)
                rhs = UU * power(du, power(dv, h, n), m) + (m * power(du, power(dv, h, n), m - 1) if m else 0)
                good = good and sp.simplify(lhs - rhs) == 0
    report("C. [d_u, d_v] = 0 and the two Leibniz identities, m, n <= 2, %d test functions" % len(tests), good)
    # D: (d - dW ^)(alpha du + beta dv) = (T_u beta - T_v alpha) du ^ dv, checked in the frame dx ^ dy
    W = UU * VV
    good = True
    for al in tests:
        for be in tests:
            p = al * sp.diff(u, x) + be * sp.diff(v, x)            # omega = p dx + q dy
            q = al * sp.diff(u, y) + be * sp.diff(v, y)
            lhs = sp.diff(q, x) - sp.diff(p, y) - (sp.diff(W, x) * q - sp.diff(W, y) * p)
            rhs = ((du(be) - VV * be) - (dv(al) - UU * al)) * J
            good = good and sp.simplify(lhs - rhs) == 0
    report("D. twisted differential in the frame du ^ dv (Proposition 2.1)", good)
    # E: operators a = u - d_v, b = v - d_u (centre 0), D_a = d_u, D_b = d_v
    A_ = lambda h: u * h - dv(h)
    B_ = lambda h: v * h - du(h)
    good = True
    for h in tests:
        good = good and sp.simplify(A_(B_(h)) - B_(A_(h))) == 0
        good = good and sp.simplify(du(A_(h)) - A_(du(h)) - h) == 0
        good = good and sp.simplify(dv(B_(h)) - B_(dv(h)) - h) == 0
        good = good and sp.simplify(du(B_(h)) - B_(du(h))) == 0
        good = good and sp.simplify(dv(A_(h)) - A_(dv(h))) == 0
    report("E. Weyl relations of (u - d_v, v - d_u, d_u, d_v)  [family module G = B]", good)

# ---------------------------------------------------------------- F
L = lambda i, j: factorial(i) if i == j else 0
bad = 0
for i in range(9):
    for j in range(9):
        # T_u (U^i V^j) = i U^(i-1) V^j - U^i V^(j+1);  T_v (U^i V^j) = j U^i V^(j-1) - U^(i+1) V^j
        tu = (i * L(i - 1, j) if i else 0) - L(i, j + 1)
        tv = (j * L(i, j - 1) if j else 0) - L(i + 1, j)
        bad += (tu != 0) + (tv != 0)
report("F. U^i V^j -> i! [i = j] kills T_u, T_v (i, j <= 8)", bad == 0)

print("ALL OK" if ok else "SOME CHECK FAILED")
