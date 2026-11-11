"""b6_prop63_special_centre.py - second audit, README 6.3 (Proposition 6.3): X = C^2 minus {1 + xy = 0},
F = inclusion, centre (0, 0), W = xy.  Own engine (eng2: charge basis x^c q^n / y^c q^n, q = 1 + xy = 1/w,
prime 2^61 - 1, filtration |c| + 2|n|; the author uses the basis x^i y^j w^k, the weight |i| + j + k, p = 2147483629).

 A. the formulas of the proof, symbolically (phi an arbitrary function of t = xy);
 B. charge 0: M_z cap R = (d/dt - 1)(t R) is spanned by g_n = -q^(n+1) + (n+2) q^n - n q^(n-1) (n in Z), of
    codimension 2 in R = C[q, 1/q] with complement <1, 1/q>; the two functionals
        lam_G(phi) = int_0^inf phi e^(-t) dt,     lam_res(phi) = Res_(t = -1) phi e^(-t) dt
    kill it; lam_G(1) = 1, lam_res(1) = 0, lam_res(1/q) = e;
 C. machine: c(N, N') at (0, 0) and (2, 3); classes of d_u^m h, d_v^m h, (d_u d_v)^m w; classes [h U^k V^(m+k)].
"""
import sys
import sympy as sp
import eng2
from eng2 import RingH, add, scal, lincomb, power, Echelon

R = RingH()
x, y, one = R.x(), R.y(), R.one()
w = R.q(-1)
wt = lambda m: abs(m[0]) + 2 * abs(m[1])
keyf = lambda m: (wt(m), m)


def partA():
    X, Y, t = sp.symbols("x y t")
    phi = sp.Function("phi")
    k = sp.symbols("k", positive=True, integer=True)
    T = X * Y
    e1 = sp.diff(X ** (k - 1) * phi(T), Y) - X * X ** (k - 1) * phi(T)                 # (d/dy - x)(x^(k-1) phi)
    target1 = X ** k * (sp.diff(phi(t), t) - phi(t)).subs(t, T)
    e2 = sp.diff(X ** (k + 1) * phi(T), X) - Y * X ** (k + 1) * phi(T)                 # (d/dx - y)(x^(k+1) phi)
    target2 = X ** k * (k * phi(t) + sp.diff(t * phi(t), t) - t * phi(t)).subs(t, T)
    e3 = sp.diff(X * phi(T), X) - Y * X * phi(T)                                       # charge 0, from x phi
    e4 = sp.diff(Y * phi(T), Y) - X * Y * phi(T)                                       # charge 0, from y phi
    target0 = (sp.diff(t * phi(t), t) - t * phi(t)).subs(t, T)
    # k <= -1 by the symmetry x <-> y
    e5 = sp.diff(Y ** (k - 1) * phi(T), X) - Y * Y ** (k - 1) * phi(T)
    target5 = Y ** k * (sp.diff(phi(t), t) - phi(t)).subs(t, T)
    e6 = sp.diff(Y ** (k + 1) * phi(T), Y) - X * Y ** (k + 1) * phi(T)
    target6 = Y ** k * (k * phi(t) + sp.diff(t * phi(t), t) - t * phi(t)).subs(t, T)
    res = [sp.simplify((a - b_).doit()) for a, b_ in ((e1, target1), (e2, target2), (e3, target0), (e4, target0), (e5, target5), (e6, target6))]
    print("A. proof formulas (k >= 1: two; charge 0: two; k <= -1: two), differences: %s" % res)


def partB():
    t = sp.symbols("t")
    q = 1 + t
    lamG = lambda f: sp.integrate(sp.expand(f) * sp.exp(-t), (t, 0, sp.oo))
    lamR = lambda f: sp.residue(f * sp.exp(-t), t, -1)
    bad = 0
    for n in range(-6, 7):
        g = -q ** (n + 1) + (n + 2) * q ** n - n * q ** (n - 1)
        chk = sp.simplify(sp.diff(t * q ** n, t) - t * q ** n - g)
        if chk != 0:
            bad += 1
        if sp.simplify(lamR(g)) != 0:
            bad += 1
        if n >= 0 and sp.simplify(lamG(g)) != 0:          # closed form available for polynomials
            bad += 1
    # lam_G on negative powers numerically (exponential integral)
    import mpmath as mp
    mp.mp.dps = 30
    num_bad = 0
    for n in range(-6, 0):
        f = lambda tt, n=n: (-(1 + tt) ** (n + 1) + (n + 2) * (1 + tt) ** n - n * (1 + tt) ** (n - 1)) * mp.e ** (-tt)
        val = mp.quad(f, [0, 1, 5, 20, mp.inf])
        if abs(val) > mp.mpf(10) ** (-20):
            num_bad += 1
    print("B. g_n = (d/dt - 1)(t q^n) = -q^(n+1) + (n+2) q^n - n q^(n-1), |n| <= 6: %d symbolic failures (identity, lam_res, lam_G for n >= 0); lam_G(g_n), n < 0, numerically (30 digits): %d nonzero"
          % (bad, num_bad))
    print("   lam_G(1) = %s, lam_res(1) = %s, lam_res(1/q) = %s, lam_G(1/q) = %s (= e E_1(1))"
          % (lamG(sp.Integer(1)), lamR(sp.Integer(1)), lamR(1 / q), mp.nstr(mp.e * mp.e1(1), 15)))
    # codimension 2 in a window of R = C[q, 1/q]: rank of g_n, |n| <= N, modulo <1, q^-1>
    rows = [{n + 1: -1, n: n + 2, n - 1: -n} for n in range(-40, 41)]
    rows = [{m: c % eng2.PRIME for m, c in r.items() if c % eng2.PRIME} for r in rows]
    rk = eng2.rank_rows(rows)
    rk2 = eng2.rank_rows(rows + [{0: 1}, {-1: 1}])
    cols = sorted({m for r in rows for m in r})
    print("   window |n| <= 40: %d generators, rank %d, columns %d (q^-41 .. q^41); rank with 1 and 1/q added: %d  => codimension 2, complement <1, w>, and 1 not in M_z"
          % (len(rows), rk, len(cols), rk2))


def build(av, bv, Np):
    U = add(x, scal(one, av), -1)
    V = add(y, scal(one, bv), -1)
    W = R.mul(U, V)
    Wx, Wy = R.dx(W), R.dy(W)
    E = Echelon(keyf)
    basis = [(c, n) for c in range(-Np, Np + 1) for n in range(-Np, Np + 1) if abs(c) + 2 * abs(n) <= Np]
    for m in basis:
        h = {m: 1}
        E.insert(add(R.dx(h), R.mul(Wx, h), -1))
        E.insert(add(R.dy(h), R.mul(Wy, h), -1))
    return E, U, V


def codim(E, N):
    dimV = sum(1 for c in range(-N, N + 1) for n in range(-N, N + 1) if abs(c) + 2 * abs(n) <= N)
    return dimV - sum(1 for m in E.piv if wt(m) <= N)


def partC():
    for (av, bv) in ((0, 0), (2, 3)):
        for Np in (28, 36):
            E, U, V = build(av, bv, Np)
            print("C. centre (%d, %d), generators of weight <= %d (%d pivots): c(N) = %s"
                  % (av, bv, Np, len(E.piv), " ".join("%d:%d" % (N, codim(E, N)) for N in (2, 4, 8, 12, 16, 20) if N <= Np - 8)))
            sys.stdout.flush()
        cls = lambda h: ("?" if (h and max(wt(m) for m in h) > Np - 8) else ("N" if not E.contains(h) else "0")) if h else "=0"
        print("   class of 1: %s, of w: %s, of w^2: %s, of x: %s, of x y: %s" % tuple(cls(h) for h in (one, w, R.mul(w, w), x, R.mul(x, y))))
        tests = (("w", w, 0), ("w^2", R.mul(w, w), 0), ("x w", R.mul(x, w), 1), ("y w", R.mul(y, w), -1), ("x^2 w", R.mul(power(R, x, 2), w), 2),
                 ("y^3 w^2", R.mul(power(R, y, 3), R.mul(w, w)), -3), ("x w + y^2 w", add(R.mul(x, w), R.mul(power(R, y, 2), w)), None))
        for name, h, ch in tests:
            outs = []
            for lab, d in (("d_u", R.dx), ("d_v", R.dy)):
                res, g = [], h
                for m in range(1, 8):
                    g = d(g)
                    res.append(cls(g))
                outs.append("%s^m: %s" % (lab, " ".join(res)))
            print("   h = %-12s %s" % (name, ";  ".join(outs)))
        res, g = [], w
        for m in range(1, 5):
            g = R.dy(R.dx(g))
            res.append(cls(g))
        print("   (d_u d_v)^m w, m = 1..4: %s" % " ".join(res))
        # cone vanishing: classes of h U^k V^(m+k) and h U^(m+k) V^k for h = w, x w, y^2 w
        for name, h in (("w", w), ("x w", R.mul(x, w)), ("y^2 w", R.mul(power(R, y, 2), w))):
            rows_ = []
            for m in range(0, 5):
                rr = []
                for k in range(0, 4):
                    e1 = R.mul(h, R.mul(power(R, U, k), power(R, V, m + k)))
                    e2 = R.mul(h, R.mul(power(R, U, m + k), power(R, V, k)))
                    rr.append(cls(e1) + cls(e2))
                rows_.append("m=%d: %s" % (m, " ".join(rr)))
            print("   cone, h = %-6s [h U^k V^(m+k)][h U^(m+k) V^k], k = 0..3:  %s" % (name, " | ".join(rows_)))
        sys.stdout.flush()


if __name__ == "__main__":
    partA()
    partB()
    partC()
