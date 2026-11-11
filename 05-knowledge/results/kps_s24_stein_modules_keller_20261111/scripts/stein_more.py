"""stein_more.py - run S6 (linear potentials on the slice) and exploratory tables (kind-pasteur-2026-11-11-S24).

S6: dim H^2(X_1, d - dW ^) for W = u, v, u + v on the slice F_1 = 1 of THM-1300 (predicted 2, 4, 6).
Exploratory (no prediction): the Hermite classes [d_u^m d_v^n h] of h = w, x on the slice (centre (2, 3)), computed as
classes of h He_(n,m)(U, V) J; other centres; automorphisms of C^2 (dimension 1, all Hermite classes eventually 0)."""
import sys, time
from math import comb, factorial
from stein import Ring, twisted, codim


def slice_data():
    R = Ring("C")
    x, y, w, one = R.mono(1, 0), R.mono(0, 1), R.mono(0, 0, 1), R.mono(0, 0)
    u = R.add(y, R.scal(R.mul(x, w), 3))
    v = R.add(R.add(R.mul(x, w), R.mul(x, R.mul(w, w))), R.mul(R.power(x, 3), R.power(w, 3)), -1)
    J = R.scal(R.power(w, 3), -2)
    return R, x, y, w, one, u, v, J


def hermite(R, U, V, n, m, one):
    """complex Hermite polynomial H_(n,m)(U, V) = sum_j (-1)^j j! C(n,j) C(m,j) U^(n-j) V^(m-j)."""
    out = {}
    for j in range(min(n, m) + 1):
        term = R.mul(R.power(U, n - j), R.power(V, m - j))
        out = R.add(out, R.scal(term, ((-1) ** j) * factorial(j) * comb(n, j) * comb(m, j)))
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "s6"
    R, x, y, w, one, u, v, J = slice_data()
    if mode == "s6":
        for name, W, pred in (("u", u, 2), ("v", v, 4), ("u + v", R.add(u, v), 6)):
            for Np in (16, 22, 28):
                E = twisted(R, W, Np)
                print("W = %-5s N' = %d: %s   predicted %d" % (name, Np, " ".join("c(%d)=%d" % (N, codim(R, E, N)) for N in (4, 8, 12, 16) if N <= Np), pred))
                sys.stdout.flush()
    elif mode == "centres":
        for (a, b) in ((2, 3), (5, -7), (1, 1), (0, 0), (-3, 4)):
            W = R.mul(R.add(u, R.scal(one, a), -1), R.add(v, R.scal(one, b), -1))
            E = twisted(R, W, 22)
            print("centre (%d, %d): %s" % (a, b, " ".join("c(%d)=%d" % (N, codim(R, E, N)) for N in (4, 8, 12, 16))))
            sys.stdout.flush()
    elif mode == "hermite":
        Np = int(sys.argv[2]) if len(sys.argv) > 2 else 40
        U = R.add(u, R.scal(one, 2), -1)
        V = R.add(v, R.scal(one, 3), -1)
        W = R.mul(U, V)
        t0 = time.time()
        E = twisted(R, W, Np)
        print("slice, centre (2, 3), N' = %d: c(30) = %d  (%.0f s)" % (Np, codim(R, E, 30), time.time() - t0))
        for name, h in (("w", w), ("x", x), ("1", one)):
            print("Hermite classes [h H_(n,m)(U, V)] for h = %s (rows n = 0.., columns m = 0..; N = nonzero, 0 = zero, ? = weight above 30):" % name)
            for n in range(0, 5):
                row = []
                for m in range(0, 5):
                    el = R.mul(R.mul(h, hermite(R, U, V, n, m, one)), J)
                    wt = max(Ring.weight(mm) for mm in el) if el else 0
                    if wt > 30:
                        row.append("?")
                    else:
                        row.append("N" if E.normal_form(el) else "0")
                print("   n = %d: %s" % (n, " ".join(row)))
            sys.stdout.flush()
    elif mode == "auto":
        A = Ring("A")
        X, Y, o = A.mono(1, 0), A.mono(0, 1), A.mono(0, 0)
        f = A.add(X, A.mul(Y, Y))                       # (x + y^2, y) then (f, y + f^2)
        g = A.add(Y, A.mul(f, f))
        U = A.add(f, A.scal(o, 2), -1)
        V = A.add(g, A.scal(o, 3), -1)
        E = twisted(A, A.mul(U, V), 30)
        print("automorphism (x + y^2, y + (x + y^2)^2), centre (2, 3): c(10) = %d, c(16) = %d" % (codim(A, E, 10), codim(A, E, 16)))
        for name, h in (("x", X), ("x y", A.mul(X, Y))):
            for n in range(0, 4):
                row = []
                for m in range(0, 4):
                    el = A.mul(h, hermite(A, U, V, n, m, o))
                    wt = max(Ring.weight(mm) for mm in el) if el else 0
                    row.append("?" if wt > 16 else ("N" if E.normal_form(el) else "0"))
                print("   h = %s, n = %d: %s" % (name, n, " ".join(row)))


def identity_check():
    """Proposition 2 (d): [d_u^m d_v^n h] = [h H_(n,m)(U, V)] in the Stein space, on the slice and on an automorphism.
    Dual fields: d_u h = (h_x v_y - h_y v_x)/J, d_v h = (u_x h_y - u_y h_x)/J, J = u_x v_y - u_y v_x."""
    R, x, y, w, one, u, v, J = slice_data()
    Jinv = R.scal(R.power(R.add(one, R.mul(x, y)), 3), -(P_INV2))      # 1/J = -(1 + xy)^3 / 2
    assert R.mul(J, Jinv) == one
    du = lambda h: R.mul(R.add(R.mul(R.dx(h), R.dy(v)), R.mul(R.dy(h), R.dx(v)), -1), Jinv)
    dv = lambda h: R.mul(R.add(R.mul(R.dx(u), R.dy(h)), R.mul(R.dy(u), R.dx(h)), -1), Jinv)
    assert du(u) == one and dv(v) == one and not du(v) and not dv(u)
    U = R.add(u, R.scal(one, 2), -1)
    V = R.add(v, R.scal(one, 3), -1)
    E = twisted(R, R.mul(U, V), 34)
    bad = tot = 0
    for name, h in (("w", w), ("x", x), ("y", y)):
        for m in range(0, 3):
            for n in range(0, 3):
                if m + n == 0 or m + n > 3:
                    continue
                d = h
                for _ in range(m):
                    d = du(d)
                for _ in range(n):
                    d = dv(d)
                lhs = R.mul(d, J)
                rhs = R.mul(R.mul(h, hermite(R, U, V, n, m, one)), J)
                diff = R.add(lhs, rhs, -1)
                wt = max(Ring.weight(mm) for mm in diff) if diff else 0
                if wt > 24:
                    continue
                tot += 1
                bad += bool(E.normal_form(diff))
    print("slice: [d_u^m d_v^n h] = [h H_(n,m)(U, V)] for h = w, x, y and 1 <= m + n <= 3 (weight <= 24): %d checked, %d failures" % (tot, bad))


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "identity":
    from stein import P
    P_INV2 = pow(2, P - 2, P)
    identity_check()
