"""stein_runs.py - the pre-registered runs S2-S5 (kind-pasteur-2026-11-11-S24).

usage: python stein_runs.py toy [N' ...]       ring B, (u, v) = (x^2, y/x), centre (2, 3)
       python stein_runs.py slice [N' ...]     ring C, slice F_1 = 1 of THM-1300, centre (2, 3)

In the frame du ^ dv the Stein kernel of (u, v) is M_F; in the frame dx ^ dy it is
M' = Im(d/dx - W_x) + Im(d/dy - W_y), and h du ^ dv = h J dx ^ dy.  So the class of h in the Stein space is the
class of h J modulo M'."""
import sys, time
from stein import Ring, twisted, codim, Echelon, P


def report(R, W, J, u0, v0, tests, Ns, Nps, label):
    print(label)
    E = None
    for Np in Nps:
        t0 = time.time()
        E = twisted(R, W, Np)
        print("   N' = %2d: %s   (%.1f s, %d pivots)" % (Np, " ".join("c(%d)=%d" % (N, codim(R, E, N)) for N in Ns if N <= Np), time.time() - t0, len(E.piv)))
        sys.stdout.flush()
    Np = Nps[-1]
    one = R.mono(0, 0)
    # Gaussian subalgebra: U^k V^l = k! [k = l]
    bad = 0
    tot = 0
    fact = [1, 1, 2, 6, 24, 120]
    for k in range(0, 5):
        for l in range(0, 5 - k):
            el = R.mul(R.power(u0, k), R.power(v0, l))
            if k == l:
                el = R.add(el, R.scal(one, fact[k]), -1)
            el = R.mul(el, J)
            wt = max(Ring.weight(m) for m in el) if el else 0
            nf = E.normal_form(el)
            tot += 1
            if nf:
                bad += 1
    print("   U^k V^l - k! [k = l] in the Stein kernel for k + l <= 4: %d of %d (at N' = %d)" % (tot - bad, tot, Np))
    nf1 = E.normal_form(R.mul(one, J))
    print("   class of 1: %s" % ("NONZERO" if nf1 else "zero"))
    for name, h, Pm, mmax in tests:
        res = []
        for m in range(1, mmax + 1):
            el = R.mul(R.mul(h, R.power(Pm, m)), J)
            wt = max(Ring.weight(mm) for mm in el)
            nf = E.normal_form(el)
            res.append("m=%d:%s(wt %d)" % (m, "NONZERO" if nf else "zero", wt))
        print("   classes of %s: %s" % (name, " ".join(res)))
    sys.stdout.flush()
    return E


if __name__ == "__main__":
    mode = sys.argv[1]
    Nps = [int(t) for t in sys.argv[2:]]
    if mode == "toy":
        R = Ring("B")
        x, y, xi, one = R.mono(1, 0), R.mono(0, 1), R.mono(-1, 0), R.mono(0, 0)
        u = R.mul(x, x)
        v = R.mul(y, xi)
        u0 = R.add(u, R.scal(one, 2), -1)
        v0 = R.add(v, R.scal(one, 3), -1)
        W = R.mul(u0, v0)
        J = R.scal(one, 2)
        tests = [("x^(-2) (v - 3)^m", R.mono(-2, 0), v0, 6), ("x^(-1) (v - 3)^m", xi, v0, 6), ("x (u - 2)^m", x, u0, 6), ("x^(-2) (u - 2)^m", R.mono(-2, 0), u0, 6)]
        report(R, W, J, u0, v0, tests, (4, 8, 12), Nps or [12, 18, 24, 30], "S2. toy: C^* x C, (u, v) = (x^2, y/x), W = (u - 2)(v - 3); predicted dimension 2")
    else:
        R = Ring("C")
        x, y, w, one = R.mono(1, 0), R.mono(0, 1), R.mono(0, 0, 1), R.mono(0, 0)
        u = R.add(y, R.scal(R.mul(x, w), 3))
        v = R.add(R.add(R.mul(x, w), R.mul(x, R.mul(w, w))), R.mul(R.power(x, 3), R.power(w, 3)), -1)
        # Jacobian check: u_x v_y - u_y v_x = -2 w^3
        Jc = R.add(R.mul(R.dx(u), R.dy(v)), R.mul(R.dy(u), R.dx(v)), -1)
        J = R.scal(R.power(w, 3), -2)
        assert Jc == J, "Jacobian mismatch"
        u0 = R.add(u, R.scal(one, 2), -1)
        v0 = R.add(v, R.scal(one, 3), -1)
        W = R.mul(u0, v0)
        tests = [("w (v - 3)^m", w, v0, 5), ("w (u - 2)^m", w, u0, 5), ("x (v - 3)^m", x, v0, 4), ("x (u - 2)^m", x, u0, 4)]
        report(R, W, J, u0, v0, tests, (4, 12, 20, 28, 34, 36), Nps or [16, 22, 28], "S3-S5. slice F_1 = 1 of THM-1300: X = C^2 minus {xy = -1}, W = (u - 2)(v - 3); predicted dimension 11")
