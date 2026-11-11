"""a6_slice.py - the slice F_1 = 1 of THM-1300 with the independent engine tw.py (torus-weight basis, two primes).

usage: python a6_slice.py dims [N' ...]        dimension tables: centre (2,3) with three filtrations; other centres;
                                               linear potentials u, v, u + v
       python a6_slice.py classes N' [mmax]    Stein classes at the centre (2, 3): [U^i V^j], [1], w V^m, w U^m, x V^m, x U^m
       python a6_slice.py timing N'

In the torus basis m(e, n):  w = m(0,-1), x = m(1,0), y = m(-1,0),
   u = y + 3 x w = m(-1,0) + 3 m(1,-1),   v = x w + x w^2 - x^3 w^3 = m(1,-1) + m(1,-2) - m(3,-3),   J = -2 w^3.
The class of h in the Stein space S_F is the class of h J modulo M' = Im(d/dx - W_x) + Im(d/dy - W_y).

The filtration (we, wpos, wneg) = (1, 2, 1) spans the same subspaces V_N as the author's weight i + j + k (a
polynomial of degree <= N, plus x^i w^k and y^j w^k with i + k, j + k <= N), so c(N, N') must agree with the author's
tables number by number; the filtrations (1, 2, 2) and (2, 3, 1) are different exhaustions of the ring.
"""
import sys, time
from math import factorial
from tw import RingC, add, mul, scal, power, deriv, table, build, codim, maxwt, P1, P2

one = {(0, 0): 1}
x, y, w = {(1, 0): 1}, {(-1, 0): 1}, {(0, -1): 1}


def data(R):
    pw = lambda el, k: power(R, el, k, one)
    u = add(y, scal(mul(R, x, w), 3))
    v = add(add(mul(R, x, w), mul(R, x, pw(w, 2))), mul(R, pw(x, 3), pw(w, 3)), -1)
    J = scal(pw(w, 3), -2)
    # Jacobian check in this basis
    Jc = add(mul(R, deriv(R, u, "x"), deriv(R, v, "y")), mul(R, deriv(R, u, "y"), deriv(R, v, "x")), -1)
    assert Jc == J, "Jacobian mismatch"
    return u, v, J, pw


def potential(R, a0, b0):
    u, v, J, pw = data(R)
    U = add(u, scal(one, a0), -1)
    V = add(v, scal(one, b0), -1)
    return U, V, mul(R, U, V), J, pw


mode = sys.argv[1]
if mode == "timing":
    Np = int(sys.argv[2])
    R = RingC(1, 2, 1)
    U, V, W, J, pw = potential(R, 2, 3)
    t0 = time.time()
    E = build(R, W, Np, P1)
    print("N' = %d: %d pivots, %.1f s; c(4) = %d, c(12) = %d" % (Np, len(E.piv), time.time() - t0, codim(R, E, 4), codim(R, E, 12)))

elif mode == "dims":
    Nps = [int(a) for a in sys.argv[2:]] or [28, 34]
    R = RingC(1, 2, 1)
    U, V, W, J, pw = potential(R, 2, 3)
    for p in (P1, P2):
        E = table(R, W, (4, 8, 12, 16, 20, 24, 28, 34, 36), Nps, p, "slice, W = (u-2)(v-3), filtration (1,2,1) [= the author's], prime %d; predicted 11" % p)
    print("     non-pivot basis elements m(e, n) of weight <= %d:" % min(Nps[-1], 20), E.nonpivots_upto(min(Nps[-1], 20)))
    for wts, Ns, Nq in (((1, 2, 2), (4, 8, 12, 16, 20, 24), [n + 6 for n in Nps]), ((2, 3, 1), (6, 12, 18, 24, 30), [n + 14 for n in Nps])):
        R2 = RingC(*wts)
        U2, V2, W2, J2, pw2 = potential(R2, 2, 3)
        table(R2, W2, Ns, Nq, P1, "slice, W = (u-2)(v-3), filtration %s, prime %d; predicted 11" % (wts, P1))
    Nc = Nps[0]
    for (a0, b0, pred) in ((5, -7, 11), (1, 1, 11), (-3, 4, 11), (7, 2, 11), (-1, -2, 11), (0, 0, 9), (2, 0, 9), (0, 3, 11)):
        Uc, Vc, Wc, Jc, _ = potential(R, a0, b0)
        table(R, Wc, (4, 8, 12, 16), (Nc,), P1, "slice, centre (%d, %d); predicted %d" % (a0, b0, pred))
    u, v, J, pw = data(R)
    for name, Wl, pred in (("u", u, 2), ("v", v, 4), ("u + v", add(u, v), 6), ("u - 2v", add(u, scal(v, -2)), 6), ("2u + 3", add(scal(u, 2), scal(one, 3)), 2)):
        for p in (P1, P2):
            table(R, Wl, (4, 8, 12, 16), (Nc,), p, "slice, W = %s, prime %d; predicted %d" % (name, p, pred))

elif mode == "classes":
    Np = int(sys.argv[2])
    mmax = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    R = RingC(1, 2, 1)
    U, V, W, J, pw = potential(R, 2, 3)
    for p in (P1, P2):
        t0 = time.time()
        E = build(R, W, Np, p)
        Ntop = Np - 10
        print("prime %d, N' = %d: %d pivots (%.0f s); c(N) for N = 4, 12, 20, ..., %d: %s" % (
            p, Np, len(E.piv), time.time() - t0, Ntop, " ".join(str(codim(R, E, N)) for N in list(range(4, Ntop, 8)) + [Ntop])))
        nf1 = E.nf(mul(R, one, J))
        print("   class of 1: %s (normal form supported on %s)" % ("NONZERO" if nf1 else "zero", sorted(nf1)))
        bad = tot = 0
        for i in range(0, 5):
            for j in range(0, 5 - i):
                el = mul(R, pw(U, i), pw(V, j))
                if i == j:
                    el = add(el, scal(one, factorial(i)), -1)
                tot += 1
                bad += bool(E.nf(mul(R, el, J)))
        print("   [U^i V^j] - i! [i = j] [1] = 0 for i + j <= 4: %d of %d" % (tot - bad, tot))
        for name, h, Pm, mm in (("w V^m", w, V, mmax), ("w U^m", w, U, mmax), ("x V^m", x, V, mmax - 1), ("x U^m", x, U, mmax - 1)):
            res = []
            for m in range(1, mm + 1):
                el = mul(R, mul(R, h, pw(Pm, m)), J)
                wt_ = maxwt(R, el)
                if wt_ > Np - 6:
                    res.append("m=%d:?(wt %d)" % (m, wt_))
                else:
                    nf = E.nf(el)
                    res.append("m=%d:%s(wt %d)" % (m, "NONZERO" if nf else "zero", wt_))
            print("   classes of %s: %s" % (name, " ".join(res)), flush=True)
