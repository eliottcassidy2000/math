"""stein_special_centre.py - a special centre where SN and CV hold for a map that is not an isomorphism
(section 6.3 of the note; added after the first audit; kind-pasteur-2026-11-11-S24).

X = C^2 minus {1 + xy = 0}, F = the inclusion (u, v) = (x, y), centre (0, 0), W = xy.  By hand: the Stein kernel is
graded by the charge deg x - deg y (w = 1/(1 + xy) has charge 0), contains every element of nonzero charge, and
has codimension 2 in charge 0.  So the classes of d_u^m h and d_v^m h vanish for m larger than the charges of h.
The same map at the centre (2, 3) is the control: there the classes do not vanish."""
import sys
from stein import Ring, twisted, codim

R = Ring("C")
x, y, w, one = R.mono(1, 0), R.mono(0, 1), R.mono(0, 0, 1), R.mono(0, 0)


def run(a, b, Np, wmax):
    U = R.add(x, R.scal(one, a), -1)
    V = R.add(y, R.scal(one, b), -1)
    E = twisted(R, R.mul(U, V), Np)
    print("centre (%d, %d), N' = %d: c(N) = %s" % (a, b, Np, " ".join("%d:%d" % (N, codim(R, E, N)) for N in (4, 8, 12, 16) if N < Np)))
    print("   class of 1: %s; class of w: %s; class of w^2: %s" % tuple("NONZERO" if E.normal_form(h) else "zero" for h in (one, w, R.mul(w, w))))
    tests = (("w", w), ("w^2", R.mul(w, w)), ("x w", R.mul(x, w)), ("y w", R.mul(y, w)), ("x^2 w", R.mul(R.mul(x, x), w)))
    for name, h in tests:
        outs = []
        for lab, d in (("d_u", R.dx), ("d_v", R.dy)):
            res, g = [], h
            for m in range(1, 7):
                g = d(g)
                if not g:
                    res.append("=0")
                elif max(Ring.weight(mm) for mm in g) > wmax:
                    res.append("?")
                else:
                    res.append("N" if E.normal_form(g) else "0")
            outs.append("%s^m: %s" % (lab, " ".join(res)))
        print("   h = %-5s  %s   (0 = zero class, N = nonzero class, =0 zero polynomial, ? weight above %d)" % (name, ";  ".join(outs), wmax))
    # mixed derivatives of w (charge zero): finite chaos must fail
    res, g = [], w
    for m in range(1, 4):
        g = R.dy(R.dx(g))
        if max(Ring.weight(mm) for mm in g) > wmax:
            res.append("?")
        else:
            res.append("N" if E.normal_form(g) else "0")
    print("   h = w, classes of (d_u d_v)^m w, m = 1..3: %s" % " ".join(res))
    sys.stdout.flush()


if __name__ == "__main__":
    run(0, 0, 34, 26)
    run(2, 3, 34, 26)
