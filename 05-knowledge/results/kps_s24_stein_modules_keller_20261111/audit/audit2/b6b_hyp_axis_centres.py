"""b6b_hyp_axis_centres.py - second audit, README 6.3: the same map (inclusion of C^2 minus {1 + xy = 0}) at centres on
the two lines a b = 0 other than (0, 0).  Dimension and the classes of d_u^m h, d_v^m h (relative to the truncation)."""
import sys
import eng2
from b6_prop63_special_centre import R, x, y, one, w, wt, build, codim
from eng2 import power

for (av, bv) in ((0, 3), (0, -7), (2, 0), (5, 0), (1, 1)):
    Np = 36
    E, U, V = build(av, bv, Np)
    cls = lambda h: ("?" if (h and max(wt(m) for m in h) > Np - 8) else ("N" if not E.contains(h) else "0")) if h else "=0"
    print("centre (%d, %d): c(N) = %s; class of 1: %s, of w: %s, of x: %s, of y: %s, of x w: %s, of y w: %s"
          % (av, bv, " ".join("%d:%d" % (N, codim(E, N)) for N in (4, 8, 12, 16, 20, 24)), cls(one), cls(w), cls(x), cls(y), cls(R.mul(x, w)), cls(R.mul(y, w))))
    for name, h in (("w", w), ("x w", R.mul(x, w)), ("y w", R.mul(y, w)), ("x^2 w", R.mul(power(R, x, 2), w)), ("y^2 w^2", R.mul(power(R, y, 2), R.mul(w, w)))):
        outs = []
        for lab, d in (("d_u", R.dx), ("d_v", R.dy)):
            res, g = [], h
            for m in range(1, 8):
                g = d(g)
                res.append(cls(g))
            outs.append("%s^m: %s" % (lab, " ".join(res)))
        print("   h = %-8s %s" % (name, ";  ".join(outs)))
    sys.stdout.flush()
