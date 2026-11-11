"""stein_centres.py - one-sided Stein nilpotence on the slice F_1 = 1 of THM-1300 at several centres (added after
the pre-audit checkpoint, for Theorem 5.2: SN must fail at centres of every open set).  For each centre (a, b):
dimension of the Stein space and the classes of w V^m, w U^m, x V^m, x U^m for m = 1, 2, 3."""
import sys
from stein import Ring, twisted, codim
from stein_more import slice_data

R, x, y, w, one, u, v, J = slice_data()
for (a, b) in ((2, 3), (5, -7), (1, 1), (-3, 4), (7, 2), (-1, -2)):
    U = R.add(u, R.scal(one, a), -1)
    V = R.add(v, R.scal(one, b), -1)
    E = twisted(R, R.mul(U, V), 30)
    out = []
    for name, h, Pm in (("w V^m", w, V), ("w U^m", w, U), ("x V^m", x, V), ("x U^m", x, U)):
        res = []
        for m in (1, 2, 3):
            el = R.mul(R.mul(h, R.power(Pm, m)), J)
            wt = max(Ring.weight(mm) for mm in el)
            res.append("?" if wt > 24 else ("N" if E.normal_form(el) else "0"))
        out.append("%s: %s" % (name, " ".join(res)))
    print("centre (%d, %d): c(24) = %d; %s   (N = nonzero class, 0 = zero, ? = weight above 24)" % (a, b, codim(R, E, 24), "; ".join(out)))
    sys.stdout.flush()
