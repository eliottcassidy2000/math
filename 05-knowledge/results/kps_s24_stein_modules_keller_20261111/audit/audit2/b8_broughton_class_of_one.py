"""b8_broughton_class_of_one.py - second audit, README section 8: W = x + x^2 y (Broughton), H^2 one-dimensional,
class of 1 not zero.  Classes of small monomials as multiples of the class of 1 (exact rationals recovered from the
residues modulo the prime), and the same for W = x^2 y (class of 1 zero, basis x)."""
from fractions import Fraction
import eng2
from eng2 import RingL, add, scal, Echelon, inv

R = RingL()
P = eng2.PRIME


def rat(v):
    """rational reconstruction of a residue modulo P (small numerator and denominator)."""
    for d in range(1, 2000):
        n = v * d % P
        if n > P // 2:
            n -= P
        if abs(n) < 10 ** 6:
            return Fraction(n, d)
    return None


def run(Wd, label, Np=30, N=10):
    Wx, Wy = R.dx(Wd), R.dy(Wd)
    E = Echelon(lambda m: (m[0] + m[1], m))
    for i in range(Np + 1):
        for j in range(Np + 1 - i):
            m = R.mono(i, j)
            E.insert(add(R.dx(m), R.mul(Wx, m), -1))
            E.insert(add(R.dy(m), R.mul(Wy, m), -1))
    dimV = (N + 1) * (N + 2) // 2
    cod = dimV - sum(1 for m in E.piv if m[0] + m[1] <= N)
    print("%s: codimension in degree <= %d: %d" % (label, N, cod))
    out = []
    for (i, j) in ((0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (2, 1), (0, 2), (3, 1), (1, 2)):
        nf, lead = E.reduce(R.mono(i, j))
        out.append("x^%d y^%d -> %s" % (i, j, " + ".join("%s*x^%d y^%d" % (rat(c), m[0], m[1]) for m, c in sorted(nf.items())) if nf else "0"))
    print("   normal forms: " + ";  ".join(out))


if __name__ == "__main__":
    run(add(R.mono(1, 0), R.mono(2, 1)), "W = x + x^2 y")
    run(R.mono(2, 1), "W = x^2 y")
    run(R.mono(1, 1), "W = x y")
