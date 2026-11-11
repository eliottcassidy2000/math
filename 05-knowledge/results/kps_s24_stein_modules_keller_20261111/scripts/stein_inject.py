"""stein_inject.py - sanity checks for section 5 of the note (added after the first audit; kind-pasteur-2026-11-11-S24).

Theorem 5.3 says: a nonzero element h of the ring lies in the Stein kernel M_c only for centres c in a countable
union of proper algebraic subsets.  The arithmetic is modulo the prime P, so "a centre in general position" is
emulated by pseudo-random residues; this is a sanity check, not a proof.

  A. the operators AA = u - d_v, BB = v - d_u commute (slice), and M_c = (AA - a) B + (BB - b) B;
  B. the dual fields of the slice are not locally nilpotent (d_u w, d_v w != 0; iterates on x, y, w);
  C. slice: classes of all monomials of weight <= 4 at pseudo-random centres and at (2, 3), (0, 0);
  D. toy on C^* x C: the same;
  E. standard pair on C^2: the same (at the centre (0, 0) the class of x^i y^j vanishes iff i != j).
"""
import random
import sys
from stein import Ring, twisted, P
from stein_more import slice_data

P_INV2 = pow(2, P - 2, P)
random.seed(20261111)


def classes(R, E, J, elems, wmax):
    """number of elements (of weight <= wmax after multiplication by J) with zero class, and the total tested."""
    zero, tot, zs = 0, 0, []
    for name, h in elems:
        el = R.mul(h, J)
        if not el:
            continue
        if max(Ring.weight(mm) for mm in el) > wmax:
            continue
        tot += 1
        if not E.normal_form(el):
            zero += 1
            zs.append(name)
    return zero, tot, zs


def run_slice():
    R, x, y, w, one, u, v, J = slice_data()
    Jinv = R.scal(R.power(R.add(one, R.mul(x, y)), 3), -P_INV2)
    assert R.mul(J, Jinv) == one
    du = lambda h: R.mul(R.add(R.mul(R.dx(h), R.dy(v)), R.mul(R.dy(h), R.dx(v)), -1), Jinv)
    dv = lambda h: R.mul(R.add(R.mul(R.dx(u), R.dy(h)), R.mul(R.dy(u), R.dx(h)), -1), Jinv)
    assert du(u) == one and dv(v) == one and not du(v) and not dv(u)
    AA = lambda h: R.add(R.mul(u, h), dv(h), -1)
    BB = lambda h: R.add(R.mul(v, h), du(h), -1)
    bad = 0
    tests = [x, y, w, R.mul(x, y), R.mul(w, w), R.mul(x, w), R.add(R.power(x, 2), R.mul(y, w))]
    for h in tests:
        bad += bool(R.add(AA(BB(h)), BB(AA(h)), -1))
    print("A. slice: (u - d_v)(v - d_u) h = (v - d_u)(u - d_v) h on %d test elements: %d failures" % (len(tests), bad))
    print("B. slice: d_u w %s 0, d_v w %s 0 (a locally nilpotent derivation of a domain kills the units)" % ("!=" if du(w) else "=", "!=" if dv(w) else "="))
    for name, h in (("x", x), ("y", y), ("w", w)):
        ru, rv = [], []
        a, b = h, h
        for m in range(1, 5):
            a, b = du(a), dv(b)
            ru.append("N" if a else "0")
            rv.append("N" if b else "0")
        print("   d_u^m %s, m = 1..4: %s;   d_v^m %s: %s   (N = nonzero polynomial)" % (name, " ".join(ru), name, " ".join(rv)))
    sys.stdout.flush()
    mons = [("x^%d y^%d w^%d" % m, R.mono(*m)) for m in R.basis(4)]
    extra = [("d_u x", du(x)), ("d_u y", du(y)), ("d_u w", du(w)), ("d_v x", dv(x)), ("d_v y", dv(y)), ("d_v w", dv(w))]
    cents = [(random.randrange(P), random.randrange(P)) for _ in range(3)] + [(2, 3), (0, 0)]
    for (a, b) in cents:
        U = R.add(u, R.scal(one, a), -1)
        V = R.add(v, R.scal(one, b), -1)
        E = twisted(R, R.mul(U, V), 30)
        z, t, zs = classes(R, E, J, mons, 24)
        z2, t2, zs2 = classes(R, E, J, extra, 24)
        print("C. slice, centre (%d, %d): monomials of weight <= 4 with zero class: %d of %d %s; first dual derivatives of x, y, w with zero class: %d of %d %s"
              % (a, b, z, t, zs if zs else "", z2, t2, zs2 if zs2 else ""))
        sys.stdout.flush()


def run_toy():
    R = Ring("B")
    x, y, xi, one = R.mono(1, 0), R.mono(0, 1), R.mono(-1, 0), R.mono(0, 0)
    u, v = R.mul(x, x), R.mul(y, xi)
    J = R.scal(one, 2)
    mons = [("x^%d y^%d" % (k, l), R.mono(k, l)) for k in range(-4, 5) for l in range(0, 4)]
    cents = [(random.randrange(P), random.randrange(P)) for _ in range(3)] + [(2, 3), (0, 3)]
    for (a, b) in cents:
        U = R.add(u, R.scal(one, a), -1)
        V = R.add(v, R.scal(one, b), -1)
        E = twisted(R, R.mul(U, V), 24)
        z, t, zs = classes(R, E, J, mons, 14)
        print("D. toy, centre (%d, %d): monomials x^k y^l (|k| <= 4, l <= 3) with zero class: %d of %d %s" % (a, b, z, t, ("(all)" if z == t else zs) if zs else ""))
        sys.stdout.flush()


def run_std():
    R = Ring("A")
    x, y, one = R.mono(1, 0), R.mono(0, 1), R.mono(0, 0)
    mons = [("x^%d y^%d" % (i, j), R.mono(i, j)) for i in range(0, 6) for j in range(0, 6 - i)]
    cents = [(random.randrange(P), random.randrange(P)) for _ in range(3)] + [(2, 3), (0, 0)]
    for (a, b) in cents:
        U = R.add(x, R.scal(one, a), -1)
        V = R.add(y, R.scal(one, b), -1)
        E = twisted(R, R.mul(U, V), 20)
        z, t, zs = classes(R, E, one, mons, 12)
        print("E. standard pair, centre (%d, %d): monomials of degree <= 5 with zero class: %d of %d %s" % (a, b, z, t, zs if zs else ""))
        sys.stdout.flush()


if __name__ == "__main__":
    run_std()
    run_toy()
    run_slice()
