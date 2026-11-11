"""a11_slice_centre00.py - the excluded centre (0, 0) of the slice (dimension 9): which classes h V^m, h U^m vanish?
At (0, 0) the potential u v is invariant under (x, y) -> (-x, -y) (u, v are odd), so the Stein space splits into an
even and an odd part and classes can vanish by parity; the table shows what else vanishes."""
import sys
from tw import RingC, add, mul, scal, power, build, codim, maxwt, P1, P2

one = {(0, 0): 1}
x, y, w = {(1, 0): 1}, {(-1, 0): 1}, {(0, -1): 1}
R = RingC(1, 2, 1)
pw = lambda el, k: power(R, el, k, one)
u = add(y, scal(mul(R, x, w), 3))
v = add(add(mul(R, x, w), mul(R, x, pw(w, 2))), mul(R, pw(x, 3), pw(w, 3)), -1)
J = scal(pw(w, 3), -2)
Np = int(sys.argv[1]) if len(sys.argv) > 1 else 70
for (a0, b0) in ((0, 0), (2, 0)):
    U = add(u, scal(one, a0), -1)
    V = add(v, scal(one, b0), -1)
    for p in (P1, P2):
        E = build(R, mul(R, U, V), Np, p)
        print("centre (%d, %d), prime %d, N' = %d: c(N) for 4 <= N <= %d: %s" % (a0, b0, p, Np, Np - 10, sorted(set(codim(R, E, N) for N in range(4, Np - 9)))))
        for name, h, Pm in (("w V^m", w, V), ("w U^m", w, U), ("x V^m", x, V), ("x U^m", x, U), ("y V^m", y, V), ("y U^m", y, U), ("w^2 V^m", pw(w, 2), V), ("x w V^m", mul(R, x, w), V)):
            res = []
            for m in range(1, 13):
                el = mul(R, mul(R, h, pw(Pm, m)), J)
                if maxwt(R, el) > Np - 10:
                    break
                res.append("N" if E.nf(el) else "0")
            print("   %-8s m = 1..%-2d: %s" % (name, len(res), " ".join(res)), flush=True)
