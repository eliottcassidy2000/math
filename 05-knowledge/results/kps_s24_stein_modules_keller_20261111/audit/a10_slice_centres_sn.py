"""a10_slice_centres_sn.py - independent replication of the (uncommitted) run runs/sn_centres.txt that accompanies the
draft Theorem 5.2: on the slice, at six centres, the classes of w V^m, w U^m, x V^m, x U^m for m = 1..4.
Also the excluded centre (2, 0) and the centre (0, 3), and a centre ON the non-properness curve S_1: (4, 0), where the
fibre has one point.
"""
from tw import RingC, add, mul, scal, power, build, codim, maxwt, P1

one = {(0, 0): 1}
x, y, w = {(1, 0): 1}, {(-1, 0): 1}, {(0, -1): 1}
R = RingC(1, 2, 1)
pw = lambda el, k: power(R, el, k, one)
u = add(y, scal(mul(R, x, w), 3))
v = add(add(mul(R, x, w), mul(R, x, pw(w, 2))), mul(R, pw(x, 3), pw(w, 3)), -1)
J = scal(pw(w, 3), -2)
Np = 46
for (a0, b0) in ((2, 3), (5, -7), (1, 1), (-3, 4), (7, 2), (-1, -2), (0, 3), (2, 0), (0, 0), (4, 0)):
    U = add(u, scal(one, a0), -1)
    V = add(v, scal(one, b0), -1)
    E = build(R, mul(R, U, V), Np, P1)
    out = []
    for name, h, Pm in (("w V^m", w, V), ("w U^m", w, U), ("x V^m", x, V), ("x U^m", x, U)):
        res = []
        for m in (1, 2, 3, 4):
            el = mul(R, mul(R, h, pw(Pm, m)), J)
            res.append("?" if maxwt(R, el) > Np - 10 else ("N" if E.nf(el) else "0"))
        out.append("%s: %s" % (name, " ".join(res)))
    one_cls = "nonzero" if E.nf(dict(J)) else "ZERO"
    print("centre (%2d, %2d): c(4..36) = %s; [1] %s; %s" % (a0, b0, sorted(set(codim(R, E, N) for N in range(4, 37))), one_cls, "; ".join(out)), flush=True)
