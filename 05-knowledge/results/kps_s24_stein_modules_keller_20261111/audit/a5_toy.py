"""a5_toy.py - the double cover (u, v) = (x^2, y/x) on C^* x C with the independent engine tw.py.

 1. dim H^2 for W = (x^2 - 2)(y/x - 3): expected 2 (two weightings, two primes)
 2. classes of x^(-2) V^m (m <= 8), x^(-1) V^m, x U^m, x^(-2) U^m: zero / nonzero
 3. QUANTITATIVE cross-check with the residue functionals of Proposition 6.1 (a2_toy_residues.py):
       lambda_+-(1) = 1/2, lambda_+-(x) = +- c/2  =>  [h] = alpha [1] + beta [x] with
       alpha = lambda_+(h) + lambda_-(h),  beta = (lambda_+(h) - lambda_-(h)) / c,   c^2 = a = 2,
    so for h = x^(-2) V^m:  alpha = (-1)^m m! / 2^(m+1), beta = 0   (exact, compared modulo the prime)
 4. [U^i V^j] = i! [i = j] [1] for i + j <= 5
 5. other centres: (5, -7) -> 2;  (0, 3) [the missed line u = 0] -> expected 0 (the fibre x(y - 3x) = t is C^*, and
    H^1(X) -> H^1(Y) is onto), so 'finite chaos' holds trivially there for a map that is not an isomorphism
 6. the open immersion C^* x C -> C^2 (u, v) = (x, y): centre (2, 3) -> 1, centre (0, 3) -> 0
"""
from math import factorial
from tw import RingB, add, mul, scal, power, table, build, codim, maxwt, P1, P2

one = {(0, 0): 1}
x, xi, y = {(1, 0): 1}, {(-1, 0): 1}, {(0, 1): 1}


def run(R, p, a0, b0, Ns, Nps, verbose=True):
    pw = lambda el, k: power(R, el, k, one)
    u = pw(x, 2)
    v = mul(R, y, xi)
    U = add(u, scal(one, a0), -1)
    V = add(v, scal(one, b0), -1)
    W = mul(R, U, V)
    E = table(R, W, Ns, Nps, p, "toy, centre (%d, %d), weights (|i| -> %d, j -> %d), prime %d" % (a0, b0, R.wx, R.wy, p))
    return E, U, V, pw


for R in (RingB(1, 1), RingB(1, 2), RingB(2, 1)):
    for p in (P1, P2):
        E, U, V, pw = run(R, p, 2, 3, (4, 8, 12), (20, 30, 40))
    print("     non-pivot monomials of weight <= 12:", E.nonpivots_upto(12))

# classes, with the (1, 1) weighting and N' = 60 (weights of the test elements are at most 2 m + 2 <= 18)
R = RingB(1, 1)
for p in (P1, P2):
    E, U, V, pw = run(R, p, 2, 3, (4, 8, 12, 20), (60,))
    J = scal(one, 2)
    nf1 = E.nf(mul(R, one, J))
    nfx = E.nf(mul(R, x, J))
    supp = sorted(set(nf1) | set(nfx))
    print("     normal forms: 1 ->", nf1, "  x ->", nfx)
    assert len(supp) == 2

    def coords(h):
        """(alpha, beta) with [h] = alpha [1] + beta [x], or None if the class is not in the span."""
        nf = E.nf(mul(R, h, J))
        m0, m1 = supp
        # solve alpha * nf1 + beta * nfx = nf on the two support monomials
        a11, a12, a21, a22 = nf1.get(m0, 0), nfx.get(m0, 0), nf1.get(m1, 0), nfx.get(m1, 0)
        det = (a11 * a22 - a12 * a21) % p
        r0, r1 = nf.get(m0, 0), nf.get(m1, 0)
        al = (r0 * a22 - a12 * r1) * pow(det, p - 2, p) % p
        be = (a11 * r1 - a21 * r0) * pow(det, p - 2, p) % p
        extra = [m for m in nf if m not in supp]
        return (al, be) if not extra else None

    good = True
    out = []
    for m in range(0, 9):
        h = mul(R, pw(xi, 2), pw(V, m))
        cd = coords(h)
        target = (-1) ** m * factorial(m) * pow(2, -(m + 1), p) % p
        good = good and cd == (target, 0)
        out.append("m=%d:%s" % (m, "match" if cd == (target, 0) else "MISMATCH %s" % (cd,)))
    print("  3. [x^-2 V^m] = ((-1)^m m!/2^(m+1)) [1] + 0 [x]  (prediction of the residue formula), m <= 8:", " ".join(out), "->", "OK" if good else "FAIL")
    res = []
    for name, h0, Pm in (("x^-1 V^m", xi, V), ("x U^m", x, U), ("x^-2 U^m", pw(xi, 2), U)):
        r = []
        for m in range(1, 7):
            cd = coords(mul(R, h0, pw(Pm, m)))
            r.append("0" if cd == (0, 0) else "N")
        res.append("%s: %s" % (name, " ".join(r)))
    print("  2. classes for m = 1..6 (N nonzero, 0 zero):", ";  ".join(res))
    bad = 0
    for i in range(0, 6):
        for j in range(0, 6 - i):
            cd = coords(mul(R, pw(U, i), pw(V, j)))
            bad += cd != ((factorial(i) if i == j else 0) % p, 0)
    print("  4. [U^i V^j] = i! [i = j] [1], i + j <= 5: %d failures of 21" % bad)

# other centres
R = RingB(1, 1)
run(R, P1, 5, -7, (4, 8, 12), (30, 40))
run(R, P1, 0, 3, (4, 8, 12), (30, 40))
run(R, P2, 0, 3, (4, 8, 12), (30, 40))

# open immersion C^* x C -> C^2
for (a0, b0) in ((2, 3), (0, 3)):
    U = add(x, scal(one, a0), -1)
    V = add(y, scal(one, b0), -1)
    table(R, mul(R, U, V), (4, 8, 12), (30, 40), P1, "open immersion (x, y) on C^* x C, centre (%d, %d): expected %d" % (a0, b0, 1 if a0 else 0))
