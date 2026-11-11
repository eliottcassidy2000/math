"""a7_slice_stability.py - how far the truncation of a6_slice.py can be trusted, and more values of m.

 1. Normal forms of the same elements at N' = 46, 58, 70 (prime P1): if M(N') cap V_N had not yet reached M cap V_N
    for the weights involved, the coordinates on the 11 non-pivot monomials would change with N'.  They are compared
    as exact vectors modulo the prime.
 2. c(N, N') for all N <= N' - 10 at N' = 70 (the plateau 11).
 3. Classes of w V^m (m <= 10), w U^m (m <= 12), x V^m, x U^m at N' = 70, and the rank of the span of
    {[w V^m]}, {[w U^m]} in the 11-dimensional quotient.

A nonzero normal form proves that the element is not in M(N'); it proves that it is not in M only if
M(N') cap V_N = M cap V_N for N = the weight of the element (see AUDIT.md, 'truncation').
"""
import sys, time
from tw import RingC, add, mul, scal, power, deriv, build, codim, maxwt, P1, P2

one = {(0, 0): 1}
x, y, w = {(1, 0): 1}, {(-1, 0): 1}, {(0, -1): 1}
R = RingC(1, 2, 1)
pw = lambda el, k: power(R, el, k, one)
u = add(y, scal(mul(R, x, w), 3))
v = add(add(mul(R, x, w), mul(R, x, pw(w, 2))), mul(R, pw(x, 3), pw(w, 3)), -1)
J = scal(pw(w, 3), -2)
U = add(u, scal(one, 2), -1)
V = add(v, scal(one, 3), -1)
W = mul(R, U, V)
p = P1

tests = {}
for m in range(1, 13):
    tests["w V^%d" % m] = mul(R, mul(R, w, pw(V, m)), J)
    tests["w U^%d" % m] = mul(R, mul(R, w, pw(U, m)), J)
    tests["x V^%d" % m] = mul(R, mul(R, x, pw(V, m)), J)
    tests["x U^%d" % m] = mul(R, mul(R, x, pw(U, m)), J)
tests["1"] = dict(J)
wts = {k: maxwt(R, el) for k, el in tests.items()}

Nps = [int(a) for a in sys.argv[1:]] or [46, 58, 70]
nfs = {}
for Np in Nps:
    t0 = time.time()
    E = build(R, W, Np, p)
    plateau = sorted(set(codim(R, E, N) for N in range(4, Np - 9)))
    print("N' = %d: %d pivots (%.0f s); values of c(N, N') for 4 <= N <= %d: %s" % (Np, len(E.piv), time.time() - t0, Np - 10, plateau), flush=True)
    nfs[Np] = {k: E.nf(el) for k, el in tests.items() if wts[k] <= Np - 10}
    lastE = E

# 1. stability of the normal forms
common = [k for k in tests if all(k in nfs[Np] for Np in Nps)]
diff = [k for k in common if any(nfs[Np][k] != nfs[Nps[0]][k] for Np in Nps[1:])]
print("1. normal forms identical at N' = %s for %d elements (weights up to %d): %s" % (
    Nps, len(common), max(wts[k] for k in common), "YES" if not diff else "NO: %s" % diff))

# 3. classes at the largest N'
big = nfs[Nps[-1]]
basis = sorted(set(m for nf in big.values() for m in nf))
print("   support of all normal forms: %d monomials %s" % (len(basis), basis))


def rank(vecs):
    rows = [[nf.get(m, 0) for m in basis] for nf in vecs]
    r = 0
    for col in range(len(basis)):
        piv = next((i for i in range(r, len(rows)) if rows[i][col] % p), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][col], p - 2, p)
        rows[r] = [a * inv % p for a in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][col] % p:
                f = rows[i][col]
                rows[i] = [(a - f * b) % p for a, b in zip(rows[i], rows[r])]
        r += 1
    return r


for fam in ("w V", "w U", "x V", "x U"):
    ks = [k for k in big if k.startswith(fam)]
    ks.sort(key=lambda s: int(s.split("^")[1]))
    print("3. %s^m at N' = %d: %s ; rank of the span: %d" % (
        fam, Nps[-1], " ".join("m=%s:%s" % (k.split("^")[1], "N" if big[k] else "0") for k in ks), rank([big[k] for k in ks])))
print("   rank of the span of all tested classes together with [1]:", rank(list(big.values())))
