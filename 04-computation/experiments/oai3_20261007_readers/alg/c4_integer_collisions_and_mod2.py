"""C4: (i) integer collisions of THM-1300's map on [-B, B]^3;
(ii) F mod 2 over F_{2^k}: image sizes (dominance test) and fibre structure.
"""
import numpy as np
from collections import defaultdict
from fractions import Fraction as Fr

def F(P):
    x, y, z = P
    u = 1 + x*y
    return (u**3*z + y**2*u*(4 + 3*x*y), y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y), 2*x - 3*x**2*y - x**3*z)

# (i) integer collisions
B = 24
r = np.arange(-B, B + 1, dtype=object)
img = defaultdict(list)
for x in range(-B, B + 1):
    for y in range(-B, B + 1):
        for z in range(-B, B + 1):
            img[F((x, y, z))].append((x, y, z))
coll = {w: Ps for w, Ps in img.items() if len(Ps) >= 2}
print("(i) integer collisions with all preimages in [-%d,%d]^3: %d targets" % (B, B, len(coll)))
def size(Ps):
    return max(max(abs(c) for c in P) for P in Ps)
for w, Ps in sorted(coll.items(), key=lambda kv: (size(kv[1]), kv[0]))[:14]:
    print("   F =", w, "<-", Ps, "  (g = 0)" if w[2] == 0 else "  (g != 0)")
print("   targets with g != 0:", sum(1 for w in coll if w[2] != 0))
# the s = 2 family: w odd
for w in range(-9, 10, 2):
    P0 = (0, 2*w, Fr(-(63*w*w + 1), 4))
    P1 = (1, Fr(w - 3, 2), Fr(13 - 3*w, 2))
    P2 = (-1, Fr(w + 3, 2), Fr(13 + 3*w, 2))
    t = (Fr(w*w - 1, 4), 2*w, 0)
    assert F(P0) == F(P1) == F(P2) == t, w
    assert all(c.denominator == 1 for P in (P0, P1, P2) for c in map(Fr, P))
print("   s=2 family F(0,2w,-(63w^2+1)/4) = F(1,(w-3)/2,(13-3w)/2) = F(-1,(w+3)/2,(13+3w)/2) = ((w^2-1)/4,2w,0): verified, w odd in [-9,9]")
# the s = 1 family: v = 7 mod 8 gives an integer double collision
for v in range(-25, 26):
    if v % 8 == 7:
        u0 = Fr(v*v - 1, 16)
        Pa = (0, v, u0 - 4*v*v)
        Pb = (2, Fr(v - 3, 4), Fr(13 - 3*v, 8))
        assert F(Pa) == F(Pb) == (u0, v, 0)
        assert all(Fr(c).denominator == 1 for c in Pa + Pb)
print("   s=1 family F(0,v,(v^2-1)/16-4v^2) = F(2,(v-3)/4,(13-3v)/8) = ((v^2-1)/16, v, 0), v = 7 mod 8: verified")

# (ii) F mod 2 over F_{2^k}
IRRED = {1: 0b11, 2: 0b111, 3: 0b1011, 4: 0b10011, 5: 0b100101, 6: 0b1000011, 7: 0b10000011}
def gf_tables(k):
    q = 1 << k
    poly = IRRED[k]
    exp = [0] * (2 * q)
    log = [0] * q
    a = 1
    # find a generator: try g = 2 (x); for these primitive polynomials x is primitive
    for i in range(q - 1):
        exp[i] = a
        log[a] = i
        a <<= 1
        if a & q:
            a ^= poly
    for i in range(q - 1, 2 * q):
        exp[i] = exp[i - (q - 1)]
    return np.array(exp, dtype=np.int64), np.array(log, dtype=np.int64)

print("\n(ii) F mod 2 over F_q, q = 2^k: image size (dominant => ~ c q^3; surface => ~ c q^2)")
for k in range(1, 8):
    q = 1 << k
    if k == 1:
        exp = np.array([1, 1, 1], dtype=np.int64); log = np.array([0, 0], dtype=np.int64)
    else:
        exp, log = gf_tables(k)
    def mul(a, b):
        res = exp[(log[a] + log[b]) % (q - 1)] if q > 2 else (a & b)
        return np.where((a == 0) | (b == 0), 0, res)
    rr = np.arange(q, dtype=np.int64)
    X, Y, Z = (A.ravel() for A in np.meshgrid(rr, rr, rr, indexing='ij'))
    XY = mul(X, Y); U = 1 ^ XY
    U2 = mul(U, U); U3 = mul(U2, U); Y2 = mul(Y, Y); Y3 = mul(Y2, Y); X2 = mul(X, X); X3 = mul(X2, X)
    # mod 2: F1 = u^3 z + y^2 u (4+3xy) = u^3 z + x y^3 u ; F2 = y + x u^2 z + x^2 y^3 ; F3 = x^3 z + x^2 y
    f1 = mul(U3, Z) ^ mul(mul(X, Y3), U)
    f2 = Y ^ mul(mul(X, U2), Z) ^ mul(X2, Y3)
    f3 = mul(X3, Z) ^ mul(X2, Y)
    code = (f1 * q + f2) * q + f3
    uniq, counts = np.unique(code, return_counts=True)
    print("   q=%4d  image=%9d  image/q^3=%.4f  image/q^2=%.3f  max fibre=%d" %
          (q, len(uniq), len(uniq) / q**3, len(uniq) / q**2, counts.max()))
