"""Lie-sphere reading: oriented circles are points of the Lie quadric Q in P^4 (form of signature (3,2)); oriented tangency is
Lie-orthogonality; the oriented solutions for a sign class e are the two points of the line orthogonal to the three oriented
data circles C_i^(e_i), i.e. a line meeting Q in two points, whose discriminant is (mod squares) disc(Q) * Gram_e, with
Gram_e = det(<C_i^(e_i), C_j^(e_j)>).  Test numerically: Delta_e * Gram_e is a rational square for each class e (so
Delta_e == Gram_e mod squares up to a common constant), and prod_e Gram_e is a square (the (Z/2)^3 relation).
Lie coordinates of the oriented circle (a, b, r): v = (a, b, r, (a^2+b^2-r^2-1)/2, (a^2+b^2-r^2+1)/2) with the form
<v,w> = a a' + b b' - r r' - (v4 w5 + v5 w4) ... we use the standard pentaspherical form: <v, w> = ( |c-c'|^2 - (r-r')^2 )/2
for oriented circles with signed radii, which vanishes iff oriented tangency (same orientation, r r' > 0 external...)."""
from fractions import Fraction as F
from math import isqrt
import random, itertools
random.seed(23)
def is_sq(fr):
    if fr == 0: return True
    if fr < 0: return False
    n, d = fr.numerator, fr.denominator
    return isqrt(n)**2 == n and isqrt(d)**2 == d
def sqfree_sign(fr):
    # crude square class: sign and squarefree part of numerator*denominator (small numbers only)
    if fr == 0: return 0
    n = abs(fr.numerator * fr.denominator); s = 1 if fr > 0 else -1
    out = 1; p = 2
    while p*p <= n:
        c = 0
        while n % p == 0: n //= p; c += 1
        if c % 2: out *= p
        p += 1
    return s * out * n
def discs(circ):
    out = []
    (a1,b1,r1),(a2,b2,r2),(a3,b3,r3) = circ
    for e in [(1,1,1),(1,1,-1),(1,-1,1),(-1,1,1)]:
        rows = []
        for (a,b,rr), ei in ((circ[1], e[1]), (circ[2], e[2])):
            rows.append((F(2*(a-a1)), F(2*(b-b1)), F(2*(ei*rr - e[0]*r1)), F(a1*a1 + b1*b1 - r1*r1 - a*a - b*b + rr*rr)))
        (m11, m12, c1, d1), (m21, m22, c2, d2) = rows
        det = m11*m22 - m12*m21
        if det == 0: return None
        xr = (-(c1*m22) + c2*m12)/det; x0 = (-(d1*m22) + d2*m12)/det
        yr = (-(c2*m11) + c1*m21)/det; y0 = (-(d2*m11) + d1*m21)/det
        A = xr*xr + yr*yr - 1; B = 2*xr*(x0 - a1) + 2*yr*(y0 - b1) - 2*e[0]*r1; Cc = (x0 - a1)**2 + (y0 - b1)**2 - r1*r1
        out.append(B*B - 4*A*Cc)
    return out
def lie(v, w):
    # oriented circles v = (a, b, r_signed): <v,w> = (|c - c'|^2 - (r - r')^2)/2 ; zero iff oriented tangency
    return (F((v[0]-w[0])**2 + (v[1]-w[1])**2) - F((v[2]-w[2])**2))/2
def gram(circ, e):
    O = [(c[0], c[1], ei*c[2]) for c, ei in zip(circ, e)]
    M = [[lie(O[i], O[j]) for j in range(3)] for i in range(3)]
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
n = 0; ok_each = 0; ok_prod = 0; classes = [(1,1,1),(1,1,-1),(1,-1,1),(-1,1,1)]; ratio_classes = set()
while n < 300:
    circ = [(random.randint(-20, 20), random.randint(-20, 20), random.randint(1, 20)) for _ in range(3)]
    D = discs(circ)
    if D is None or any(d == 0 for d in D): continue
    G = [gram(circ, e) for e in classes]
    if any(g == 0 for g in G): continue
    n += 1
    # test Delta_e * Gram_e square for each e (up to a universal constant: record the square class of the ratio)
    cls = tuple(sqfree_sign(D[k]*G[k]) for k in range(4))
    ratio_classes.add(cls)
    if all(is_sq(abs(D[k]*G[k])) for k in range(4)): ok_each += 1
    if is_sq(G[0]*G[1]*G[2]*G[3]): ok_prod += 1
print(f"{n} triples: Delta_e * Gram_e a square (all four classes) in {ok_each}; square classes of Delta_e*Gram_e seen: {sorted(ratio_classes)[:6]} ({len(ratio_classes)} distinct)")
print(f"product of the four Gram determinants a square in {ok_prod} of {n}")
# Ford / Descartes: Ford circle at p/q has curvature 2q^2; with the line (curvature 0) and neighbours a/b, c/d (ad - bc = +-1):
from math import gcd
checks = 0
for q in range(2, 40):
    for p in range(1, q):
        if gcd(p,q) != 1: continue
        for b in range(1, q):
            a = (p*b + 1)//q if (p*b + 1) % q == 0 else None
            if a is None: continue
            d = q - b; c = p - a
            if c < 0 or d <= 0 or a*d - b*c != 1: continue
            k1, k2, k3 = 2*b*b, 2*d*d, 0
            # Descartes: k4 = k1 + k2 + k3 +- 2 sqrt(k1 k2 + k2 k3 + k3 k1)
            s = isqrt(k1*k2); assert s*s == k1*k2
            k4 = {k1 + k2 + 2*s, k1 + k2 - 2*s}
            assert 2*q*q in k4 and 2*(b-d)**2 in k4, (p, q, b, d)
            checks += 1; break
print(f"Ford circles: Descartes' rule with the line gives exactly the mediant (b+d)^2 and the parent (b-d)^2 curvatures in all {checks} cases, q < 40")
