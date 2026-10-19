"""Apollonius over F_q: for random generic triples of circles (a_i, b_i, r_i), count the circles (x, y, r) with r != 0 tangent to
all three in the sense (x-a_i)^2 + (y-b_i)^2 = (r + e_i r_i)^2 for some signs e_i (solutions (x,y,r) and (x,y,-r) with all
signs flipped are the same circle; we count unordered circles = solutions with r taken up to sign).  Generic over C: 8.
Chebotarev: the distribution of the F_q-count over random triples reflects the monodromy group acting on the 8 solutions."""
import random, sys
from collections import Counter
qs = [int(v) for v in sys.argv[1:]] or [11, 13, 17, 19, 23]
random.seed(7)
def count(q, circ):
    sq = {}
    for t in range(q): sq.setdefault(t*t % q, []).append(t)
    # for each sign class (4), solve: two linear equations in x,y,r from differences, then a quadratic in r.
    sols = set()
    (a1,b1,r1),(a2,b2,r2),(a3,b3,r3) = circ
    for e in [(1,1,1),(1,1,-1),(1,-1,1),(-1,1,1)]:
        # equation i: x^2 - 2 a_i x + a_i^2 + y^2 - 2 b_i y + b_i^2 - r^2 - 2 e_i r_i r - r_i^2 = 0
        # difference (1)-(i): 2(a_i - a_1) x + 2(b_i - b_1) y + 2(e_i r_i - e_1 r_1) r + (a_1^2 + b_1^2 - r_1^2 - a_i^2 - b_i^2 + r_i^2) = 0
        rows = []
        for (a,b,rr), ei in ((circ[1], e[1]), (circ[2], e[2])):
            rows.append(((2*(a-a1)) % q, (2*(b-b1)) % q, (2*(ei*rr - e[0]*r1)) % q, (a1*a1 + b1*b1 - r1*r1 - a*a - b*b + rr*rr) % q))
        # solve for x, y as affine functions of r: M [x y]^T = -(c r + d)
        (m11, m12, c1, d1), (m21, m22, c2, d2) = rows
        det = (m11*m22 - m12*m21) % q
        if det == 0: return None     # degenerate configuration; skip
        inv = pow(det, -1, q)
        # x = (-(c1 r + d1) m22 + (c2 r + d2) m12) / det ; y = (-(c2 r + d2) m11 + (c1 r + d1) m21) / det
        xr = (-(c1*m22) + c2*m12) * inv % q; x0 = (-(d1*m22) + d2*m12) * inv % q
        yr = (-(c2*m11) + c1*m21) * inv % q; y0 = (-(d2*m11) + d1*m21) * inv % q
        # substitute into equation 1: (x - a1)^2 + (y - b1)^2 - (r + e1 r1)^2 = 0, quadratic in r
        A = (xr*xr + yr*yr - 1) % q
        B = (2*xr*(x0 - a1) + 2*yr*(y0 - b1) - 2*e[0]*r1) % q
        Cc = ((x0 - a1)**2 + (y0 - b1)**2 - r1*r1) % q
        if A == 0:
            if B == 0: return None
            rs = [(-Cc * pow(B, -1, q)) % q]
        else:
            disc = (B*B - 4*A*Cc) % q
            if disc not in sq: rs = []
            else: rs = [((-B + s) * pow(2*A, -1, q)) % q for s in sq[disc]]
        for rv in rs:
            if rv == 0: continue
            xv = (xr*rv + x0) % q; yv = (yr*rv + y0) % q
            sols.add((xv, yv, min(rv, q - rv)))
    return len(sols)
for q in qs:
    hist = Counter(); n = 0; degenerate = 0
    while n < 3000:
        circ = [(random.randrange(q), random.randrange(q), random.randrange(1, q)) for _ in range(3)]
        if len({c[:2] for c in circ}) < 3: continue
        c = count(q, circ)
        if c is None: degenerate += 1; continue
        hist[c] += 1; n += 1
    mean = sum(k*v for k, v in hist.items())/n
    print(f"q = {q}: counts of rational tangent circles over {n} generic triples: {dict(sorted(hist.items()))}; mean {mean:.3f}; degenerate skipped {degenerate}")
