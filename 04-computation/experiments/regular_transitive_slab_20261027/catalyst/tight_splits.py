"""Window-free exceptionality (THM-4820 (C)).  For x in Z[omega] with Q(x) = 1 (mod 3), a merge x = y + z is tight
(zero slack against Phi = floor(Q/3)) iff y is perpendicular to z (beta = 0) or beta(y, z) = -1; any merge tree of
cost floor(Q(x)/3) has a tight split of x with y, z != 0.  So if x has neither, g(x) >= floor(Q(x)/3) + 1 whatever the
tree.  Splits are found exactly: with y = a + b omega, u = 2a - b, the conditions read (u - U)^2 + 3(b - B)^2-type norm
equations; here we simply scan y with Q(y) <= Q(x) + 1 (both parts of a tight split have norm <= Q(x) + 1)."""
import sys
def Q(a, b): return a * a - a * b + b * b
def tight_split(xa, xb):
    q = Q(xa, xb); R = int((q + 2) ** 0.5) + 2
    for a in range(-2 * R, 2 * R + 1):
        for b in range(-2 * R, 2 * R + 1):
            if (a, b) in ((0, 0), (xa, xb)): continue
            qy = Q(a, b)
            if qy > q + 1: continue
            za, zb = xa - a, xb - b; beta = q - qy - Q(za, zb)
            if beta in (0, -1): return (a, b, za, zb, beta)
    return None
pts = [tuple(map(int, l.split()[1:])) for l in open(sys.argv[1]) if l.startswith('P ')]
seen = set(); proved = []; window = []
for a, b, q, g in sorted(pts, key=lambda t: t[2]):
    if q > 2000 or q <= 1 or g - q // 3 != 1 or q in seen: continue
    seen.add(q); s = tight_split(a, b)
    (window if s else proved).append((q, (a, b), s))
print("exceptional norms with NO tight split (exceptional for every tree, window-free):", [p[0] for p in proved])
print("exceptional norms with a tight split whose parts are not tight-buildable in the window:", [(w[0], w[2]) for w in window])
