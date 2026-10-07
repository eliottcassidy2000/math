import sympy as sp, itertools
from collections import defaultdict
x, y, z, t = sp.symbols('x y z t')
def F(P):
    X, Y, Z = P; U = 1 + X*Y
    return tuple(sp.expand(e) for e in (U**3*Z + Y**2*U*(4 + 3*X*Y), Y + 3*X*U**2*Z + 3*X*Y**2*(4 + 3*X*Y), 2*X - 3*X**2*Y - X**3*Z))
w = 2*t + 1   # general odd w
P0 = (0, 2*w, -(63*w**2 + 1)/4); P1 = (1, (w - 3)/2, (13 - 3*w)/2); P2 = (-1, (w + 3)/2, (13 + 3*w)/2)
T = ((w**2 - 1)/4, 2*w, 0)
print("family identity (symbolic in t, w = 2t+1):", all(sp.simplify(a - b) == 0 for P in (P0, P1, P2) for a, b in zip(F(P), T)))
print("coordinates integral for integer t:", [sp.expand(c) for c in P0 + P1 + P2])
print("example w=1:", F((1, -1, 5)), F((-1, 2, 8)), F((0, 2, -16)))
# integer collisions in a box, sorted by max-norm of the pair
B = 6
img = defaultdict(list)
for P in itertools.product(range(-B, B+1), repeat=3):
    X, Y, Z = P; U = 1 + X*Y
    img[(U**3*Z + Y**2*U*(4 + 3*X*Y), Y + 3*X*U**2*Z + 3*X*Y**2*(4 + 3*X*Y), 2*X - 3*X**2*Y - X**3*Z)].append(P)
pairs = sorted((max(max(map(abs, P)) for P in Ps), w0, Ps) for w0, Ps in img.items() if len(Ps) > 1)
for p in pairs[:6]: print("collision (max-norm %d):" % p[0], p[1], "<-", p[2])
