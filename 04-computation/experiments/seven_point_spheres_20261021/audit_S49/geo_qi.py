"""Independent check of the Q(i) example of THM-4742 (4): V0 = 0, V1 = e1, V2 = e2, V3 = (0, i, i). Own Gaussian rationals."""
from fractions import Fraction as Fr


class Gi:
    __slots__ = ('r', 'i')
    def __init__(self, r, i=0): self.r, self.i = Fr(r), Fr(i)
    def __add__(s, o): o = o if isinstance(o, Gi) else Gi(o); return Gi(s.r + o.r, s.i + o.i)
    __radd__ = __add__
    def __sub__(s, o): o = o if isinstance(o, Gi) else Gi(o); return Gi(s.r - o.r, s.i - o.i)
    def __rsub__(s, o): return Gi(o) - s
    def __mul__(s, o): o = o if isinstance(o, Gi) else Gi(o); return Gi(s.r * o.r - s.i * o.i, s.r * o.i + s.i * o.r)
    __rmul__ = __mul__
    def inv(s):
        n = s.r * s.r + s.i * s.i
        if n == 0: raise ZeroDivisionError
        return Gi(s.r / n, -s.i / n)
    def __truediv__(s, o): o = o if isinstance(o, Gi) else Gi(o); return s * o.inv()
    def __rtruediv__(s, o): return Gi(o) * s.inv()
    def __neg__(s): return Gi(-s.r, -s.i)
    def __eq__(s, o): o = o if isinstance(o, Gi) else Gi(o); return s.r == o.r and s.i == o.i
    def __hash__(s): return hash((s.r, s.i))
    def __repr__(s): return f'({s.r}{"+" if s.i >= 0 else ""}{s.i}i)'


def dot(x, y): return sum((a * b for a, b in zip(x, y)), Gi(0))
def sub(x, y): return [a - b for a, b in zip(x, y)]
def add(x, y): return [a + b for a, b in zip(x, y)]
def sc(x, t): return [a * t for a in x]


def solve(M, rhs):
    n = len(M); A = [list(M[i]) + [rhs[i]] for i in range(n)]
    for col in range(n):
        pr = next((i for i in range(col, n) if not A[i][col] == 0), None)
        if pr is None: raise ZeroDivisionError
        A[col], A[pr] = A[pr], A[col]
        iv = A[col][col].inv(); A[col] = [x * iv for x in A[col]]
        for i in range(n):
            if i != col and not A[i][col] == 0:
                f = A[i][col]; A[i] = [x - f * y for x, y in zip(A[i], A[col])]
    return [A[i][n] for i in range(n)]


def rank(M):
    M = [list(r) for r in M]; rows, cols = len(M), len(M[0]); r = 0
    for col in range(cols):
        pr = next((i for i in range(r, rows) if not M[i][col] == 0), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]
        iv = M[r][col].inv(); M[r] = [x * iv for x in M[r]]
        for i in range(rows):
            if i != r and not M[i][col] == 0:
                f = M[i][col]; M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        r += 1
        if r == rows: break
    return r


I = Gi(0, 1); Z = Gi(0); One = Gi(1)
V = [[Z, Z, Z], [One, Z, Z], [Z, One, Z], [Z, I, I]]
O = solve([sc(sub(V[i], V[0]), 2) for i in (1, 2, 3)], [dot(V[i], V[i]) - dot(V[0], V[0]) for i in (1, 2, 3)])
R2 = dot(sub(V[0], O), sub(V[0], O))
G = [sum((Vi[j] for Vi in V), Gi(0)) / 4 for j in range(3)]
V1, V2, V3 = V[1:]
e1, e2 = sub(V2, V1), sub(V3, V1)
xy = solve([[2 * dot(e1, e1), 2 * dot(e1, e2)], [2 * dot(e2, e1), 2 * dot(e2, e2)]], [dot(e1, e1), dot(e2, e2)])
O1 = add(V1, add(sc(e1, xy[0]), sc(e2, xy[1])))
f = sub(V3, V2)
xy2 = solve([[dot(e1, f), dot(e2, f)], [dot(e1, e2), dot(e2, e2)]], [Z, dot(e1, e2)])
H = add(V1, add(sc(e1, xy2[0]), sc(e2, xy2[1])))
B = add(O, solve([sub(V[i], O) for i in (1, 2, 3)], [R2] * 3))
g = sub(G, O); E = add(O, sc(g, R2 / dot(g, g)))
dv = sub(G, V[0]); lam = (-2 * dot(sub(V[0], O), dv)) / dot(dv, dv); A = add(V[0], sc(dv, lam))
D = add(O, sc(g, 4)); w = sub(O1, O); F = add(O, sc(w, 4))
Sig = sum((dot(sub(V[0], V[i]), sub(V[0], V[i])) for i in (1, 2, 3)), Gi(0))
cr = lambda a, b: [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
delta = dot(sub(V[0], O), cr(sub(add(add(V1, V2), V3), sc(O, 3)), w))
L = dict(A=A, B=B, C=V[0], D=D, E=E, F=F, H=H)
print('O', O, 'R2', R2, 'Sigma', Sig, 'delta', delta, '|G-O|^2', dot(g, g), '|G-V0|^2', dot(dv, dv), '|w|^2', dot(w, w))
print('labels', L)
print('D != E', not all(x == y for x, y in zip(D, E)), '; A == C', all(x == y for x, y in zip(A, V[0])))
P = [L[x] for x in 'ABCDEFH']
M = [[-2 * x for x in pt] + [One, -dot(pt, pt)] for pt in P]
print('seven cospherical:', rank([r[:-1] for r in M]) == rank(M), '; rank', rank([r[:-1] for r in M]))
print('distinct labels:', len({tuple(L[x]) for x in 'ABCDEFH'}))
print('tangent: (V0-O).(G-V0) =', dot(sub(V[0], O), sub(G, V[0])))
# orthocentre altitude check
print('altitudes:', dot(sub(H, V1), sub(V2, V3)), dot(sub(H, V2), sub(V1, V3)), dot(sub(H, V3), sub(V1, V2)))
