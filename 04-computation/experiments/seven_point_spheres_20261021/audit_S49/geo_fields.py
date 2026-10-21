"""Independent audit of THM-4742 statement 4 over F_p (exhaustive p = 5 with V0 = 0; random for larger p) and Q(i).
A configuration is 'defined' if every label exists: nondegenerate tetrahedron, face circumcentre and orthocentre unique,
B unique (tangent planes), |G-O|^2 != 0, |G-V0|^2 != 0.  'Regime' adds delta != 0 and D != E.
Cosphericity = consistency of the linear system |x|^2 - 2u.x + c = 0 (rank test, independent of a choice of 4 points)."""
import sys, random, itertools
from fractions import Fraction as Fr


def make_field(p):
    def inv(x):
        x %= p
        if x == 0:
            raise ZeroDivisionError
        return pow(x, p - 2, p)
    return inv


def rank_mod(M, p):
    M = [[x % p for x in r] for r in M]
    rows, cols = len(M), len(M[0]); r = 0
    for col in range(cols):
        pr = next((i for i in range(r, rows) if M[i][col]), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        iv = pow(M[r][col], p - 2, p)
        M[r] = [x * iv % p for x in M[r]]
        for i in range(rows):
            if i != r and M[i][col]:
                f = M[i][col]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        r += 1
        if r == rows:
            break
    return r


def solve_mod(M, rhs, p):
    n = len(M)
    A = [list(M[i]) + [rhs[i]] for i in range(n)]
    for col in range(n):
        pr = next((i for i in range(col, n) if A[i][col] % p), None)
        if pr is None:
            raise ZeroDivisionError
        A[col], A[pr] = A[pr], A[col]
        iv = pow(A[col][col] % p, p - 2, p)
        A[col] = [x * iv % p for x in A[col]]
        for i in range(n):
            if i != col and A[i][col] % p:
                f = A[i][col]
                A[i] = [(x - f * y) % p for x, y in zip(A[i], A[col])]
    return [A[i][n] % p for i in range(n)]


def dot(x, y, p): return sum(a * b for a, b in zip(x, y)) % p
def sub(x, y, p): return [(a - b) % p for a, b in zip(x, y)]
def add(x, y, p): return [(a + b) % p for a, b in zip(x, y)]
def sc(x, t, p): return [a * t % p for a in x]
def cross(a, b, p): return [(a[1] * b[2] - a[2] * b[1]) % p, (a[2] * b[0] - a[0] * b[2]) % p, (a[0] * b[1] - a[1] * b[0]) % p]


def labels_mod(V, p):
    inv = make_field(p)
    O = solve_mod([sc(sub(V[i], V[0], p), 2, p) for i in (1, 2, 3)], [(dot(V[i], V[i], p) - dot(V[0], V[0], p)) % p for i in (1, 2, 3)], p)
    R2 = dot(sub(V[0], O, p), sub(V[0], O, p), p)
    i4 = inv(4)
    G = [sum(Vi[j] for Vi in V) * i4 % p for j in range(3)]
    V1, V2, V3 = V[1], V[2], V[3]
    e1, e2 = sub(V2, V1, p), sub(V3, V1, p)
    # face circumcentre in plane: O' = V1 + x e1 + y e2, 2 ei.(O'-V1) = |ei|^2
    xy = solve_mod([[2 * dot(e1, e1, p), 2 * dot(e1, e2, p)], [2 * dot(e2, e1, p), 2 * dot(e2, e2, p)]], [dot(e1, e1, p), dot(e2, e2, p)], p)
    O1 = add(V1, add(sc(e1, xy[0], p), sc(e2, xy[1], p), p), p)
    # orthocentre in plane: H = V1 + x e1 + y e2, (H - V1).(V3 - V2) = 0, (H - V2).(V3 - V1) = 0
    f = sub(V3, V2, p)
    # (x e1 + y e2).f = 0 ; (x e1 + y e2 - e1).e2 = 0
    xy2 = solve_mod([[dot(e1, f, p), dot(e2, f, p)], [dot(e1, e2, p), dot(e2, e2, p)]], [0, dot(e1, e2, p)], p)
    H = add(V1, add(sc(e1, xy2[0], p), sc(e2, xy2[1], p), p), p)
    # B: (B - O).(Vi - O) = R2
    Bm = solve_mod([sub(V[i], O, p) for i in (1, 2, 3)], [R2] * 3, p)
    B = add(O, Bm, p)
    g = sub(G, O, p); gg = dot(g, g, p)
    E = add(O, sc(g, R2 * inv(gg), p), p)
    dv = sub(G, V[0], p); ddv = dot(dv, dv, p)
    lam = (-2 * dot(sub(V[0], O, p), dv, p)) * inv(ddv) % p
    A = add(V[0], sc(dv, lam, p), p)
    D = add(O, sc(g, 4, p), p)
    w = sub(O1, O, p)
    F = add(O, sc(w, 4, p), p)
    Sig = sum(dot(sub(V[0], V[i], p), sub(V[0], V[i], p), p) for i in (1, 2, 3)) % p
    sv = sub(add(add(V1, V2, p), V3, p), sc(O, 3, p), p)
    delta = dot(sub(V[0], O, p), cross(sv, w, p), p)
    return dict(A=A, B=B, C=list(V[0]), D=D, E=E, F=F, H=H, O=O, R2=R2, Sig=Sig, delta=delta)


def cosph(P, p):
    M = [[(-2 * x) % p for x in pt] + [1, (-dot(pt, pt, p)) % p] for pt in P]
    return rank_mod([r[:-1] for r in M], p) == rank_mod(M, p)


def run(p, configs):
    st = dict(defined=0, regime=0, six_fail=0, five_fail=0, seven=0, sig0=0, mism_seven=0, mism_AC=0, R2zero_in_regime=0,
              BeqF_in_regime=0, distinct7_on_sphere=0)
    out_regime = {}
    for V in configs:
        V = [list(x) for x in V]
        M = [sub(V[i], V[0], p) for i in (1, 2, 3)]
        if rank_mod(M, p) < 3:
            continue
        try:
            L = labels_mod(V, p)
        except ZeroDivisionError:
            continue
        st['defined'] += 1
        sig0 = L['Sig'] == 0
        seven = cosph([L[x] for x in 'ABCDEFH'], p)
        if L['delta'] == 0 or L['D'] == L['E']:
            key = ('delta0' if L['delta'] == 0 else 'D=E', seven, sig0, L['A'] == L['C'],
                   len({tuple(L[x]) for x in 'ABCDEFH'}) == 7)
            out_regime[key] = out_regime.get(key, 0) + 1
            continue
        st['regime'] += 1
        st['R2zero_in_regime'] += L['R2'] == 0
        st['BeqF_in_regime'] += L['B'] == L['F']
        six = cosph([L[x] for x in 'ABDEFH'], p)
        five = cosph([L[x] for x in 'BCDEF'], p)
        st['six_fail'] += not six; st['five_fail'] += not five
        st['seven'] += seven; st['sig0'] += sig0
        st['mism_seven'] += seven != sig0
        st['mism_AC'] += sig0 != (L['A'] == L['C'])
        st['distinct7_on_sphere'] += seven and len({tuple(L[x]) for x in 'ABCDEFH'}) == 7
    return st, out_regime


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'exh5':
        p = 5
        rng_pts = list(itertools.product(range(p), repeat=3))
        def gen():
            for V1 in rng_pts:
                for V2 in rng_pts:
                    for V3 in rng_pts:
                        yield ((0, 0, 0), V1, V2, V3)
        st, outr = run(p, gen())
        print('p = 5 exhaustive (V0 = 0):', st)
        print('   outside regime (reason, seven-cospherical, Sigma=0, A=C, 7 distinct): count', outr)
    else:
        for p in (7, 11, 13, 17, 19, 23):
            rng = random.Random(1000 + p)
            def gen():
                for _ in range(40000):
                    yield [tuple(rng.randrange(p) for _ in range(3)) for _ in range(4)]
            st, outr = run(p, gen())
            print(f'p = {p} random 40000:', st)
            print('   outside regime:', outr)
