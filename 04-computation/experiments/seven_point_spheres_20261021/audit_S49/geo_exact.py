"""Independent audit of THM-4742: actual-coordinate exact checks in dimension n = 2..6 (own constructions).
Labels: A second point of V0G on circumsphere; B pole of face hyperplane (B-O).(Vi-O) = R^2; C = V0;
D = O + (n+1)(G-O); E = inverse of G; F = O + (n+1)(O'-O); H' = O' + n(G'-O') (the note's analogue; = orthocentre
for n = 3, checked against altitude equations).  Spheres |x|^2 - 2u.x + c = 0 fitted by exact elimination."""
import random, itertools
from fractions import Fraction as Fr


def dot(x, y): return sum(a * b for a, b in zip(x, y))
def sub(x, y): return [a - b for a, b in zip(x, y)]
def add(x, y): return [a + b for a, b in zip(x, y)]
def sc(x, t): return [a * t for a in x]


def rref(M):
    """Exact RREF; returns (rank, rows, pivots)."""
    M = [list(r) for r in M]
    rows, cols = len(M), len(M[0])
    piv = []; r = 0
    for col in range(cols):
        pr = next((i for i in range(r, rows) if M[i][col] != 0), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        inv = 1 / Fr(M[r][col])
        M[r] = [x * inv for x in M[r]]
        for i in range(rows):
            if i != r and M[i][col] != 0:
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv.append(col); r += 1
        if r == rows:
            break
    return r, M, piv


def solve_unique(M, rhs):
    n = len(M[0])
    rk, R, piv = rref([list(M[i]) + [rhs[i]] for i in range(len(M))])
    if rk != n or piv[-1] == n:
        raise ValueError('not unique')
    return [R[i][n] for i in range(n)]


def labels(V):
    n = len(V) - 1
    O = solve_unique([sc(sub(V[i], V[0]), 2) for i in range(1, n + 1)], [dot(V[i], V[i]) - dot(V[0], V[0]) for i in range(1, n + 1)])
    R2 = dot(sub(V[0], O), sub(V[0], O))
    G = [sum(Vi[j] for Vi in V) / Fr(n + 1) for j in range(n)]
    F_ = V[1:]
    # face circumcentre: O' = V1 + sum c_j (Vj - V1), 2 (Vi - V1).(O' - V1) = |Vi - V1|^2
    E_ = [sub(F_[j], F_[0]) for j in range(1, n)]
    cc = solve_unique([[2 * dot(E_[i], E_[j]) for j in range(n - 1)] for i in range(n - 1)], [dot(E_[i], E_[i]) for i in range(n - 1)])
    O1 = list(F_[0])
    for j in range(n - 1):
        O1 = add(O1, sc(E_[j], cc[j]))
    G1 = [sum(Vi[j] for Vi in F_) / Fr(n) for j in range(n)]
    H = add(O1, sc(sub(G1, O1), n))
    B = add(O, solve_unique([sub(F_[i], O) for i in range(n)], [R2] * n))
    g = sub(G, O); w = sub(O1, O); dvec = sub(G, V[0])
    E = add(O, sc(g, R2 / dot(g, g)))
    lam = -2 * dot(sub(V[0], O), dvec) / dot(dvec, dvec)
    A = add(V[0], sc(dvec, lam))
    D = add(O, sc(g, n + 1)); F = add(O, sc(w, n + 1))
    Sig = sum(dot(sub(V[0], V[i]), sub(V[0], V[i])) for i in range(1, n + 1))
    return dict(A=A, B=B, C=list(V[0]), D=D, E=E, F=F, H=H, O=O, R2=R2, G=G, O1=O1, G1=G1, Sig=Sig)


def sphere_system(P):
    """rows for unknowns (u_1..u_n, c): -2 u.x + c = -|x|^2."""
    return [[-2 * x for x in pt] + [Fr(1), -dot(pt, pt)] for pt in P]


def cospherical(P):
    n = len(P[0])
    M = sphere_system(P)
    rk_aug, R, piv = rref(M)
    rk = rref([r[:-1] for r in M])[0]
    return rk == rk_aug, n + 1 - rk                     # consistent?, dimension of solution family


def affine_rank(P):
    return rref([sub(p, P[0]) for p in P[1:]])[0]


rng = random.Random(4742)


def rand_simplex(n):
    while True:
        V = [[Fr(rng.randint(-12, 12), rng.randint(1, 4)) for _ in range(n)] for _ in range(n + 1)]
        if rref([sub(V[i], V[0]) for i in range(1, n + 1)])[0] == n:
            return V


# ---- n = 3: full THM-4742 statements 1-3 in actual coordinates
bad = 0; cnt = 0
for _ in range(300):
    V = rand_simplex(3)
    try:
        L = labels(V)
    except (ValueError, ZeroDivisionError):
        continue
    if L['D'] == L['E']:
        continue
    # orthocentre check (altitudes)
    H = L['H']; V1, V2, V3 = V[1], V[2], V[3]
    ok = dot(sub(H, V1), sub(V2, V3)) == 0 and dot(sub(H, V2), sub(V1, V3)) == 0
    # S2 through A,B,D,E,F,H'
    good, dim = cospherical([L[x] for x in 'ABDEFH'])
    ok &= good and dim == 0
    good1, dim1 = cospherical([L[x] for x in 'BCDEF'])
    ok &= good1 and dim1 == 0
    # S2 explicit, power at O and at C
    M = sphere_system([L[x] for x in 'DFHB'])
    sol = solve_unique([r[:-1] for r in M], [r[-1] for r in M]); u, c = sol[:3], sol[3]
    pw = lambda X: dot(X, X) - 2 * dot(u, X) + c
    ok &= pw(L['O']) == 4 * L['R2'] and pw(L['A']) == 0 and pw(L['E']) == 0 and pw(L['C']) == L['Sig']
    ok &= L['Sig'] == -8 * dot(sub(V[0], L['O']), sub(L['G'], V[0]))
    M1 = sphere_system([L[x] for x in 'CDFB'])
    sol1 = solve_unique([r[:-1] for r in M1], [r[-1] for r in M1]); u1, c1 = sol1[:3], sol1[3]
    pw1 = lambda X: dot(X, X) - 2 * dot(u1, X) + c1
    ok &= pw1(L['H']) == L['Sig'] and pw1(L['O']) == 4 * L['R2']
    # radical plane {x : pw1 - pw2 = 0} = {x: 2(u - u1).x + (c1 - c) = 0} contains O, G, O'
    rad = lambda X: pw1(X) - pw(X)
    ok &= rad(L['O']) == 0 and rad(L['G']) == 0 and rad(L['O1']) == 0
    ok &= all(rad(L[x]) == 0 for x in 'BDEF')
    ok &= affine_rank([L[x] for x in 'ABCDEFH'] + [L['O']]) == 3
    cnt += 1; bad += not ok
print(f'n=3: {cnt} rational tetrahedra, failures {bad}')

# ---- other dimensions: affine rank of the labels and cospherical sets
for n in (2, 3, 4, 5, 6):
    stats = {}
    trials = 0
    while trials < 12:
        V = rand_simplex(n)
        try:
            L = labels(V)
        except (ValueError, ZeroDivisionError):
            continue
        if L['D'] == L['E']:
            continue
        trials += 1
        rk = affine_rank([L[x] for x in 'ABCDEFH'] + [L['O']])
        six, dim6 = cospherical([L[x] for x in 'ABDEFH'])
        five, dim5 = cospherical([L[x] for x in 'BCDEF'])
        seven, _ = cospherical([L[x] for x in 'ABCDEFH'])
        # power of C on any sphere of the six-label family (unique within the 3-flat): use min-norm solution? use one solution
        key = (rk, six, dim6, five, dim5, seven)
        stats[key] = stats.get(key, 0) + 1
    print(f'n={n}: (affine rank of labels+O, ABDEFH cospherical, family dim, BCDEF cospherical, family dim, all 7): {stats}')
    if n == 2:
        cnt2 = {}
        for _ in range(12):
            V = rand_simplex(2)
            L = labels(V)
            if L['D'] == L['E']:
                continue
            for S in itertools.combinations('ABCDEFH', 4):
                pts = [L[x] for x in S]
                if len({tuple(p) for p in pts}) < 4:
                    continue
                if cospherical(pts)[0]:
                    cnt2[''.join(S)] = cnt2.get(''.join(S), 0) + 1
        print('   n=2 concyclic 4-sets of distinct points (count over 12 triangles):', cnt2)

# ---- n = 4: with C, does the pencil persist? power of C on the six-label spheres = Sigma within the 3-flat
for n in (4, 5):
    V = rand_simplex(n)
    L = labels(V)
    M = sphere_system([L[x] for x in 'ABDEFH'])
    rk, R, piv = rref(M)
    # pick particular solution (free vars = 0)
    sol = [Fr(0)] * (n + 1)
    for i, col in enumerate(piv):
        sol[col] = R[i][-1]
    u, c = sol[:n], sol[n]
    pw = lambda X: dot(X, X) - 2 * dot(u, X) + c
    print(f'n={n}: one sphere through ABDEFH: residuals {[pw(L[x]) for x in "ABDEFH"]}, Pow(C) = {pw(L["C"])}, '
          f'Sigma = {L["Sig"]}, Pow(O) = {pw(L["O"])}, (n+1)R^2 = {(n + 1) * L["R2"]}')
