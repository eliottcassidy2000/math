# Audit check H (results note section 3 and THM-4780 (7)).
#  (a) Monsky matrix for the congruent-number curve E_n: y^2 = x^3 - n^2 x, n = p_1...p_k odd squarefree (Monsky's
#      appendix to Heath-Brown, Invent. Math. 118 (1994); as stated in Watkins-Donnelly-Elkies-Fisher-Granville-Rogers,
#      "Ranks of quadratic twists of elliptic curves", section 4.1):
#        A_ij = [p_j/p_i] (i != j), A_ii = sum_(j != i) A_ij, [a/p] = (1 - (a/p))/2 in F_2, D_u = diag([u/p_i]),
#        M = [[A + D_2, D_2], [D_2, A + D_-2]],  dim Sel_2(E_n) - 2 = nullity(M).
#      Sanity: parity of nullity(M) = root-number parity (odd iff n = 5, 6, 7 mod 8) for all odd squarefree n < 3000.
#  (b) the four-prime witnesses (all p = 3 mod 8), their Legendre-tournament classes (p_i -> p_j iff (p_j/p_i) = 1),
#      and the class -> nullity table over all 64 labelled 4-tournaments, in BOTH orientation conventions; the |Pf| = 3
#      classes (THM-1455).
#  (c) THM-4780 (7): TT4 >= n(n-1)(n-3)(n-7)/64 on regular tournaments, equality iff doubly regular, and the class-count
#      relations, exhaustively on all regular tournaments of orders 7, 9, 11 (gentourng) and on P_q, q <= 43.
#  (d) The Paley-graph K4 count vs Evans-Pulham-Sheehan, with eta(4z)^6 coefficients from PARI.   Audit 2026-10-25.
import itertools, subprocess, sys
from math import comb
import numpy as np
from sympy import legendre_symbol, factorint, primerange

def br(a, p): return (1 - legendre_symbol(a % p, p)) // 2
def nullity_F2(M):
    M = [[x % 2 for x in row] for row in M]; n = len(M); m = len(M[0]) if M else 0; r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if M[i][c]), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(n):
            if i != r and M[i][c]: M[i] = [a ^ b for a, b in zip(M[i], M[r])]
        r += 1
    return m - r
def monsky(ps):
    k = len(ps)
    A = [[br(ps[j], ps[i]) if i != j else 0 for j in range(k)] for i in range(k)]
    for i in range(k): A[i][i] = sum(A[i][j] for j in range(k) if j != i) % 2
    D2 = [br(2, p) for p in ps]; Dm2 = [br(-2, p) for p in ps]
    M = [[0] * (2 * k) for _ in range(2 * k)]
    for i in range(k):
        for j in range(k):
            M[i][j] = (A[i][j] + (D2[i] if i == j else 0)) % 2
            M[k + i][k + j] = (A[i][j] + (Dm2[i] if i == j else 0)) % 2
        M[i][k + i] = D2[i]; M[k + i][i] = D2[i]
    return nullity_F2(M) if k else 0
out = []
bad = 0; cnt = 0
for n in range(1, 3000, 2):
    f = factorint(n)
    if any(e > 1 for e in f.values()): continue
    cnt += 1
    s = monsky(sorted(f))
    if (s % 2 == 1) != (n % 8 in (5, 7)): bad += 1
out.append(f'(a) Monsky nullity parity = root-number parity (odd iff n = 5, 7 mod 8) for all {cnt} odd squarefree n < 3000: failures {bad}')

def cls(adj):
    sc = tuple(sorted(sum(adj[i][j] for j in range(4) if j != i) for i in range(4)))
    return {(0, 1, 2, 3): 'TT4', (1, 1, 1, 3): 'source+C3', (0, 2, 2, 2): 'sink+C3', (1, 1, 2, 2): 'strong'}[sc]
W = [(3, 11, 59, 83), (3, 11, 19, 107), (3, 11, 19, 163), (3, 11, 19, 43)]
for ps in W:
    assert all(p % 8 == 3 for p in ps)
    adj = [[1 if i != j and legendre_symbol(ps[j], ps[i]) == 1 else 0 for j in range(4)] for i in range(4)]
    out.append(f'(b) {ps}: class (p_i -> p_j iff (p_j/p_i) = 1) = {cls(adj)}; converse convention = '
               f'{cls([[adj[j][i] for j in range(4)] for i in range(4)])}; Monsky nullity = {monsky(list(ps))}')
# abstract table: Monsky A for a labelled tournament T (T_ij = 1 iff i -> j, i.e. (p_j/p_i) = 1): A_ij = 1 - T_ij
def table(conv):
    res = {}
    for bits in itertools.product((0, 1), repeat=6):
        T = [[0] * 4 for _ in range(4)]; t = 0
        for i in range(4):
            for j in range(i + 1, 4):
                T[i][j] = bits[t]; T[j][i] = 1 - bits[t]; t += 1
        A = [[(1 - T[i][j] if conv == 'note' else T[i][j]) if i != j else 0 for j in range(4)] for i in range(4)]
        for i in range(4): A[i][i] = sum(A[i]) % 2
        M = [[A[i][j] + (1 if i == j else 0) for j in range(4)] + [1 if i == j else 0 for j in range(4)] for i in range(4)] + \
            [[1 if i == j else 0 for j in range(4)] + [A[i][j] for j in range(4)] for i in range(4)]
        A2 = [[sum(A[i][s] * A[s][j] for s in range(4)) % 2 for j in range(4)] for i in range(4)]
        Q = [[(A2[i][j] + A[i][j] + (i == j)) % 2 for j in range(4)] for i in range(4)]
        S = [[0 if i == j else (1 if T[i][j] else -1) for j in range(4)] for i in range(4)]
        pf = S[0][1] * S[2][3] - S[0][2] * S[1][3] + S[0][3] * S[1][2]
        res.setdefault(cls(T), set()).add((nullity_F2(M), nullity_F2(Q), abs(pf)))
    return {c: sorted(v) for c, v in res.items()}
out.append(f'(b) class -> {{(nullity M, nullity(A^2+A+I), |Pf|)}}, A from (p_j/p_i) as in the note: {table("note")}')
out.append(f'    same with the converse orientation (bsd_selmer_growth_20261013a.md section 3 convention): {table("converse")}')
# PARI 2-descent on the witnesses (ellrank: [r, R, s, points])
gp_in = ''.join(f'print({n}, " ", ellrank(ellinit([0,0,0,-{n}^2,0]))[1..3]);' for n in [int(np.prod(w)) for w in W]) + 'quit;'
res = subprocess.run(['gp', '-q', '-D', 'parisize=200000000'], input=gp_in, capture_output=True, text=True).stdout.strip().split('\n')
out.append('(b) PARI ellrank [r, R, s] on E_n for the witnesses n = ' + '; '.join(res))

# (c) TT4 lemma
def tt4_and_classes(A):
    n = len(A); counts = {'TT4': 0, 'source+C3': 0, 'sink+C3': 0, 'strong': 0}
    for S in itertools.combinations(range(n), 4):
        sub = [[A[a][b] for b in S] for a in S]
        counts[cls(sub)] += 1
    return counts
def bound(n): return n * (n - 1) * (n - 3) * (n - 7) / 64
def doubly_regular(A):
    A = np.array(A); M = A @ A.T; n = len(A)
    return bool(np.all(M[~np.eye(n, dtype=bool)] == (n - 3) / 4)) if n % 4 == 3 else False
for n in (7, 9, 11):
    m = (n - 1) // 2
    codes = subprocess.run(['gentourng', '-q', f'-d{m}', f'-D{m}', str(n)], capture_output=True, text=True).stdout.split()
    viol = eqbad = relbad = 0; vals = []
    iu = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for c in codes:
        A = [[0] * n for _ in range(n)]
        for (i, j), ch in zip(iu, c):
            if ch == '1': A[i][j] = 1
            else: A[j][i] = 1
        cc = tt4_and_classes(A); t = cc['TT4']; vals.append(t)
        if t < bound(n): viol += 1
        if (t == bound(n)) != doubly_regular(A): eqbad += 1
        e = n * comb(m, 3) - t
        if not (cc['source+C3'] == e and cc['sink+C3'] == e and cc['strong'] == comb(n, 4) - t - 2 * e): relbad += 1
    out.append(f'(c) n = {n}: {len(codes)} regular tournaments; bound {bound(n)}; TT4 in [{min(vals)}, {max(vals)}]; below bound {viol}; '
               f'equality != doubly regular {eqbad}; class-count relation failures {relbad}')
for q in (7, 11, 19, 23, 31, 43):
    sq = {(x * x) % q for x in range(1, q)}
    A = [[1 if (j - i) % q in sq else 0 for j in range(q)] for i in range(q)]
    cc = tt4_and_classes(A); k = (q - 1) // 2
    out.append(f'(c) P_{q}: classes {cc}; TT4 = bound {cc["TT4"] == bound(q)}; source+C3 = sink+C3 = q C(k,3) - TT4 = {q * comb(k, 3) - cc["TT4"]}')

# (d) K4 in Paley graphs
gp_in = 'mf=mfinit([16,3,Mod(15,16)],0); B=mfeigenbasis(mf); print(#B); c=mfcoefs(B[1],700); forprime(p=5,700, if(p%4==1, print(p," ",c[p+1])));quit;'
res = subprocess.run(['gp', '-q'], input=gp_in, capture_output=True, text=True).stdout.split('\n')
cm = {int(l.split()[0]): int(l.split()[1]) for l in res[1:] if l.strip()}
bad_formula = bad_eps = 0; nchk = 0
for p in sorted(cm):
    if p > 200: break
    sq = {(x * x) % p for x in range(1, p)}
    Q = sorted(sq)
    adj = np.zeros((p, p), dtype=np.int64)
    for x in range(p):
        for s in Q: adj[x, (x + s) % p] = 1
    # K4 count: sum over edges {u,v} of edges inside the common neighbourhood, / 6
    cn = adj @ adj
    tri_edges = 0
    for u in range(p):
        for v in range(u + 1, p):
            if adj[u, v]:
                C = np.nonzero(adj[u] & adj[v])[0]
                tri_edges += int(adj[np.ix_(C, C)].sum()) // 2
    k4 = tri_edges // 6
    nchk += 1
    if 1536 * k4 != p * (p - 1) * (p * p - 20 * p + 81 + cm[p]): bad_formula += 1
    y = next(y for y in range(0, p, 2) if int(round((p - y * y) ** 0.5)) ** 2 == p - y * y)   # q = x^2 + y^2, y even
    if 1536 * k4 != p * (p - 1) * ((p - 9) ** 2 - 4 * y * y): bad_eps += 1
out.append(f'(d) Paley graph K4 counts (brute force) at the {nchk} primes p = 1 mod 4, 5 <= p <= 200: '
           f'q(q-1)(q^2-20q+81+a_q(eta(4z)^6, PARI))/1536 failures {bad_formula}; EPS q(q-1)((q-9)^2-4y^2)/1536 failures {bad_eps}')
print('\n'.join(out))
with open('audit_selmer_tt4.out', 'w') as fh: fh.write('\n'.join(out) + '\n')

# (e) results note section 3.2: the normalised K4 sum equals a_q(eta(4z)^6) at q = 1 mod 4 and vanishes at q = 3 mod 4
from audit_graph_sums import S_norm
K4 = [(i, j) for i in range(4) for j in range(i + 1, 4)]
bad1 = [p for p in sorted(cm) if p <= 200 and S_norm(p, 4, K4) != cm[p]]
bad3 = [p for p in primerange(7, 200) if p % 4 == 3 and S_norm(p, 4, K4) != 0]
line = (f'(e) S_K4 (x0 = 0, x1 = 1, order orientation) = a_q(eta(4z)^6, PARI 16.3.c.a) at all primes q = 1 mod 4 below 200: '
        f'failures {bad1}; S_K4 = 0 at all q = 3 mod 4 in 7..199: failures {bad3}')
print(line)
with open('audit_selmer_tt4.out', 'a') as fh: fh.write(line + '\n')
