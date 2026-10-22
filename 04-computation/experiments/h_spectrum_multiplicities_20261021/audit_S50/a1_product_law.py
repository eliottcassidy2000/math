"""Audit S50 - a1: THM-4744 (1), the product law H(IC(p)) = prod (2^(p_i-1)+1), and the proof's path-shape claims.
  (1a) every composition of N <= 15 into parts >= 2 (n <= 16), by subset DP; a random subset re-checked by
       inclusion-exclusion; extra compositions with large parts (up to a single part 16, n = 17);
  (1b) H(NT_s) = 2^(s-2) + 1 for s = 3..17;
  (1c) the proof, checked path by path for every composition with N <= 9 and k >= 2 (all Hamiltonian paths enumerated by
       DFS): with A = first block, v = x_(p_1), B = the rest, every Hamiltonian path has the shape X_a X_b v Y_a Y_b,
       neither bad shape (i) nor (ii) occurs, every path restricts to Hamiltonian paths of A and B, and the restriction
       map is a bijection onto pairs (so every pair has exactly one merge);
  (1d) optimal parts: splitting a part p >= 5 into (2, p-2) strictly lowers the product; block costs per vertex;
       for m <= 40 the compositions (parts >= 2) attaining min 3^x 5^y are exactly the {2,3,4} ones counted by
       q, 1, q(q+3)/2 (m <= 26 over ALL compositions with parts >= 2; beyond, over parts {2,...,8})."""
import random, time, math
import numpy as np
import alib as L

T0 = time.time()
FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


rng = random.Random(7)
# (1a)
bad = []; tested = 0; ie_tested = 0
for N in range(2, 16):
    for p in L.compositions(N):
        adj = L.chain(p)
        h = L.H(adj)
        tested += 1
        if h != L.prod_law(p):
            bad.append(p)
        if N <= 12 and rng.random() < 0.15:
            ie_tested += 1
            if int(L.ham_ie(len(adj), np.array(adj, dtype=np.int64))) != h:
                bad.append(('ie', p))
extra = [(16,), (2, 14), (14, 2), (5, 5, 6), (7, 9), (8, 8), (3, 13), (6, 2, 8), (2, 2, 2, 2, 2, 2, 2, 2), (9, 2, 5),
         (4, 4, 4, 4), (5, 11), (10, 6), (12, 4), (2, 12, 2)]
for p in extra:
    tested += 1
    if L.H(L.chain(p)) != L.prod_law(p):
        bad.append(p)
check(not bad, f'(1a) H(IC(p)) = prod (2^(p_i - 1) + 1) for {tested} compositions (all with N <= 15, parts >= 2; plus '
               f'{len(extra)} with large parts, n = 17); {ie_tested} re-checked by inclusion-exclusion' + (f' BAD {bad[:5]}' if bad else ''))
# (1b)
check(all(L.H(L.chain((s - 1,))) == 2 ** (s - 2) + 1 for s in range(3, 18)), '(1b) H(NT_s) = 2^(s-2) + 1 for s = 3..17')


# (1c)
def all_paths(adj):
    n = len(adj)
    out = []
    full = (1 << n) - 1

    def dfs(path, used):
        if used == full:
            out.append(tuple(path))
            return
        u = path[-1]
        nxt = adj[u] & ~used
        while nxt:
            low = nxt & -nxt
            w = low.bit_length() - 1
            path.append(w)
            dfs(path, used | low)
            path.pop()
            nxt ^= low

    for s in range(n):
        dfs([s], 1 << s)
    return out


def arc(adj, u, w):
    return (adj[u] >> w) & 1 == 1


ok_shape = ok_i = ok_ii = ok_restr = ok_bij = True
nchains = 0; npaths = 0
for N in range(4, 10):
    for p in L.compositions(N):
        if len(p) < 2:
            continue
        nchains += 1
        adj = L.chain(p)
        n = len(adj)
        p1 = p[0]
        v = p1
        Aset = set(range(0, p1 + 1)); Ap = Aset - {v}
        Bset = set(range(p1, n)); Bp = Bset - {v}
        paths = all_paths(adj)
        npaths += len(paths)
        pairs = set()
        for P in paths:
            iv = P.index(v)
            X, Y = P[:iv], P[iv + 1:]

            def split(seq):
                # seq must be (A' vertices)(B' vertices)
                k = 0
                while k < len(seq) and seq[k] in Ap:
                    k += 1
                return seq[:k], seq[k:], all(x in Bp for x in seq[k:])
            Xa, Xb, okx = split(X)
            Ya, Yb, oky = split(Y)
            ok_shape &= okx and oky
            # bad shape (i): Xa, Xb nonempty and v -> last(Xa)
            if Xa and Xb and arc(adj, v, Xa[-1]):
                ok_i = False
            # bad shape (ii): Ya, Yb nonempty and first(Yb) -> v
            if Ya and Yb and arc(adj, Yb[0], v):
                ok_ii = False
            PA = tuple(x for x in P if x in Aset); PB = tuple(x for x in P if x in Bset)
            rA = all(arc(adj, PA[t], PA[t + 1]) for t in range(len(PA) - 1))
            rB = all(arc(adj, PB[t], PB[t + 1]) for t in range(len(PB) - 1))
            ok_restr &= rA and rB
            pairs.add((PA, PB))
        hA = 2 ** (p1 - 1) + 1
        hB = L.prod_law(p[1:])
        ok_bij &= (len(pairs) == len(paths) == hA * hB)
check(ok_shape, f'(1c) every Hamiltonian path of every chain with N = 4..9, k >= 2 ({nchains} chains, {npaths} paths) has the shape X_a X_b v Y_a Y_b')
check(ok_i and ok_ii, '(1c) the bad shapes (i) [X_a, X_b nonempty, v -> last(X_a)] and (ii) [Y_a, Y_b nonempty, first(Y_b) -> v] never occur')
check(ok_restr and ok_bij, '(1c) every path restricts to Hamiltonian paths of A and B, and restriction is a bijection onto the H(A)H(B) pairs')

# (1d)
check(all(2 ** (p - 1) + 1 > 3 * (2 ** (p - 3) + 1) for p in range(5, 60)), '(1d) 2^(p-1)+1 > 3(2^(p-3)+1) for p >= 5: a part >= 5 is never optimal')
costs = [round(math.log(2 ** (p - 1) + 1) / p, 3) for p in (2, 3, 4, 5, 6)]
check(costs[:4] == [0.549, 0.536, 0.549, 0.567], f'(1d) per-vertex block costs log(2^(p-1)+1)/p, p = 2..6: {costs}')
ok = True; ok_closed = True; rows = []
for m in range(3, 25):                       # ALL compositions with parts >= 2
    f = L.busch_min(m)
    ok_closed &= f == L.busch(m)
    allc = list(L.compositions(m - 1, 2))
    opt = [c for c in allc if L.prod_law(c) == f]
    below = [c for c in allc if L.prod_law(c) < f]
    q, r = divmod(m, 3)
    pred = {0: q, 1: 1, 2: q * (q + 3) // 2}[r]
    ok &= (len(opt) == pred) and all(max(c) <= 4 for c in opt) and not below
    if m <= 17:
        rows.append(f'{m}:{len(opt)}')
check(ok_closed, '(1d) the closed form f(3q)=3*5^(q-1), f(3q+1)=5^q, f(3q+2)=9*5^(q-1) equals min 3^x 5^y (2x+3y=m-1), m <= 24')
check(ok, '(1d) over ALL compositions with parts >= 2: optimal ones have parts in {2,3,4}, none is below f(m), and their number is '
          f'q, 1, q(q+3)/2 for m = 3..24: ' + ', '.join(rows) + ' ...')
# beyond: by the exchange inequality only parts {2,3,4} matter; count multinomially
ok = True
for m in range(3, 121):
    f = L.busch_min(m); ok &= f == L.busch(m)
    cnt = 0
    for n4 in range(0, (m - 1) // 4 + 1):
        for n3 in range(0, (m - 1) // 3 + 1):
            rest = m - 1 - 4 * n4 - 3 * n3
            if rest < 0 or rest % 2:
                continue
            n2 = rest // 2
            if 3 ** n2 * 5 ** n3 * 9 ** n4 == f:
                cnt += math.factorial(n2 + n3 + n4) // (math.factorial(n2) * math.factorial(n3) * math.factorial(n4))
    q, r = divmod(m, 3)
    ok &= cnt == {0: q, 1: 1, 2: q * (q + 3) // 2}[r]
check(ok, '(1d) multinomial count of optimal {2,3,4}-compositions equals q, 1, q(q+3)/2 and the closed form of f holds, m = 3..120')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
