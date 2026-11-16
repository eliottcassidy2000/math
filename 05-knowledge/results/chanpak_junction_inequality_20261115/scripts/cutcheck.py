"""Interior cuts are bottlenecks (sanity check of the lemma used for the non-strict growth bound).
For random span-3 factors W (untruncated letters, n vertices): if the induced tournament T[W] is an ordinal sum
A => B with 7 <= max(A) <= n - 3, then the operator of W has rank <= 1 in the minimal representation.
Also: for words of rank >= 2 the largest strong component of T[W] has at least n - 12 vertices."""
import pickle, random, numpy as np
B = pickle.load(open('lit_base.pkl', 'rb'))
A, O, Rb = B['A'], B['O'], B['Rb']
P1 = 2147483629
def rank_mod(M, p=P1):
    M = [[int(x) % p for x in row] for row in M]; r = 0; rows = len(M); cols = len(M[0])
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i][c]), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]; inv = pow(M[r][c], p - 2, p); M[r] = [x * inv % p for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        r += 1
    return r
p = P1
OR = (O @ Rb)
def pivots(M):
    M = [[int(x) % p for x in r] for r in M]; rows = len(M); cols = len(M[0]); piv = []; r = 0
    for c in range(cols):
        pr = next((i for i in range(r, rows) if M[i][c]), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]; iv = pow(M[r][c], p - 2, p); M[r] = [x * iv % p for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == rows: break
    return piv
def inv_mod(M):
    n = M.shape[0]; M = [[int(x) % p for x in r] + [1 if i == j else 0 for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        piv = next(i for i in range(c, n) if M[i][c]); M[c], M[piv] = M[piv], M[c]
        iv = pow(M[c][c], p - 2, p); M[c] = [x * iv % p for x in M[c]]
        for i in range(n):
            if i != c and M[i][c]:
                f = M[i][c]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[c])]
    return np.array([r[n:] for r in M], dtype=np.int64)
def modmul(X, Y):
    Xh, Xl = X >> 16, X & 0xFFFF
    return (((Xh @ Y) % p) * 65536 + (Xl @ Y)) % p
J = pivots(OR.tolist()); Iv = inv_mod(OR[:, J])
Mm = [modmul((O @ A[l] @ Rb)[:, J] % p, Iv) for l in range(8)]
def oprank(word):
    M = Mm[word[0]]
    for l in word[1:]: M = modmul(Mm[l], M)
    return rank_mod(M.tolist())
def components(word):
    n = len(word); beats = [[False] * n for _ in range(n)]
    for q in range(n):
        for s in range(1, 4):
            if q - s >= 0:
                if (word[q] >> (s - 1)) & 1: beats[q][q - s] = True
                else: beats[q - s][q] = True
        for r in range(q - 3): beats[r][q] = True
    reach = [[beats[i][j] or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    if reach[k][j]: reach[i][j] = True
    comp = {}
    for i in range(n):
        comp.setdefault(frozenset(j for j in range(n) if reach[i][j] and reach[j][i]), None)
    comps = sorted(comp, key=lambda c: -sum(1 for i in range(n) if any(reach[x][i] for x in c)))   # by reach size: first = source side
    return comps
random.seed(11)
ncut = nok = 0; nbf = 0; worst = 99
for trial in range(1500):
    n = random.randint(11, 16)
    word = [random.choice([0, 0, 0, 1, 2, 4, 3, 5, 6, 7]) for _ in range(n)]
    comps = components(word)
    # interior cut: a union of initial components with 7 <= max <= n-3 (1-based)
    interior = False; acc = set()
    for c in comps[:-1]:
        acc |= c
        if 7 <= max(acc) + 1 <= n - 3: interior = True
    r = oprank(word)
    if interior:
        ncut += 1
        if r <= 1: nok += 1
        else: print("COUNTEREXAMPLE", word, r)
    if r >= 2:
        nbf += 1; worst = min(worst, max(len(c) for c in comps) - n)
print("words with an interior cut: %d ; of these with operator rank <= 1: %d" % (ncut, nok))
print("words of operator rank >= 2: %d ; smallest (largest strong component size) - n over them: %d" % (nbf, worst))
