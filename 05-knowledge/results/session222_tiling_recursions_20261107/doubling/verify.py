import itertools, json
from pathlib import Path

CHECKS=0
def require(condition, message):
 global CHECKS
 CHECKS += 1
 if not condition: raise RuntimeError(message)


def tournament_word(w):
    n = len(w)
    A = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(n):
        for j in range(i + 1, n):
            A[i][j] = 1
        A[i][n] = int(w[i])
        A[n][i] = 1 - int(w[i])
    for i in range(n):
        for j in range(i):
            A[i][j] = 1 - A[j][i]
    return A

def H(A):
    n = len(A)
    out = [sum((a << j for j, a in enumerate(row))) for row in A]
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for s in range(1, 1 << n):
        for i, c in enumerate(dp[s]):
            if not c:
                continue
            rest = out[i] & ~s
            while rest:
                b = rest & -rest
                rest -= b
                dp[s | b][b.bit_length() - 1] += c
    return sum(dp[-1])

def C(w):
    return sum(((1 - int(w[i])) * int(w[j]) * (1 if j - i == 1 else 2 ** (j - i - 2)) for i in range(len(w)) for j in range(i + 1, len(w))))

def double(A):
    n = len(A)
    M = [[0 if i == j else 2 * A[i][j] - 1 for j in range(n)] for i in range(n)]
    return [[int((M[i % n][j % n] if i < n and j < n else -M[i % n][j % n] if i >= n and j >= n else M[i % n][j % n] + (int(i == j % n) if i < n else -int(i % n == j))) > 0) for j in range(2 * n)] for i in range(2 * n)]

def weight(A, p):
    return sum((A[p[i]][p[j]] for i in range(len(p)) for j in range(i + 2, len(p))))

def matrix_square(M):
    return [[sum((x * y for x, y in zip(row, col))) for col in zip(*M)] for row in M]
rows = []
for N in range(1, 33):
    M = N - 1
    d = '' if M == 0 else ''.join((str(1 - int(b)) for b in bin(M)[2:]))
    for strong in [False, True]:
        tail = '1' if strong else ''
        x = d + '110' + tail
        y = d + '101' + tail
        xx = d + '0110' + tail
        yy = d + '0101' + tail
        if strong and N == 1:
            x, y, xx, yy = ('0101', '0011', '00101', '00011')
            d = ''
        vals = [H(tournament_word(z)) for z in [x, y, xx, yy]]
        require(vals[1] - vals[0] == 2 * N and vals[3] - vals[2] == 4 * N, 'check at original line 47')
        require(C(xx) == (vals[2] - 1) // 2, 'check at original line 48')
        for old, new in [(x, xx), (y, yy)]:
            A = tournament_word(old)
            B = tournament_word(new)
            emb = [j if j < len(d) else j + 1 for j in range(len(A))]
            require(all((A[i][j] == B[emb[i]][emb[j]] for i in range(len(A)) for j in range(len(A)))), 'check at original line 52')
        rows.append(dict(N=N, strong=strong, words=[x, y, xx, yy], counts=vals))
words = 0
for k in range(8):
    for t in itertools.product('01', repeat=k):
        w = ''.join(t)
        require(H(tournament_word(w)) == 1 + 2 * C(w), 'check at original line 58')
        words += 1
nt = []
for n in range(3, 13):
    A = [[int(i < j) for j in range(n)] for i in range(n)]
    A[0][-1] = 0
    A[-1][0] = 1
    h = H(A)
    require(h == 1 + 2 ** (n - 2), 'check at original line 63')
    nt.append([n, h, (h - 1) // 2])
dcases = 0
for n in range(1, 5):
    edges = list(itertools.combinations(range(n), 2))
    for bits in itertools.product([0, 1], repeat=len(edges)):
        A = [[0] * n for _ in range(n)]
        for (i, j), b in zip(edges, bits):
            A[i][j] = b
            A[j][i] = 1 - b
        B = double(A)
        require(all((B[i][j] + B[j][i] == int(i != j) for i in range(2 * n) for j in range(2 * n))), 'check at original line 71')
        M = [[0 if i == j else 2 * A[i][j] - 1 for j in range(n)] for i in range(n)]
        MM = [[0 if i == j else 2 * B[i][j] - 1 for j in range(2 * n)] for i in range(2 * n)]
        sq = matrix_square(M)
        sq2 = matrix_square(MM)
        require(all((sq2[i][j] == (2 * sq[i % n][j % n] - int(i % n == j % n) if i // n == j // n else 0) for i in range(2 * n) for j in range(2 * n))), 'check at original line 75')
        require(H(A) % 2 == H(B) % 2 == 1, 'check at original line 76')
        dcases += 1
A = [[0]]
p = [0]
tower = []
for k in range(5):
    n = len(A)
    w = weight(A, p)
    tower.append(dict(n=n, weight=w, H=H(A) if n <= 8 else None))
    if k < 4:
        B = double(A)
        q = p + [n + i for i in reversed(p)]
        require(all((B[q[i]][q[i + 1]] for i in range(len(q) - 1))), 'check at original line 82')
        require(weight(B, q) == 2 * w + n - 1 + n * (n - 1) // 2, 'check at original line 83')
        A, p = (B, q)
blue = []
for n in range(2, 101):
    cells = lambda k: {(x, y) for x in range(1, k + 1) for y in range(1, x - 1)}
    sig = lambda k, t: (k + 1 - t[1], k + 1 - t[0])
    old, new = (cells(n), cells(n + 2))
    embedded = {(x + 1, y + 1) for x, y in old}
    border = new - embedded
    require(len(border) == 2 * n - 1, 'check at original line 91')
    require(all((sig(n + 2, t) in border for t in border)), 'check at original line 92')
    require(len([t for t in border if sig(n + 2, t) == t]) == 1, 'check at original line 93')
    require(all((sig(n + 2, (x + 1, y + 1)) == tuple((a + 1 for a in sig(n, (x, y)))) for x, y in old)), 'check at original line 94')
    f = sum((sig(n, t) == t for t in old))
    ff = sum((sig(n + 2, t) == t for t in new))
    b = (len(old) + f) // 2
    bb = (len(new) + ff) // 2
    require(ff - f == 1 and bb - b == n, 'check at original line 97')
    blue.append(dict(n=n, blue_dimension=b, new_dimension=bb, glue_dimension=f, new_glue_dimension=ff))
framed = 0
for n in range(1, 7):
    free = [(i, j) for i in range(n) for j in range(i + 2, n)]
    cells = [(i, j) for i in range(2 * n) for j in range(i + 2, 2 * n)]
    reflection = lambda t: (2 * n - 1 - t[1], 2 * n - 1 - t[0])
    path = list(range(n)) + list(range(2 * n - 1, n - 1, -1))
    c = {t: int(t[0] < n <= t[1] and path[t[0]] != path[t[1]] - n) for t in cells}
    require(all((c[t] == c[reflection(t)] for t in cells)), 'defect is blue')
    require(any(c.values()) == (n >= 2), 'nonblue guard')
    images = set()
    t0 = None
    for bits in itertools.product([0, 1], repeat=len(free)):
        A = [[int(i < j) for j in range(n)] for i in range(n)]
        for (i, j), b in zip(free, bits):
            A[i][j] = b
            A[j][i] = 1 - b
        D = double(A)
        t = {z: D[path[z[0]]][path[z[1]]] for z in cells}
        if t0 is None:
            t0 = t
        require(all((t[z] ^ t[reflection(z)] == c[z] for z in cells)), 'D image defect')
        require(all((t[z] ^ t0[z] == t[reflection(z)] ^ t0[reflection(z)] for z in cells)), 'representative coset')
        require(all((1 - t[z] ^ 1 - t[reflection(z)] == c[z] for z in cells)), 'complement gauge defect')
        images.add(tuple((t[z] for z in cells)))
        framed += 1
    require(len(images) == 2 ** len(free), 'D tiling injective')
out = dict(framed_defect_cases=framed, active_checks=CHECKS, one_apex_words=words, doubling_tournaments=dcases, NT=nt, tower=tower, blue=blue, selector_examples=rows)
Path(__file__).with_name('RESULT.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({k: v for k, v in out.items() if k != 'selector_examples'}, indent=2))
print('selector checks', len(rows), 'PASS')
