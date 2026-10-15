"""Cut-vertex transfer chain for uniform Hamiltonian paths (HYP-9279 (Z)).

T strong, z a cut vertex: T - z has strong components D_1 => D_2 => ... => D_k (k >= 2).  A Hamiltonian path is
(A-path) z (B-path) with A = vertices before z; each of T[A], T[B] visits the blocks in order, each block contiguously.
Scanning blocks left to right with states pre (B not started, A continues), gap (A ended at the previous block, B not
started), mid (B started, A continues), post (A ended), block j contributes (h = #Ham paths, h+ = starting in N+(z),
h- = ending in N-(z), sums over X = A cap D_j):
   pre->pre  h(D)                     pre->gap  h-(D)                pre->mid  sum_{X != D} h(X) h+(D\\X)
   pre->post sum_{X != 0, D} h-(X) h+(D\\X)    gap->post h+(D)      mid->mid  sum_X h(X) h(D\\X)
   mid->post sum_{X != 0} h-(X) h(D\\X)        post->post h(D)
start = pre + gap (gap = empty A), accept = gap (empty B) + post.  For v in D_i the 'v in A' parts restrict X to X containing v
(pre->pre, pre->gap count fully; gap->post, post->post count 0).  mid->mid is symmetric under X <-> D\\X (fair-coin lemma).
Validated against the full subset DP (tbal-style) for small n; then used on random cut-vertex tournaments up to n ~ 40.
"""
import random, sys
from fractions import Fraction as Fr

def ham_counts(A, verts, zout):
    """For every subset X of verts (bitmask over local indices): h(X), h+(X) (start in zout), h-(X) (end not in zout)."""
    k = len(verts); N = 1 << k
    f = [[0] * k for _ in range(N)]      # f[S][v]: paths on S ending at v
    for v in range(k): f[1 << v][v] = 1
    for S in range(1, N):
        if S & (S - 1) == 0: continue
        for v in range(k):
            if S >> v & 1:
                R = S ^ (1 << v)
                f[S][v] = sum(f[R][u] for u in range(k) if R >> u & 1 and A[verts[u]][verts[v]])
    g = [[0] * k for _ in range(N)]      # g[S][v]: paths on S starting at v
    for v in range(k): g[1 << v][v] = 1
    for S in range(1, N):
        if S & (S - 1) == 0: continue
        for v in range(k):
            if S >> v & 1:
                R = S ^ (1 << v)
                g[S][v] = sum(g[R][u] for u in range(k) if R >> u & 1 and A[verts[v]][verts[u]])
    h = [sum(f[S]) if S else 1 for S in range(N)]
    hplus = [sum(g[S][v] for v in range(k) if zout[verts[v]]) if S else 1 for S in range(N)]
    hminus = [sum(f[S][v] for v in range(k) if not zout[verts[v]]) if S else 1 for S in range(N)]
    return h, hplus, hminus

def block_matrices(A, verts, zout, vloc=None):
    h, hp, hm = ham_counts(A, verts, zout); k = len(verts); full = (1 << k) - 1
    has = (lambda X: True) if vloc is None else (lambda X: X >> vloc & 1)
    M = {}
    M['pp'] = h[full] if vloc is None or True else 0
    M['pg'] = hm[full]
    M['pm'] = sum(h[X] * hp[full ^ X] for X in range(full) if has(X))
    M['pq'] = sum(hm[X] * hp[full ^ X] for X in range(1, full) if has(X))
    M['gq'] = hp[full] if vloc is None else 0
    M['mm'] = sum(h[X] * h[full ^ X] for X in range(full + 1) if has(X))
    M['mq'] = sum(hm[X] * h[full ^ X] for X in range(1, full + 1) if has(X))
    M['qq'] = h[full] if vloc is None else 0
    return M

def chain(mats):
    # state vector over (pre, gap, mid, post)
    pre, gap, mid, post = 1, 1, 0, 0
    for M in mats:
        pre, gap, mid, post = (pre * M['pp'],
                               pre * M['pg'],
                               pre * M['pm'] + mid * M['mm'],
                               pre * M['pq'] + gap * M['gq'] + mid * M['mq'] + post * M['qq'])
    return gap + post

def components(A, verts):
    """strong components of T[verts], in condensation order."""
    vs = list(verts); reach = {}
    for v in vs:
        seen = {v}; st = [v]
        while st:
            u = st.pop()
            for w in vs:
                if w not in seen and A[u][w]: seen.add(w); st.append(w)
        reach[v] = seen
    comps = []; done = set()
    for v in vs:
        if v in done: continue
        c = [u for u in vs if u in reach[v] and v in reach[u]]; done |= set(c); comps.append(c)
    comps.sort(key=lambda c: -len(reach[c[0]]))       # earlier components reach more vertices
    return comps

def probs_cut(A, z):
    n = len(A); zout = [A[z][x] for x in range(n)]
    blocks = components(A, [x for x in range(n) if x != z])
    base = [block_matrices(A, b, zout) for b in blocks]
    H = chain(base); res = {}
    for i, b in enumerate(blocks):
        for loc, v in enumerate(b):
            mats = base[:i] + [block_matrices(A, b, zout, loc)] + base[i + 1:]
            res[v] = Fr(chain(mats), H)            # P(v before z)
    return H, blocks, res

def brute(A):
    n = len(A); N = 1 << n; full = N - 1
    f = [[0] * n for _ in range(N)]; g = [[0] * n for _ in range(N)]
    for v in range(n): f[1 << v][v] = g[1 << v][v] = 1
    for S in range(1, N):
        if S & (S - 1) == 0: continue
        for v in range(n):
            if S >> v & 1:
                R = S ^ (1 << v)
                f[S][v] = sum(f[R][u] for u in range(n) if R >> u & 1 and A[u][v])
                g[S][v] = sum(g[R][u] for u in range(n) if R >> u & 1 and A[v][u])
    H = sum(f[full]); cnt = [[0] * n for _ in range(n)]
    for x in range(n):
        for S in range(N):
            if S >> x & 1 and f[S][x]:
                pr = f[S][x] * g[(full ^ S) | (1 << x)][x]
                for y in range(n):
                    if (full ^ S) >> y & 1: cnt[x][y] += pr
    return H, cnt

def random_strong_block(s, rng):
    while True:
        B = [[False] * s for _ in range(s)]
        for i in range(s):
            for j in range(i + 1, s):
                if rng.random() < 0.5: B[i][j] = True
                else: B[j][i] = True
        if s == 1 or len(components(B, range(s))) == 1: return B

def random_cut_tournament(sizes, rng, pz=0.5):
    """blocks D_1 => ... => D_k of the given sizes plus z (index 0); z's arcs random, forced strong."""
    n = 1 + sum(sizes); A = [[False] * n for _ in range(n)]; idx = []; c = 1
    for s in sizes:
        B = random_strong_block(s, rng); ids = list(range(c, c + s)); c += s; idx.append(ids)
        for i in range(s):
            for j in range(s):
                if B[i][j]: A[ids[i]][ids[j]] = True
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            for x in idx[a]:
                for y in idx[b]: A[x][y] = True
    while True:
        for x in range(1, n):
            out = rng.random() < pz
            A[0][x] = out; A[x][0] = not out
        if any(A[0][x] for x in idx[0]) and any(A[x][0] for x in idx[-1]): return A, idx

if __name__ == '__main__':
    rng = random.Random(20261014)
    # (1) validation against brute force
    bad = 0
    for trial in range(300):
        k = rng.randint(2, 4); sizes = [rng.choice([1, 1, 3, 4]) for _ in range(k)]
        if sum(sizes) + 1 > 12: continue
        A, idx = random_cut_tournament(sizes, rng)
        H, blocks, res = probs_cut(A, 0); Hb, cnt = brute(A)
        if H != Hb or any(res[v] != Fr(cnt[v][0], Hb) for v in res): bad += 1
    print("validation vs brute force (n <= 12, 300 random cut-vertex tournaments): mismatches", bad)
    # (2) I1 / I3 on random large instances
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    worst_i1 = Fr(1, 2); worst_i3 = Fr(1, 2); f1 = f3 = 0; tested1 = 0
    for trial in range(trials):
        k = rng.randint(2, 9); sizes = [rng.choice([1, 1, 1, 3, 3, 4, 5]) for _ in range(k)]
        pz = rng.choice([0.2, 0.5, 0.8])
        A, idx = random_cut_tournament(sizes, rng, pz)
        H, blocks, res = probs_cut(A, 0)
        bal = {v: min(p, 1 - p) for v, p in res.items()}
        best = max(bal.values()); worst_i3 = min(worst_i3, best)
        if best < Fr(1, 3): f3 += 1; print("I3 FAILS", sizes, pz)
        if len(blocks) >= 3:
            tested1 += 1
            inner = min(bal[v] for b in blocks[1:-1] for v in b); worst_i1 = min(worst_i1, inner)
            if inner < Fr(1, 3): f1 += 1; print("I1 FAILS", sizes, pz, float(inner))
    print(f"random cut-vertex tournaments: {trials}; I3 failures {f3} (worst best-partner balance {float(worst_i3):.4f}); "
          f"I1 tested {tested1}, failures {f1} (worst interior balance {float(worst_i1):.4f})")
