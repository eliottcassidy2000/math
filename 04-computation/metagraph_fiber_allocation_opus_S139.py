"""
metagraph_fiber_allocation_opus_S139.py

THE PER-CLASS FIBER ALLOCATION of the tiling model (owner directive; differentiated from
mac-mini-S46/THM-643 + klein-S161, who hold the LINE-side parity theorems / node-type
formulas / n<=7 line censuses — this script is the "which particular tilings map to which
class" layer).

STRICT DEFINITIONS (CLAUDE.md / tournament-tiling-explorer):
  * tiles (x,y), x >= y+2; m = C(n-1,2); base path n -> n-1 -> ... -> 1; bit 0 = x->y.
  * grid involution sigma: (x,y) -> (n-y+1, n-x+1); fixed tiles have x+y = n+1,
    f = floor((n-1)/2).  isGridSym(t): bit(sigma(tile)) = bit(tile) for every tile.
  * LINE = pair {t, flip(t)} (flip = complement every tile bit); BLUE iff t gridsym
    (equivalently both endpoints gridsym), else BLACK.
  * #gridsym tilings = 2^{(m+f)/2}; #blue lines = 2^{(m+f)/2 - 1}; #black = 2^{m-1} - blue.
    (Verified below.)
  * PURE BLUE class: all tilings gridsym; PURE BLACK: none; MIXED: both.

NEW OBJECTS (the fiber allocation):
  * N(C) = #tilings in class C  (sanity: N(C)*|Aut(C)| = H(C), LEM-003)
  * g(C) = #gridsym tilings in C, b(C) = N(C) - g(C)
  * the FLIP MULTI-MAP: for each class C, the multiset {class(flip(t)) : t in C} — the
    line-multigraph on classes with multiplicities; blue sub-multi-map from gridsym t.
  * cross-tabs: node type x transpose-self(T^op) x parities of (N, g, b).
  * OWNER'S PARITY PROBE: "blues contribute odd amounts and blacks even amounts" — tested
    against per-class g(C) / b(C) parities and per-class line-in-degree parities.

n = 4..6 exact by full canonicalization; n = 7 by numpy-vectorized canonical min over S_7.
"""
import sys, time
from math import comb
from itertools import permutations
from collections import Counter, defaultdict

try:
    import numpy as np
except ImportError:
    np = None

def tiles_of(n):
    return [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2]

def sigma_map(n, T):
    idx = {t: i for i, t in enumerate(T)}
    return [idx[(n - y + 1, n - x + 1)] for (x, y) in T]

def arcs_of_tiling(n, T, bits):
    """Return adjacency matrix a[i][j]=1 if i->j (players 1..n as 0..n-1)."""
    a = [[0] * n for _ in range(n)]
    for k in range(n, 1, -1):               # base path k -> k-1
        a[k - 1][k - 2] = 1
    for i, (x, y) in enumerate(T):
        if (bits >> i) & 1:
            a[y - 1][x - 1] = 1             # bit 1: y -> x (backward)
        else:
            a[x - 1][y - 1] = 1             # bit 0: x -> y (forward)
    return a

def pack(a, n):
    v = 0; idx = 0
    for i in range(n):
        for j in range(n):
            if i != j:
                if a[i][j]: v |= (1 << idx)
                idx += 1
    return v

def canon(a, n, perms):
    best = None
    for p in perms:
        v = 0; idx = 0
        for i in range(n):
            for j in range(n):
                if i != j:
                    if a[p[i]][p[j]]: v |= (1 << idx)
                    idx += 1
        if best is None or v < best: best = v
    return best

def census(n, verbose=True):
    T = tiles_of(n)
    m = len(T)
    assert m == comb(n - 1, 2)
    sig = sigma_map(n, T)
    f = sum(1 for i, s in enumerate(sig) if s == i)
    perms = list(permutations(range(n)))
    t0 = time.time()
    cls_of = {}
    canon_cache = {}
    for bits in range(1 << m):
        a = arcs_of_tiling(n, T, bits)
        c = canon(a, n, perms)
        cls_of[bits] = c
    # gridsym test
    def gridsym(bits):
        for i, s in enumerate(sig):
            if ((bits >> i) & 1) != ((bits >> s) & 1):
                return False
        return True
    G = [bits for bits in range(1 << m) if gridsym(bits)]
    gs = set(G)
    full = (1 << m) - 1
    # per-class data
    N = Counter(cls_of.values())
    g = Counter(cls_of[b] for b in G)
    flipmap = defaultdict(Counter)          # C -> Counter of class(flip(t)) over t in C
    bluemap = defaultdict(Counter)
    for bits, c in cls_of.items():
        c2 = cls_of[bits ^ full]
        flipmap[c][c2] += 1
        if bits in gs:
            bluemap[c][c2] += 1
    # transpose-self: T^op
    selfop = {}
    for bits, c in cls_of.items():
        if c in selfop: continue
        a = arcs_of_tiling(n, T, bits)
        aop = [[a[j][i] for j in range(n)] for i in range(n)]
        selfop[c] = (canon(aop, n, perms) == c)
    # node types
    types = {}
    for c in N:
        gc = g.get(c, 0)
        types[c] = "pureblue" if gc == N[c] else ("pureblack" if gc == 0 else "mixed")
    # ---- report ----
    nblue = len(G) // 2
    nblack = (1 << (m - 1)) - nblue
    print(f"\n===== n={n}: m={m}, f={f}; tilings {1<<m}; gridsym {len(G)} "
          f"(= 2^((m+f)/2) {'OK' if len(G) == 1 << ((m+f)//2) else '***'}); "
          f"blue lines {nblue}, black lines {nblack}; classes {len(N)} [{time.time()-t0:.0f}s]")
    tc = Counter(types.values())
    st = Counter((types[c], selfop[c]) for c in N)
    print(f"  node types: {dict(tc)}")
    print(f"  type x transpose-self: {dict(st)}")
    # parity probe
    par = Counter((types[c], N[c] % 2, g.get(c, 0) % 2) for c in N)
    print(f"  (type, N mod 2, g mod 2) distribution: {dict(par)}")
    # g-spectrum among mixed/pureblue
    gspec = Counter(g.get(c, 0) for c in N if types[c] != "pureblack")
    print(f"  g(C) spectrum on non-pure-black classes: {dict(sorted(gspec.items()))}")
    # flip multi-map structure: self-paired vs cross-paired lines at class level
    selfl = crossl = 0; blueself = bluecross = 0
    for c, ctr in flipmap.items():
        for c2, mult in ctr.items():
            if c == c2: selfl += mult
            else: crossl += mult
    for c, ctr in bluemap.items():
        for c2, mult in ctr.items():
            if c == c2: blueself += mult
            else: bluecross += mult
    # each line counted twice (once per endpoint) except... t and flip(t) both enumerated
    print(f"  line endpoints: class-self {selfl//2} lines, class-cross {crossl//2} lines "
          f"(blue: self {blueself//2}, cross {bluecross//2})")
    # the flip multi-map on classes: is class(flip(t)) CONSTANT on each class?
    nonconst = sum(1 for c, ctr in flipmap.items() if len(ctr) > 1)
    print(f"  classes where flip-partner class is NOT constant: {nonconst} / {len(N)}")
    # N*|Aut| = H sanity on a few classes is standard; skip (LEM-003 proved).
    return dict(n=n, N=N, g=g, types=types, selfop=selfop, flipmap=flipmap)

def census7_numpy():
    n = 7
    T = tiles_of(n); m = len(T)
    sig = sigma_map(n, T)
    f = sum(1 for i, s in enumerate(sig) if s == i)
    if np is None:
        print("numpy unavailable; skipping n=7"); return None
    t0 = time.time()
    NT = 1 << m
    bits = np.arange(NT, dtype=np.int64)
    tilebits = np.zeros((NT, m), dtype=bool)
    for i in range(m):
        tilebits[:, i] = (bits >> i) & 1
    # arc list: (i, j) directed pairs among 0..6; build per-tiling arc matrix as bit array
    # order arcs canonically: all ordered pairs (i, j), i != j
    pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
    pidx = {p: k for k, p in enumerate(pairs)}
    A = np.zeros((NT, len(pairs)), dtype=bool)
    for k in range(n, 1, -1):
        A[:, pidx[(k - 1, k - 2)]] = True
    for i, (x, y) in enumerate(T):
        A[tilebits[:, i], pidx[(y - 1, x - 1)]] = True
        A[~tilebits[:, i], pidx[(x - 1, y - 1)]] = True
    # canonical form: min over perms of packed big-int — pack 42 arc bits into int64
    weights = (1 << np.arange(len(pairs), dtype=np.int64))
    best = None
    for p in permutations(range(n)):
        perm_cols = np.array([pidx[(p[i], p[j])] for (i, j) in pairs], dtype=np.int64)
        v = (A[:, perm_cols] * weights).sum(axis=1)
        best = v if best is None else np.minimum(best, v)
    canon_arr = best
    grids = np.ones(NT, dtype=bool)
    for i, s in enumerate(sig):
        if s != i:
            grids &= (tilebits[:, i] == tilebits[:, s])
    full = NT - 1
    flip_arr = canon_arr[full - bits]      # class of flip(t) (flip = full XOR bits = full - bits since complement)
    print(f"\n===== n=7 (numpy): m={m}, f={f}; tilings {NT}; gridsym {int(grids.sum())} "
          f"(= 2^((m+f)/2) {'OK' if int(grids.sum()) == 1 << ((m+f)//2) else '***'}) "
          f"[{time.time()-t0:.0f}s]")
    # per-class tallies
    classes, inv = np.unique(canon_arr, return_inverse=True)
    Nc = np.bincount(inv)
    gc = np.bincount(inv, weights=grids.astype(np.int64)).astype(np.int64)
    types = np.where(gc == Nc, 0, np.where(gc == 0, 2, 1))   # 0 pureblue 1 mixed 2 pureblack
    tc = Counter(["pureblue", "mixed", "pureblack"][t] for t in types)
    print(f"  classes {len(classes)}; node types: {dict(tc)}")
    gspec = Counter(int(x) for x in gc[types != 2])
    print(f"  g(C) spectrum on non-pure-black: {dict(sorted(gspec.items()))}")
    par = Counter((["pb","mx","pk"][types[i]], int(Nc[i]) % 2, int(gc[i]) % 2) for i in range(len(classes)))
    print(f"  (type, N mod 2, g mod 2): {dict(par)}")
    # flip-map constancy + self/cross lines
    cls_flip = flip_arr
    nonconst = 0; selfl = 0; crossl = 0
    order = np.argsort(canon_arr, kind="stable")
    # count lines: pair t with full-t; enumerate t < full-t
    tmask = bits < (full - bits)
    same = canon_arr[tmask] == cls_flip[tmask]
    selfl = int(same.sum()); crossl = int((~same).sum())
    blue_t = grids & tmask
    bsame = canon_arr[blue_t] == cls_flip[blue_t]
    print(f"  lines: class-self {selfl}, class-cross {crossl} (blue: self {int(bsame.sum())}, cross {int((~bsame).sum())})")
    # flip-partner constancy per class
    dd = defaultdict(set)
    ca = canon_arr; fa = cls_flip
    for t in range(NT):
        dd[int(ca[t])].add(int(fa[t]))
    nonconst = sum(1 for s in dd.values() if len(s) > 1)
    print(f"  classes where flip-partner class NOT constant: {nonconst} / {len(classes)}")
    return None

if __name__ == "__main__":
    for n in (4, 5, 6):
        census(n)
    census7_numpy()
