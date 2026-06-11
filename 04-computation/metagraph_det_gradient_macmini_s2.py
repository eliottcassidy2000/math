#!/usr/bin/env python3
"""
metagraph_det_gradient_macmini_s2.py — macmini-2026-06-10-S2

HYP-2386: Attach the determinant invariant d(T) = det(I+S)/2^(n-1)
(S = A - A^T, skew +-1 matrix) to the nodes of the iso-class metagraph G_n
for n = 5, 6, 7 and measure:

  1. Edge-orientation stats: d-gradient vs the known H-gradient
  2. Correlation d vs H (Pearson, Spearman, log-log, within-score-multiset)
  3. Second-eigenvector test (is d, like H, ~ the 2nd adjacency eigenvector?)
     on G_n and on the merged G_n/Z_2
  4. SC/NS decomposition (spine / ribs / sea edge |delta d|)
  5. Location of argmax(d) vs argmax(H): BFS distances, distance from transitive

Construction reuses:
  - 04-computation/iso_class_graph_fast.py (arc-flip metagraph via canonical
    forms; here applied to the tiling enumeration, which hits every iso class
    since every tournament has a Hamiltonian path)
  - 04-computation/fast_metagraph_n7.py (tiling model + rich invariant
    (scores, c3, H, sorted vertex-deleted H) which separates all 456 classes
    at n=7)

Two edge sets are built and compared:
  - ARC-FLIP graph: flip each of the C(n,2) arcs of one representative per
    class (exact: arc flips commute with isomorphism, so one rep suffices).
    This is THE metagraph G_n of iso_class_graph_fast.py.
  - TILE-FLIP (wiggly, d=1) graph: flip each of the C(n-1,2) tiles over ALL
    tilings (this is what fast_metagraph_n7.py builds). Tile flips are arc
    flips, so tile-flip edges are a subset of arc-flip edges.

All determinants are EXACT (fraction-free Bareiss on Python ints).
Verification: det(I+S) = 2^(n-1) * d with d a positive odd integer for every
class; E[det(I+S)] over all 2^C(n,2) labeled tournaments equals the involution
number I(n) (class sizes from orbit-stabilizer: |class| = n! * #tilings / H).
"""

import sys
import time
import math
from itertools import permutations, product as iproduct
from collections import defaultdict, Counter, deque

import numpy as np

INVOLUTION = {3: 4, 4: 10, 5: 26, 6: 76, 7: 232}
NCLASS_EXPECT = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456}
SC_EXPECT = {5: 8, 6: 12, 7: 88}   # from (V - SC)/2 = 2, 22, 184 (CLAUDE.md)
H_LEVEL_CANON = {5: 1, 6: 15, 7: 136}  # canon: same-H cross-class edges


def sgn(z):
    return (z > 0) - (z < 0)


# ---------------------------------------------------------------- tilings
def make_tiles(n):
    tiles = []
    for y in range(1, n - 1):
        for x in range(n, y + 1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles


def tiling_to_adj(n, tiles, bits):
    A = [[0] * n for _ in range(n)]
    for k in range(n - 1):
        A[k + 1][k] = 1
    for i, (x, y) in enumerate(tiles):
        a, b = x - 1, y - 1
        if bits & (1 << i):
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A


def transpose(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]


def out_bits(A, n):
    return [sum(1 << j for j in range(n) if A[i][j]) for i in range(n)]


# ---------------------------------------------------------------- H count
def count_hp(out, n):
    """Hamiltonian path count, DP over subsets (out = out-neighbor bitmasks)."""
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(full + 1)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, full):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c:
                continue
            avail = out[v] & ~mask & full
            while avail:
                ub = avail & (-avail)
                avail ^= ub
                dp[mask | ub][ub.bit_length() - 1] += c
    return sum(dp[full])


# ---------------------------------------------------------------- exact det
def det_bareiss(Min):
    """Exact integer determinant, fraction-free Bareiss."""
    M = [list(map(int, row)) for row in Min]
    nn = len(M)
    sign = 1
    prev = 1
    for k in range(nn - 1):
        if M[k][k] == 0:
            piv = None
            for r in range(k + 1, nn):
                if M[r][k] != 0:
                    piv = r
                    break
            if piv is None:
                return 0
            M[k], M[piv] = M[piv], M[k]
            sign = -sign
        pk = M[k][k]
        Mk = M[k]
        for i in range(k + 1, nn):
            Mi = M[i]
            mik = Mi[k]
            for j in range(k + 1, nn):
                Mi[j] = (Mi[j] * pk - mik * Mk[j]) // prev
        prev = pk
    return sign * M[-1][-1]


def det_IpS(A, n):
    """Exact det(I+S), S = A - A^T (off-diag +1 if i->j else -1)."""
    M = [[1 if i == j else (1 if A[i][j] else -1) for j in range(n)]
         for i in range(n)]
    return det_bareiss(M)


# ---------------------------------------------------------------- stats
def pearson(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.std() == 0 or y.std() == 0:
        return float('nan')
    return float(np.corrcoef(x, y)[0, 1])


def ranks(x):
    order = sorted(range(len(x)), key=lambda i: x[i])
    r = [0.0] * len(x)
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]:
            j += 1
        avg = (i + j) / 2 + 1
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def spearman(x, y):
    return pearson(ranks(list(x)), ranks(list(y)))


def eta2(groups, x):
    """Between-group variance share (1 - SSwithin/SStotal)."""
    x = np.asarray(x, dtype=float)
    tot = float(((x - x.mean()) ** 2).sum())
    if tot == 0:
        return float('nan')
    ssw = 0.0
    for g in groups.values():
        xs = x[list(g)]
        ssw += float(((xs - xs.mean()) ** 2).sum())
    return 1 - ssw / tot


def pooled_within_corr(groups, x, y, use_ranks=False):
    xs, ys = [], []
    for g in groups.values():
        if len(g) < 2:
            continue
        gx = [x[i] for i in g]
        gy = [y[i] for i in g]
        if use_ranks:
            gx = ranks(gx)
            gy = ranks(gy)
        mx = sum(gx) / len(gx)
        my = sum(gy) / len(gy)
        xs.extend(v - mx for v in gx)
        ys.extend(v - my for v in gy)
    if len(xs) < 3:
        return float('nan')
    return pearson(xs, ys)


def regress_resid(y, v):
    """Residual of y after least-squares fit on [1, v]."""
    v = np.asarray(v, dtype=float)
    y = np.asarray(y, dtype=float)
    X = np.vstack([np.ones(len(v)), v]).T
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def bfs_dist(adjlist, sources):
    dist = {}
    dq = deque()
    for s in sources:
        dist[s] = 0
        dq.append(s)
    while dq:
        u = dq.popleft()
        for v in adjlist[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                dq.append(v)
    return dist


def eig_second(edges, nnodes, named_vecs):
    """Adjacency spectrum stats: R^2 of each named vector vs the 2nd
    (algebraically second-largest) eigenvector, plus best-matching eigvec."""
    A = np.zeros((nnodes, nnodes))
    for a, b in edges:
        A[a, b] = A[b, a] = 1.0
    w, V = np.linalg.eigh(A)
    v2 = V[:, -2]
    out = {'lam0': float(w[-1]), 'lam1': float(w[-2]), 'lam_min': float(w[0]),
           'v2': v2}
    for name, x in named_vecs.items():
        r = pearson(x, v2)
        best_k, best_r = -1, -1.0
        for k in range(len(w)):
            rr = abs(pearson(x, V[:, k]))
            if rr > best_r:
                best_r, best_k = rr, k
        out[name] = {'r': r, 'R2': r * r,
                     'best_idx_from_top': len(w) - 1 - best_k,
                     'best_lam': float(w[best_k]), 'best_R2': best_r ** 2}
    return out


# ---------------------------------------------------------------- build
def build(n):
    tiles = make_tiles(n)
    m = len(tiles)
    ntil = 1 << m
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    delh_cache = {}

    def del_h(A, drop):
        idxs = [u for u in range(n) if u != drop]
        key = 0
        b = 0
        for ii in range(n - 1):
            ri = A[idxs[ii]]
            for jj in range(ii + 1, n - 1):
                if ri[idxs[jj]]:
                    key |= (1 << b)
                b += 1
        v = delh_cache.get(key)
        if v is None:
            B = [[A[idxs[i]][idxs[j]] for j in range(n - 1)]
                 for i in range(n - 1)]
            v = count_hp(out_bits(B, n - 1), n - 1)
            delh_cache[key] = v
        return v

    # EXACT canonical form with vertex-invariant pruning.
    # NOTE: the rich invariant of fast_metagraph_n7.py
    # (scores, c3, H, sorted vertex-deleted H) is NOT complete at n=7:
    # it yields only 388 groups instead of 456 (and c3 is determined by
    # the score multiset, so it adds nothing). We therefore canonicalize
    # exactly: minimize the adjacency string over all labelings in which
    # vertices are sorted by the vertex key (score, H(T-v)); ties only
    # allow permutations within equal-key blocks. This min is a genuine
    # adjacency string of the class, hence a COMPLETE iso invariant.
    def key_of(A):
        keys = []
        for v in range(n):
            keys.append((sum(A[v]), del_h(A, v)))
        order = sorted(range(n), key=lambda v: keys[v])
        blocks = []
        i = 0
        while i < n:
            j = i
            while j + 1 < n and keys[order[j + 1]] == keys[order[i]]:
                j += 1
            blocks.append(order[i:j + 1])
            i = j + 1
        pools = [list(permutations(b)) for b in blocks]
        best = None
        for combo in iproduct(*pools):
            p = []
            for c in combo:
                p.extend(c)
            form = tuple(A[p[i]][p[j]] for (i, j) in pairs)
            if best is None or form < best:
                best = form
        return best

    # --- classify all tilings ---
    t0 = time.time()
    key2id = {}
    reps = []       # representative adjacency matrix per class
    ngroup = []     # tilings per class
    cm = [0] * ntil  # tiling bits -> class id
    for bits in range(ntil):
        A = tiling_to_adj(n, tiles, bits)
        k = key_of(A)
        cid = key2id.get(k)
        if cid is None:
            cid = len(reps)
            key2id[k] = cid
            reps.append(A)
            ngroup.append(0)
        ngroup[cid] += 1
        cm[bits] = cid
    nc = len(reps)
    t_classify = time.time() - t0

    # --- per-class invariants ---
    H = [count_hp(out_bits(reps[c], n), n) for c in range(nc)]
    dets = [det_IpS(reps[c], n) for c in range(nc)]
    dval = []
    n_odd_d = 0
    for c in range(nc):
        q, r = divmod(dets[c], 1 << (n - 1))
        assert r == 0 and q > 0, f"det divisibility fails: class {c}, det={dets[c]}"
        if q % 2 == 1:
            n_odd_d += 1
        dval.append(q)
    score = [tuple(sorted(sum(row) for row in reps[c])) for c in range(nc)]
    rev = []
    for c in range(nc):
        kT = key_of(transpose(reps[c]))
        assert kT in key2id, f"reversal of class {c} not classified"
        rev.append(key2id[kT])
    sc = [rev[c] == c for c in range(nc)]
    for c in range(nc):
        assert dval[rev[c]] == dval[c], "d not reversal-invariant!"
        assert H[rev[c]] == H[c], "H not reversal-invariant!"

    # --- class sizes via orbit-stabilizer: |class| = n! * #tilings / H ---
    nfact = math.factorial(n)
    sizes = []
    for c in range(nc):
        num = nfact * ngroup[c]
        assert num % H[c] == 0, f"orbit-stabilizer fails at class {c}"
        sizes.append(num // H[c])
    total_T = 1 << (n * (n - 1) // 2)
    assert sum(sizes) == total_T, "class sizes do not sum to 2^C(n,2)"
    # E[det(I+S)] exact check vs involution number
    sum_det = sum(sizes[c] * dets[c] for c in range(nc))
    Edet_ok = (sum_det == INVOLUTION[n] * total_T)

    # --- ARC-FLIP edges (one rep per class is exact) ---
    t1 = time.time()
    edges_arc = set()
    for cid in range(nc):
        A = [row[:] for row in reps[cid]]
        for (i, j) in pairs:
            A[i][j], A[j][i] = A[j][i], A[i][j]
            cj = key2id[key_of(A)]
            if cj != cid:
                e = (cid, cj) if cid < cj else (cj, cid)
                edges_arc.add(e)
            A[i][j], A[j][i] = A[j][i], A[i][j]
    t_arc = time.time() - t1

    # --- TILE-FLIP edges + oriented flip events (0->1 = away from transitive) ---
    t2 = time.time()
    edges_tile = set()
    ev_H = Counter()
    ev_d = Counter()
    ev_joint = Counter()
    edge_sH = defaultdict(set)   # signs of dH along 0->1, per edge
    edge_sd = defaultdict(set)
    n_events = 0
    for bits in range(ntil):
        c1 = cm[bits]
        for t in range(m):
            if bits & (1 << t):
                continue
            c2 = cm[bits | (1 << t)]
            if c1 == c2:
                continue
            e = (c1, c2) if c1 < c2 else (c2, c1)
            edges_tile.add(e)
            sH = sgn(H[c2] - H[c1])
            sd = sgn(dval[c2] - dval[c1])
            ev_H[sH] += 1
            ev_d[sd] += 1
            ev_joint[(sH, sd)] += 1
            edge_sH[e].add(sH)
            edge_sd[e].add(sd)
            n_events += 1
    t_tile = time.time() - t2

    return dict(n=n, nc=nc, m=m, ntil=ntil, reps=reps, H=H, d=dval, dets=dets,
                score=score, rev=rev, sc=sc, sizes=sizes, ngroup=ngroup,
                edges_arc=edges_arc, edges_tile=edges_tile,
                ev_H=ev_H, ev_d=ev_d, ev_joint=ev_joint, n_events=n_events,
                edge_sH=edge_sH, edge_sd=edge_sd,
                n_odd_d=n_odd_d, Edet_ok=Edet_ok,
                sum_det=sum_det, total_T=total_T,
                t_classify=t_classify, t_arc=t_arc, t_tile=t_tile)


# ---------------------------------------------------------------- analysis
def edge_dir_class(sign_set):
    """Classify an edge by its observed 0->1 flip signs."""
    if sign_set == {0}:
        return 'level'
    if 1 in sign_set and -1 not in sign_set:
        return 'always-up'
    if -1 in sign_set and 1 not in sign_set:
        return 'always-down'
    return 'mixed'


def analyze(B):
    n = B['n']
    nc = B['nc']
    H = B['H']
    d = B['d']
    score = B['score']
    sc = B['sc']
    edges_arc = sorted(B['edges_arc'])
    edges_tile = B['edges_tile']

    print("\n" + "#" * 76)
    print(f"#  n = {n}   ({nc} iso classes)")
    print("#" * 76)

    # ---------- build verification ----------
    nsc = sum(sc)
    print(f"\n[build] classes = {nc} (expect {NCLASS_EXPECT[n]}): "
          f"{'OK' if nc == NCLASS_EXPECT[n] else '*** MISMATCH ***'}")
    print(f"[build] SC classes = {nsc} (expect {SC_EXPECT[n]}): "
          f"{'OK' if nsc == SC_EXPECT[n] else '*** MISMATCH ***'}")
    print(f"[build] arc-flip edges = {len(edges_arc)}, "
          f"tile-flip (wiggly d=1) edges = {len(edges_tile)}, "
          f"tile subset of arc: {edges_tile <= set(edges_arc)}, "
          f"arc-only edges = {len(set(edges_arc) - edges_tile)}")
    print(f"[verify] det(I+S) = 2^{n-1} * d with d a positive integer: OK "
          f"(all {nc} classes); d parity: {B['n_odd_d']} odd, "
          f"{nc - B['n_odd_d']} even")
    print(f"[verify] E[det(I+S)] = {B['sum_det']}/{B['total_T']} "
          f"= {B['sum_det'] // math.gcd(B['sum_det'], B['total_T'])}"
          f"/{B['total_T'] // math.gcd(B['sum_det'], B['total_T'])} "
          f"vs involution I({n}) = {INVOLUTION[n]}: "
          f"{'OK (exact)' if B['Edet_ok'] else '*** MISMATCH ***'}")
    dd = sorted(set(d))
    print(f"[d] distinct d values = {len(dd)}, min = {dd[0]}, max = {dd[-1]}")
    if len(dd) <= 25:
        cnt = Counter(d)
        print(f"[d] value -> #classes: " +
              ", ".join(f"{v}:{cnt[v]}" for v in dd))
    n_d1 = sum(1 for v in d if v == 1)
    print(f"[d] d=1 (locally transitive) classes: {n_d1}")
    print(f"[timing] classify {B['t_classify']:.1f}s, arc {B['t_arc']:.1f}s, "
          f"tile {B['t_tile']:.1f}s")

    # full node table for small n
    if n <= 5:
        print(f"\n  {'cid':>3} {'score':>15} {'H':>4} {'d':>4} {'det':>5} SC")
        for c in sorted(range(nc), key=lambda c: (H[c], d[c])):
            print(f"  {c:>3} {str(score[c]):>15} {H[c]:>4} {d[c]:>4} "
                  f"{B['dets'][c]:>5} {'Y' if sc[c] else ''}")

    # ---------------- [1] edge orientation -------------------
    print(f"\n[1] EDGE ORIENTATION (is there a d-gradient like the H-gradient?)")
    EA = len(edges_arc)
    h_lvl = sum(1 for (a, b) in edges_arc if H[a] == H[b])
    d_lvl = sum(1 for (a, b) in edges_arc if d[a] == d[b])
    both_lvl = sum(1 for (a, b) in edges_arc if H[a] == H[b] and d[a] == d[b])
    print(f"  ARC graph: {EA} edges")
    print(f"    H-level edges: {h_lvl} ({h_lvl/EA:.4f})   "
          f"[canon H-level for n={n}: {H_LEVEL_CANON[n]}]")
    print(f"    d-level edges: {d_lvl} ({d_lvl/EA:.4f})")
    print(f"    both-level   : {both_lvl}")
    et = len(edges_tile)
    h_lvl_t = sum(1 for (a, b) in edges_tile if H[a] == H[b])
    d_lvl_t = sum(1 for (a, b) in edges_tile if d[a] == d[b])
    print(f"  TILE graph: {et} edges; H-level {h_lvl_t} ({h_lvl_t/et:.4f}), "
          f"d-level {d_lvl_t} ({d_lvl_t/et:.4f})")
    # orient by H, ask d
    up = lvl = dn = 0
    for (a, b) in edges_arc:
        if H[a] == H[b]:
            continue
        lo, hi = (a, b) if H[a] < H[b] else (b, a)
        s = sgn(d[hi] - d[lo])
        if s > 0:
            up += 1
        elif s == 0:
            lvl += 1
        else:
            dn += 1
    tot = up + lvl + dn
    print(f"  Orient arc edges low-H -> high-H ({tot} non-H-level edges):")
    print(f"    d increases {up} ({up/tot:.4f}), level {lvl} ({lvl/tot:.4f}), "
          f"decreases {dn} ({dn/tot:.4f})")
    # orient by d, ask H
    up2 = lvl2 = dn2 = 0
    for (a, b) in edges_arc:
        if d[a] == d[b]:
            continue
        lo, hi = (a, b) if d[a] < d[b] else (b, a)
        s = sgn(H[hi] - H[lo])
        if s > 0:
            up2 += 1
        elif s == 0:
            lvl2 += 1
        else:
            dn2 += 1
    tot2 = up2 + lvl2 + dn2
    print(f"  Orient arc edges low-d -> high-d ({tot2} non-d-level edges):")
    print(f"    H increases {up2} ({up2/tot2:.4f}), level {lvl2} "
          f"({lvl2/tot2:.4f}), decreases {dn2} ({dn2/tot2:.4f})")
    nboth = sum(1 for (a, b) in edges_arc if H[a] != H[b] and d[a] != d[b])
    agree = sum(1 for (a, b) in edges_arc
                if H[a] != H[b] and d[a] != d[b]
                and sgn(H[a] - H[b]) == sgn(d[a] - d[b]))
    print(f"  Sign agreement (edges with both dH,dd != 0): {agree}/{nboth} "
          f"= {agree/nboth:.4f}")
    # flip events, 0->1 = away from transitive
    ne = B['n_events']
    evH, evd, evj = B['ev_H'], B['ev_d'], B['ev_joint']
    print(f"  FLIP EVENTS (tile 0->1 = away from transitive; {ne} cross-class events):")
    print(f"    H: up {evH[1]} ({evH[1]/ne:.4f}), level {evH[0]} "
          f"({evH[0]/ne:.4f}), down {evH[-1]} ({evH[-1]/ne:.4f})")
    print(f"    d: up {evd[1]} ({evd[1]/ne:.4f}), level {evd[0]} "
          f"({evd[0]/ne:.4f}), down {evd[-1]} ({evd[-1]/ne:.4f})")
    ja = evj[(1, 1)] + evj[(-1, -1)]
    jb = evj[(1, -1)] + evj[(-1, 1)]
    print(f"    joint (both nonzero): same-sign {ja}, opposite {jb}, "
          f"agreement {ja/(ja+jb):.4f}")
    # per-edge flip direction classes
    cH = Counter(edge_dir_class(s) for s in B['edge_sH'].values())
    cd = Counter(edge_dir_class(s) for s in B['edge_sd'].values())
    print(f"  Per-edge 0->1 direction (tile graph, {et} edges):")
    print(f"    H: always-up {cH['always-up']}, always-down {cH['always-down']}, "
          f"mixed {cH['mixed']}, level {cH['level']}")
    print(f"    d: always-up {cd['always-up']}, always-down {cd['always-down']}, "
          f"mixed {cd['mixed']}, level {cd['level']}")

    # ---------------- [2] correlation -------------------
    print(f"\n[2] CORRELATION d vs H (across {nc} nodes)")
    logd = [math.log(v) for v in d]
    logH = [math.log(v) for v in H]
    print(f"  Pearson(d, H)        = {pearson(d, H):+.4f}")
    print(f"  Spearman(d, H)       = {spearman(d, H):+.4f}")
    print(f"  Pearson(log d, log H)= {pearson(logd, logH):+.4f}")
    groups = defaultdict(list)
    for c in range(nc):
        groups[score[c]].append(c)
    g2 = {k: v for k, v in groups.items() if len(v) >= 2}
    print(f"  Score multisets: {len(groups)} total, {len(g2)} with >=2 classes")
    print(f"  Between-score-group variance share (eta^2): "
          f"d {eta2(groups, d):.4f}, H {eta2(groups, H):.4f}, "
          f"log d {eta2(groups, logd):.4f}, log H {eta2(groups, logH):.4f}")
    pw = pooled_within_corr(groups, d, H)
    pwr = pooled_within_corr(groups, d, H, use_ranks=True)
    print(f"  Pooled WITHIN-score-group Pearson(d,H)  = {pw:+.4f}")
    print(f"  Pooled WITHIN-score-group Spearman(d,H) = {pwr:+.4f}")
    big = sorted((k for k in g2 if len(g2[k]) >= 3),
                 key=lambda k: -len(g2[k]))
    show = [k for k in big if len(g2[k]) >= (10 if n == 7 else 3)]
    if show:
        print(f"  Within-group r(d,H) for groups with >= "
              f"{10 if n == 7 else 3} classes:")
        for k in show:
            g = g2[k]
            r = pearson([d[i] for i in g], [H[i] for i in g])
            rs = spearman([d[i] for i in g], [H[i] for i in g])
            print(f"    score {k}: {len(g)} classes, "
                  f"Pearson {r:+.4f}, Spearman {rs:+.4f}")

    # ---------------- [3] second eigenvector -------------------
    print(f"\n[3] SECOND-EIGENVECTOR TEST (adjacency of arc-flip G_{n})")
    es = eig_second(edges_arc, nc, {'H': H, 'd': d, 'logd': logd})
    print(f"  spectrum: lam0 = {es['lam0']:.4f}, lam1 = {es['lam1']:.4f}, "
          f"lam_min = {es['lam_min']:.4f}")
    for name in ('H', 'd', 'logd'):
        e = es[name]
        print(f"  {name:>4} vs v2: r = {e['r']:+.4f}, R^2 = {e['R2']:.4f} | "
              f"best-matching eigvec: #{e['best_idx_from_top']} from top "
              f"(lam = {e['best_lam']:.4f}), R^2 = {e['best_R2']:.4f}")
    v2 = es['v2']
    rd = regress_resid(d, v2)
    rh = regress_resid(H, v2)
    print(f"  corr(resid_d, resid_H) after projecting out v2: "
          f"{pearson(rd, rh):+.4f}   [raw corr(d,H) = {pearson(d, H):+.4f}]")
    rld = regress_resid(logd, v2)
    print(f"  corr(resid_logd, resid_H) = {pearson(rld, rh):+.4f}")

    # merged graph G_n/Z_2
    rev = B['rev']
    pair_id = {}
    merged_nodes = []
    for c in range(nc):
        p = min(c, rev[c])
        if p not in pair_id:
            pair_id[p] = len(merged_nodes)
            merged_nodes.append(p)
    mid = [pair_id[min(c, rev[c])] for c in range(nc)]
    medges = set()
    for (a, b) in edges_arc:
        ma, mb = mid[a], mid[b]
        if ma != mb:
            medges.add((min(ma, mb), max(ma, mb)))
    mH = [H[p] for p in merged_nodes]
    md = [d[p] for p in merged_nodes]
    mlogd = [math.log(v) for v in md]
    print(f"  MERGED G_{n}/Z_2: V = {len(merged_nodes)}, E = {len(medges)}")
    esm = eig_second(sorted(medges), len(merged_nodes),
                     {'H': mH, 'd': md, 'logd': mlogd})
    print(f"  merged spectrum: lam0 = {esm['lam0']:.4f}, "
          f"lam1 = {esm['lam1']:.4f}, lam_min = {esm['lam_min']:.4f}")
    for name in ('H', 'd', 'logd'):
        e = esm[name]
        print(f"  merged {name:>4} vs v2: r = {e['r']:+.4f}, "
              f"R^2 = {e['R2']:.4f} | best eigvec #{e['best_idx_from_top']} "
              f"from top, R^2 = {e['best_R2']:.4f}")
    print(f"  [S268 reference for merged H R^2: n=5: 0.723, n=6: 0.786, "
          f"n=7: 0.733]")
    v2m = esm['v2']
    rdm = regress_resid(md, v2m)
    rhm = regress_resid(mH, v2m)
    print(f"  merged corr(resid_d, resid_H) = {pearson(rdm, rhm):+.4f}   "
          f"[raw merged corr(d,H) = {pearson(md, mH):+.4f}]")

    # ---------------- [4] SC/NS decomposition -------------------
    print(f"\n[4] SC/NS DECOMPOSITION (spine / ribs / sea)")
    dsc = [d[c] for c in range(nc) if sc[c]]
    dns = [d[c] for c in range(nc) if not sc[c]]
    print(f"  SC nodes ({len(dsc)}): mean d = {np.mean(dsc):.4f}, "
          f"min = {min(dsc)}, max = {max(dsc)}")
    print(f"  NS nodes ({len(dns)}): mean d = {np.mean(dns):.4f}, "
          f"min = {min(dns)}, max = {max(dns)}")
    hsc = [H[c] for c in range(nc) if sc[c]]
    hns = [H[c] for c in range(nc) if not sc[c]]
    print(f"  (H for context: SC mean {np.mean(hsc):.2f} "
          f"[{min(hsc)},{max(hsc)}], NS mean {np.mean(hns):.2f} "
          f"[{min(hns)},{max(hns)}])")
    etype = defaultdict(list)
    for (a, b) in edges_arc:
        t = ('SC-SC' if sc[a] and sc[b] else
             'NS-NS' if not sc[a] and not sc[b] else 'SC-NS')
        etype[t].append((a, b))
    for t in ('SC-SC', 'SC-NS', 'NS-NS'):
        es_ = etype[t]
        if not es_:
            print(f"  {t} (0 edges)")
            continue
        add = [abs(d[a] - d[b]) for (a, b) in es_]
        adh = [abs(H[a] - H[b]) for (a, b) in es_]
        print(f"  {t} ({len(es_)} edges): mean|dd| = {np.mean(add):.4f} "
              f"(max {max(add)}), mean|dH| = {np.mean(adh):.4f} "
              f"(max {max(adh)})")

    # ---------------- [5] argmax location -------------------
    print(f"\n[5] WHERE d-MAXIMIZERS AND H-MAXIMIZERS SIT")
    adjl = defaultdict(list)
    for (a, b) in edges_arc:
        adjl[a].append(b)
        adjl[b].append(a)
    trans = [c for c in range(nc) if H[c] == 1]
    assert len(trans) == 1
    t_id = trans[0]
    dist_t = bfs_dist(adjl, [t_id])
    dmax = max(d)
    hmax = max(H)
    Dset = [c for c in range(nc) if d[c] == dmax]
    Hset = [c for c in range(nc) if H[c] == hmax]
    print(f"  transitive class: cid {t_id} (H=1, d={d[t_id]})")
    print(f"  d-max = {dmax}: {len(Dset)} class(es)")
    for c in Dset:
        print(f"    cid {c}: H = {H[c]}, score {score[c]}, "
              f"SC = {'Y' if sc[c] else 'N'}, dist-from-transitive = "
              f"{dist_t.get(c, 'inf')}")
    detmax = dmax << (n - 1)
    if n % 2 == 1:
        ceil = (n + 1) ** ((n - 1) // 2)
        print(f"    [odd-n ceiling on det(I+S): (n+1)^((n-1)/2) = {ceil}; "
              f"max det = {detmax}; attained iff DRT switching class: "
              f"{'ATTAINED' if detmax == ceil else 'not attained'}]")
    else:
        ceil = n ** (n // 2)
        print(f"    [even-n ceiling on det(I+S): n^(n/2) = {ceil}; "
              f"max det = {detmax}; attained iff skew-conference "
              f"SS^T=(n-1)I: {'ATTAINED' if detmax == ceil else 'not attained'}]")
    print(f"  H-max = {hmax}: {len(Hset)} class(es)")
    for c in Hset:
        print(f"    cid {c}: d = {d[c]}, score {score[c]}, "
              f"SC = {'Y' if sc[c] else 'N'}, dist-from-transitive = "
              f"{dist_t.get(c, 'inf')}")
    inter = set(Dset) & set(Hset)
    distDH = 0 if inter else min(
        bfs_dist(adjl, Dset).get(c, float('inf')) for c in Hset)
    print(f"  overlap argmax(d) & argmax(H): {sorted(inter) if inter else 'none'}")
    print(f"  graph distance between argmax sets: {distDH}")
    dist_t_vals = [v for v in dist_t.values()]
    print(f"  eccentricity of transitive node: {max(dist_t_vals)}; "
          f"graph diameter sample (from transitive) covers {len(dist_t)}/{nc} nodes")
    # top-10 by d for medium n
    if n >= 6:
        print(f"  Top 10 classes by d:")
        for c in sorted(range(nc), key=lambda c: -d[c])[:10]:
            print(f"    cid {c}: d = {d[c]}, H = {H[c]}, score {score[c]}, "
                  f"SC = {'Y' if sc[c] else 'N'}, dist-from-trans = "
                  f"{dist_t.get(c, 'inf')}")

    return dict(nc=nc)


def main():
    ns = [int(a) for a in sys.argv[1:]] or [5, 6, 7]
    print("=" * 76)
    print("METAGRAPH DET GRADIENT — d(T) = det(I+S)/2^(n-1) on iso-class graph G_n")
    print("macmini-2026-06-10-S2, HYP-2386")
    print("=" * 76)
    for n in ns:
        B = build(n)
        analyze(B)
    print("\nDONE.")


if __name__ == "__main__":
    main()
