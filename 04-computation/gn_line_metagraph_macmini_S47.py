"""
mac-mini-2026-07-07-S47 (HYP-5017, THM-646) -- THE LINE METAGRAPH: the explicit flip
matching on iso classes, its deterministic skeleton, and its generating pattern.

Builds on the triple-proved parity layer (THM-643 / THM-644 / klein-S161).

DERIVED LAW TO VERIFY (THM-646 target, proof in canon file):
  For any tiling t with labeled score vector s and partner t' = flip(t) with s':
      s(v) + s'(v) = c(v)  for every vertex v,  c = (n-2, n-1, n-1, ..., n-1, n)
  (in the base-path labeling; c(v) = deg_tiles(v) + 2*[v>=2], deg_tiles(v) =
   (v-2)_+ + (n-1-v)_+).  The flip is a SCORE-COMPLEMENT: the matching's
  deterministic skeleton.  Corollaries: the transitive class's unique line partner is
  the class with score multiset {1,2,...,n-3, n-2, n-2} (conjecturally the principal-
  line SC neighbor at H = 2^{n-2}+1); N(C,D) = 0 unless the score multisets of C and D
  are c-compatible (a labeled matching between them summing to c).

CENSUS (n = 4..7):
  1. verify the law on every tiling;
  2. the LINE KERNEL N(C,D) = #lines between classes (with color split);
     support size vs the c-compatibility bound; klein-S161 C5 quasi-randomness tested
     ON the compatible support: N vs fiber(C)*fiber(D)/2^{m-1};
  3. the transitive-partner identification (H value = 2^{n-2}+1?);
  4. THE STRATA QUOTIENT (the metagraph of the line metagraph): nodes = score
     multisets, edges = induced line-pairs; is the quotient near-deterministic
     (per-stratum out-degree distribution)?
  5. even-graph projection: K* = cycle-space address of the all-ones tiling; lines
     project to the +K* translation (verify: address(flip t) = address(t) XOR K*).
"""
import numpy as np
from itertools import permutations, combinations
from math import comb
from collections import defaultdict

def run(n):
    pairs = list(combinations(range(1, n+1), 2))
    pidx = {p: i for i, p in enumerate(pairs)}
    P = len(pairs)
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2: tiles.append((x, y))
    m = len(tiles)
    tile_pair_idx = [pidx[(y, x)] for (x, y) in tiles]
    T = 1 << m
    tb = ((np.arange(T)[:, None] >> np.arange(m)[None, :]) & 1).astype(np.uint8)
    tourn = np.zeros((T, P), dtype=np.uint8)
    for t_i, p_i in enumerate(tile_pair_idx):
        tourn[:, p_i] = tb[:, t_i]     # bit=1 iff lower->higher (i->j)
    # scores (labeled): s(v) = #wins of v
    S = np.zeros((T, n), dtype=np.int16)
    for a, (i, j) in enumerate(pairs):
        S[:, i-1] += tourn[:, a]           # i->j when bit 1
        S[:, j-1] += (1 - tourn[:, a])     # j->i when bit 0
    # canonicalization
    perms = list(permutations(range(1, n+1)))
    weights = (1 << np.arange(P)[::-1]).astype(object)
    canon = None
    for perm in perms:
        colmap = np.empty(P, dtype=np.int64); flip = np.empty(P, dtype=np.uint8)
        for a, (i, j) in enumerate(pairs):
            pi, pj = perm[i-1], perm[j-1]
            if pi < pj: colmap[a] = pidx[(pi, pj)]; flip[a] = 0
            else:       colmap[a] = pidx[(pj, pi)]; flip[a] = 1
        v = ((tourn[:, colmap] ^ flip[None, :]).astype(object) * weights[None, :]).sum(axis=1)
        canon = v if canon is None else np.minimum(canon, v)
    uniq, cls = np.unique(canon, return_inverse=True)
    nC = len(uniq)
    fiber = np.bincount(cls, minlength=nC)
    # gridsym
    gs_map = [tiles.index((n-y+1, n-x+1)) for (x, y) in tiles]
    gs = np.ones(T, dtype=bool)
    for i, j in enumerate(gs_map):
        if i < j: gs &= (tb[:, i] == tb[:, j])
    # H per class (for the transitive-partner check)
    def ham_count(bits):
        adj = np.zeros((n, n), dtype=bool)
        for a, (i, j) in enumerate(pairs):
            if bits[a]: adj[i-1, j-1] = True
            else:       adj[j-1, i-1] = True
        full = 1 << n
        dp = np.zeros((full, n), dtype=np.int64)
        for v in range(n): dp[1 << v, v] = 1
        for Sm in range(full):
            for v in range(n):
                if not dp[Sm, v] or not (Sm >> v) & 1: continue
                for w in range(n):
                    if (Sm >> w) & 1: continue
                    if adj[v, w]: dp[Sm | (1 << w), w] += dp[Sm, v]
        return int(dp[full-1].sum())
    reps = [int(np.argmax(cls == c)) for c in range(nC)]
    H = np.array([ham_count(tourn[r]) for r in reps])

    print(f"\n===== n={n}: m={m}, classes={nC}, tilings={T} =====")
    # 1. THE SCORE-COMPLEMENT LAW
    c_vec = np.array([(max(v-2,0) + max(n-1-v,0) + (2 if v >= 2 else 0)) for v in range(1, n+1)])
    partner = np.arange(T) ^ (T-1)
    law_ok = np.all(S + S[partner] == c_vec[None, :])
    print(f"  SCORE-COMPLEMENT LAW s + s' = c, c = {tuple(c_vec)}: "
          f"{'HOLDS on all ' + str(T) + ' tilings' if law_ok else 'VIOLATED'}")

    # 2. line kernel
    Nk = defaultdict(int); NkBlue = defaultdict(int)
    seen = np.arange(T) < partner
    for t in np.nonzero(seen)[0]:
        key = tuple(sorted((int(cls[t]), int(cls[partner[t]]))))
        Nk[key] += 1
        if gs[t]: NkBlue[key] += 1
    support = len(Nk)
    # c-compatibility bound on the support: score multisets must admit a labeled
    # c-complement matching; cheap necessary form: multiset(s') = multiset(c - s) for
    # SOME labeling in the class -- we use the observed per-class score-multiset SETS
    cls_scoresets = defaultdict(set)
    for t in range(T):
        cls_scoresets[int(cls[t])].add(tuple(sorted(S[t])))
    compat = 0
    for a in range(nC):
        for b in range(a, nC):
            ok = False
            for t in np.nonzero(cls == a)[0][:200]:
                if tuple(sorted(c_vec - S[t])) in cls_scoresets[b]: ok = True; break
            if ok: compat += 1
    print(f"  line-kernel support: {support} class-pairs (of {nC*(nC+1)//2} possible; "
          f"c-compatible pairs (sampled nec. cond.): {compat})")
    # quasi-randomness on the support
    tot_lines = T // 2
    obs, exp = [], []
    for (a, b), v in Nk.items():
        e = fiber[a] * fiber[b] / (2 * tot_lines) * (1 if a != b else 0.5)
        obs.append(v); exp.append(e)
    obs = np.array(obs, float); exp = np.array(exp, float)
    r = np.corrcoef(obs, exp)[0, 1] if len(obs) > 2 else float('nan')
    ratio = obs.sum() / exp.sum()
    print(f"  quasi-randomness on support: corr(N, fiber-product model) = {r:.4f}; "
          f"sum ratio = {ratio:.3f}; max |N/E| dev = {float(np.max(np.abs(obs/np.maximum(exp,1e-9) - 1))):.2f}")

    # 3. transitive partner
    trans_c = int(cls[int(np.argmin(S.sum(axis=1) - S.max(axis=1)))])  # find class with scores 0..n-1
    # robust: find tiling with score multiset {0,...,n-1}
    target = tuple(range(n))
    tr_t = next(t for t in range(T) if tuple(sorted(S[t])) == target)
    tr_partner_cls = int(cls[partner[tr_t]])
    print(f"  transitive partner: H = {H[tr_partner_cls]} (2^(n-2)+1 = {2**(n-2)+1}), "
          f"score multiset {tuple(sorted(S[partner[tr_t]]))} "
          f"(predicted {tuple(sorted(c_vec - S[tr_t]))})")

    # 4. strata quotient
    strata = defaultdict(int)
    for t in np.nonzero(seen)[0]:
        key = (tuple(sorted(S[t])), tuple(sorted(S[partner[t]])))
        key = tuple(sorted(key))
        strata[key] += 1
    n_strata_nodes = len({k for pair in strata for k in pair})
    outdeg = defaultdict(set)
    for (sa, sb) in strata:
        outdeg[sa].add(sb); outdeg[sb].add(sa)
    dist = sorted(len(v) for v in outdeg.values())
    print(f"  STRATA QUOTIENT: {n_strata_nodes} score-multiset nodes, {len(strata)} stratum-edges; "
          f"per-stratum partner-multiset degree distribution: {dist}")
    det = sum(1 for v in outdeg.values() if len(v) == 1)
    print(f"    strata with a UNIQUE partner multiset: {det}/{len(outdeg)} "
          f"({'DETERMINISTIC quotient' if det == len(outdeg) else 'non-deterministic'})")

    # 5. even-graph projection: address = XOR of fundamental cycles of chosen tiles
    # fundamental cycle of tile (x,y) wrt base path = edge set {(x,y)} xor path edges x..y
    # address as a vector over EDGE set of K_n (undirected pairs): cycle space element
    edge_idx = pidx
    def tile_cycle(x, y):
        vec = np.zeros(P, dtype=np.uint8)
        vec[edge_idx[(y, x)]] ^= 1
        for k in range(y, x):
            vec[edge_idx[(k, k+1)]] ^= 1
        return vec
    cyc = np.array([tile_cycle(x, y) for (x, y) in tiles])
    addr = (tb @ cyc) % 2
    Kstar = cyc.sum(axis=0) % 2
    proj_ok = np.all((addr[partner] - addr) % 2 == Kstar[None, :] % 2)
    print(f"  even-graph projection: address(flip t) = address(t) XOR K* : {bool(proj_ok)}; "
          f"K* edge count = {int(Kstar.sum())} (of {P})")
    return dict(n=n, support=support, r=r)

for n in (4, 5, 6, 7):
    run(n)
