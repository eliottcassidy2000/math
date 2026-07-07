"""
mac-mini-2026-07-07-S46 (HYP-4977, THM-643) -- the BLUE/BLACK LINE structure of the
tiling fibration, STRICT explorer definitions (CLAUDE.md), n = 3..7 exhaustive.

OBJECTS (strict):
  m = C(n-1,2) tiles = non-consecutive pairs (x,y), x-y >= 2; base path n -> ... -> 1 fixed.
  LINE = (t, flip(t)) with flip = complement tiling (XOR all-ones)  [d = m layer].
  BLUE line: t grid-symmetric (bit(x,y) = bit(n-y+1, n-x+1) for all tiles); else BLACK.
  Node = iso class; fiber(C) = tilings whose tournament is in C.
  PURE BLUE node: all fiber tilings gridsym; PURE BLACK: none; MIXED: both.

THEOREM TARGETS (verified here numerically, proofs in THM-643):
  T1 fiber(C) = H(C)/|Aut(C)| and is ODD for every class (Redei x odd Aut).
  T2 gridsym tiling => its class is SC (self-converse); SC class => #gridsym in fiber
     is ODD (>= 1); non-SC => 0.  Hence PURE BLACK = NON-SC exactly, all n.
  T3 blue lines connect only SC nodes (endpoint classes both SC).
  T4 per-node allocation parity: SC: (#blue tilings odd, #black even); NS: (0, odd).
     Cross-class line-ENDPOINT counts: every node odd total; SC: odd blue-cross,
     even black-cross; NS: odd black-cross.
  T5 global: #gridsym = 2^{(m+f)/2}, f = floor((n-1)/2); #blue lines = 2^{(m+f)/2 - 1};
     #black lines = 2^{m-1} - 2^{(m+f)/2 - 1}.
CENSUS OUTPUTS: node-type counts (pure-blue / mixed / pure-black) vs SC/NS; line
allocation (blue/black x self-loop/cross x endpoint types); per-SC-node blue degrees;
pure-blue characterization data; conjecture mining.
"""
import numpy as np
from itertools import permutations, combinations
from math import comb

def run(n, verbose=True):
    pairs = list(combinations(range(1, n+1), 2))          # (i<j) vertex pairs
    pidx = {p: i for i, p in enumerate(pairs)}
    P = len(pairs)
    # base path arcs k -> k-1 : pair (k-1,k), direction j->i => bit 0 (bit=1 iff i->j)
    base_bits = np.zeros(P, dtype=np.uint8)
    base_mask = np.zeros(P, dtype=bool)
    for k in range(2, n+1):
        base_mask[pidx[(k-1, k)]] = True
    # tiles: (x,y) x-y>=2, in explorer order (for y=1..n-2: for x=n..y+2)
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2: tiles.append((x, y))
    m = len(tiles)
    assert m == comb(n-1, 2)
    tile_pair_idx = [pidx[(y, x)] for (x, y) in tiles]     # pair (y<x)
    # grid symmetry on tile positions: (x,y) -> (n-y+1, n-x+1)
    gs_map = [tiles.index((n-y+1, n-x+1)) for (x, y) in tiles]
    f = sum(1 for i, j in enumerate(gs_map) if i == j)
    # enumerate tilings as bit rows
    T = 1 << m
    tb = ((np.arange(T)[:, None] >> np.arange(m)[None, :]) & 1).astype(np.uint8)
    # tournament pair-bit matrix: bit(i<j) = 1 iff i -> j
    tourn = np.zeros((T, P), dtype=np.uint8)
    for t_i, p_i in enumerate(tile_pair_idx):
        # tile bit b=1 => y -> x (i->j, since y<x) => bit 1 ; b=0 => x->y => bit 0
        tourn[:, p_i] = tb[:, t_i]
    # base path bits already 0 (k -> k-1 = j -> i)
    # canonicalization over S_n
    perms = list(permutations(range(1, n+1)))
    weights = (1 << np.arange(P)[::-1]).astype(object)
    def pack(bits):   # bits (T,P) -> object ints
        return (bits.astype(object) * weights[None, :]).sum(axis=1)
    canon = None
    for perm in perms:
        colmap = np.empty(P, dtype=np.int64); flip = np.empty(P, dtype=np.uint8)
        for a, (i, j) in enumerate(pairs):
            pi, pj = perm[i-1], perm[j-1]
            if pi < pj: colmap[a] = pidx[(pi, pj)]; flip[a] = 0
            else:       colmap[a] = pidx[(pj, pi)]; flip[a] = 1
        permuted = tourn[:, colmap] ^ flip[None, :]
        v = pack(permuted)
        canon = v if canon is None else np.minimum(canon, v)
    # classes
    uniq, cls = np.unique(canon, return_inverse=True)
    nC = len(uniq)
    # converse canon (reverse all arcs = complement all pair bits)
    conv_canon = None
    conv_t = tourn ^ 1
    for perm in perms:
        colmap = np.empty(P, dtype=np.int64); flip = np.empty(P, dtype=np.uint8)
        for a, (i, j) in enumerate(pairs):
            pi, pj = perm[i-1], perm[j-1]
            if pi < pj: colmap[a] = pidx[(pi, pj)]; flip[a] = 0
            else:       colmap[a] = pidx[(pj, pi)]; flip[a] = 1
        permuted = conv_t[:, colmap] ^ flip[None, :]
        v = pack(permuted)
        conv_canon = v if conv_canon is None else np.minimum(conv_canon, v)
    is_sc_tiling = (canon == conv_canon)
    # SC per class (constant on classes)
    sc_class = np.zeros(nC, dtype=bool)
    for c in range(nC):
        sc_class[c] = is_sc_tiling[np.argmax(cls == c)]
    # gridsym per tiling
    gs = np.ones(T, dtype=bool)
    for i, j in enumerate(gs_map):
        if i < j: gs &= (tb[:, i] == tb[:, j])
    n_gs = int(gs.sum())
    # H and |Aut| per class via a representative tiling
    def ham_count(bits):
        adj = np.zeros((n, n), dtype=bool)
        for a, (i, j) in enumerate(pairs):
            if bits[a]: adj[i-1, j-1] = True
            else:       adj[j-1, i-1] = True
        full = 1 << n
        dp = np.zeros((full, n), dtype=np.int64)
        for v in range(n): dp[1 << v, v] = 1
        for S in range(full):
            for v in range(n):
                if not dp[S, v]: continue
                if not (S >> v) & 1: continue
                for w in range(n):
                    if (S >> w) & 1: continue
                    if adj[v, w]: dp[S | (1 << w), w] += dp[S, v]
        return int(dp[full-1].sum())
    def aut_count(bits):
        cnt = 0
        for perm in perms:
            ok = True
            for a, (i, j) in enumerate(pairs):
                pi, pj = perm[i-1], perm[j-1]
                if pi < pj: b2 = bits[pidx[(pi, pj)]]
                else:       b2 = bits[pidx[(pj, pi)]] ^ 1
                if b2 != bits[a]: ok = False; break
            if ok: cnt += 1
        return cnt
    reps = [np.argmax(cls == c) for c in range(nC)]
    H = np.array([ham_count(tourn[r]) for r in reps])
    Aut = np.array([aut_count(tourn[r]) for r in reps])
    fiber = np.bincount(cls, minlength=nC)
    gs_in_fiber = np.bincount(cls[gs], minlength=nC)
    # lines: t <-> t ^ full
    full_mask = T - 1
    partner_cls = cls[np.arange(T) ^ full_mask]
    # per-tiling line data (each line counted from both endpoints; divide by 2 for lines)
    self_line = (cls == partner_cls)
    # node types
    pure_blue = (gs_in_fiber == fiber)
    pure_black = (gs_in_fiber == 0)
    mixed = ~pure_blue & ~pure_black
    # T-checks
    ok_T1 = np.all(fiber * Aut == H) and np.all(fiber % 2 == 1)
    ok_T2 = np.all((gs_in_fiber > 0) == sc_class) and np.all((gs_in_fiber % 2 == 1) == sc_class)
    # blue lines endpoints SC
    blue_ok = np.all(sc_class[cls[gs]])
    # per-node endpoint parities
    blue_tilings = gs_in_fiber
    black_tilings = fiber - gs_in_fiber
    cross_end = np.bincount(cls[~self_line], minlength=nC)
    blue_cross_end = np.bincount(cls[gs & ~self_line], minlength=nC)
    black_cross_end = np.bincount(cls[~gs & ~self_line], minlength=nC)
    ok_T4a = np.all((blue_tilings % 2 == 1) == sc_class) and np.all(black_tilings[sc_class] % 2 == 0)
    ok_T4b = np.all(cross_end % 2 == 1)
    ok_T4c = np.all(blue_cross_end[sc_class] % 2 == 1) and np.all(black_cross_end[sc_class] % 2 == 0) \
             and np.all(black_cross_end[~sc_class] % 2 == 1)
    n_blue_lines = n_gs // 2
    n_black_lines = (T - n_gs) // 2
    ok_T5 = (n_gs == 2**((m+f)//2)) and ((m+f) % 2 == 0)
    # line-type allocation by endpoint node types (count each LINE once)
    seen = np.arange(T) < (np.arange(T) ^ full_mask)
    def node_type(c):
        return 'PB' if pure_blue[c] else ('PK' if pure_black[c] else 'MX')
    alloc = {}
    for t in np.nonzero(seen)[0]:
        u = t ^ full_mask
        col = 'blue' if gs[t] else 'black'
        kinds = tuple(sorted((node_type(cls[t]), node_type(cls[u]))))
        loop = 'self' if cls[t] == cls[u] else 'cross'
        alloc[(col, loop, kinds)] = alloc.get((col, loop, kinds), 0) + 1
    res = dict(n=n, m=m, f=f, nC=nC, SC=int(sc_class.sum()),
               pure_blue=int(pure_blue.sum()), mixed=int(mixed.sum()),
               pure_black=int(pure_black.sum()), n_gs=n_gs,
               blue_lines=n_blue_lines, black_lines=n_black_lines,
               T1=bool(ok_T1), T2=bool(ok_T2), T3=bool(blue_ok),
               T4a=bool(ok_T4a), T4b=bool(ok_T4b), T4c=bool(ok_T4c), T5=bool(ok_T5),
               alloc=alloc)
    if verbose:
        print(f"\n===== n={n}: m={m}, f={f}, classes={nC} (SC={res['SC']}), tilings={T} =====")
        print(f"  node types: PURE-BLUE={res['pure_blue']}  MIXED={res['mixed']}  PURE-BLACK={res['pure_black']}"
              f"   [pure-black == non-SC? {res['pure_black'] == nC - res['SC']}]")
        print(f"  gridsym tilings = {n_gs} (formula 2^((m+f)/2) = {2**((m+f)//2)});"
              f" blue lines = {n_blue_lines}, black lines = {n_black_lines} (total {T//2})")
        print(f"  THEOREM CHECKS: T1(fiber=H/Aut odd):{ok_T1}  T2(pureblack=nonSC,parity):{ok_T2}"
              f"  T3(blue on SC):{blue_ok}")
        print(f"                  T4a(SC blue-odd/black-even):{ok_T4a}  T4b(all cross-ends odd):{ok_T4b}"
              f"  T4c(colored cross parities):{ok_T4c}  T5(gs count):{ok_T5}")
        print("  line allocation (color, self/cross, endpoint-types): count")
        for k in sorted(alloc): print(f"    {k}: {alloc[k]}")
        # pure-blue characterization + SC blue-degree distribution
        pb = np.nonzero(pure_blue)[0]
        print(f"  PURE-BLUE classes: {[(int(H[c]), int(Aut[c]), int(fiber[c])) for c in pb]}  (H,|Aut|,fiber)")
        scn = np.nonzero(sc_class)[0]
        bd = sorted(int(blue_tilings[c]) for c in scn)
        print(f"  SC-node #blue-tilings distribution: {bd}")
        mx = np.nonzero(mixed)[0]
        print(f"  MIXED classes (H,|Aut|,fiber,#gs): {[(int(H[c]),int(Aut[c]),int(fiber[c]),int(gs_in_fiber[c])) for c in mx][:12]}")
    return res

results = [run(n) for n in (3, 4, 5, 6, 7)]

print("\n===== GLOBAL TABLE =====")
print(f"{'n':>2} {'m':>3} {'f':>2} {'classes':>8} {'SC':>4} {'PB':>4} {'MX':>4} {'PK':>6} "
      f"{'gridsym':>8} {'blueL':>7} {'blackL':>9}")
for r in results:
    print(f"{r['n']:2d} {r['m']:3d} {r['f']:2d} {r['nC']:8d} {r['SC']:4d} {r['pure_blue']:4d} "
          f"{r['mixed']:4d} {r['pure_black']:6d} {r['n_gs']:8d} {r['blue_lines']:7d} {r['black_lines']:9d}")
print("\nn=8 by formula: m=21, f=3 -> gridsym = 2^12 = 4096, blue lines = 2048, "
      "black lines = 2^20 - 2^11 = 1046528; pure-black = 6880 - SC(8) (SC(8) from census lit).")
