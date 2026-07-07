#!/usr/bin/env python3
"""
monad-explorer-2026-07-07-S6 -- THE LINE-INCIDENCE ALGEBRA of blue/black lines (HYP-4967).

Owner directive: beyond counting blue/black lines and node types (mac-mini-S46 THM-643,
klein-S161 own those), determine HOW the lines structure the metagraph:
  Q1  Is the flip-induced class map CLASS-FUNCTIONAL (does class(flip t) depend only on
      class(t))?  If yes it is a metagraph involution; if no, measure its concentration.
      Identify it against known maps (converse? identity? H-preserving?).
  Q2  Cross-line node-type incidence: which node types (pure-blue / mixed / pure-black=NS)
      do blue and black CROSS lines connect, with counts, per n.
  Q3  Allocation vs n: totals (lines, self-loops, cross) per color; growth shape.
  Q4  Reconstruction bookkeeping: per class C the pair (B(C), K(C)) with the σ-parity
      facts (B odd on SC, 0 on NS; K even on SC, odd on NS -- the owner's 'blues odd /
      blacks even', proved via the anti-transpose involution + Redei).

Conventions (CLAUDE.md strict): base path n -> n-1 -> ... -> 1 (arcs k -> k-1 fixed);
tiles = pairs (x,y), x-y >= 2, bit=1 means x -> y; anti-transpose sigma: (x,y) ->
(n-y+1, n-x+1) (verified: #gridsym = 2^{f+(m-f)/2} matches canon exponents);
flip = complement tiling = XOR all tile bits (path arcs kept); LINE = {t, flip t};
BLUE iff t gridsym.  Class = tournament isomorphism class (canonical min over S_n).
"""
import itertools
import numpy as np
from collections import Counter, defaultdict

def analyze(n, verbose=True):
    verts = list(range(1, n + 1))
    pairs = [(i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1)]  # i<j slots
    pidx = {p: k for k, p in enumerate(pairs)}
    npairs = len(pairs)
    tiles = [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2]
    m = len(tiles)
    ntil = 1 << m

    # full arc code: bit per pair slot (i<j): 1 iff arc i -> j
    # base path arcs: k -> k-1, i.e. pair (k-1, k) oriented high->low => bit 0 for slot (k-1,k)
    base_bits = np.zeros(npairs, dtype=bool)  # all path arcs high->low: arc j->i for (i,j): bit 0
    # tile (x,y): x>y: bit=1 in tiling means x->y => slot (y,x) arc from higher to lower => code bit 0
    #             bit=0 means y->x => code bit 1 (arc from lower i to higher j... define code bit = 1 iff i->j)
    tile_slots = np.array([pidx[(y, x)] for (x, y) in tiles])

    tilings = np.arange(ntil, dtype=np.int64)
    tilebits = ((tilings[:, None] >> np.arange(m)[None, :]) & 1).astype(bool)   # ntil x m

    codes = np.zeros((ntil, npairs), dtype=bool)
    codes[:, tile_slots] = ~tilebits          # bit=1 in tiling => x->y => slot bit (i->j) = 0

    # canonicalization over S_n
    perms = list(itertools.permutations(range(1, n + 1)))
    canon = np.full(ntil, (1 << npairs) - 1, dtype=np.int64)
    weights = (1 << np.arange(npairs, dtype=np.int64))
    for perm in perms:
        pm = {v: perm[v - 1] for v in verts}
        colmap = np.empty(npairs, dtype=np.int64)
        flip = np.zeros(npairs, dtype=bool)
        for k, (i, j) in enumerate(pairs):
            a, b = pm[i], pm[j]
            if a < b:
                colmap[k] = pidx[(a, b)]
            else:
                colmap[k] = pidx[(b, a)]
                flip[k] = True
        transformed = np.empty((ntil, npairs), dtype=bool)
        transformed[:, colmap] = codes ^ flip[None, :] if False else codes ^ flip[None, :]
        # careful: bit at new slot colmap[k] comes from old slot k (value possibly flipped)
        # rebuild properly:
        tr = np.empty((ntil, npairs), dtype=bool)
        tr[:, colmap] = codes ^ flip[None, :]
        vals = tr @ weights
        np.minimum(canon, vals, out=canon)

    classes = {c: idx for idx, c in enumerate(sorted(set(canon.tolist())))}
    cls = np.array([classes[c] for c in canon.tolist()])
    ncls = len(classes)

    # gridsym: sigma on tile positions
    sig = {}
    for k, (x, y) in enumerate(tiles):
        sx, sy = n - y + 1, n - x + 1
        sig[k] = tiles.index((sx, sy))
    fixed = [k for k in sig if sig[k] == k]
    porb = [(k, sig[k]) for k in sig if k < sig[k]]
    gs = np.ones(ntil, dtype=bool)
    for a, b in porb:
        gs &= (tilebits[:, a] == tilebits[:, b])

    # converse class map (flip ALL npairs bits)
    full_mask = (1 << npairs) - 1
    codes_int = codes @ weights
    # converse of canonical rep: conv canon by transforming codes ^ all-ones
    conv_canon = np.full(ntil, (1 << npairs) - 1, dtype=np.int64)
    codes_conv = ~codes
    for perm in perms:
        pm = {v: perm[v - 1] for v in verts}
        colmap = np.empty(npairs, dtype=np.int64)
        flip = np.zeros(npairs, dtype=bool)
        for k, (i, j) in enumerate(pairs):
            a, b = pm[i], pm[j]
            if a < b:
                colmap[k] = pidx[(a, b)]
            else:
                colmap[k] = pidx[(b, a)]
                flip[k] = True
        tr = np.empty((ntil, npairs), dtype=bool)
        tr[:, colmap] = codes_conv ^ flip[None, :]
        vals = tr @ weights
        np.minimum(conv_canon, vals, out=conv_canon)
    cls_conv = np.array([classes.get(c, -1) for c in conv_canon.tolist()])
    # class-level converse map (well-defined by iso-invariance)
    conv_map = {}
    for t in range(ntil):
        conv_map.setdefault(cls[t], set()).add(cls_conv[t])
    assert all(len(v) == 1 for v in conv_map.values())
    conv_map = {k: v.pop() for k, v in conv_map.items()}
    SC = {c for c, cc in conv_map.items() if cc == c}

    # per-class tallies
    T = np.bincount(cls, minlength=ncls)
    B = np.bincount(cls[gs], minlength=ncls)
    K = T - B
    ntype = {}
    for c in range(ncls):
        if B[c] == 0:
            ntype[c] = "pureblack"
        elif K[c] == 0:
            ntype[c] = "pureblue"
        else:
            ntype[c] = "mixed"

    # parity facts
    p1 = all(T[c] % 2 == 1 for c in range(ncls))
    p2 = all((B[c] % 2 == 1) == (c in SC) for c in range(ncls))
    p3 = all((K[c] % 2 == 0) == (c in SC) for c in range(ncls))
    p4 = all((ntype[c] == "pureblack") == (c not in SC) for c in range(ncls))

    # flip lines
    flip_t = tilings ^ (ntil - 1)
    lines_blue = []
    lines_black = []
    for t in range(ntil):
        ft = int(flip_t[t])
        if t < ft:
            (lines_blue if gs[t] else lines_black).append((cls[t], cls[ft]))
    # Q1: class-functionality of flip
    img = defaultdict(Counter)
    for t in range(ntil):
        img[cls[t]][cls[int(flip_t[t])]] += 1
    functional = sum(1 for c in img if len(img[c]) == 1)
    conc = np.mean([max(img[c].values()) / sum(img[c].values()) for c in img])
    # is single-image class = converse? or identity?
    single_id = sum(1 for c in img if len(img[c]) == 1 and next(iter(img[c])) == c)
    single_conv = sum(1 for c in img if len(img[c]) == 1 and next(iter(img[c])) == conv_map[c])

    def line_stats(lines):
        self_l = Counter(ntype[a] for a, b in lines if a == b)
        cross = [(a, b) for a, b in lines if a != b]
        tp = Counter(tuple(sorted((ntype[a], ntype[b]))) for a, b in cross)
        return len(lines), dict(self_l), len(cross), dict(tp)

    nb = line_stats(lines_blue)
    nk = line_stats(lines_black)
    func_types = Counter(ntype[c] for c in img if len(img[c]) == 1)
    res = dict(n=n, m=m, ntil=ntil, ncls=ncls, SC=len(SC),
               types=Counter(ntype.values()), parity=(p1, p2, p3, p4),
               blue=nb, black=nk, functional=functional, conc=conc,
               single_id=single_id, single_conv=single_conv,
               TBK=[(int(T[c]), int(B[c]), int(K[c]), ntype[c], c in SC)
                    for c in range(ncls)])
    if verbose:
        print(f"n={n}: m={m} tilings={ntil} classes={ncls} SC={len(SC)} types={dict(res['types'])}")
        print(f"  PARITY: T odd all: {p1}; B odd<=>SC: {p2}; K even<=>SC: {p3}; pureblack<=>NS: {p4}")
        print(f"  BLUE  lines: total={nb[0]} self={nb[1]} cross={nb[2]} cross-node-types={nb[3]}")
        print(f"  BLACK lines: total={nk[0]} self={nk[1]} cross={nk[2]} cross-node-types={nk[3]}")
        print(f"  FLIP class-map: functional on {functional}/{ncls} classes "
              f"(avg concentration {conc:.3f}); single-image = identity: {single_id}, "
              f"= converse: {single_conv}; functional-class types: {dict(func_types)}")
    return res

if __name__ == "__main__":
    all_res = []
    for n in (4, 5, 6, 7):
        all_res.append(analyze(n))
        print()
    print("ALLOCATION vs n (blue total / black total / blue cross / black cross):")
    for r in all_res:
        print(f"  n={r['n']}: blue {r['blue'][0]} (cross {r['blue'][2]}), "
              f"black {r['black'][0]} (cross {r['black'][2]}); "
              f"nodes pb/mix/pk = {r['types'].get('pureblue',0)}/"
              f"{r['types'].get('mixed',0)}/{r['types'].get('pureblack',0)}")
    print("\nDONE.")
