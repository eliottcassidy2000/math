"""
mac-mini-2026-07-07-S48 (HYP-5037, THM-648) -- PROOF VERIFICATION: blue self-loop lines
exist only at EVEN n (THM-643 C2 / klein-S161 C1, now a theorem).

PROOF (canon THM-648; verified step-by-step here):
  Let t be gridsym with tournament T, labeled scores s.
  (1) THM-644 (opus): gridsym <=> rho(i) = n+1-i is an anti-automorphism
      => s(rho(v)) = n-1-s(v) for all v; in particular s(n) = n-1-s(1).
  (2) THM-646 (mac-mini): scores of flip(t): s'(v) = c(v)-s(v), c = (n-2, n-1,...,n-1, n).
  (3) Blue self-loop at class C <=> t, flip(t) both in fiber(C) => multiset{s'} = multiset{s}.
      multiset{s'} = multiset{n-1-s(v)} - {n-1-s(1), n-1-s(n)} + {n-2-s(1), n-s(n)}
                   = multiset{s}       - {n-1-a, a}            + {n-2-a, a+1}   (a = s(1))
      (using (1) twice).  Equality <=> {n-1-a, a} = {n-2-a, a+1}
      <=> a = (n-2)/2  =>  n EVEN.  QED (odd-n impossibility).
  Even n: s(1) = (n-2)/2 is NECESSARY for a blue self-loop (not claimed sufficient).

VERIFICATION (n=4..7): (a) step (1)'s score relation on every gridsym tiling;
(b) the multiset bookkeeping identity on every gridsym tiling; (c) the census of blue
self-loops (0,1,0,2,0 at n=3..7) with the s(1) = (n-2)/2 condition checked on each;
(d) at even n: which classes carry them + does every gridsym tiling with s(1)=(n-2)/2
give a self-loop? (sufficiency probe, reported honestly).
"""
import numpy as np
from itertools import permutations, combinations
from math import comb

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
        tourn[:, p_i] = tb[:, t_i]
    S = np.zeros((T, n), dtype=np.int16)
    for a, (i, j) in enumerate(pairs):
        S[:, i-1] += tourn[:, a]
        S[:, j-1] += (1 - tourn[:, a])
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
    gs_map = [tiles.index((n-y+1, n-x+1)) for (x, y) in tiles]
    gs = np.ones(T, dtype=bool)
    for i, j in enumerate(gs_map):
        if i < j: gs &= (tb[:, i] == tb[:, j])
    partner = np.arange(T) ^ (T-1)

    print(f"\n===== n={n} =====")
    # (a) rho-score relation on gridsym tilings
    idx_gs = np.nonzero(gs)[0]
    rho = np.arange(n)[::-1]   # rho(v)=n+1-v in 0-index: v -> n-1-v
    rel_ok = np.all(S[idx_gs][:, rho] == (n-1) - S[idx_gs])
    print(f"  (a) s(rho v) = n-1-s(v) on all {len(idx_gs)} gridsym tilings: {bool(rel_ok)}")
    # (b) multiset bookkeeping: multiset(c - s) vs modified multiset
    c_vec = np.array([(max(v-2,0) + max(n-1-v,0) + (2 if v >= 2 else 0)) for v in range(1, n+1)])
    ok_b = True
    for t in idx_gs:
        s = S[t]; a = int(s[0])
        ms  = sorted(c_vec - s)
        base = sorted(s)
        # remove {n-1-a, a}, add {n-2-a, a+1}
        tmp = list(base)
        try:
            tmp.remove(n-1-a); tmp.remove(a)
        except ValueError:
            ok_b = False; break
        tmp += [n-2-a, a+1]
        if sorted(tmp) != ms: ok_b = False; break
    print(f"  (b) bookkeeping identity multiset(c-s) = multiset(s) -{{n-1-a,a}}+{{n-2-a,a+1}}: {ok_b}")
    # (c) blue self-loops census + condition
    loops = [t for t in idx_gs if cls[t] == cls[partner[t]] and t < partner[t]]
    print(f"  (c) blue self-loop lines: {len(loops)}")
    for t in loops:
        a = int(S[t][0])
        print(f"      loop at class {int(cls[t])}: s = {tuple(S[t])}, s(1) = {a} "
              f"((n-2)/2 = {(n-2)/2}) {'OK' if a == (n-2)//2 and n % 2 == 0 else 'VIOLATION'}")
    if n % 2 == 1:
        print(f"      odd n: theorem predicts 0 loops -> {'CONFIRMED' if not loops else 'VIOLATION'}")
    # (d) sufficiency probe at even n
    if n % 2 == 0:
        cand = [t for t in idx_gs if int(S[t][0]) == (n-2)//2]
        suff = sum(1 for t in cand if cls[t] == cls[partner[t]])
        print(f"  (d) gridsym tilings with s(1)=(n-2)/2: {len(cand)}; of these, self-loops: {suff} "
              f"(condition necessary{'/NOT sufficient' if suff < len(cand) else ' AND sufficient here'})")

for n in (3, 4, 5, 6, 7):
    run(n)
