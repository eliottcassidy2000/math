"""aud_ttwords.py -- audit S11: H-values of seam words built only from two-port transitive gadgets TT_k (k <= K).
Own implementation of the THM-4879 (B) closed form (checked independently for k <= 9 by aud_tt.c).
Single-gadget states: (TT_k, z) has v = (1, [z = 0], 2^z).  States with R <= X are expanded.
Two ways of recording attained H <= X:
  (a) 'kept only'  : H of a state is recorded only when the state itself has H <= X and R <= X
                     (the procedure of the audited transitive_words_bfs.py);
  (b) 'all images' : H of every image with H <= X is recorded, whatever its R (complete: H'' >= R' for every gadget,
                     so states with R > X have no descendant with H <= X, but their own H can be <= X).
Usage: python3 aud_ttwords.py X K1 K2 ..."""
import sys
X = int(sys.argv[1]); Ks = [int(a) for a in sys.argv[2:]]

def closed(k, u, z):
    Hrow = (1, -(1 if u == k - 1 else 0), 2 ** (k - 1 - u))
    if k == 1: Srow = (1, -1, 0)
    elif z == 0: Srow = (1, 0, 0) if u == 0 else Hrow
    else: Srow = (0, 0, 2 ** (k - 1 - z)) if u == 0 else (0, 0, 0)
    if z == k - 1 and k - 1 > u: mu = 2 ** u
    elif u == k - 1 and k - 1 > z: mu = 2 ** z
    else: mu = 0
    c = 2 ** (k - 1 - u + z) if z <= u else 2 ** (k - z + u) * 3 ** (z - u - 1)
    return (Hrow, Srow, (2 ** z, -mu, c))

for K in Ks:
    maps = sorted({closed(k, u, z) for k in range(1, K + 1) for u in range(k) for z in range(k)})
    seeds = {(1, 1, 1)} | {(1, 0, 2 ** z) for z in range(1, K) if 2 ** z <= X}
    V = set(seeds); fr = list(V)
    hit_all = {s[0] for s in seeds}
    while fr:
        nx = []
        for (H, S, R) in fr:
            for (a, b, c) in maps:
                h = a[0] * H + a[1] * S + a[2] * R
                if h > X: continue
                hit_all.add(h)
                r = c[0] * H + c[1] * S + c[2] * R
                if r > X: continue
                p = (h, b[0] * H + b[1] * S + b[2] * R, r)
                if p not in V: V.add(p); nx.append(p)
        fr = nx
    hit_kept = {h for h, s, r in V}
    miss_a = [n for n in range(1, X + 1, 2) if n not in hit_kept]
    miss_b = [n for n in range(1, X + 1, 2) if n not in hit_all]
    print(f"X={X} K={K}: distinct maps {len(maps)}, kept states {len(V)}; "
          f"missed odd H (a) kept-only {len(miss_a)} first {miss_a[:14]}; (b) all images {len(miss_b)} first {miss_b[:14]}")
