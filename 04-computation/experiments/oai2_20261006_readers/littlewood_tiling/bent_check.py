import numpy as np, itertools
from tower_merit import tower
from switched_tower import rs_seq
for k in [2, 4]:
    H = tower(k).astype(np.int64); N = H.shape[0]
    # all d in {+-1}^N with d_0 = 1
    cnt = 0; best = 0
    D = np.array(list(itertools.product([1, -1], repeat=N - 1)), dtype=np.int64)
    D = np.concatenate([np.ones((D.shape[0], 1), dtype=np.int64), D], axis=1)
    Y = D @ H.T  # rows: (H d)
    m = np.abs(Y).min(axis=1)
    bent = np.all(np.abs(Y) == int(np.sqrt(N)), axis=1).sum()
    print(f"k={k} N={N}: #d (d_0=1) with H d in sqrt(N){{+-1}}^N (tower-bent): {bent};  max over d of min_i |(Hd)_i| = {m.max()}")
for k in range(2, 13):
    H = tower(k).astype(np.int64); a = rs_seq(k).astype(np.int64); N = H.shape[0]
    y = H @ a
    print(f"k={k} N={N}: RS vector: (H a)_i / sqrtN  min |.| = {np.abs(y).min()/np.sqrt(N):.3f}, max |.| = {np.abs(y).max()/np.sqrt(N):.3f}, #zeros = {(y==0).sum()}")
