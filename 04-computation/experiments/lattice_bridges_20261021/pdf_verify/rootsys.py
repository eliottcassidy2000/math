"""Root systems as integer vectors with an integer Gram scaling.
All returned root sets R are numpy int arrays; inner products are computed as R @ R.T // scale
so that roots have norm 2."""
import numpy as np
from itertools import product, combinations

def roots_A(k):
    R = []
    for i in range(k+1):
        for j in range(k+1):
            if i != j:
                v = [0]*(k+1); v[i] = 1; v[j] = -1; R.append(v)
    return np.array(R, dtype=np.int64), 1

def roots_D(k):
    R = []
    for i, j in combinations(range(k), 2):
        for si, sj in product((1, -1), repeat=2):
            v = [0]*k; v[i] = si; v[j] = sj; R.append(v)
    return np.array(R, dtype=np.int64), 1

def roots_E8():
    # doubled coordinates: D8 roots (2*(+-e_i +- e_j)) and (+-1)^8 with even # of minus signs
    R = []
    for i, j in combinations(range(8), 2):
        for si, sj in product((2, -2), repeat=2):
            v = [0]*8; v[i] = si; v[j] = sj; R.append(v)
    for s in product((1, -1), repeat=8):
        if s.count(-1) % 2 == 0:
            R.append(list(s))
    return np.array(R, dtype=np.int64), 4

def roots_E7():
    R, sc = roots_E8()
    a = R[np.where((R == 1).all(axis=1))[0][0]]       # (1,...,1)/2 is a root
    keep = (R @ a) == 0
    return R[keep], sc

def roots_E6():
    R, sc = roots_E8()
    IP = (R @ R.T) // sc
    a = 0
    b = int(np.where(IP[a] == -1)[0][0])
    keep = (IP[a] == 0) & (IP[b] == 0)
    return R[keep], sc

def roots(typ, k=None):
    if typ == 'A': return roots_A(k)
    if typ == 'D': return roots_D(k)
    if typ == 'E':
        return {6: roots_E6, 7: roots_E7, 8: roots_E8}[k]()
    raise ValueError

def gram(R, sc):
    return (R @ R.T) // sc

def direct_sum(parts):
    """parts: list of (R, sc). Returns block-diagonal root set with common scale and component labels."""
    L = 1
    for _, sc in parts:
        L = int(np.lcm(L, sc))
    dims = [R.shape[1] for R, _ in parts]
    tot = sum(dims)
    blocks, labels = [], []
    off = 0
    for ci, (R, sc) in enumerate(parts):
        f = int(round((L / sc) ** 0.5))
        assert f*f*sc == L, "scales must differ by squares"
        B = np.zeros((R.shape[0], tot), dtype=np.int64)
        B[:, off:off+R.shape[1]] = R * f
        blocks.append(B); labels += [ci]*R.shape[0]
        off += R.shape[1]
    return np.vstack(blocks), L, np.array(labels)

def count_tuples(IP, G, first_candidates=None):
    """Number of ordered tuples (v_1..v_m) of rows with IP[v_i,v_j] == G[i][j] (i != j)
    and IP[v_i,v_i] == G[i][i]. If first_candidates is given, v_1 ranges only over it.
    Last two levels done by a matrix-vector product."""
    G = np.asarray(G); m = G.shape[0]
    N = IP.shape[0]
    diag = np.diag(IP)
    base = [(diag == G[i, i]) for i in range(m)]
    eq = {}
    def E(val):
        if val not in eq:
            eq[val] = (IP == val)
        return eq[val]
    if m == 1:
        c = base[0].copy()
        if first_candidates is not None:
            mask = np.zeros(N, bool); mask[first_candidates] = True; c &= mask
        return int(c.sum())
    def rec(chosen, masks):
        t = len(chosen)
        if t == m - 2:
            c_a = masks[m-2].copy(); c_b = masks[m-1].copy()
            if t == 0 and first_candidates is not None:
                fm = np.zeros(N, bool); fm[first_candidates] = True; c_a &= fm
            for i, v in enumerate(chosen):
                c_a &= E(int(G[i, m-2]))[v]
                c_b &= E(int(G[i, m-1]))[v]
            M = E(int(G[m-2, m-1]))
            return int(c_a.astype(np.int64) @ (M.astype(np.int64) @ c_b.astype(np.int64)))
        cand = masks[t].copy()
        for i, v in enumerate(chosen):
            cand &= E(int(G[i, t]))[v]
        if t == 0 and first_candidates is not None:
            mask = np.zeros(N, bool); mask[first_candidates] = True; cand &= mask
        total = 0
        for v in np.nonzero(cand)[0]:
            total += rec(chosen + [int(v)], masks)
        return total
    return rec([], base)

def cartan(typ, k):
    C = np.zeros((k, k), dtype=np.int64)
    for i in range(k): C[i, i] = 2
    if typ == 'A':
        for i in range(k-1): C[i, i+1] = C[i+1, i] = -1
    elif typ == 'D':
        for i in range(k-2): C[i, i+1] = C[i+1, i] = -1
        C[k-3, k-1] = C[k-1, k-3] = -1
    elif typ == 'E':
        # Bourbaki: 1-3-4-5-6(-7-8), 2 attached to 4
        edges = [(1, 3), (3, 4), (4, 5), (5, 6), (2, 4)] + ([(6, 7)] if k >= 7 else []) + ([(7, 8)] if k >= 8 else [])
        for a, b in edges:
            C[a-1, b-1] = C[b-1, a-1] = -1
    return C
