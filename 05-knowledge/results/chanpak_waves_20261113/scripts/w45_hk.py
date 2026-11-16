"""wave 45 helpers: Hamiltonian path counts with start/end splits (bitmask DP)."""
def ham_split(n, out):
    """out[v] = bitmask of out-neighbours. returns (H, start[v], end[v])"""
    if n == 0:
        return 1, [], []
    full = (1 << n) - 1
    # forward dp: f[S][v] = paths covering S ending at v
    f = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        f[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = f[S]
        for v in range(n):
            c = row[v]
            if c:
                t = out[v] & ~S
                while t:
                    u = (t & -t).bit_length() - 1
                    t &= t - 1
                    f[S | 1 << u][u] += c
    end = f[full][:]
    H = sum(end)
    # start counts: reverse graph
    inn = [0] * n
    for v in range(n):
        t = out[v]
        while t:
            u = (t & -t).bit_length() - 1
            t &= t - 1
            inn[u] |= 1 << v
    g = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        g[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = g[S]
        for v in range(n):
            c = row[v]
            if c:
                t = inn[v] & ~S
                while t:
                    u = (t & -t).bit_length() - 1
                    t &= t - 1
                    g[S | 1 << u][u] += c
    start = g[full][:]
    return H, start, end

def H_of(n, out):
    return ham_split(n, out)[0]

def sub(out, keep):
    idx = {v: i for i, v in enumerate(keep)}
    o = [0] * len(keep)
    for v in keep:
        t = out[v]
        m = 0
        for u in keep:
            if t >> u & 1:
                m |= 1 << idx[u]
        o[idx[v]] = m
    return len(keep), o
