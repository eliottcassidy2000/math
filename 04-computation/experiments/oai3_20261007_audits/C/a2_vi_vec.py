# Vectorized independent VI for the Terras-clock chain on B(J,C); N = c*3^max(0,-j) integer representation.
# Formulas re-derived from the c-table: for j>=0 c=N; for j<0 c=N/3^d, d=-j.
import numpy as np, sys, time
J, C, SW = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
t0 = time.time()
off = {}; tot = 0
for j in range(-J, J + 1):
    R = C * 3**abs(j); off[j] = (tot, R); tot += 2*R + 1
n = tot
nxt = np.full((2, n), n, dtype=np.int64)   # n = OUT (value 0)
def index_of(jn, Nn):
    # jn, Nn arrays; returns index or n if outside
    res = np.full(len(Nn), n, dtype=np.int64)
    for jj in np.unique(jn):
        if abs(jj) > J: continue
        b, R = off[int(jj)]
        m = (jn == jj) & (np.abs(Nn) <= R)
        res[m] = b + Nn[m] + R
    return res
for j in range(-J, J + 1):
    b, R = off[j]
    N = np.arange(-R, R + 1, dtype=np.int64)
    e = N & 1                      # numpy & on negative int64 is two's complement: parity correct
    d = max(0, -j)
    for p in (0, 1):
        jn = np.full(len(N), j, dtype=np.int64); Nn = np.empty(len(N), dtype=np.int64)
        if p == 0:
            # e=0: c/2 -> N/2 ; e=1: j+1, c'=(3c+1)/2
            Nn = np.where(e == 0, N // 2, 0)
            if j >= 0:
                alt = (3*N + 1) // 2
            else:  # d>=1: N' = 3^(d-1) (3c+1)/2 = (N + 3^(d-1))/2
                alt = (N + 3**(d-1)) // 2
            Nn = np.where(e == 1, alt, Nn); jn = np.where(e == 1, j + 1, j)
        else:
            # e=0: c'=(3c+1-3^j)/2 ; e=1: j-1, c'=(c-3^(j-1))/2
            if j >= 0:
                a0 = (3*N + 1 - 3**j) // 2
            else:
                a0 = (3*N + 3**d - 1) // 2
            if j >= 1:
                a1 = (N - 3**(j-1)) // 2
            else:   # j<=0: d' = d+1, N' = 3^(d+1)(c - 3^(j-1))/2 = (3N - 1)/2
                a1 = (3*N - 1) // 2
            Nn = np.where(e == 0, a0, a1); jn = np.where(e == 0, j, j - 1)
        nxt[p, b:b + 2*R + 1] = index_of(jn, Nn)
z = off[0][0] + off[0][1]
print(f"B({J},{C}): {n} states built in {time.time()-t0:.1f}s", flush=True)
V = np.zeros(n + 1)
def at(j, N): b, R = off[j]; return b + N + R
small = np.array([at(j, N) for j in (-1, 0, 1) for N in range(-10*3**abs(j), 10*3**abs(j) + 1) if (j, N) != (0, 0)])
rep = {'(0,1)': at(0, 1), '(2,1)': at(2, 1), '(1,1)': at(1, 1)}
chk = sorted(set([SW // 4, SW]))
n0, n1 = nxt[0], nxt[1]
for it in range(1, SW + 1):
    Vn = V[n0]; Vn += V[n1]; Vn *= 0.5; Vn[z] = 1.0
    V[:n] = Vn
    if it in chk:
        print(f"  sweeps {it}: " + "; ".join(f"{k} >= {V[v]:.6f}" for k, v in rep.items()) + f"; min B(1,10) >= {V[small].min():.6f}   ({time.time()-t0:.0f}s)", flush=True)
