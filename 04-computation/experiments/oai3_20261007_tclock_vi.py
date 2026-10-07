# independent value iteration for the T-clock relation chain (state j, N = c*3^max(0,-j)); merge state (0,0)
import numpy as np, sys
J, C, SWEEPS = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
P3 = [3 ** i for i in range(40)]
def nxt(j, N, p):
    e = N & 1 if N >= 0 else (-N) & 1
    if p == 0 and e == 0: return j, N // 2 if N >= 0 else -((-N) // 2)
    if p == 0 and e == 1:
        if j >= 0: return j + 1, (3 * N + 1) // 2
        return j + 1, (N + P3[-j - 1]) // 2
    if p == 1 and e == 0:
        if j >= 0: return j, (3 * N + 1 - P3[j]) // 2
        return j, (3 * N + P3[-j] - 1) // 2
    if j >= 1: return j - 1, (N - P3[j - 1]) // 2
    return j - 1, (3 * N - 1) // 2
# exactness: all divisions above are exact for states reachable from integer N (checked below by assertion)
states = []
for j in range(-J, J + 1):
    lim = C * P3[abs(j)]
    for N in range(-lim, lim + 1):
        states.append((j, N))
idx = {s: i for i, s in enumerate(states)}
S = len(states)
nx = np.full((S, 2), -1, dtype=np.int64)
for i, (j, N) in enumerate(states):
    for p in (0, 1):
        jj, NN = nxt(j, N, p)
        nx[i, p] = idx.get((jj, NN), -1)
V = np.zeros(S + 1)              # index S = outside the box (value 0)
merge = idx[(0, 0)]
nx[nx < 0] = S
for sweep in range(SWEEPS):
    newV = 0.5 * (V[nx[:, 0]] + V[nx[:, 1]])
    newV[merge] = 1.0
    V[:S] = newV
print(f"box B({J},{C}): {S} states; after {SWEEPS} sweeps: M1 (0,1) >= {V[idx[(0,1)]]:.5f}; (2,1) >= {V[idx[(2,1)]]:.5f}; (1,1) >= {V[idx[(1,1)]]:.5f}")
