"""Two new integer sequences from the digit-1 census (not in OEIS on 2026-11-04).

(1) F = empty: the digit-1 ready tournament whose entry attachment Y = Phi(T) is the pure backward cycle
    (v_i -> v_{i+1}, v_j -> v_i for j >= i+2).  Claim: T is the tower over a single vertex of height n-2 (the full tiling
    tournament T*_{n-1} of THM-4731 (6)) plus one more sole-exit vertex pointing at the BASE vertex x_0, so
    H = h_{n-2} - e_{x_0}(T_{n-2}) + B_{x_0}(T_{n-2}) by the one-exit law.  Both constructions are computed and compared.
(2) the digit-1 ready floors 5, 9, 17, 33, 61, 115, 219, 423, 823, 1615, 3183 (orders 4-14) are read off the census files
    of the sweep packet (data/gen1_order*.txt, data/extgen_order*.txt, data/floor_*.out) and are only listed here."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tower_check import ham_paths_ending, h_of, B_of, tower

def backcycle_phi_inverse(n):
    # Y: positions 0..n-1, i -> i+1 (and n-1 -> 0), j -> i for j >= i+2 except the pair (0, n-1) which is the cycle arc
    Y = [0] * n
    for i in range(n - 1): Y[i] |= 1 << (i + 1)
    Y[n - 1] |= 1
    for i in range(n):
        for j in range(i + 2, n):
            if not (i == 0 and j == n - 1): Y[j] |= 1 << i
    # T = reverse every arc not incident to a = 0
    T = [0] * n
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if i == 0 or j == 0:
                if Y[i] >> j & 1: T[i] |= 1 << j
            else:
                if Y[j] >> i & 1: T[i] |= 1 << j
    return T

seq_direct = []; seq_tower = []
for n in range(4, 14):
    T = backcycle_phi_inverse(n)
    seq_direct.append(h_of(T, range(n)))
    # tower over a point to height n-2 (order n-1), then re-root at the base vertex 0
    rows = [0]; target = 0
    for k in range(1, n - 1):
        x = len(rows); rows.append(0)
        for v in range(x):
            if v == target: rows[x] |= 1 << v
            else: rows[v] |= 1 << x
        target = x
    m = len(rows)
    h = h_of(rows, range(m)); e0 = ham_paths_ending(rows, range(m))[0]
    B0 = B_of(rows, range(m), 0) if m <= 10 else None
    if B0 is None:
        # compute H directly by adding the vertex
        rows2 = list(rows) + [0]
        for v in range(m):
            if v == 0: rows2[m] |= 1 << v
            else: rows2[v] |= 1 << m
        seq_tower.append(h_of(rows2, range(m + 1)))
    else:
        seq_tower.append(h - e0 + B0)
print("F = empty family, orders 4-13, direct:", seq_direct)
print("tower over a point re-rooted at the base, orders 4-13:", seq_tower)
print("agree:", seq_direct == seq_tower)
print("digit-1 ready floors, orders 4-14 (from the sweep packet's censuses): [5, 9, 17, 33, 61, 115, 219, 423, 823, 1615, 3183]")
