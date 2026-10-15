"""Six-vertex cores of the Paley run sums at 12 steps, and the K6-minus-P6 terminal.

U(G) = unrestricted character sum over all vertex maps (vertex 0 = 0, vertex 1 = 1, times q(q-1); even edge counts).
Claims checked: cores a, c, d vanish (odd-reversal automorphisms, THM-4652); the octahedron satisfies
U(K_{2,2,2}) = q(q-1)((5-q) a_q + 3q) with a_q the coefficients of eta(2z)^4 eta(4z)^4 (all odd primes 5..107);
K6-minus-P6 (tiling-modular-atoms terminal V_q, relabelled to keep edge (0,1)) is tabulated, not identified.
Session opus-2026-10-14-S35.
"""
import numpy as np, numba as nb, sys
CORES = {
 'K6-P6':        [(0,1),(0,3),(0,4),(0,5),(1,4),(1,5),(2,3),(2,4),(2,5),(3,4)],
 'a_K6-(P3+K2)': [(0,1),(0,2),(0,3),(0,4),(0,5),(1,2),(1,3),(1,4),(2,3),(2,4),(3,5),(4,5)],
 'b_octahedron': [(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,5),(2,4),(2,5),(3,4),(3,5),(4,5)],
 'c_K24+2K2':    [(0,1),(0,2),(0,3),(0,4),(1,2),(1,5),(2,5),(3,4),(3,5),(4,5)],
 'd_K33+e':      [(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,5),(3,5),(4,5)],
}

@nb.njit(cache=True)
def U6(q, chi, eu, ev):
    # vertex 0 = 0, vertex 1 = 1 (edge (0,1) present, even edge count): U = q (q-1) * S
    tot = 0
    v = np.zeros(6, dtype=np.int64)
    v[0] = 0; v[1] = 1
    for a in range(q):
        v[2] = a
        for b in range(q):
            v[3] = b
            for c in range(q):
                v[4] = c
                for d in range(q):
                    v[5] = d
                    p = 1
                    for k in range(eu.shape[0]):
                        p *= chi[(v[ev[k]] - v[eu[k]]) % q]
                        if p == 0:
                            break
                    tot += p
    return tot * q * (q - 1)

def chi_table(q):
    chi = np.zeros(q, dtype=np.int64); sq = set((x * x) % q for x in range(1, q))
    for x in range(1, q): chi[x] = 1 if x in sq else -1
    return chi


N = 220
P = [0] * (N + 1); P[0] = 1
for n in range(1, N + 1):
    for step in (2 * n, 4 * n):
        if step > N: continue
        for _ in range(4):
            for i in range(N, step - 1, -1):
                P[i] -= P[i - step]
a = [0] * (N + 2)
for i in range(N + 1): a[i + 1] = P[i]

if __name__ == '__main__':
    QS = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107]
    for name, E in CORES.items():
        eu = np.array([e[0] for e in E], dtype=np.int64); ev = np.array([e[1] for e in E], dtype=np.int64)
        line = []
        for q in QS:
            U = int(U6(q, chi_table(q), eu, ev))
            y = U // (q * (q - 1))
            assert U % (q * (q - 1)) == 0
            tag = ''
            if name == 'b_octahedron':
                tag = 'OK' if y == (5 - q) * a[q] + 3 * q else 'FAIL'
            line.append(f'{q}:{y}{"" if not tag else "/" + tag}')
        print(name, 'U/(q(q-1)) at primes:', ' '.join(line), flush=True)
