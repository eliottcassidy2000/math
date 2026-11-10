#!/usr/bin/env python3
"""Check: the gate Gamma(R) = (a -> b -> R -> a) (Session 210; THM-4892; THM-4900) satisfies
H(Gamma R) = H(R) + P_R(2), P_R(2) = sum over ALL subsets S of H(R[S]) H(R[V - S]) (H(empty) = 1), and, by the
parity law (4.5), P_R(2) = 0 mod 4, i.e. H(Gamma R) = H(R) mod 4, for every tournament R on m >= 2 vertices.
All classes of orders 2..7 (gentourng)."""
import subprocess
def ham_all(n, out):
    """hS[S] = number of Hamiltonian paths of the subtournament induced on S (hS[0] = 1)"""
    inn = [sum(1 << w for w in range(n) if out[w] >> v & 1) for v in range(n)]
    f = [[0] * n for _ in range(1 << n)]; hS = [0] * (1 << n); hS[0] = 1
    for S in range(1, 1 << n):
        tot = 0
        for v in range(n):
            if not S >> v & 1: continue
            R = S & ~(1 << v)
            if R == 0: f[S][v] = 1; tot += 1; continue
            a = 0; mm = R & inn[v]
            while mm: w = (mm & -mm).bit_length() - 1; mm &= mm - 1; a += f[R][w]
            f[S][v] = a; tot += a
        hS[S] = tot
    return hS
bad = 0; tot = 0
for m in range(2, 8):
    res = subprocess.run(['/opt/homebrew/bin/gentourng', '-q', str(m)], capture_output=True, text=True).stdout.split()
    cnt = 0
    for line in res:
        out = [0] * m; k = 0
        for i in range(m):
            for j in range(i + 1, m):
                if line[k] == '1': out[i] |= 1 << j
                else: out[j] |= 1 << i
                k += 1
        hS = ham_all(m, out); full = (1 << m) - 1; H = hS[full]
        P2 = sum(hS[S] * hS[full ^ S] for S in range(1 << m))
        # gate: vertices 0..m-1 = R, a = m, b = m+1 ; a -> b, b -> R, R -> a
        og = [out[v] | (1 << m) for v in range(m)] + [1 << (m + 1), full]
        HG = ham_all(m + 2, og)[(1 << (m + 2)) - 1]
        ok = (HG == H + P2) and (P2 % 4 == 0) and ((HG - H) % 4 == 0)
        bad += (not ok); cnt += 1
    tot += cnt
    print("order %d: %d classes checked" % (m, cnt))
print("total %d classes; failures %d" % (tot, bad))
