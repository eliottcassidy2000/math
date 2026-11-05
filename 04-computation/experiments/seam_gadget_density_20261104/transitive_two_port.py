"""Two-port transitive gadget transfers (the session's THM-4878 (A)), checked against SEAM's own matrix().

G = TT_k with vertex 0 the source (i -> j iff i < j), entry u, exit z, SEAM convention v = (H, S, R).
Closed form (proved by slot counting in THM-4878):
  H-row: (1, -[u = k-1], 2^(k-1-u))
  S-row: k = 1: (1, -1, 0); z = 0: (1, 0, 0) if u = 0 (k >= 2) else the H-row; z >= 1: (0, 0, 2^(k-1-z)) if u = 0 else 0
  R-row: (2^z, -mu, c) with mu = 2^u if z = k-1 > u, 2^z if u = k-1 > z, else 0;
         c = 2^(k-1-u+z) if z <= u, and 2^(k-z+u) 3^(z-u-1) if z > u.
Run from anywhere: python3 transitive_two_port.py [kmax]
"""
import sys, os
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', '..', '..', '05-knowledge', 'results', 'seven_return_native_recursions_20261029'))
from seam_transfer import matrix  # SEAM.md's universal seam transfer (independent construction)

def closed(k, u, z):
    H = [1, -1 if u == k - 1 else 0, 2 ** (k - 1 - u)]
    if k == 1:
        S = [1, -1, 0]
    elif z == 0:
        S = [1, 0, 0] if u == 0 else list(H)
    else:
        S = [0, 0, 2 ** (k - 1 - z)] if u == 0 else [0, 0, 0]
    if z == k - 1 and u < k - 1: mu = 2 ** u
    elif u == k - 1 and z < k - 1: mu = 2 ** z
    else: mu = 0
    c = 2 ** (k - 1 - u + z) if z <= u else 2 ** (k - z + u) * 3 ** (z - u - 1)
    return [H, S, [2 ** z, -mu, c]]

if __name__ == '__main__':
    def TT(k): return [[1 if i < j else 0 for j in range(k)] for i in range(k)]
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    checks = bad = 0
    for k in range(1, kmax + 1):
        for u in range(k):
            for z in range(k):
                checks += 1
                if matrix(TT(k), u, z) != closed(k, u, z):
                    bad += 1; print('MISMATCH', k, u, z, matrix(TT(k), u, z), closed(k, u, z))
    print(f'two-port transitive transfers: {checks} (k, u, z) cases for k <= {kmax}, mismatches {bad}')
