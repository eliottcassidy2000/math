"""Independent small controls and deterministic replay of the three research lanes."""
from fractions import Fraction
from itertools import product, permutations
from pathlib import Path
import subprocess
import sys


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def main():
    root = Path(__file__).resolve().parents[2]
    for lane in ('graph', 'arithmetic', 'flow'):
        path = root / '04-computation/experiments' / ('crossroads10_20260926_' + lane + '.py')
        a = subprocess.check_output([sys.executable, '-B', str(path)], cwd=root).replace(b'\r\n', b'\n')
        b = subprocess.check_output([sys.executable, '-B', '-O', str(path)], cwd=root).replace(b'\r\n', b'\n')
        retained = (root / '05-knowledge/results' / (path.stem + '.out')).read_bytes().replace(b'\r\n', b'\n')
        require(a == b == retained, lane + ' replay mismatch')
        print(lane + ': normal = optimized = retained output')

    # Independent literal permutation path count: no imported lane code or DP.
    adj = [252,169,242,164,202,144,42,64]
    corners = []
    for e, f in ((0,0), (1,0), (0,1), (1,1)):
        x = list(adj)
        for bit, i, j in ((e,0,3), (f,2,3)):
            if bit:
                x[i] ^= 1 << j
                x[j] ^= 1 << i
        corners.append(sum(all((x[u] >> v) & 1 for u,v in zip(p,p[1:]))
                           for p in permutations(range(8))))
    require(corners == [233,291,123,223], 'response square')
    print('independent response square:', corners)

    # Prime rows 5,7,13 of (91,175,325), independently evaluated determinant.
    M = [[0,2,2], [1,1,0], [1,0,1]]
    determinant = sum((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
                      * M[0][p[0]]*M[1][p[1]]*M[2][p[2]]
                      for p in permutations(range(3)))
    require(determinant == -4, 'exclusive-prime triangle')
    print('orbit-27 prime triangle determinant:', determinant)

    # Full ten-bit legal universe and the rounding-only false gain.
    legal = []
    for bits in product((0,1), repeat=10):
        x = [0,0] + list(bits)
        if not all(x[i]+x[j]==1 for i,j in ((2,3),(4,6),(6,9))):
            continue
        if not all(x[i]<=x[j] for i,j in ((4,3),(3,5),(7,5),(5,8),(10,7),(7,11))):
            continue
        A = [sum(x[:i+1]) for i in range(12)]
        legal.append(A)
    vals = [min(A[11]+A[4] for A in legal), min(A[4]+A[1] for A in legal),
            min(A[11]+2*A[4]+A[1] for A in legal),
            min(A[11]+2*A[4]+A[2] for A in legal)]
    require(len(legal)==16 and vals==[3,1,4,5], 'ten-bit common-cutoff control')
    print('ten-bit legal assignments:', len(legal), 'costs:', vals)
    print('all independent controls PASS')


if __name__ == '__main__':
    main()
