# Two-block cut-vertex family: T = z + (B1 => B2), B1, B2 strong extremal blocks, z attached by every pattern (strong only).
# Writes gentourng-style ascii lines (vertex 0 = z, then B1, then B2) to stdout for ./tbal n.
import sys
from check_minimiser import dec
def lines(b1, b2):
    A1, A2 = dec(b1), dec(b2); n1, n2 = len(A1), len(A2); n = 1 + n1 + n2
    for pat in range(1 << (n1 + n2)):
        A = [[False] * n for _ in range(n)]
        for i in range(n1):
            for j in range(n1): A[1 + i][1 + j] = A1[i][j]
        for i in range(n2):
            for j in range(n2): A[1 + n1 + i][1 + n1 + j] = A2[i][j]
        for i in range(n1):
            for j in range(n2): A[1 + i][1 + n1 + j] = True
        for x in range(1, n):
            if pat >> (x - 1) & 1: A[0][x] = True
            else: A[x][0] = True
        if not any(A[0][x] for x in range(1, 1 + n1)) or not any(A[x][0] for x in range(1 + n1, n)): continue
        yield ''.join('1' if A[i][j] else '0' for i in range(n) for j in range(i + 1, n))
if __name__ == '__main__':
    for L in lines(sys.argv[1], sys.argv[2]): print(L)
