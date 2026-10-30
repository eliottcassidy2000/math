"""Exploratory exact signed-permutation census; not a proof certificate."""
from collections import Counter

def table(d):
    if not d:
        return [[1]]
    old = table(d-1)
    h = len(old)
    return [[(old[i % h][j % h] * (-1 if i % h else 1) ** (j // h)
              if i < h else old[j % h][i % h] * (-(-1 if i % h else 1)) ** (j // h))
             for j in range(2*h)] for i in range(2*h)]

def positive(t, i, j):
    step = i ^ j
    return sum(t[j][b]*t[i][b^step]*t[j][b^step]*t[i][b] == 1
               for b in range(len(t))) // 2

def fast(d, i, j):
    if i == j:
        raise ValueError('distinct indices required')
    if not i or not j:
        return 0
    h = 1 << (d-1)
    u, v = i % h, j % h
    if (i < h) == (j < h):
        return 2 * fast(d-1, u, v)
    if not u or not v or u == v:
        return 0
    return h - 4 - 2 * fast(d-1, u, v)

if __name__ == '__main__':
    for d in range(1, 8):
        t = table(d)
        hist = Counter(positive(t,i,j) for i in range(len(t)) for j in range(i+1,len(t)))
        for i in range(len(t)):
            for j in range(i+1,len(t)):
                if positive(t,i,j) != fast(d,i,j):
                    raise RuntimeError((d,i,j))
        print(d, sorted(hist.items()), flush=True)
        if d == 4:
            print('first defect', next((i,j) for i in range(16) for j in range(i+1,16) if positive(t,i,j)), flush=True)
