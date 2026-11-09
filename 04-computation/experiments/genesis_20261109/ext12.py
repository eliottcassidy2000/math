"""Order-12 test of THM-4942's five-operation grammar: extend each of the seven 3-arc-strong order-11
NON-creators (kind-pasteur S17, lambda3_classes_order11.txt, status inS10) by a 12th vertex in all 2^11 ways,
and classify every extension with grammar.classify.  Prints escapers (if any) and a summary.
S17 code convention: bit b over pairs (i, j), i < j, row-major; bit = 1 iff i -> j.
"""
import os
import sys
from grammar import profile, classify

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '..', '..', '05-knowledge', 'results', 'first_creator_generation_20261108',
                   'lambda3_classes_order11.txt')


def decode(c, n):
    adj = [[False] * n for _ in range(n)]
    b = 0
    for i in range(n):
        for j in range(i + 1, n):
            if c >> b & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
            b += 1
    return adj


def code_str(adj):
    n = len(adj)
    return ''.join('1' if adj[i][j] else '0' for i in range(n) for j in range(i + 1, n))


def main():
    parents = []
    with open(SRC) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.split()
            if len(p) >= 3 and p[2] == 'inS10':
                parents.append((int(p[0]), int(p[1])))
    print('# parents (3-arc-strong order-11 non-creators):', len(parents))
    total = 0
    esc = []
    stat = {}
    for code, h in parents:
        a11 = decode(code, 11)
        items = []
        for x in range(1 << 11):
            adj = [row[:] + [False] for row in a11] + [[False] * 12]
            for a in range(11):
                if x >> a & 1:
                    adj[a][11] = True
                else:
                    adj[11][a] = True
            items.append((12, code_str(adj)))
        for p in profile(items):
            labels, st = classify(p)
            total += 1
            key = (st, ','.join(labels))
            stat[key] = stat.get(key, 0) + 1
            if st != 'generated':
                esc.append((code, h, p, st))
    print('# extensions classified:', total)
    for k in sorted(stat):
        print('#', k, stat[k])
    for code, h, p, st in esc:
        print(st, 'parent', code, h, 'H', p['H'], 'lam', p['lam'], 'prime', p['prime'], 'dels', p['dels'], p['code'])


if __name__ == '__main__':
    main()
