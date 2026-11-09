# Standalone independent audit; all checks remain active under python -O.
# No producer imports.
from itertools import combinations

def decode(c, n):
    a = [[0] * n for _ in range(n)]
    for b, (i, j) in enumerate(combinations(range(n), 2)):
        a[i][j] = c >> b & 1
        a[j][i] = 1 - a[i][j]
    return a

def require(condition, message):
    if not condition:
        raise ValueError(message)

def main():
    cases = [(135528708, [(6,), (5, 7), (2, 4), (0, 1), (3,)]), (36175876, [(5, 7), (4, 6), (1, 2), (0,), (3,)])]
    for c, blocks in cases:
        a = decode(c, 8)
        for i, B in enumerate(blocks):
            for j, C in enumerate(blocks):
                if i == j:
                    continue
                expected = i < j if i < 4 and j < 4 else j % 2 == 0 if i == 4 else i % 2 == 1
                require(all((a[b][d] == expected for b in B for d in C)), 'Check failed: all((a[b][d] == expected for b in B for d in C))')
        print(c, blocks, 'W5 sizes', list(map(len, blocks)))
    print('EXPLICIT FACTORS PASS')
if __name__ == '__main__':
    main()
