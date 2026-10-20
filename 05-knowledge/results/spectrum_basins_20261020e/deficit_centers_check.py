"""Two independent positive controls for the nonmonomial defect-six series."""
from itertools import combinations, permutations
from pathlib import Path
from hashlib import sha256
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def inspect(n, mask, expected):
    pairs = list(combinations(range(n), 2))
    adj = [[False] * n for _ in range(n)]
    for k, (i, j) in enumerate(pairs):
        adj[i][j] = bool(mask >> k & 1)
        adj[j][i] = not adj[i][j]
    for s in range(n):
        seen, todo = {s}, [s]
        while todo:
            u = todo.pop()
            for v in range(n):
                if adj[u][v] and v not in seen:
                    seen.add(v)
                    todo.append(v)
        need(len(seen) == n, 'not strong')
    aut, anti, paths, blue = [], [], [], []
    for p in permutations(range(n)):
        if all(adj[p[i]][p[j]] == adj[i][j] for i, j in pairs):
            aut.append(p)
        if all(adj[p[i]][p[j]] != adj[i][j] for i, j in pairs):
            anti.append(p)
        if all(adj[p[i]][p[i + 1]] for i in range(n - 1)):
            paths.append(p)
            if all(adj[p[i]][p[j]] == adj[p[n - 1 - j]][p[n - 1 - i]]
                   for i, j in pairs):
                blue.append(p)
    a, h_raw, b_raw = len(aut), len(paths), len(blue)
    need(bool(anti), 'not self-converse')
    need(h_raw % a == b_raw % a == 0, 'nonintegral orbit count')
    h, b = h_raw // a, b_raw // a
    need((h_raw, a, b, h - b) == expected, 'wrong invariant ledger')
    return dict(n=n, mask=mask,
                outmasks=[sum(int(adj[i][j]) << j for j in range(n))
                          for i in range(n)],
                H=h_raw, aut=a, raw_blue=b_raw, h=h, b=b, delta=h-b,
                anti_witness=anti[0], path_witness=paths[0])


if __name__ == '__main__':
    here = Path(__file__).resolve()
    result = {'status': 'PASS',
              'universe': 'two specified labelled tournaments; positive controls only',
              'source_sha256': sha256(here.read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
              'witnesses': [inspect(5, 10, (9, 1, 3, 6)),
                            inspect(6, 408, (45, 3, 9, 6))]}
    payload = json.dumps(result, indent=2) + '\n'
    here.with_suffix('.json').write_text(payload, encoding='utf-8', newline='\n')
    print(payload, end='')
