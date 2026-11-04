"""Independent small controls for the specific growing zero-channel cap."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent


def need(ok, why):
    if not ok:
        raise ValueError(why)


def birth(rows, target):
    n = len(rows)
    return [r | (0 if i == target else 1 << n) for i, r in enumerate(rows)] + [1 << target]


def paths(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask, ends in enumerate(dp):
        for i, count in enumerate(ends):
            rest = rows[i] & ~mask
            if not count:
                continue
            while rest:
                bit = rest & -rest
                rest ^= bit
                dp[mask | bit][bit.bit_length() - 1] += count
    return sum(dp[-1]), dp[-1][0]


def transitive(rows, mask):
    scores = [(rows[i] & mask).bit_count() for i in range(len(rows)) if mask >> i & 1]
    return sorted(scores) == list(range(len(scores)))


def main():
    born, born_names, aux = [0], ['x0'], 0
    records = []
    for j in range(7):
        if j:
            n = len(born)
            born = birth(birth(born, aux), n)
            born_names += [f'x{j}', f'y{j}']
            aux = n
        names = ['x0'] + [s for i in range(1, j + 1) for s in (f'y{i}', f'x{i}')]
        old = [born_names.index(s) for s in names]
        rows = [sum(1 << k for k, b in enumerate(old) if born[a] >> b & 1) for a in old]
        expected = [sum(1 << k for k in range(i + 1, len(names))) for i in range(len(names))]
        for i in range(1, j + 1):
            a, b = names.index(f'x{i-1}'), names.index(f'x{i}')
            expected[a] ^= 1 << b
            expected[b] ^= 1 << a
        need(rows == expected, 'birth graph versus ordered sparse arcs')
        h, endpoint = paths(rows)
        need((h, endpoint) == (3 ** j, 1), 'path/anchor invariants')
        n = len(rows)
        best = n - max(mask.bit_count() for mask in range(1 << n) if transitive(rows, mask))
        need(best == (j + 1) // 2, 'complete feedback minimum')
        removed = [names.index(f'x{i}') for i in range(1, j + 1, 2)]
        keep = ((1 << n) - 1) ^ sum(1 << v for v in removed)
        need(transitive(rows, keep), 'claimed achieving set')
        records.append(dict(depth=j, names=names, rows=rows, H=h, endpoint_x0=endpoint,
                            minimum_feedback=best, achieving_set=[names[i] for i in removed]))
    data = dict(status='FINITE-EXACT independent controls of proved cap-width law',
                universe='all deletion subsets at each depth0..6, orders1..13', records=records,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest())
    (HERE / 'cap_width.json').write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8', newline='\n')
    print('PASS: seven cap matrices, endpoint1, H=3^j, exact feedback ceil(j/2)')


if __name__ == '__main__':
    main()
