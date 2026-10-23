"""Exhaustive fixed-path coordinates, with independent involution orbits."""
from pathlib import Path
from hashlib import sha256
import argparse
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(why)


def positions(n):
    return [(i, j) for j in range(2, n) for i in range(j-1)]


def source(n, mask):
    a = [[0]*n for _ in range(n)]
    for i in range(n-1):
        a[i][i+1], a[i+1][i] = 1, -1
    for b, (i, j) in enumerate(positions(n)):
        a[i][j] = 1-2*((mask >> b) & 1)
        a[j][i] = -a[i][j]
    return a


def curls(a):
    return {(i, j): 1+a[i+1][j]-a[i][j] for i, j in positions(len(a))}


def decode(n, d):
    need(type(n) is int and n >= 3, 'literal order')
    need(type(d) is dict and set(d) == set(positions(n)), 'complete labelled array')
    need(all(type(p) is tuple and len(p) == 2 and all(type(i) is int for i in p)
             for p in d), 'literal position labels')
    need(all(type(v) is int for v in d.values()), 'literal integer curls')
    a = source(n, 0)
    for j in range(2, n):
        total = 0
        for i in range(j-2, -1, -1):
            total += d[i, j]
            v = j-i-total
            need(v in (-1, 1), 'column capacity exceeded')
            a[i][j], a[j][i] = v, -v
    return a


def rejected(n, d):
    try:
        decode(n, d)
    except ValueError:
        return
    raise AssertionError('invalid array accepted')


def main():
    rows = []
    for n in range(3, 8):
        pos = positions(n); m = len(pos); full = (1 << m)-1
        idx = {p:b for b, p in enumerate(pos)}
        perm = [idx[n-1-j, n-1-i] for i, j in pos]
        fixed = sum(b == c for b, c in enumerate(perm))
        need(fixed == (n-1)//2, 'antidiagonal count')
        blue = 0; orbits = set(); sigs = set()
        for mask in range(1 << m):
            a = source(n, mask); d = curls(a)
            need(decode(n, d) == a, 'suffix inverse')
            sig = tuple(d[p] for p in pos)
            need(sig not in sigs, 'injective coordinates'); sigs.add(sig)
            other = source(n, full ^ mask)
            transport = {p:(4 if p[1] == p[0]+2 else 2)-d[p] for p in pos}
            need(curls(other) == transport, 'direct affine complement')
            need(decode(n, transport) == other, 'transport inversion')
            reflected = sum(((mask >> b) & 1) << c for b, c in enumerate(perm))
            direct = [[-a[n-1-i][n-1-j] for j in range(n)] for i in range(n)]
            need(source(n, reflected) == direct, 'path reversal plus converse')
            is_blue = mask == reflected
            blue += is_blue
            orbit = tuple(sorted({mask, full ^ mask, reflected, full ^ reflected}))
            need(len(orbit) == (2 if is_blue else 4), 'exact orbit size')
            need(mask != (full ^ reflected), 'no complement-reflection fixed point')
            orbits.add(orbit)
        k = (m+fixed)//2
        b2 = sum(len(o) == 2 for o in orbits)
        b4 = sum(len(o) == 4 for o in orbits)
        need(blue == 2**k and b2 == 2**(k-1) and b4 == ((1 << m)-blue)//4,
             'independent orbit totals')
        rows.append({'n':n, 'tilings':1 << m, 'blue_tilings':blue,
                     'fundamental_blue_lines':b2, 'fundamental_black_lines':b4})
    invalid = {(0,2):1, (0,3):3, (1,3):3}
    rejected(4, invalid)
    good = curls(source(4, 0))
    for value in (True, 1.0, '1'):
        bad = dict(good); bad[0,2] = value; rejected(4, bad)
    bad = dict(good); del bad[0,2]; rejected(4, bad)
    bad = dict(good); bad[0,4] = 1; rejected(4, bad)
    for label in ((False, 2), (0.0, 2)):
        bad = dict(good); del bad[0,2]; bad[label] = 1; rejected(4, bad)
    result = {'scope':'all fixed-path masks n=3..7', 'rows':rows,
              'hostile':{'d02':1,'d03':3,'d13':3,'impossible_a03':-3}, 'checks':CHECKS}
    raw = (json.dumps(result, indent=2, sort_keys=True)+'\n').encode()
    args = argparse.ArgumentParser(); args.add_argument('--write', action='store_true')
    if args.parse_args().write:
        (HERE/'receipt.json').write_bytes(raw)
    else:
        need((HERE/'receipt.json').read_bytes().replace(b'\r\n', b'\n') == raw, 'saved receipt')
    print(json.dumps({'checks':CHECKS, 'sha256_lf':sha256(raw).hexdigest(), 'rows':rows}))


if __name__ == '__main__':
    main()
