"""Exact controls for neutral cycle-opening and delayed seam activation."""
from itertools import combinations, permutations
from functools import lru_cache
import json, pathlib

def tournaments(n):
    pairs = list(combinations(range(n), 2))
    for k in range(1 << len(pairs)):
        r = [0] * n
        for i, (a, b) in enumerate(pairs):
            if k >> i & 1:
                r[a] |= 1 << b
            else:
                r[b] |= 1 << a
        yield tuple(r)

def paths(r, S=None):
    if S is None:
        S = range(len(r))
    for p in permutations(S):
        if all((r[a] >> b & 1 for a, b in zip(p, p[1:]))):
            yield p

@lru_cache(None)
def H(r, mask=None):
    if mask is None:
        mask = (1 << len(r)) - 1
    if not mask:
        return 1
    d = {(1 << v, v): 1 for v in range(len(r)) if mask >> v & 1}
    for s in range(1, mask + 1):
        if s & ~mask:
            continue
        for v in range(len(r)):
            c = d.get((s, v), 0)
            if c:
                rem = r[v] & mask & ~s
                while rem:
                    bit = rem & -rem
                    rem -= bit
                    w = bit.bit_length() - 1
                    d[s | bit, w] = d.get((s | bit, w), 0) + c
    return sum((d.get((mask, v), 0) for v in range(len(r))))

def join(*blocks):
    out = []
    total = sum(map(len, blocks))
    off = 0
    for block in blocks:
        tail = (1 << total) - 1 ^ (1 << off + len(block)) - 1
        out += [x << off | tail for x in block]
        off += len(block)
    return tuple(out)

def tt(n):
    return tuple(((1 << n) - 1 ^ (1 << i + 1) - 1 for i in range(n)))

def extend(r, cut):
    n = len(r)
    return tuple((x | (0 if cut >> i & 1 else 1) << n for i, x in enumerate(r))) + (cut,)

def first_scc(r):
    n = len(r)
    reach = [x | 1 << i for i, x in enumerate(r)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return {v for v in range(n) if reach[v] == (1 << n) - 1}

def endpoint(r, x, start):
    result = 0
    n = len(r)
    for S in range(1, 1 << n):
        if not S >> x & 1:
            continue
        pcount = sum((1 for p in paths(r, [i for i in range(n) if S >> i & 1]) if p[0 if start else -1] == x))
        result += pcount * H(r, (1 << n) - 1 ^ S)
    return result

def main():
    neutral = positive = sectors = fixed_name_second = 0
    for n in range(1, 6):
        for r in tournaments(n):
            h = H(r)
            A = first_scc(r)
            for cut in range(1 << n):
                u = extend(r, cut)
                if H(u) != h:
                    continue
                neutral += 1
                found = {}
                for p in paths(u):
                    j = p.index(n)
                    if 0 < j < n and p[0] != n and (p[-1] != n):
                        v, t, a, b = (p[0], p[-1], p[j - 1], p[j + 1])
                        if len({v, t, a, b}) == 4 and r[a] >> v & 1:
                            key = (a, v, b, t)
                            found[key] = found.get(key, 0) + 1
                expected = len(A) >= 3 and n - len(A) >= 2 and (cut == sum((1 << i for i in range(n) if i not in A)))
                if not bool(found) == expected:
                    raise AssertionError((r, cut, found, A))
                if found:
                    positive += 1
                for (a, v, b, t), count in found.items():
                    if not ({a, v} <= A and b not in A and (t not in A)):
                        raise AssertionError('check failed')
                    ca = sum((p[0] == v and p[-1] == a for p in paths(r, A)))
                    pb = sum((p[0] == b and p[-1] == t for p in paths(r, set(range(n)) - A)))
                    if not count == ca * pb:
                        raise AssertionError('check failed')
                    sectors += 1
                if found:
                    for cut2 in range(1 << n + 1):
                        u2 = extend(u, cut2)
                        if H(u2) != h:
                            continue
                        for p in paths(u2):
                            j = p.index(n + 1)
                            if 0 < j < n + 1:
                                key = (p[j - 1], p[0], p[j + 1], p[-1])
                                if not key not in found:
                                    raise AssertionError((r, cut, cut2, key))
                        fixed_name_second += 1
    activation = 0
    for na in range(1, 4):
        for nb in range(1, 4):
            for A in tournaments(na):
                for B in tournaments(nb):
                    for a in range(na):
                        for b in range(nb):
                            L = endpoint(A, a, True)
                            R = endpoint(B, b, False)
                            for m in range(0, 4):
                                u = list(join(A, tt(m), B))
                                bb = na + m + b
                                u[a] ^= 1 << bb
                                u[bb] ^= 1 << a
                                u = tuple(u)
                                correction = 0 if m else sum((p[-1] == a for p in paths(A))) * sum((p[0] == b for p in paths(B)))
                                if not H(u) == H(A) * H(B) + (1 << m) * L * R - correction:
                                    raise AssertionError((A, B, a, b, m, H(u), L, R))
                                activation += 1
    A = (2, 4, 1)
    B = tt(2)
    rows = []
    for m in range(1, 7):
        u = join(A, tt(m), B)
        a = 0
        bb = len(u) - 1
        v = list(u)
        v[a] ^= 1 << bb
        v[bb] ^= 1 << a
        rows.append([m, H(u), H(tuple(v))])
        if not rows[-1][1:] == [3, 3 + 6 * (1 << m)]:
            raise AssertionError('check failed')
    result = {'neutral_births_orders_1_to_5': neutral, 'positive_neutral_openings': positive, 'marked_positive_sectors': sectors, 'fixed_name_second_neutral_births': fixed_name_second, 'two_sided_activation_cases': activation, 'C3_TT2_rows': rows, 'status': 'PASS'}
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
