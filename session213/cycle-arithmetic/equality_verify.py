#!/usr/bin/env python3
"""Exact tests of endpoint-bound equality with a structural classifier."""
import itertools
import json


def check(ok, why):
    if not ok:
        raise RuntimeError(why)


def endpoints(a):
    n = len(a)
    d = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        d[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask >> v) & 1:
                continue
            old = mask ^ (1 << v)
            if old:
                d[mask][v] = sum(d[old][u] for u in range(n)
                                  if (old >> u) & 1 and (a[u] >> v) & 1)
    return d[-1]


def transpose(a):
    n = len(a)
    return [sum(((a[j] >> i) & 1) << j for j in range(n)) for i in range(n)]


def chain(a, vertices):
    mask = sum(1 << v for v in vertices)
    order = sorted(vertices, key=lambda v: -(a[v] & mask).bit_count())
    if all((a[u] >> v) & 1 for i, u in enumerate(order) for v in order[i + 1:]):
        return order
    return None


def reachability(a):
    r = [row | (1 << i) for i, row in enumerate(a)]
    for k in range(len(a)):
        for i in range(len(a)):
            if (r[i] >> k) & 1:
                r[i] |= r[k]
    return r


def recognized(a, v, reach):
    n = len(a)
    c = [w for w in range(n) if (reach[v] >> w) & 1 and (reach[w] >> v) & 1]
    if len(c) < 3:
        return False
    pre = [w for w in range(n) if w not in c]
    if chain(a, pre) is None or not all((a[p] >> w) & 1 for p in pre for w in c):
        return False
    out = [w for w in c if (a[v] >> w) & 1]
    if len(out) != 1:
        return False
    r = out[0]
    x = chain(a, [w for w in c if w not in (v, r)])
    if x is None or not x:
        return False
    return [w for w in x if (a[r] >> w) & 1] == [x[0]]


def tournaments(n):
    pairs = list(itertools.combinations(range(n), 2))
    for code in range(1 << len(pairs)):
        a = [0] * n
        for k, (v, w) in enumerate(pairs):
            if (code >> k) & 1:
                a[v] |= 1 << w
            else:
                a[w] |= 1 << v
        yield a


def near_transitive(m, j):
    n = m + j
    a = [sum(1 << w for w in range(v + 1, n)) for v in range(n)]
    a[j] ^= 1 << (n - 1)
    a[n - 1] |= 1 << j
    return a, n - 2


def is_path(a, p):
    return len(p) == len(a) and len(set(p)) == len(a) and all((a[u] >> v) & 1 for u, v in zip(p, p[1:]))


def main():
    report = {"status": "PASS", "scope": "FINITE-EXACT controls; all-order proof in EQUALITY.md",
              "census": []}
    for n in range(1, 7):
        ts = roots = positive_equalities = zero_equalities = 0
        for a in tournaments(n):
            ts += 1
            for orientation in (a, transpose(a)):
                end = endpoints(orientation)
                h = sum(end)
                reach = reachability(orientation)
                transitive = chain(orientation, list(range(n))) is not None
                check((h == 1) == transitive, ("unique path", n, orientation))
                for v in range(n):
                    if not orientation[v]:
                        continue
                    roots += 1
                    equality = h == 2 * end[v] + 1
                    structural = recognized(orientation, v, reach)
                    check((equality and end[v] > 0) == structural,
                          ("classification", n, orientation, v, h, end))
                    if end[v] == 0:
                        check(equality == transitive, ("zero boundary", n, orientation, v))
                        zero_equalities += equality
                    elif equality:
                        positive_equalities += 1
                        check(end[v] & (end[v] - 1) == 0, ("power of two", n, orientation, v))
                    elif end[v] & (end[v] - 1):
                        check(h >= 2 * end[v] + 3, ("nonpower floor", n, orientation, v))
        report["census"].append({"n": n, "labeled_tournaments": ts,
                                 "nonsink_roots_in_both_orientations": roots,
                                 "positive_equalities_both_orientations": positive_equalities,
                                 "zero_equalities_both_orientations": zero_equalities})
    orphan_cases = 0
    for ell in range(2, 9):
        for outmask in range(1 << ell):
            ys = [i for i in range(ell) if (outmask >> i) & 1]
            if len(ys) < 2:
                continue
            n = ell + 2
            aidx, ridx = ell, ell + 1
            a = [0] * n
            for i in range(ell):
                a[i] = sum(1 << k for k in range(i + 1, ell)) | (1 << aidx)
                if not (outmask >> i) & 1:
                    a[i] |= 1 << ridx
            a[aidx] = 1 << ridx
            a[ridx] = outmask
            i, j = ys[:2]
            p = [k for k in range(j + 1) if k != i] + [aidx, ridx, i] + list(range(j + 1, ell))
            check(is_path(a, p), ("orphan is not path", ell, outmask, p))
            reduced = [x for x in p if x != aidx]
            check(any(not((a[u] >> v) & 1) for u, v in zip(reduced, reduced[1:])),
                  ("orphan deletion remained path", ell, outmask))
            orphan_cases += 1
    family_cases = 0
    for m in range(3, 11):
        for j in range(3):
            a, v = near_transitive(m, j)
            end = endpoints(a)
            check(sum(end) == (1 << (m - 2)) + 1 and end[v] == 1 << (m - 3),
                  ("family counts", m, j))
            check(recognized(a, v, reachability(a)), ("family classification", m, j))
            family_cases += 1
    report["orphan_word_cases"] = orphan_cases
    report["near_transitive_family_cases"] = family_cases
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
