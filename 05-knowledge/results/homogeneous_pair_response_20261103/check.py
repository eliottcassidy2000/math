"""Exact homogeneous-pair identities and actual two-port separation controls."""
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import argparse
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def tournament(n, code):
    a = [0]*n
    for bit, (u, v) in enumerate(combinations(range(n), 2)):
        if code >> bit & 1:
            a[u] |= 1 << v
        else:
            a[v] |= 1 << u
    return a


def opposite(a):
    return [((1 << len(a))-1) ^ (1 << u) ^ row for u, row in enumerate(a)]


def ends(a):
    n = len(a)
    dp = [[0]*n for _ in range(1 << n)]
    for u in range(n):
        dp[1 << u][u] = 1
    for mask in range(1, 1 << n):
        for u, count in enumerate(dp[mask]):
            if count:
                todo = a[u] & ~mask
                while todo:
                    bit = todo & -todo
                    todo -= bit
                    dp[mask | bit][bit.bit_length()-1] += count
    return dp


def h(a):
    return sum(ends(a)[-1]) if a else 1


def literal_h(a):
    return sum(all(a[u] >> v & 1 for u, v in zip(word, word[1:]))
               for word in permutations(range(len(a))))


def sidecar(a):
    n, full = len(a), (1 << len(a))-1
    e, s = ends(a), ends(opposite(a))
    c = [[0]*n for _ in a]
    for mask in range(1, full):
        for u, eu in enumerate(e[mask]):
            if eu:
                for v, sv in enumerate(s[full ^ mask]):
                    c[u][v] += eu*sv
    return sum(e[-1]), s[-1], e[-1], c


def literal_sidecar(a):
    n = len(a)
    s, e, c = [0]*n, [0]*n, [[0]*n for _ in a]
    for word in permutations(range(n)):
        bad = [i for i in range(n-1) if not a[word[i]] >> word[i+1] & 1]
        if not bad:
            s[word[0]] += 1
            e[word[-1]] += 1
        for cut in range(1, n):
            if all(i == cut-1 for i in bad):
                c[word[cut-1]][word[cut]] += 1
    return sum(s), s, e, c


def inflate(q, module, root=0):
    """Module names first, then quotient vertices in their original order, omitting root."""
    exterior = [u for u in range(len(q)) if u != root]
    groups = {root: list(range(len(module)))}
    groups.update({u: [len(module)+i] for i, u in enumerate(exterior)})
    out = module[:] + [0]*len(exterior)
    for u in range(len(q)):
        for v in range(len(q)):
            if q[u] >> v & 1:
                for i in groups[u]:
                    for j in groups[v]:
                        out[i] |= 1 << j
    return out


def extend(a, x, y):
    n = len(a)
    out = a[:] + [1 << (n+1), 0]
    for j, mask in enumerate((x, y)):
        for u in range(n):
            if mask >> u & 1:
                out[u] |= 1 << (n+j)
            else:
                out[n+j] |= 1 << u
    return out


def calibrated(a, v):
    return h(extend(a, 1 << v, 0))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    source = Path(__file__).read_bytes()
    counts = dict(quotients=0, rooted_expansions=0, cross_identities=0,
                  calibrated_queries=0, literal_sidecars=0, literal_queries=0,
                  ambient_module_pairs=0)
    for n in range(1, 6):
        for code in range(1 << (n*(n-1)//2)):
            q = tournament(n, code)
            hq, sq, eq, cq = sidecar(q)
            counts['quotients'] += 1
            for root in range(n):
                t = inflate(q, [2, 0], root)
                ht, st, et, ct = sidecar(t)
                counts['rooted_expansions'] += 1
                need(st[0]-st[1] == sq[root], 'start difference')
                need(et[1]-et[0] == eq[root], 'end difference')
                need(ct[0][1] == ct[1][0] == hq, 'internal pair cover')
                for w, u in enumerate([v for v in range(n) if v != root], 2):
                    need(ct[1][w]-ct[0][w] == cq[root][u], 'outgoing cover difference')
                    need(ct[w][0]-ct[w][1] == cq[u][root], 'incoming cover difference')
                    counts['cross_identities'] += 2
                r = [et[v]+sum(ct[v]) for v in (0, 1)]
                expected_gap = eq[root]+sum(cq[root])
                need(r[1]-r[0] == expected_gap >= hq > 0, 'strict exact dominance')
                # All expanded orders <=5, and every rooted six-vertex expansion with root0.
                if n <= 4 or root == 0:
                    actual = [calibrated(t, v) for v in (0, 1)]
                    need(actual == [ht+2*x for x in r], 'actual calibrated two-port formula')
                    counts['calibrated_queries'] += 2
                    if n <= 3:
                        need(actual == [literal_h(extend(t, 1 << v, 0)) for v in (0, 1)],
                             'literal calibrated queries')
                        counts['literal_queries'] += 2
                if n <= 4 or (root == 0 and code < 8):
                    need((ht, st, et, ct) == literal_sidecar(t), 'literal cover/path state')
                    counts['literal_sidecars'] += 1
    # The reversed module has a reversed homogeneous pair inside every exterior context.
    for qn in range(1, 4):
        for qc in range(1 << (qn*(qn-1)//2)):
            m = inflate(tournament(qn, qc), [2, 0])
            for rn in range(1, 4):
                for rc in range(1 << (rn*(rn-1)//2)):
                    r = tournament(rn, rc)
                    a, b = inflate(r, m), inflate(r, opposite(m))
                    f = [calibrated(a, v) for v in (0, 1)]
                    g = [calibrated(b, v) for v in (0, 1)]
                    need(f[1] > f[0] and g[1] < g[0] and f != g, 'ambient separation')
                    counts['ambient_module_pairs'] += 1
    # Endpoints alone vanish on a pair inside a middle ordered component.
    q = [2, 4, 0]  # 0->1->2 and 0->2 added below, with root1.
    q[0] |= 4
    t = inflate(q, [2, 0], 1)
    u = t[:]
    u[0] ^= 2
    u[1] ^= 1
    dt, du = sidecar(t), sidecar(u)
    need(dt[:3] == du[:3], 'endpoint-only hostile equality')
    f, g = [calibrated(t, v) for v in (0, 1)], [calibrated(u, v) for v in (0, 1)]
    need(f != g, 'endpoint-only hostile separation')
    need(Path(__file__).read_bytes() == source, 'source changed')
    result = dict(status='FINITE-EXACT controls for all-order pair separation', counts=counts,
                  universe='All1099 labelled quotients orders1..5, every distinguished root; '
                           'every labelled quotient/module context orders1..3 in the declared nested construction',
                  endpoint_only_hostile=dict(rows=t, changed_rows=u, common_HSE=dt[:3],
                                            pair_queries=f, reversed_pair_queries=g),
                  source_sha256=sha256(source).hexdigest())
    text = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text, encoding='utf8', newline='\n')
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
