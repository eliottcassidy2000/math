"""Exact retained-mark reader and a bounded integer inverse.

Vertex rows are outgoing masks. A 1 appends a global sink; a 0 appends a
vertex beating only the retained mark. No graph is identified by its fields.
"""
from itertools import combinations


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def birth(rows, out):
    n = len(rows)
    return tuple(x | (0 if out >> i & 1 else 1 << n)
                 for i, x in enumerate(rows)) + (out,)


def strong(rows):
    n = len(rows)
    for start in range(n):
        seen = 1 << start
        while True:
            nxt = seen
            for i in range(n):
                if seen >> i & 1:
                    nxt |= rows[i]
            if nxt == seen:
                break
            seen = nxt
        if seen != (1 << n) - 1:
            return False
    return True


def tables(rows):
    n = len(rows)
    full = (1 << n) - 1
    ends = [[0] * n for _ in range(full + 1)]
    starts = [[0] * n for _ in range(full + 1)]
    hs = [0] * (full + 1)
    hs[0] = 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask >> v & 1):
                continue
            rest = mask ^ (1 << v)
            if not rest:
                ends[mask][v] = starts[mask][v] = 1
            else:
                ends[mask][v] = sum(ends[rest][u] for u in range(n)
                                    if rest >> u & 1 and rows[u] >> v & 1)
                starts[mask][v] = sum(starts[rest][u] for u in range(n)
                                      if rest >> u & 1 and rows[v] >> u & 1)
        hs[mask] = sum(ends[mask])
    return hs, ends, starts


def all_fields(rows):
    """Small literal subset DP, independent of the reader's recurrence."""
    n = len(rows)
    full = (1 << n) - 1
    hs, ends, starts = tables(rows)
    p2 = []
    for mask in range(full + 1):
        total = 0
        sub = mask
        while True:
            total += hs[sub] * hs[mask ^ sub]
            if not sub:
                break
            sub = (sub - 1) & mask
        p2.append(total)
    answer = []
    for w in range(n):
        S2 = sum(starts[m][w] * hs[full ^ m] for m in range(full + 1))
        F = sum(ends[m][w] * hs[full ^ m] for m in range(full + 1))
        S3 = sum(starts[m][w] * p2[full ^ m] for m in range(full + 1))
        answer.append((hs[full], ends[full][w], S2, F,
                       p2[full ^ (1 << w)], hs[full ^ (1 << w)],
                       p2[full], S3))
    return answer


def step(x, bit):
    h, e, t, F, p, q, P2, S3 = x
    if bit:
        return h, 0, 2*t, F, 2*p, q, 2*P2, 3*S3
    return h+t-e, q, 2*t-q, F+p, 2*p, q, 2*(P2-F+S3), 3*S3-p


def run(x, word):
    for bit in word:
        x = step(x, int(bit))
    return x


def graph(rows, root, word):
    for bit in word:
        rows = birth(rows, 0 if bit == '1' else 1 << root)
    return rows


def word_h(word):
    h, binary, last = 1, 0, 0
    for bit in map(int, word):
        h += (1-bit) * (binary+last)
        binary, last = 2*binary+bit, bit
    return h


def decode_r(r):
    need(r >= 1, 'Positive reversed zero mask')
    return ''.join('0' if r >> j & 1 else '1' for j in range(r.bit_length()))


def inverse(seed, target, include_empty=True):
    """All strong-ending words and all nonempty transitive prefix clocks.

    For strong outputs the seed must have at least two vertices and its
    mark must reach every seed vertex. Empty words are included only when
    the caller has a strong seed. Singleton inputs still obey the arithmetic
    formulas, but a zero-ending word alone is not their strongness test.
    Explicit finite bounds are proved in NOTE.
    """
    h, e, t, F, p, q, _, _ = seed
    need(0 <= e <= q and t >= q and p > 0 and F > 0, 'Seed bounds')
    hits = []
    tried = 0
    M, s = 2, 1
    while h + M*F <= target:
        residual = target-h-M*F
        if include_empty and residual == 0:
            hits.append(dict(word='', prefix=s, fields=list(seed)))
        K, d = t-q+M*p, q-e
        lo = max(1, -(-(residual-d) // (K+q)))
        hi = residual // K
        for r in range(lo, hi+1):
            tried += 1
            word = decode_r(r)
            value = K*r + q*(word_h(word)-1) + d*(r & 1)
            if value == residual:
                output = run(seed, word)
                need(output[0]+M*output[3] == target, 'Independent state replay')
                hits.append(dict(word=word, prefix=s, fields=list(output)))
        M, s = 2*M, s+1
    return hits, tried


def small_bank(maximum=5):
    """Every labelled strong seed through maximum, reduced by eight fields.

    A full graph representative is kept for each key. This quotient is
    valid for these reader operations and wrapper observations only.
    """
    bank = {}
    labelled = 0
    for n in range(3, maximum+1):
        pairs = list(combinations(range(n), 2))
        for code in range(1 << len(pairs)):
            rows = [0]*n
            for k, (a, b) in enumerate(pairs):
                u, v = (a, b) if code >> k & 1 else (b, a)
                rows[u] |= 1 << v
            rows = tuple(rows)
            if not strong(rows):
                continue
            labelled += 1
            for root, fields in enumerate(all_fields(rows)):
                bank.setdefault(fields, dict(rows=rows, root=root, fields=fields))
    return list(bank.values()), labelled


if __name__ == '__main__':
    import json
    targets = [170807, 512423, 512425, 1537271, 1537273, 1537277, 1537279]
    bank, labelled = small_bank()
    print('Small bank:', len(bank), 'eight-field keys;', labelled, 'labelled strong seeds')
    for target in targets:
        hits, tried = [], 0
        for seed in bank:
            found, count = inverse(seed['fields'], target)
            tried += count
            hits += [dict(seed=seed, **item) for item in found]
        hits.sort(key=lambda z: (len(z['seed']['rows'])+len(z['word'])+z['prefix']+1,
                                 len(z['word']), z['word'], z['prefix']))
        print(target, 'hits', len(hits), 'candidate integers', tried,
              'first', json.dumps(hits[:2]))
