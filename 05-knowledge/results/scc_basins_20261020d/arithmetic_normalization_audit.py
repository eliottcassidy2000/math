"""Independent end-ledger and golden-normalization hostile controls.

No producer module is imported. The matrix control checks all standard-path
sink contexts through order seven. The word control checks every word in
{1,3,4} of vertex weight at most 24. Word steps are certified rewrites whose
ambient realizations are proved in arithmetic_interior.md and moves.md;
the word experiment does not itself certify every intermediate tournament.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def need(value, message):
    if not value:
        raise RuntimeError(message)


def arc(t, u, v):
    return bool(t[u] >> v & 1)


def join(a, b):
    m, n = len(a), len(b)
    return (tuple(x | (((1 << n)-1) << m) for x in a)
            + tuple(x << m for x in b))


def blue(t, p):
    n = len(t)
    return all(arc(t, p[i], p[j]) == arc(t, p[n-1-j], p[n-1-i])
               for i, j in combinations(range(n), 2))


def complement_nonpath(t, p):
    need(sorted(p) == list(range(len(t))), 'not a permutation')
    need(all(arc(t, u, v) for u, v in zip(p, p[1:])), 'absent path arc')
    retained = {frozenset((u, v)) for u, v in zip(p, p[1:])}
    result = list(t)
    for u, v in combinations(range(len(t)), 2):
        if frozenset((u, v)) not in retained:
            result[u] ^= 1 << v
            result[v] ^= 1 << u
    return tuple(result)


def end_ledger():
    b, d = (6, 4, 8, 3), (10, 0, 3, 6)
    counts = []
    for k in range(1, 8):
        free = [(u, v) for u, v in combinations(range(k-1), 2)
                if v != u+1]
        count = 0
        for mask in range(1 << len(free)):
            w = [0]*k
            for u in range(k-1):
                w[u] |= 1 << (k-1)
            for u in range(k-2):
                w[u] |= 1 << (u+1)
            for j, (u, v) in enumerate(free):
                if mask >> j & 1:
                    w[u] |= 1 << v
                else:
                    w[v] |= 1 << u
            state = join(b, w)
            wp = tuple(range(4, 4+k))
            words = [(0, 1, 2, 3)+wp, (3,)+wp+(2, 0, 1),
                     (3, 1)+wp+(2, 0), (2, 0, 3, 1)+wp]
            excess = []
            for p in words:
                excess.append(state[p[0]].bit_count()
                              + state[p[-1]].bit_count()-(k+3))
                need(not blue(state, p), 'end ledger used a blue path')
                state = complement_nonpath(state, p)
            need(excess == [-1, -k, k+1, 1], 'endpoint excess mismatch')
            need(state == join(d, w), 'end ledger failed literal restoration')
            count += 1
        counts.append(dict(context_order=k, presentations=count))
    direct = []
    for p in permutations(range(4)):
        if all(arc(b, u, v) for u, v in zip(p, p[1:])) and not blue(b, p):
            target = complement_nonpath(b, p)
            scores = sorted(x.bit_count() for x in target)
            need(scores in ([0, 2, 2, 2], [1, 1, 1, 3]), 'wrong R4 target')
            direct.append(dict(path=p, target=target))
    need(len(direct) == 2, 'R4 does not have exactly two black paths')
    return dict(census=counts, total=sum(x['presentations'] for x in counts),
                unpadded_R4_black_paths=direct)


def golden_words(bound, prefix=()):
    yield prefix
    for letter in (1, 3, 4):
        if letter <= bound:
            yield from golden_words(bound-letter, prefix+(letter,))


def active(w):
    return 4 in w or w != w[::-1]


def inventory(w):
    return sum(w), sum(x != 1 for x in w)


def normalize(w, branches):
    original = inventory(w)
    need(active(w), 'normalization starts pure blue')
    while 4 in w:
        i = w.index(4)
        a, b = w[:i], w[i+1:]
        right = a+(3, 1)+b
        left = a+(1, 3)+b
        need(active(right) or active(left), 'both replacements are pure blue')
        if a and b:
            target = right if active(right) else left
            branches['interior'] += 1
            branches['interior_alternate'] += int(not active(right))
        elif a:
            if active(right):
                target = right
                branches['right_ordinary'] += 1
            else:
                need(a[0] == 1, 'pure right target lacks singleton prefix')
                target = left
                branches['right_boundary_repair'] += 1
        elif b:
            if active(left):
                target = left
                branches['left_ordinary'] += 1
            else:
                need(b[-1] == 1, 'pure left target lacks singleton suffix')
                target = right
                branches['left_boundary_repair'] += 1
        else:
            target = left
            branches['unpadded_R4'] += 1
        need(active(target), 'dissolution reached pure blue')
        need(inventory(target) == original, 'rewrite changed n or f')
        need(target.count(4) == w.count(4)-1, 'R4 count did not decrease')
        w = target
    r, t = w.count(1), w.count(3)
    need(r > 0 and t > 0, 'active rigid inventory has an empty letter type')
    inv = sum(w[i] == 3 and w[j] == 1
              for i in range(len(w)) for j in range(i+1, len(w)))
    if 2*inv == r*t:
        branches['midpoint_nonpalindrome'] += 1
    down = 2*inv <= r*t
    pair = (3, 1) if down else (1, 3)
    while inv != (0 if down else r*t):
        i = next(i for i in range(len(w)-1) if w[i:i+2] == pair)
        w = w[:i]+pair[::-1]+w[i+2:]
        inv += -1 if down else 1
        need(active(w), 'inversion normalization hit a palindrome')
        need(inventory(w) == original, 'adjacent swap changed inventory')
        branches['rigid_swaps'] += 1
    canonical = min(w, w[::-1])
    need(canonical == (1,)*r+(3,)*t, 'wrong converse-merged normal form')
    return canonical


def word_control():
    branches = Counter()
    rows = [Counter() for _ in range(25)]
    strata = {}
    for w in golden_words(24):
        n, f = inventory(w)
        rows[n]['words'] += 1
        if active(w):
            rows[n]['active'] += 1
            need(1 <= f <= (n-1)//3, 'active inventory outside exact bounds')
            result = normalize(w, branches)
            if (n, f) in strata:
                need(strata[n, f] == result, 'normal form depends on starting word')
            strata[n, f] = result
        else:
            rows[n]['pure_blue'] += 1
    for n in range(25):
        need({f for nn, f in strata if nn == n} == set(range(1, (n-1)//3+1)),
             'a predicted active stratum is missing')
    need(all(branches[name] > 0 for name in (
        'interior_alternate', 'right_boundary_repair', 'left_boundary_repair',
        'midpoint_nonpalindrome', 'unpadded_R4')), 'hostile branch unexercised')
    return dict(vertex_weight_bound=24, ordinary_words=sum(x['words'] for x in rows),
                active_words=sum(x['active'] for x in rows),
                active_strata=len(strata), branches=dict(sorted(branches.items())),
                census=[dict(order=n, **row) for n, row in enumerate(rows)])


def main():
    report = dict(status='FINITE-EXACT independent controls; all-order proof audited separately',
                  end_ledger=end_ledger(), symbolic_normalization=word_control(),
                  source_sha256_lf=sha256(Path(__file__).read_text().encode()).hexdigest())
    output = json.dumps(report, indent=2)+'\n'
    (HERE/'arithmetic_normalization_audit.json').write_text(output, encoding='utf-8')
    print(output, end='')


if __name__ == '__main__':
    main()
