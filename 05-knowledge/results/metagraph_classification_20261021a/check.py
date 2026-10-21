"""Source-bound checks for the all-order three-presentation classification.

The analytic proof imports the audited bypass, pure-blue, and deficit-two
theorems. Finite controls below do not replace those inputs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools as it
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ValueError(why)


def integer(x, lo=0):
    if type(x) is not int or x < lo:
        raise ValueError('exact integer domain')
    return x


def tournament(t):
    if type(t) is not tuple:
        raise ValueError('immutable labelled source required')
    n = len(t)
    for i, r in enumerate(t):
        integer(r)
        if r >= 1 << n or r >> i & 1:
            raise ValueError('row range/diagonal')
        for j in range(i):
            if (r >> j & 1) + (t[j] >> i & 1) != 1:
                raise ValueError('not a tournament')
    return t


def arc(t, i, j):
    return bool(t[i] >> j & 1)


def rows(code, n):
    integer(n); integer(code)
    if code >= 1 << (n * (n - 1) // 2):
        raise ValueError('code range')
    a = [0] * n
    for bit, (i, j) in enumerate(it.combinations(range(n), 2)):
        x, y = (i, j) if code >> bit & 1 else (j, i)
        a[x] |= 1 << y
    return tuple(a)


def census_rows(code, n):
    """Inherited TSV uses most-significant-first triangular bits."""
    width = n * (n-1) // 2
    reversed_code = sum(((code >> (width-1-i)) & 1) << i for i in range(width))
    return rows(reversed_code, n)


def local(t, vertices):
    return tuple(sum(int(arc(t, x, y)) << j for j, y in enumerate(vertices))
                 for x in vertices)


def sccs(t):
    tournament(t)
    n = len(t)
    reach = [t[i] | 1 << i for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    todo, result = set(range(n)), []
    while todo:
        x = min(todo)
        c = tuple(y for y in sorted(todo) if reach[x] >> y & 1 and reach[y] >> x & 1)
        result.append(c); todo.difference_update(c)
    # sorted creates a separate list; list.sort would temporarily empty the
    # very container needed by this key, losing the ordered SCC sidecar.
    return tuple(sorted(result, key=lambda c: -sum(arc(t, c[0], d[0]) for d in result if c != d)))


def cyclic(n):
    integer(n, 1)
    if not n % 2:
        raise ValueError('odd carousel required')
    return tuple(sum(1 << ((i + j) % n) for j in range(1, (n + 1) // 2))
                 for i in range(n))


def join(parts):
    out, offset = [], 0
    size = sum(map(len, parts))
    for p in parts:
        tournament(p)
        higher = ((1 << size) - 1) ^ ((1 << (offset + len(p))) - 1)
        out.extend((r << offset) | higher for r in p)
        offset += len(p)
    return tuple(out)


def stack(word):
    if type(word) is not tuple or any(type(a) is not int or a not in (1, 3, 5) for a in word):
        raise ValueError('exact atom word required')
    return join(tuple(cyclic(a) for a in word))


def relabel(t, p):
    tournament(t)
    if type(p) is not tuple or any(type(x) is not int for x in p) or sorted(p) != list(range(len(t))):
        raise ValueError('literal permutation required')
    return local(t, p)


def path_ok(t, p):
    return (type(p) is tuple and all(type(x) is int for x in p)
            and sorted(p) == list(range(len(t)))
            and all(arc(t, a, b) for a, b in zip(p, p[1:])))


def is_blue(t, p):
    if not path_ok(t, p):
        raise ValueError('not a source Hamilton word')
    rev = dict(zip(p, reversed(p)))
    return all(arc(t, rev[a], rev[b]) != arc(t, a, b)
               for a, b in it.combinations(range(len(t)), 2))


def black_target(t, p):
    tournament(t)
    if not path_ok(t, p):
        raise ValueError('wrong source word')
    kept = {frozenset((a, b)) for a, b in zip(p, p[1:])}
    a = [0] * len(t)
    for x, y in it.combinations(range(len(t)), 2):
        forward = arc(t, x, y) ^ (frozenset((x, y)) not in kept)
        i, j = (x, y) if forward else (y, x)
        a[i] |= 1 << j
    return tuple(a)


def hamilton_paths(t):
    n = len(t)
    def grow(p, used):
        if len(p) == n:
            yield p
        else:
            for x in range(n):
                if not used >> x & 1 and (not p or arc(t, p[-1], x)):
                    yield from grow(p + (x,), used | 1 << x)
    return grow((), 0)


def invariants(t):
    """Independent labelled paths and score-compatible automorphism enumeration."""
    tournament(t)
    groups = {}
    for i, r in enumerate(t):
        groups.setdefault(r.bit_count(), []).append(i)
    ordered = list(groups.values())
    aut = 0
    for images in it.product(*(it.permutations(g) for g in ordered)):
        p = list(range(len(t)))
        for g, image in zip(ordered, images):
            for x, y in zip(g, image):
                p[x] = y
        aut += all(arc(t, p[i], p[j]) == arc(t, i, j)
                   for i, j in it.combinations(range(len(t)), 2))
    h = b = 0
    for p in hamilton_paths(t):
        h += 1; b += is_blue(t, p)
    need(h % aut == b % aut == 0, 'free action divisions')
    return h, aut, h // aut, b // aut


def atom_frame(t, c):
    n = len(c)
    if n not in (1, 3, 5):
        return None
    for p in it.permutations(c):
        if local(t, p) == cyclic(n):
            return p
    return None


def source_digest(t):
    return hashlib.sha256(json.dumps(t, separators=(',', ':')).encode()).hexdigest()


def decode_three(t):
    """Exact iff compiler, justified by the analytic theorem, not a census oracle."""
    tournament(t)
    components = sccs(t)
    if sum(len(c) == 5 for c in components) != 1:
        return None
    frames = tuple(atom_frame(t, c) for c in components)
    if any(p is None for p in frames):
        return None
    return {'source': t, 'digest': source_digest(t), 'frames': frames}


def audit_three(t, rec):
    tournament(t)
    if type(rec) is not dict or set(rec) != {'source', 'digest', 'frames'}:
        raise ValueError('receipt schema')
    tournament(rec['source'])
    if rec['source'] != t or rec['digest'] != source_digest(t):
        raise ValueError('source substitution')
    fs = rec['frames']
    if type(fs) is not tuple or any(type(p) is not tuple or any(type(x) is not int for x in p) for p in fs):
        raise ValueError('exact labelled frames required')
    flat = tuple(x for p in fs for x in p)
    if sorted(flat) != list(range(len(t))) or sum(len(p) == 5 for p in fs) != 1:
        raise ValueError('frame partition')
    if any(len(p) not in (1, 3, 5) or local(t, p) != cyclic(len(p)) for p in fs):
        raise ValueError('atom incidence')
    if any(not arc(t, x, y) for i, p in enumerate(fs) for q in fs[i+1:] for x in p for y in q):
        raise ValueError('ordered cross incidence')
    return tuple(map(len, fs))


def port_receipts(t, rec):
    """The three actual orbit representatives; keeps source/word/target/ports."""
    audit_three(t, rec)
    local_reps = tuple(p for p in hamilton_paths(cyclic(5)) if p[0] == 0)
    ans = []
    for rp in local_reps:
        p = tuple(x for f in rec['frames'] for x in
                  (tuple(f[j] for j in rp) if len(f) == 5 else f))
        ans.append({'source': t, 'word': p, 'ports': (p[0], p[-1]),
                    'target': black_target(t, p), 'blue': is_blue(t, p)})
    return tuple(ans)


def audit_port(t, rec):
    tournament(t)
    if type(rec) is not dict or set(rec) != {'source', 'word', 'ports', 'target', 'blue'}:
        raise ValueError('port schema')
    tournament(rec['source']); tournament(rec['target'])
    if rec['source'] != t or not path_ok(t, rec['word']):
        raise ValueError('literal source/word mismatch')
    p = rec['word']
    if type(rec['ports']) is not tuple or any(type(x) is not int for x in rec['ports']):
        raise ValueError('port type')
    if rec['ports'] != (p[0], p[-1]) or type(rec['blue']) is not bool:
        raise ValueError('port/colour sidecar')
    if rec['target'] != black_target(t, p) or rec['blue'] != is_blue(t, p):
        raise ValueError('target/colour replay')
    return True


@lru_cache(None)
def _cycle(t):
    # Private cache, all callers supply an authenticated tuple.
    for p in hamilton_paths(t):
        if p and p[0] == 0 and arc(t, p[-1], p[0]):
            return p
    raise ValueError('no cycle')


def cut_cycle(c, first=None, last=None):
    if (first is None) == (last is None):
        raise ValueError('one marked endpoint')
    j = c.index(first) if first is not None else (c.index(last) + 1) % len(c)
    return c[j:] + c[:j]


def four_patterns(t):
    """Certificate for a Hamilton frame whose step-two rotation is an automorphism."""
    tournament(t)
    n = len(t)
    if n < 6 or n % 4 != 2:
        raise ValueError('two odd orbits of size at least three required')
    if not all(arc(t, i, (i + 1) % n) for i in range(n)):
        raise ValueError('not a marked Hamilton cycle')
    if not all(arc(t, i, j) == arc(t, (i + 2) % n, (j + 2) % n)
               for i, j in it.combinations(range(n), 2)):
        raise ValueError('missing step-two symmetry')
    a, b = tuple(range(0, n, 2)), tuple(range(1, n, 2))
    ca = tuple(a[j] for j in _cycle(local(t, a)))
    cb = tuple(b[j] for j in _cycle(local(t, b)))
    ab = next((x, y) for x in a for y in b if arc(t, x, y))
    ba = next((y, x) for y in b for x in a if arc(t, y, x))
    words = (tuple(range(n)), tuple(range(1, n)) + (0,),
             cut_cycle(ca, last=ab[0]) + cut_cycle(cb, first=ab[1]),
             cut_cycle(cb, last=ba[0]) + cut_cycle(ca, first=ba[1]))
    return {'source': t, 'orbits': (a, b), 'cycles': (ca, cb), 'words': words}


def audit_four(t, rec):
    tournament(t)
    if type(rec) is not dict or set(rec) != {'source', 'orbits', 'cycles', 'words'}:
        raise ValueError('four-pattern schema')
    tournament(rec['source'])
    if rec['source'] != t:
        raise ValueError('four-pattern source substitution')
    for key in ('orbits', 'cycles', 'words'):
        if type(rec[key]) is not tuple or any(type(p) is not tuple or any(type(x) is not int for x in p)
                                              for p in rec[key]):
            raise ValueError('exact literal words required')
    if tuple(map(len, (rec['orbits'], rec['cycles'], rec['words']))) != (2, 2, 4):
        raise ValueError('four-pattern arity')
    n = len(t)
    a, b = rec['orbits']; ca, cb = rec['cycles']
    need(n >= 6 and n % 4 == 2, 'odd-orbit boundary')
    need(a == tuple(range(0, n, 2)) and b == tuple(range(1, n, 2)), 'marked orbit ownership')
    need(sorted(ca) == list(a) and sorted(cb) == list(b), 'literal simple internal cycles')
    need(all(arc(t, c[i], c[(i+1) % len(c)]) for c in (ca, cb) for i in range(len(c))), 'internal cycles')
    sa = {t[x].bit_count() for x in a}; sb = {t[x].bit_count() for x in b}
    need(len(sa) == len(sb) == 1 and sa != sb, 'intrinsic distinct score classes')
    need(all(path_ok(t, p) for p in rec['words']), 'four actual Hamilton words')
    patterns = [tuple(t[x].bit_count() for x in p) for p in rec['words']]
    need(len(set(patterns)) == 4, 'automorphism-invariant separation')
    return True


def circulant(m, code):
    a = [0] * m
    for d in range(1, (m + 1) // 2):
        e = d if code >> (d - 1) & 1 else m - d
        for i in range(m):
            a[i] |= 1 << ((i + e) % m)
    return tuple(a)


def double_frame(m, ia, ib, cross):
    a, b = circulant(m, ia), circulant(m, ib)
    t = [0] * (2*m)
    for i in range(m):
        for j in range(m):
            if arc(a, i, j): t[2*i] |= 1 << (2*j)
            if arc(b, i, j): t[2*i+1] |= 1 << (2*j+1)
            if cross >> ((j-i) % m) & 1: t[2*i] |= 1 << (2*j+1)
            else: t[2*j+1] |= 1 << (2*i)
    return tuple(t)


def series(N):
    integer(N)
    a = [0] * (N+1); a[0] = 1
    for n in range(1, N+1):
        a[n] = a[n-1] + (a[n-3] if n >= 3 else 0)
    f = [sum(a[i]*a[n-5-i] for i in range(n-4)) if n >= 5 else 0 for n in range(N+1)]
    sc = [a[(n-5)//2] if n >= 5 and n % 2 else 0 for n in range(N+1)]
    return a, f, sc


def atom_words(n, alphabet=(1, 3, 5)):
    if n == 0:
        yield ()
    for a in alphabet:
        if a <= n:
            for w in atom_words(n-a, alphabet):
                yield (a,) + w


def expect_reject(fn):
    try:
        fn()
    except ValueError:
        need(True, 'expected rejection'); return
    need(False, 'bad input accepted')


def verify():
    global GATES
    GATES = 0
    labelled = []
    for n in range(1, 6):
        hits = strong_hits = 0
        for code in range(1 << (n*(n-1)//2)):
            t = rows(code, n)
            h, aut, tc, blue = invariants(t)
            r = decode_three(t)
            need((r is not None) == (tc == 3), 'independent small iff')
            if r:
                need(audit_three(t, r) == (5,), 'small unique strong atom')
                hits += 1; strong_hits += len(sccs(t)) == 1
        labelled.append({'n': n, 'tournaments': 1 << (n*(n-1)//2), 'tc3': hits, 'strong_tc3': strong_hits})

    frames = []
    sample = None
    for m in (3, 5, 7):
        total = 0
        for ia in range(1 << ((m-1)//2)):
            for ib in range(1 << ((m-1)//2)):
                for free in range(1 << (m-2)):
                    t = double_frame(m, ia, ib, 1 | (free << 1))
                    rec = four_patterns(t); audit_four(t, rec)
                    if sample is None: sample = rec
                    total += 1
        need(total == 1 << (2*m-3), 'complete marked two-orbit universe')
        frames.append({'orbit_size': m, 'n': 2*m, 'frames': total})

    # Existing ordinary-isomorphism census is a finite cross-check, not a proof premise.
    census = []
    expected = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}
    _, f, sc = series(40)
    for n in range(3, 9):
        p = ROOT / f'04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n{n}.tsv'
        stream = p.open(encoding='utf-8')
        count = hits = fixed = 0
        with stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                count += 1
                tc = int(row['tc'])
                if tc != 3:
                    continue
                t = census_rows(int(row['canon']), n)
                r = decode_three(t)
                need(r is not None, 'source census tc3 decoder')
                w = audit_three(t, r)
                iv = invariants(t)
                need(iv[2:] == (3, 3 if w == w[::-1] else 0), 'independent tc3/blue replay')
                hits += 1; fixed += w == w[::-1]
        need((count, hits, fixed) == (expected[n], f[n], sc[n]), 'complete inherited census and new formula')
        census.append({'n': n, 'classes': count, 'tc3': hits, 'self_converse_tc3': fixed})

    port_samples = []
    for n in range(5, 16):
        ws = [w for w in atom_words(n) if w.count(5) == 1]
        need(len(ws) == f[n], 'ordered atom enumeration versus rational GF')
        need(sum(w == w[::-1] for w in ws) == sc[n], 'SC center enumeration')
        need(len({min(w, w[::-1]) for w in ws}) == (f[n]+sc[n])//2, 'converse quotient')
        for w in ws:
            t = stack(w); r = decode_three(t)
            need(audit_three(t, r) == w, 'source reconstruction')
            ports = port_receipts(t, r)
            need(len(ports) == 3, 'complete three presentations')
            for p in ports:
                need(audit_port(t, p), 'exact target/port consumer')
                need(p['blue'] == (w == w[::-1]), 'SC iff pure on tc3 stratum')
            if w in ((5,), (1, 5), (1, 5, 1), (3, 5, 1)):
                port_samples.append({'word': w, 'receipt': r, 'presentations': ports})
    for n in range(6, len(f)):
        at = lambda k: f[k] if k >= 0 else 0
        need(f[n] == 2*at(n-1)-at(n-2)+2*at(n-3)-2*at(n-4)-at(n-6), 'six-lag recurrence')

    # Hostiles to source/target compression and to dropping hypotheses.
    t5 = cyclic(5); r = decode_three(t5)
    p = port_receipts(t5, r)[0]
    expect_reject(lambda: audit_three(tuple(float(x) for x in t5), r))
    expect_reject(lambda: audit_three(t5, {**r, 'frames': ((False, 1, 2, 3, 4),)}))
    expect_reject(lambda: audit_port(t5, {**p, 'ports': (p['ports'][1], p['ports'][0])}))
    expect_reject(lambda: audit_port(t5, {**p, 'target': t5}))
    changed = relabel(t5, (1, 0, 2, 3, 4))
    expect_reject(lambda: audit_three(changed, r))
    expect_reject(lambda: four_patterns(cyclic(3)))
    expect_reject(lambda: audit_four(sample['source'], {**sample, 'words':
        (tuple(float(x) for x in sample['words'][0]),) + sample['words'][1:]}))
    expect_reject(lambda: audit_four(sample['source'], {**sample, 'cycles':
        (sample['cycles'][0] * 2, sample['cycles'][1])}))
    need(sccs((0, 1)) == ((1,), (0,)), 'ordered SCC key cannot inspect its own in-place sort')
    left, right = stack((1, 5)), stack((5, 1))
    need(left != right and decode_three(left) is not None and decode_three(right) is not None,
         'unordered atoms lose source order')
    # Two interlaced C3 orbits give a genuine frame, not a three-presentation class.
    small = invariants(tuple(sample['source']))
    need(small[2] >= 4, 'two-orbit lower bound independently counted')
    s4 = (2, 4, 9, 3)
    need(invariants(s4)[2:] == (5, 3), 'deficit-two exception is not tc3')
    need(invariants(t5) == (15, 5, 3, 3), 'regular C5 positive atom')
    strange5 = tuple(sum(1 << j for i0, j in
        ((0,1),(1,2),(2,3),(3,0),(0,4),(2,4),(4,1),(4,3),(0,2),(1,3)) if i0 == i)
        for i in range(5))
    inv = invariants(strange5)
    need(inv == (15, 3, 5, 1) and all(arc(strange5, p[-1], p[0]) for p in hamilton_paths(strange5)),
         'small-order bypass exception differs from regular C5')

    return {'status': 'FINITE-EXACT controls; all-order proof in README', 'checks': GATES,
            'labelled_small': labelled, 'two_orbit_frames': frames, 'census': census,
            'coefficient_rows': [{'n': n, 'ordinary_tc3': f[n], 'self_converse_tc3': sc[n],
                                  'merged_tc3': (f[n]+sc[n])//2,
                                  'nonself_trivalent': (f[n]-sc[n])//2} for n in range(31)],
            'four_pattern_example': sample, 'port_samples': port_samples,
            'hostiles': {'S4': [5, 3], 'nonregular_order5': list(inv),
                         'two_orbit_n6': list(small)},
            'dependencies': ['audited all-order pure-blue classification',
                             'audited strong deficit-two classification',
                             'CITED Benhocine-Wojda 1983 n>=7 bypass theorem',
                             'unique SCC factorization and inherited merged degree law']}


def serialize(x):
    return json.dumps(x, sort_keys=True, indent=2) + '\n'


def main():
    p = argparse.ArgumentParser(); p.add_argument('--write', action='store_true'); p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    text = serialize(verify())
    target = HERE / 'receipt.json'
    if args.write: target.write_text(text, encoding='utf-8', newline='\n')
    if args.verify: need(target.read_text(encoding='utf-8') == text, 'saved receipt equality')
    print(text if not (args.write or args.verify) else
          f'PASS {json.loads(text)["checks"]} checks; LF SHA256 {hashlib.sha256(text.encode()).hexdigest()}')


if __name__ == '__main__':
    main()
