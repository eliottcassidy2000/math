"""Quadratic order bound and literal deficit-six center/port controls.

The all-order argument is in README.md. The inherited class universes are
filtered, not independently re-enumerated; each of their three surviving
sources is independently replayed by literal permutations.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools as it
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ValueError(why)


def exact(x, lo=0):
    if type(x) is not int or x < lo:
        raise ValueError('exact integer domain')
    return x


def source(t):
    if type(t) is not tuple:
        raise ValueError('immutable source rows')
    n = len(t)
    for i, r in enumerate(t):
        exact(r)
        if r >= 1 << n or r >> i & 1:
            raise ValueError('row range or diagonal')
        for j in range(i):
            if (r >> j & 1) + (t[j] >> i & 1) != 1:
                raise ValueError('tournament incidence')
    return t


def arc(t, a, b):
    return bool(t[a] >> b & 1)


def word(p, universe):
    if type(p) is not tuple or any(type(x) is not int for x in p) or sorted(p) != sorted(universe):
        raise ValueError('literal permutation of declared vertices')
    return p


def path(t, p, universe=None):
    word(p, range(len(t)) if universe is None else universe)
    if not all(arc(t, a, b) for a, b in zip(p, p[1:])):
        raise ValueError('non-native Hamilton path')


def local(t, p):
    return tuple(sum(int(arc(t, a, b)) << j for j, b in enumerate(p)) for a in p)


def anti(t, a):
    word(a, range(len(t)))
    return all(arc(t, a[x], a[y]) != arc(t, x, y) for x, y in it.combinations(range(len(t)), 2))


def blue(t, p):
    path(t, p)
    r = dict(zip(p, p[::-1]))
    return anti(t, tuple(r[i] for i in range(len(t))))


def target(t, p):
    path(t, p)
    kept = {frozenset(e) for e in zip(p, p[1:])}
    rows = [0] * len(t)
    for a, b in it.combinations(range(len(t)), 2):
        i, j = (a, b) if arc(t, a, b) ^ (frozenset((a, b)) not in kept) else (b, a)
        rows[i] |= 1 << j
    return tuple(rows)


def strong(t):
    source(t)
    n = len(t)
    reach = [r | 1 << i for i, r in enumerate(t)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(r == (1 << n) - 1 for r in reach)


def order(p):
    word(p, range(len(p)))
    remaining = set(range(len(p)))
    o = 1
    while remaining:
        a = min(remaining); k = 0
        while a in remaining:
            remaining.remove(a); a = p[a]; k += 1
        o = math.lcm(o, k)
    return o


def cutoff(d):
    exact(d)
    if d % 2:
        raise ValueError('self-converse defect budget is even')
    r = 1 << (d + 3).bit_length()
    return {'d': d, 'r': r, 'q_max': d + 2,
            'q1_m_max': 3*r-1, 'q_ge2_m_max': 2*r-1,
            'order_max': max(3*r-1, (d+2)*(2*r-1)),
            'forced_defect': r-2}


def residual_orders(d, after=0, nonautomorphic_cycle=False):
    exact(after)
    if type(nonautomorphic_cycle) is not bool:
        raise ValueError('literal choice of inherited cycle theorem')
    if nonautomorphic_cycle and after < 6:
        raise ValueError('closing-prefix theorem requires n>=7')
    c = cutoff(d)
    out = set() if nonautomorphic_cycle else set(range(1, c['q1_m_max']+1, 2))
    for q in range(2, c['q_max']+1):
        out.update(q*m for m in range(1, c['q_ge2_m_max']+1, 2))
    return sorted(n for n in out if n > after and n >= 3)


def digest(t):
    return hashlib.sha256(json.dumps(t, separators=(',', ':')).encode()).hexdigest()


def audit_prefix(t, rec):
    """Read only: authenticates supplied words, never searches for a suffix."""
    source(t)
    if type(rec) is not dict or set(rec) != {'source', 'digest', 'anti', 'cycle', 'suffix', 'exits', 'words', 'targets'}:
        raise ValueError('receipt schema')
    source(rec['source'])
    if rec['source'] != t or rec['digest'] != digest(t):
        raise ValueError('source substitution')
    a, s, w, exits = rec['anti'], rec['cycle'], rec['suffix'], rec['exits']
    if not anti(t, a):
        raise ValueError('source is not authenticated self-converse')
    if type(s) is not tuple or len(s) < 4 or len(s) & (len(s)-1):
        raise ValueError('power-of-two prefix cycle')
    word(s, set(s))
    if any(x < 0 or x >= len(t) for x in s):
        raise ValueError('prefix source labels')
    word(w, set(range(len(t))) - set(s))
    if len(t) < 2*len(s) or not w:
        raise ValueError('half-support or suffix boundary')
    path(t, s, set(s)); path(t, w, set(w))
    if not arc(t, s[-1], s[0]):
        raise ValueError('prefix does not close')
    if type(exits) is not tuple or any(type(x) is not int for x in exits) or len(set(exits)) != len(exits):
        raise ValueError('distinct exact exits')
    if not exits or any(x not in s or not arc(t, x, w[0]) for x in exits):
        raise ValueError('actual source joins required')
    if type(rec['words']) is not tuple or type(rec['targets']) is not tuple or len(rec['words']) != len(exits) or len(rec['targets']) != len(exits):
        raise ValueError('word/target multiplicity')
    for x, p, dest in zip(exits, rec['words'], rec['targets']):
        j = s.index(x)+1
        if p != s[j:]+s[:j]+w:
            raise ValueError('literal cut or suffix mismatch')
        path(t, p); source(dest)
        if target(t, p) != dest:
            raise ValueError('wrong complement target')
    for p, z in it.combinations(rec['words'], 2):
        perm = dict(zip(p, z)); v = tuple(perm[i] for i in range(len(t)))
        support = {i for i in range(len(t)) if v[i] != i}
        reflection = dict(zip(p, reversed(p)))
        need(order(v) % 2 == 0, 'even relative order')
        need(not support.intersection(reflection[x] for x in support), 'disjoint reflection support')
    need(sum(blue(t, p) for p in rec['words']) <= 1, 'actual blue control')
    return 2*(len(exits)-1)


def carousel(n):
    exact(n, 3)
    if not n % 2:
        raise ValueError('odd carousel')
    return tuple(sum(1 << ((i+j) % n) for j in range(1, (n+1)//2)) for i in range(n))


def two_shores(m):
    """Interleaved carousel shores, a_i -> b_j iff i=j."""
    c = carousel(m); n = 2*m; t = [0]*n
    for a, b in it.combinations(range(n), 2):
        if a % 2 == b % 2:
            forward = arc(c, a//2, b//2)
        elif a % 2 == 0:
            forward = a//2 == b//2
        else:
            forward = a//2 != b//2
        x, y = (a, b) if forward else (b, a)
        t[x] |= 1 << y
    return source(tuple(t))


def prefix_example(m, q):
    t = carousel(m) if q == 1 else two_shores(m)
    n = len(t); r = 16
    s = tuple(q*(j*m//r) for j in range(r))
    outside_shore = [q*i for i in range(m) if q*i not in s]
    v = max(outside_shore, key=lambda x: (sum(arc(t, a, x) for a in s), -x))
    rem = tuple(x for x in range(n) if x not in s)
    j = rem.index(v); w = rem[j:]+rem[:j]
    exits = tuple(a for a in s if arc(t, a, v))
    words = tuple(s[(s.index(a)+1):]+s[:(s.index(a)+1)]+w for a in exits)
    a = tuple((-i) % n for i in range(n)) if q == 1 else tuple(n-1-i for i in range(n))
    rec = {'source': t, 'digest': digest(t), 'anti': a, 'cycle': s, 'suffix': w,
           'exits': exits, 'words': words, 'targets': tuple(target(t, p) for p in words)}
    need(strong(t), 'example strong source')
    need(len(exits) >= r//2, 'regular-shore average attained')
    need(audit_prefix(t, rec) >= 14, 'source-bound prefix obstruction')
    return rec


def hamilton(t):
    for p in it.permutations(range(len(t))):
        if all(arc(t, a, b) for a, b in zip(p, p[1:])):
            yield p


def center_record(t):
    """Complete literal path orbit and black-target receipt for one small source."""
    source(t); need(strong(t), 'center strong')
    n = len(t); perms = tuple(it.permutations(range(n)))
    aut = tuple(p for p in perms if local(t, p) == t)
    aa = next(p for p in perms if anti(t, p))
    classes = {}
    for p in hamilton(t):
        classes.setdefault(local(t, p), []).append(p)
    records = []
    keys = sorted(classes)
    for k in keys:
        ps = classes[k]; p = min(ps)
        need(len(ps) == len(aut), 'free literal orbit mass')
        bp = blue(t, p)
        need(all(blue(t, z) == bp for z in ps), 'orbit blue invariance')
        mate = tuple(aa[x] for x in p[::-1]); path(t, mate)
        mi = keys.index(local(t, mate))
        need((mi == keys.index(k)) == bp, 'converse mate fixed iff blue')
        records.append({'word': p, 'blue': bp, 'mate_orbit': mi,
                        'target': target(t, p), 'path_ports': (p[0], p[-1])})
    b = sum(r['blue'] for r in records)
    need(len(records)-b == 6, 'literal defect six')
    for i, rec in enumerate(records):
        need(records[rec['mate_orbit']]['mate_orbit'] == i, 'mate involution')
    return {'source': t, 'source_digest': digest(t), 'anti': aa, 'automorphisms': aut,
            'H': sum(map(len, classes.values())), 'Aut': len(aut), 'h': len(records),
            'b': b, 'black_pairs': [(i, r['mate_orbit']) for i, r in enumerate(records)
                                   if not r['blue'] and i < r['mate_orbit']], 'orbits': records}


def normalized_two_shore(m, mask):
    exact(m, 3); exact(mask)
    if m % 2 == 0 or mask >= 1 << (2*m-3):
        raise ValueError('odd-shore mask domain')
    sa = {j if mask >> (j-1) & 1 else m-j for j in range(1,(m+1)//2)}
    sb = {j if mask >> ((m-1)//2+j-1) & 1 else m-j for j in range(1,(m+1)//2)}
    f = {0} | {j for j in range(1,m-1) if mask >> (m+j-2) & 1}
    rows = [0]*(2*m)
    for a,b in it.combinations(range(2*m),2):
        if a % 2 == b % 2:
            yes = (b//2-a//2) % m in (sa if a % 2 == 0 else sb)
        elif a % 2 == 0:
            yes = (b//2-a//2) % m in f
        else:
            yes = (a//2-b//2) % m not in f
        x,y = (a,b) if yes else (b,a); rows[x] |= 1 << y
    return source(tuple(rows))


def subset_H(t):
    """Exact Hamilton words, preserving their final vertex in each subset."""
    source(t); n = len(t)
    dp = [[0]*n for _ in range(1<<n)]
    for a in range(n):
        dp[1<<a][a] = 1
    for mask in range(1,1<<n):
        for a,x in enumerate(dp[mask]):
            if x:
                rem = t[a] & ~mask
                while rem:
                    bit = rem & -rem; b = bit.bit_length()-1
                    dp[mask|bit][b] += x; rem -= bit
    return sum(dp[-1])


def frame_counts(t):
    """All score-compatible maps, not merely rotation centralizers."""
    source(t); n = len(t); groups = {}
    for i,r in enumerate(t):
        groups.setdefault(r.bit_count(),[]).append(i)
    aut = []; antis = []
    for flip in (False,True):
        keys = list(groups)
        targets = [groups.get(n-1-k if flip else k,[]) for k in keys]
        if any(len(groups[k]) != len(g) for k,g in zip(keys,targets)):
            continue
        for ims in it.product(*(it.permutations(g) for g in targets)):
            p = [0]*n
            for k,im in zip(keys,ims):
                for a,b in zip(groups[k],im):
                    p[a] = b
            if all(arc(t,p[a],p[b]) == (arc(t,a,b)^flip) for a,b in it.combinations(range(n),2)):
                (antis if flip else aut).append(tuple(p))
    involutions = [a for a in antis if all(a[a[x]] == x for x in range(n))]
    need(bool(antis) == bool(involutions), 'literal finite anti-involution existence')
    bh = 0; blue_port = None
    for a in involutions:
        pairs = [(x,a[x]) for x in range(n) if x < a[x]]
        need(len(pairs)*2 == n, 'even source anti has no fixed vertex')
        for pp in it.permutations(pairs):
            for signs in it.product(range(2), repeat=n//2):
                half = tuple(p[s] for p,s in zip(pp,signs))
                if all(arc(t,x,y) for x,y in zip(half,half[1:])) and arc(t,half[-1],a[half[-1]]):
                    bh += 1
                    if blue_port is None:
                        blue_port = half+tuple(a[x] for x in half[::-1])
    if blue_port is not None:
        need(blue(t,blue_port), 'actual blue sidecar')
    H = subset_H(t)
    need(H % len(aut) == 0 and bh % len(aut) == 0, 'full free orbit divisions')
    black_port = next((p for p in hamilton(t) if not blue(t,p)),None)
    need(black_port is not None, 'actual black sidecar')
    return {'H': H, 'Aut': len(aut), 'h': H//len(aut), 'blue_H': bh, 'b': bh//len(aut),
            'self_converse': bool(antis), 'automorphisms': aut, 'anti_count': len(antis),
            'anti_involutions': involutions, 'blue_word': blue_port,
            'black_word': black_port, 'black_target': target(t,black_port)}


def order_ten_census():
    records = []; controls = []; histogram = Counter()
    for m in (3,5):
        sources = set()
        for mask in range(1 << (2*m-3)):
            t = normalized_two_shore(m,mask); n = len(t); sources.add(t)
            p = tuple(range(n)); path(t,p)
            need(arc(t,n-1,0), 'literal normalizing Hamilton cycle')
            need(local(t,tuple((i+2)%n for i in range(n))) == t, 'rotation by two is native')
            c = frame_counts(t)
            if m == 3:
                # Independent literal permutation path count and blue count.
                hp = tuple(hamilton(t))
                need(c['H'] == len(hp) and c['blue_H'] == sum(blue(t,z) for z in hp), 'half-word/DP independent control')
                controls.append({'mask': mask, 'source': t, 'H': c['H'], 'Aut': c['Aut'], 'h': c['h'], 'b': c['b']})
            else:
                if c['self_converse']:
                    histogram[c['h']-c['b']] += 1
                need(not c['self_converse'] or c['h']-c['b'] != 6, 'complete order10 center exclusion')
                records.append({'mask': mask, 'source': t, 'digest': digest(t), 'hamilton_cycle': p, **c})
        need(len(sources) == 1 << (2*m-3), 'literal normalized mask injectivity')
    need(sum(x['h']-x['b'] == 6 for x in controls) == 4, 'known six-center occurs in four normalized frames')
    need(histogram == {1528:8,1628:8,2908:12,2926:12,3064:12,3114:12}, 'complete SC10 deficit histogram')
    return {'frame_count': len(records), 'self_converse_frame_count': sum(histogram.values()),
            'self_converse_deficit_histogram': dict(sorted(histogram.items())),
            'six_vertex_positive_controls': controls, 'frames': records,
            'scope': 'complete necessary normalized universe, not 128 isomorphism classes'}


def census_rows(code, n):
    out = [0]*n; width = n*(n-1)//2
    for k, (a, b) in enumerate(it.combinations(range(n), 2)):
        i, j = (a, b) if code >> (width-1-k) & 1 else (b, a)
        out[i] |= 1 << j
    return source(tuple(out))


def inherited_census():
    expected_classes = [2, 4, 12, 56, 456, 6880, 191536]
    expected_sc = [2, 2, 8, 12, 88, 176, 2752]
    answer = []
    for n in range(3, 10):
        p = (ROOT / f'04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n{n}.tsv') if n < 9 else HERE.parent / 'blue_quotient_20261020b/n9_classes.tsv.gz'
        op = gzip.open if n == 9 else open
        with op(p, 'rt', encoding='utf-8', newline='') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        sc = hits = 0; survivors = []
        for row in rows:
            is_sc = int(row['SC']) if n < 9 else int(row['code']) == int(row['converse_code'])
            if not is_sc:
                continue
            sc += 1
            if n == 9:
                need(int(row['SC_H']) % int(row['Aut']) == 0, 'inherited SC path-orbit division')
            h = int(row['tc']) if n < 9 else int(row['SC_H']) // int(row['Aut'])
            b = int(row['bm']) if n < 9 else int(row['blue_mass'])
            if h-b != 6:
                continue
            code = int(row['canon']) if n < 9 else int(row['code'])
            t = census_rows(code, n)
            if strong(t):
                hits += 1; survivors.append({'code': code, 'source': t})
        need(len(rows) == expected_classes[n-3] and sc == expected_sc[n-3], 'inherited universe counts')
        need(hits == int(n in (5, 6, 7)), 'inherited R6 initial coefficients')
        answer.append({'n': n, 'class_count': len(rows), 'self_converse_count': sc,
                       'strong_delta6': survivors, 'path': p.relative_to(ROOT).as_posix(),
                       'sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
                       'scope': 'inherited complete universe; filters replayed here'})
    return answer


def rejected(fn, *args):
    try:
        fn(*args)
    except (ValueError, TypeError):
        need(True, 'rejected hostile')
    else:
        need(False, 'accepted hostile')


def plus(a, b):
    c = [0]*max(len(a), len(b))
    for i, x in enumerate(a):
        c[i] += x
    for i, x in enumerate(b):
        c[i] += x
    while len(c) > 1 and c[-1] == 0:
        c.pop()
    return c


def times(a, b):
    c = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return c


def poly_gcd(a,b):
    a,b = list(map(Fraction,a)),list(map(Fraction,b))
    def trim(p):
        while p and not p[-1]:
            p.pop()
        return p
    trim(a); trim(b)
    while b:
        while len(a) >= len(b):
            k = len(a)-len(b); c = a[-1]/b[-1]
            for j,x in enumerate(b):
                a[k+j] -= c*x
            trim(a)
        a,b = b,a
    return [x/a[-1] for x in a] if a else []


def recurrence_control():
    # Inherited trivalent SCC formula. The top coefficient is synthetic;
    # its only use is authenticating the new safe recurrence onset.
    p = [1,-1,0,-1]; u = [1,0,-1,0,0,0,-1]
    pp, uu = times(p,p), times(u,u)
    delta = plus(uu, [-x for x in times(pp,u)])
    need(all(x % 2 == 0 for x in delta), 'integral half-difference numerator')
    r = [0]*249
    for k in (5,6,7,248):
        r[k] = 1
    num = plus([0]*5+[x//2 for x in delta], times(times(pp,u),r))
    num = plus(num, [0]*10+times(pp,[1,1,0,1]))
    den = times(pp,uu)
    need(len(den)-1 == 18 and len(num)-1 == 260 and num[-1] == -1, 'sharp formal degree bill')
    need(poly_gcd(p,u) == [1], 'coprime inherited factors')
    need(poly_gcd(p,[i*p[i] for i in range(1,len(p))]) == [1], 'squarefree p')
    need(poly_gcd(u,[i*u[i] for i in range(1,len(u))]) == [1], 'squarefree p(z2)')
    need(poly_gcd(u,[1,1,0,1]) == [1], 'no last-summand cancellation')
    need(poly_gcd(num,den) == [1], 'independent synthetic reduced denominator')
    seq = []
    for n in range(281):
        seq.append((num[n] if n < len(num) else 0)-sum(den[j]*seq[n-j] for j in range(1,min(n,len(den)-1)+1)))
    need(sum(den[j]*seq[260-j] for j in range(len(den))) == -1, 'synthetic pre-onset hostile')
    for n in range(261,281):
        need(sum(den[j]*seq[n-j] for j in range(len(den))) == 0, 'new safe recurrence onset')
    return {'denominator_coefficients': den, 'safe_onset': 261,
            'minimal_recurrence_order': 18,
            'synthetic_R6_degrees': [5,6,7,248], 'synthetic_numerator_degree': 260,
            'scope': 'formal identity only; synthetic top coefficient is not a classified center'}


def verify():
    global GATES
    GATES = 0
    for d in range(0, 100, 2):
        c = cutoff(d)
        need(c['r'] >= d+4 and c['r']//2 < d+4, 'least power two')
        need(c['forced_defect'] > d, 'strict budget contradiction')
        need(c['order_max'] < 4*(d+2)*(d+4), 'quadratic asymptotic bound')
    need(cutoff(6)['order_max'] == 248, 'R6 new cutoff')
    need(8-2 == 6, 'eight-cycle average equality is not a contradiction')
    for d in (True, 6.0, -2, 3):
        rejected(cutoff, d)
    examples = [prefix_example(49, 1), prefix_example(33, 2)]
    e = examples[0]; t = e['source']
    bad = dict(e); bad['cycle'] = e['cycle']*2; rejected(audit_prefix, t, bad)
    bad = dict(e); bad['exits'] = e['exits']+(e['exits'][0],); rejected(audit_prefix, t, bad)
    bad = dict(e); bad['anti'] = tuple(float(x) for x in e['anti']); rejected(audit_prefix, t, bad)
    bad = dict(e); bad['source'] = tuple(float(x) for x in t); rejected(audit_prefix, t, bad)
    bad = dict(e); bad['targets'] = (t,)+e['targets'][1:]; rejected(audit_prefix, t, bad)
    bad = dict(e); bad['suffix'] = e['suffix'][::-1]; rejected(audit_prefix, t, bad)
    # Independent regular-orbit averaging on every subset of small carousels.
    average_controls = 0
    for n in (3, 5, 7, 9, 11):
        t = carousel(n)
        for r in range(2, n, 2):
            for s in it.combinations(range(n), r):
                outside = set(range(n))-set(s)
                counts = [sum(arc(t, a, b) for a in s) for b in outside]
                need(2*sum(counts) == r*(n-r), 'exact regular subset exit sum')
                need(max(counts) >= r//2, 'attained average')
                average_controls += 1
    # Half-support cannot be removed: carousel5 is pure blue, yet has
    # two actual paths with an even-order relative permutation.
    c5 = carousel(5); paths5 = tuple(hamilton(c5)); hostile = None
    for p, z in it.combinations(paths5, 2):
        v = dict(zip(p, z)); v = tuple(v[i] for i in range(5))
        if order(v) % 2 == 0 and blue(c5, p) and blue(c5, z):
            hostile = {'source': c5, 'paths': (p, z), 'relative_order': order(v)}
            break
    need(hostile is not None, 'overlapping-support pure-blue hostile')
    census = inherited_census()
    centers = [center_record(tuple(s['source'])) for row in census for s in row['strong_delta6']]
    need([(x['H'], x['Aut'], x['h'], x['b']) for x in centers] == [(9,1,9,3),(45,3,15,9),(189,21,9,3)], 'three exact centers')
    # C3[K2] identifies the six-vertex source through actual homogeneous pairs.
    t6 = centers[1]['source']; pairs = ((0,4),(1,3),(2,5))
    for a, b in pairs:
        need(all(arc(t6,a,x) == arc(t6,b,x) for x in range(6) if x not in (a,b)), 'literal twin block')
    need(strong(local(t6, tuple(p[0] for p in pairs))), 'cyclic three-block quotient')
    # Paley7 identification retains an explicit source isomorphism.
    p7 = tuple(sum(1 << ((i+j)%7) for j in (1,2,4)) for i in range(7))
    paley_map = next(p for p in it.permutations(range(7)) if local(centers[2]['source'],p) == p7)
    need(local(centers[2]['source'],paley_map) == p7, 'literal Paley7 map')
    coarse = residual_orders(6, 9)
    pre_ten = residual_orders(6, 9, True)
    remaining = [n for n in pre_ten if n != 10]
    need(len(coarse) == 103 and len(pre_ten) == 92 and len(remaining) == 91 and min(remaining) == 12 and max(remaining) == 248, 'declared OPEN order support')
    need(all(n % 16 for n in remaining), 'cycle-cut odd orbit excludes16')
    ten = order_ten_census()
    rec = recurrence_control()
    return {'status': 'PROVED quadratic cutoff; FINITE-EXACT small source replays; remaining census OPEN',
            'checks': GATES, 'bound': cutoff(6), 'remaining_orders_after9': remaining,
            'remaining_order_count': len(remaining), 'prefix_certificates': examples,
            'regular_average_controls': average_controls, 'half_support_hostile': hostile,
            'inherited_census': census, 'centers': centers,
            'order_ten_complete_normalized_census': ten,
            'six_center_pairs': pairs, 'seven_center_paley_map': paley_map,
            'coarse_remaining_order_count': len(coarse),
            'trivalent_recurrence': rec,
            'trivalent_recurrence_scope': 'inherited SCC formula, arbitrary R6 degree<=248; no unknown coefficient supplied'}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--write', action='store_true'); ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = verify(); data = json.dumps(result, indent=2, sort_keys=True)+'\n'
    receipt = HERE/'receipt.json'
    if args.write:
        receipt.write_text(data, encoding='utf-8', newline='\n')
    if args.verify:
        need(receipt.read_text(encoding='utf-8') == data, 'saved receipt equality')
    print(json.dumps({'checks': result['checks'], 'sha256': hashlib.sha256(data.encode()).hexdigest(),
                      'order_cutoff': result['bound']['order_max'], 'open_orders': len(result['remaining_orders_after9'])}, sort_keys=True))


if __name__ == '__main__':
    main()
