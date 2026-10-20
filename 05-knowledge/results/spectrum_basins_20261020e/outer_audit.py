"""Independent local-ledger and endpoint-recursion audit; no producer imports.

Exact universe: every standard-HP presentation of orders 4,5,6;
all strong members, all applicable triangle/ear sites, and all endpoint-
rigid recursive choices. Four fixed-blue centers test empty and nonempty
middle blocks. This is a control of the separately written universal proof.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def arc(t, u, v):
    return bool(t[u] & (1 << v))


def induced(t, p):
    return tuple(sum(1 << j for j, v in enumerate(p) if arc(t, u, v))
                 for u in p)


def presentations(n):
    free = [(i, j) for i, j in combinations(range(n), 2) if j != i+1]
    for mask in range(1 << len(free)):
        t = [1 << (i+1) if i+1 < n else 0 for i in range(n)]
        for k, (u, v) in enumerate(free):
            x, y = (u, v) if mask & (1 << k) else (v, u)
            t[x] |= 1 << y
        yield tuple(t)


def scc(t):
    n = len(t)
    reach = [t[v] | (1 << v) for v in range(n)]
    for v in range(n):
        for u in range(n):
            if reach[u] & (1 << v):
                reach[u] |= reach[v]
    remaining = set(range(n))
    parts = []
    while remaining:
        v = min(remaining)
        part = {u for u in remaining if reach[v] & (1 << u)
                and reach[u] & (1 << v)}
        remaining -= part
        parts.append(tuple(sorted(part)))
    return sorted(parts, key=lambda p: -sum(arc(t, p[0], q[0])
                                          for q in parts if q != p))


def paths_to(t, z):
    def walk(word, unused):
        if not unused:
            yield word
        for u in sorted(unused):
            if arc(t, u, word[0]):
                yield from walk((u,)+word, unused-{u})
    yield from walk((z,), set(range(len(t)))-{z})


def blue(t, p):
    return all(arc(t, p[i], p[j]) == arc(t, p[-1-j], p[-1-i])
               for i, j in combinations(range(len(t)), 2))


def k_move(t, p):
    need(sorted(p) == list(range(len(t))), 'not a vertex permutation')
    need(all(arc(t, u, v) for u, v in zip(p, p[1:])), 'not a directed HP')
    keep = {frozenset((u, v)) for u, v in zip(p, p[1:])}
    out = list(t)
    for u, v in combinations(range(len(t)), 2):
        if frozenset((u, v)) not in keep:
            out[u] ^= 1 << v
            out[v] ^= 1 << u
    return tuple(out)


def padded(s, center):
    r, c = len(s), len(center)
    left = tuple(sum(1 << v for v in range(r)
                     if arc(s, r-1-v, r-1-u)) for u in range(r))
    blocks = (left, center, s)
    out = []
    offset = 0
    n = 2*r+c
    for block in blocks:
        after = sum(1 << v for v in range(offset+len(block), n))
        out.extend((row << offset) | after for row in block)
        offset += len(block)
    t = tuple(out)
    need(blue(t, tuple(range(n))), 'padding not fixed-blue')
    return t, tuple(range(r+c)), tuple(range(r+c, n))


def triangle_ledger(a, q, i):
    first, left = q[0], q[1:i]
    b, c, d = q[i:i+3]
    right = q[i+3:]
    return (a+left+(b,c,d)+right+(first,),
            (c,)+a+left+(b,d)+right+(first,),
            (c,)+a+(first,)+left+(b,d)+right)


def ear_ledger(f, d):
    a, b, z = d[0], d[1], d[-1]
    middle = d[2:-1]
    return (f+(b,)+middle+(z,a), middle+f+(b,z,a),
            middle+f+(a,b,z))


def near_ledger(f, d):
    a, c, z = d[0], d[-2], d[-1]
    v = d[1:-2]
    need(bool(v), 'empty near-tail middle')
    return (f+v+(c,z,a), (z,)+f+v+(c,a), (a,z)+f+v+(c,),
            (a,c,z)+f+v, (c,z)+f+(a,)+v)


def check_ledger(t, ledger, expected_excess):
    target = k_move(t, tuple(range(len(t))))
    current = t
    excess = []
    for p in ledger:
        need(not blue(current, p), 'blue step in proposed black ledger')
        excess.append(current[p[0]].bit_count()+current[p[-1]].bit_count()
                      -(len(t)-1))
        current = k_move(current, p)
    need(current == target, 'wrong labelled final target')
    need(excess[1:] == expected_excess, ('score identity', excess,
                                        expected_excess))


def recursive_choice(s):
    """Return a local ledger site, checking each claimed recursion invariant."""
    n = len(s)
    for i in range(1, n-2):
        if arc(s, i+2, i):
            return 'triangle', i
    d = tuple(range(n))
    depth = 0
    while len(d) > 1:
        a, b, z = d[0], d[1], d[-1]
        need(arc(s,z,a), 'endpoint path not closed')
        if arc(s,z,b):
            need(len(d) >= 4, 'missed consecutive C3')
            return 'ear', (d, depth)
        need(all(arc(s,a,u) for u in d[1:-1]), 'prefix lemma fails')
        remainder = d[1:]
        smaller = induced(s, remainder)
        parts = scc(smaller)
        need(len(parts) >= 2, 'deletion remains strong')
        for part in parts[:-1]:
            need(len(part) in (1,3), 'nonrigid preceding SCC')
            need(len(part) == 1, 'missed preceding consecutive C3')
        next_d = tuple(remainder[i] for i in parts[-1])
        dt = induced(s,next_d)
        need(all(induced(dt,p) == dt for p in paths_to(dt,len(dt)-1)),
             'endpoint rigidity not inherited')
        depth += 1
        if len(next_d) == 1:
            need(len(d) >= 4, 'missed terminal consecutive C3')
            return 'near', (d,depth)
        d = next_d
    raise RuntimeError('unclassified endpoint-rigid presentation')


def main():
    centers = ((), (0,), (2,4,1), (6,12,24,17,3))
    counts = Counter()
    by_order = {}
    for r in range(4,7):
        local = Counter()
        for s in presentations(r):
            local['presentations'] += 1
            if len(scc(s)) != 1:
                continue
            local['strong'] += 1
            different = next((p for p in paths_to(s,r-1)
                              if induced(s,p) != s), None)
            rigid = different is None
            if rigid:
                local['endpoint_rigid'] += 1
                choice = recursive_choice(s)
                local['recursive_'+choice[0]] += 1
            for center in centers:
                t, a, q = padded(s,center)
                n = len(t)
                scores = [row.bit_count() for row in s]
                if not rigid:
                    changed = a+tuple(q[u] for u in different)
                    original = tuple(range(n))
                    new_rotation = changed[1:]+changed[:1]
                    old_rotation = original[1:]+original[:1]
                    total = t[0].bit_count()+t[1].bit_count()
                    check_ledger(t,(changed,new_rotation,old_rotation),
                                 [n-total,total-n])
                    counts['flexible_common_hubs'] += 1
                rotation_nonauto = induced(s,tuple(range(1,r))+(0,)) != s
                if arc(s,r-1,0) and rotation_nonauto:
                    for i in range(1,r-2):
                        if arc(s,i+2,i):
                            ledger = triangle_ledger(a,q,i)
                            expected = [n-2-scores[i+1]-scores[0],
                                        scores[i+1]+scores[-1]+2-n]
                            check_ledger(t,ledger,expected)
                            counts['triangle_ledgers'] += 1
                if rigid:
                    kind, data = choice
                    if kind == 'triangle':
                        ledger = triangle_ledger(a,q,data)
                        expected = [n-2-scores[data+1]-scores[0],
                                    scores[data+1]+scores[-1]+2-n]
                    else:
                        d, depth = data
                        offset = d[0]
                        labelled_d = tuple(q[u] for u in d)
                        f = a+q[:offset]
                        if kind == 'ear':
                            ledger = ear_ledger(f,labelled_d)
                            c = d[2]
                            # These are the literal endpoint-score formulas;
                            # C may contain more than one vertex.
                            expected = [n-2-scores[c]-scores[d[0]],
                                        scores[c]+scores[-1]+2-n]
                        else:
                            ledger = near_ledger(f,labelled_d)
                            k = len(d)-3
                            vscore = induced(s,d[1:-2])[-1].bit_count()
                            expected = [n-k-scores[-1]-3, k+3-n,
                                        n-k-vscore-3, vscore+4-n]
                        counts['recursion_depth_'+str(depth)] += 1
                    check_ledger(t,ledger,expected)
                    counts['recursive_ledgers'] += 1
        by_order[r] = dict(sorted(local.items()))
    result = dict(status='PASS', universe='all standard-HP presentations n=4..6',
                  center_orders=[len(c) for c in centers], by_order=by_order,
                  checks=dict(sorted(counts.items())), producer_imports=False)
    payload = json.dumps(result,indent=2,sort_keys=True)+'\n'
    (HERE/'outer_audit.json').write_text(payload,encoding='utf-8',newline='\n')
    print(payload,end='')
    print('sha256',hashlib.sha256(payload.encode()).hexdigest())


if __name__ == '__main__':
    main()
