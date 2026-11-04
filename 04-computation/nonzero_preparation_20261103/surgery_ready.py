#!/usr/bin/env python3
"""Exact finite-bank inverse for actual one/two-ready funnel targets.

The bank is the 34 explicitly listed rooted-entry alternatives below.  The
search is complete over every finite tail word for this bank and the requested
finite target set, by its proved lower bound; it is NOT a complete ready-class
inverse.  No imported research counter or response receipt is used.
"""
import argparse
from array import array
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

GATES = 0


def require(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def converse(rows):
    full = (1 << len(rows)) - 1
    return tuple(full ^ (1 << i) ^ row for i, row in enumerate(rows))


def strong(rows):
    reach = [row | (1 << i) for i, row in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(x == (1 << len(rows)) - 1 for x in reach)


def tables(rows):
    n = len(rows)
    require(1 <= n <= 19, 'declared exact uint64 table domain, n! < 2^64')
    size = 1 << n
    dp = array('Q', [0]) * (n * size)
    H = array('Q', [0]) * size
    H[0] = 1
    for i in range(n):
        dp[(1 << i) * n + i] = 1
    for mask in range(1, size):
        vertices = mask
        while vertices:
            bit = vertices & -vertices
            vertices -= bit
            i = bit.bit_length() - 1
            value = dp[mask*n+i]
            H[mask] += value
            choices = rows[i] & (size - 1 - mask)
            while choices and value:
                bit2 = choices & -choices
                choices -= bit2
                j = bit2.bit_length() - 1
                dp[(mask | bit2)*n+j] += value
    return H, dp


def counts(rows, roots=()):
    H, E = tables(rows)
    n, full = len(rows), (1 << len(rows)) - 1
    return {'H': int(H[full]),
            'F': {r: sum(int(E[m*n+r])*int(H[full^m])
                         for m in range(1, full+1)) for r in roots},
            'ends': [int(E[full*n+r]) for r in range(n)]}


@lru_cache(None)
def core_data(rows, ports):
    H, E = tables(rows)
    _, S = tables(converse(rows))
    n, full = len(rows), (1 << len(rows)) - 1
    e = [sum(int(E[m*n+r]) for r in ports) for m in range(full+1)]
    s = [sum(int(S[m*n+r]) for r in ports) for m in range(full+1)]
    return {'h': int(H[full]),
            'p2': sum(int(H[m])*int(H[full^m]) for m in range(full+1)),
            'e': e[full], 'b': s[full],
            'F': sum(e[m]*int(H[full^m]) for m in range(full+1)),
            'G': sum(s[m]*int(H[full^m]) for m in range(full+1)),
            'C': sum(e[m]*s[full^m] for m in range(full+1))}


def tail_coefficients(word, m):
    """Literal inherited weighted word recurrence; 1 means tail -> root."""
    t = len(word)
    L = (m-1)**t + sum(m**i*(m-1)**(t-i-1)
                       for i,b in enumerate(word) if b)
    D = (m-1)**t + sum((m-1)**i*m**(t-i-1)
                       for i,b in enumerate(word) if not b)
    states = {(0,0): 1}
    for b in word:
        nxt = {}
        def add(k,v):
            nxt[k] = nxt.get(k,0) + v
        for (last,started), value in states.items():
            add((last,started), (m-1)*value)
            add((1 if b else 2, started), value)
            if started or not b:
                add((last,1), value)
        states = nxt
    return L, D, sum(v for (last,_),v in states.items() if last != 2)


def ready_count(data, word):
    L = tail_coefficients(word,2)[0]
    _, D, B = tail_coefficients(word,1)
    return (L*(data['G']+data['p2']) + B*data['h'] - D*data['e']
            - data['C'] - data['F'] - data['b'])


def funnel(rows, ports, word):
    n, t = len(rows), len(word)
    root = n+t
    tail = ((1 << t)-1) << n
    out = [row | tail | (0 if i in ports else 1 << root)
           for i,row in enumerate(rows)]
    out += [sum(1 << (n+j) for j in range(i+1,t))
            | ((1 << root) if b else 0) for i,b in enumerate(word)]
    out.append(sum(1 << i for i in ports)
               + sum(1 << (n+i) for i,b in enumerate(word) if not b))
    return tuple(out), root


def birth(rows, incoming):
    n = len(rows)
    return tuple(row | (((incoming >> i)&1) << n)
                 for i,row in enumerate(rows)) + (((1 << n)-1)^incoming,)


def ready_target(rows, ports, word):
    parent, root = funnel(rows,ports,word)
    # Newborn beaten only by root, followed by converse: auxiliary has
    # the retained root as its SOLE out-neighbour.
    target = converse(birth(parent,1 << root))
    return target, root, len(target)-1


BANK = (
    (1,(2,4,1),(0,)), (1,(2,4,1),(1,)), (1,(2,4,1),(2,)),
    (1,(2,4,9,3),(0,)), (1,(2,4,9,3),(1,)), (1,(2,4,9,3),(2,)),
    (1,(2,4,9,19,7),(0,)), (1,(2,4,9,19,7),(1,)), (1,(2,4,9,19,7),(2,)),
    (1,(2,4,25,19,3),(0,)), (1,(2,4,25,19,3),(1,)), (1,(2,4,25,19,3),(2,)),
    (1,(2,12,25,17,3),(0,)), (1,(2,12,25,17,3),(1,)), (1,(0,),(0,)),
    (2,(2,4,1),(0,1)), (2,(2,4,1),(0,2)), (2,(2,4,1),(1,2)),
    (2,(2,4,9,3),(0,1)), (2,(2,4,9,3),(0,2)), (2,(2,4,9,3),(1,2)),
    (2,(6,4,8,3),(0,)),
    (2,(2,4,9,19,7),(0,1)), (2,(2,4,9,19,7),(0,2)), (2,(2,4,9,19,7),(1,2)),
    (2,(2,4,25,19,3),(0,1)), (2,(2,4,25,19,3),(0,2)), (2,(2,4,25,19,3),(1,2)),
    (2,(2,12,25,17,3),(0,1)), (2,(6,20,24,19,1),(0,)),
    (2,(6,4,24,19,3),(0,)), (2,(10,4,9,18,7),(0,)),
    (2,(6,4,8,19,7),(0,)), (2,(2,12,9,17,7),(0,)))


def search(targets, bank=BANK):
    targets = set(targets)
    bound, tested, hits, cuts = max(targets), 0, [], []
    for gap,rows,ports in bank:
        d = core_data(rows,ports)
        require(d['b'] == gap and (len(rows)==1 or strong(rows)), 'actual bank entry count')
        C = d['G']+d['p2']
        loss = d['e']+d['C']+d['F']+d['b']
        n = 0
        while n == 0 or (1+2**(n-1))*C-loss <= bound:
            for prefix in product((0,1),repeat=max(n-1,0)):
                word = prefix+(1,) if n else ()
                if not word and len(ports)==len(rows):
                    continue
                tested += 1
                H = ready_count(d,word)
                if H in targets:
                    target,r,a = ready_target(rows,ports,word)
                    hits.append({'H':H,'gap':gap,'core':rows,'ports':ports,
                                 'word':word,'rows':target,'root':r,'aux':a,
                                 'order':len(target)})
            n += 1
        cuts.append(n)
    return {'targets':sorted(targets),'bank_entries':len(bank),'words_tested':tested,
            'first_excluded_lengths':cuts,'hits':hits}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    # All bank words of length <=3: a literal adjacency control of the
    # inherited transfer, including the allowed empty-tail boundary.
    controls = 0
    for gap,rows,ports in BANK:
        d = core_data(rows,ports)
        for n in range(4):
            for prefix in product((0,1),repeat=max(n-1,0)):
                word = prefix+(1,) if n else ()
                if not word and len(ports)==len(rows): continue
                target,r,a = ready_target(rows,ports,word)
                data = counts(target,(r,))
                dual = counts(converse(target),(a,))
                require(strong(target), 'strong literal ready target')
                require(data['H']==ready_count(d,word), 'word transfer against actual paths')
                require(dual['F'][a]-data['H']==gap, 'actual start-ready count')
                require(target[a]==1 << r, 'sole-out-neighbour joint sidecar')
                controls += 1
    result = search((1429,2941,4289,4291,12869,12871,12875,12877))
    result['scope'] = 'all finite tails for the 34 literal bank entries, not all ready families'
    result['literal_short_word_controls'] = controls
    result['gates'] = GATES
    raw = (json.dumps(result,indent=2,sort_keys=True)+'\n').encode()
    args.output.write_bytes(raw)
    print(json.dumps({'gates':GATES,'bank':len(BANK),'literal_controls':controls,
                      'words':result['words_tested'],'hits':len(result['hits']),
                      'receipt_bytes':len(raw),'receipt_sha256':sha256(raw).hexdigest()},sort_keys=True))


if __name__ == '__main__':
    main()
