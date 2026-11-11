#!/usr/bin/env python3
"""Independent F2 module orbits and actual central-switch graph controls."""
import itertools
import json
from pathlib import Path


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def matrices(r, c):
    mask = (1 << c) - 1
    for code in range(1 << (r*c)):
        yield tuple((code >> (i*c)) & mask for i in range(r))


def rank(rows):
    pivots = {}
    for x in rows:
        while x:
            p = x.bit_length()-1
            if p in pivots:
                x ^= pivots[p]
            else:
                pivots[p] = x
                break
    return len(pivots)


def mul(a, b):
    out = []
    for row in a:
        s = 0
        for j, br in enumerate(b):
            if row >> j & 1:
                s ^= br
        out.append(s)
    return tuple(out)


def inverse(a):
    n = len(a)
    rows = [x | (1 << (n+i)) for i, x in enumerate(a)]
    for j in range(n):
        p = next(i for i in range(j, n) if rows[i] >> j & 1)
        rows[p], rows[j] = rows[j], rows[p]
        for i in range(n):
            if i != j and rows[i] >> j & 1:
                rows[i] ^= rows[j]
    return tuple(x >> n for x in rows)


def normal(m):
    # Ordered multiplicities: S1,S2,L12,L21,L121.
    a, b = [], []
    for kind, copies in enumerate(m):
        for _ in range(copies):
            v, w = len(a), len(b)
            if kind == 0:
                a.append(0)
            elif kind == 1:
                b.append(0)
            elif kind == 2:
                a.append(0)
                b.append(1 << v)
            elif kind == 3:
                a.append(1 << w)
                b.append(0)
            else:
                a.extend((0, 1 << w))
                b.append(1 << v)
    return tuple(a), tuple(b)


def module_census():
    gl = {n: [(p, inverse(p)) for p in matrices(n,n) if rank(p) == n]
          for n in range(4)}
    report = []
    total = 0
    for v in range(4):
        for w in range(4):
            if v+w > 5:
                continue
            lookup = {}
            classes = 0
            for m in itertools.product(range(4), repeat=5):
                if m[0]+m[2]+m[3]+2*m[4] != v:
                    continue
                if m[1]+m[2]+m[3]+m[4] != w:
                    continue
                a,b = normal(m)
                classes += 1
                for p,pi in gl[v]:
                    for q,qi in gl[w]:
                        key = (mul(mul(p,a),qi), mul(mul(q,b),pi))
                        require(key not in lookup or lookup[key] == m,
                                ('overlapping normal-form orbits',v,w,m))
                        lookup[key] = m
            valid = 0
            for a in matrices(v,w):
                for b in matrices(w,v):
                    if any(mul(b,a)):
                        continue
                    valid += 1
                    ra,rb,r = rank(a),rank(b),rank(mul(a,b))
                    expected = (v-ra-rb,w-ra-rb+r,rb-r,ra-r,r)
                    require(lookup.get((a,b)) == expected,
                            ('rank/orbit disagreement',v,w,a,b,expected))
            require(valid == len(lookup),('orbit completeness',v,w))
            total += valid
            report.append({'dimensions':[v,w], 'valid_pairs':valid,
                           'isomorphism_classes':classes})
    return {'field':2,'total_valid_pairs':total,'by_dimension':report}


def hp(rows):
    n = len(rows)
    full = (1 << n)-1
    dp = [[0]*n for _ in range(full+1)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1,full+1):
        for v,ways in enumerate(dp[mask]):
            nxt = rows[v] & (full ^ mask)
            while ways and nxt:
                bit = nxt & -nxt
                nxt -= bit
                dp[mask | bit][bit.bit_length()-1] += ways
    vals = [sum(x) for x in dp]
    vals[0] = 1
    return vals,dp[-1]


def positive_signature(rows, cut):
    n = len(rows)
    full = (1 << n)-1
    residual, ends = hp(rows)
    rev = [sum((1 << u) for u in range(n) if rows[u] >> v & 1)
           for v in range(n)]
    _,starts = hp(rev)
    kp,ip,op,jp,jm = 0,0,0,0,0
    polynomials = {name:[0]*(n+1) for name in ('I','O','J')}
    def walk(mask,first,last):
        nonlocal kp,ip,op,jp,jm
        weight = residual[full ^ mask]
        degree = mask.bit_count()
        kp += weight
        if not (cut >> last & 1):
            ip += weight
            polynomials['I'][degree] += weight
        if cut >> first & 1:
            op += weight
            polynomials['O'][degree] += weight
            if not (cut >> last & 1):
                jp += weight
                jm += (-1 if mask.bit_count() & 1 else 1)*weight
                polynomials['J'][degree] += weight
        nxt = rows[last] & (full ^ mask)
        while nxt:
            bit = nxt & -nxt
            nxt -= bit
            walk(mask | bit,first,bit.bit_length()-1)
    for v in range(n):
        walk(1 << v,v,v)
    d = residual[-1]
    s = sum(starts[v] for v in range(n) if not (cut >> v & 1))
    e = sum(ends[v] for v in range(n) if cut >> v & 1)
    return (d+jp+jm,d,d+kp,ip,op,s,e), (jp,jm), polynomials


def g(a):
    answer = 1
    prefix = 0
    for bit in bin(a)[2:]:
        if bit == '0':
            answer += prefix+(prefix & 1)
        prefix = 2*prefix+int(bit)
    return answer


def literal(rows,cut,a,b):
    left = [int(x) for x in bin(a)[2:]] if a else []
    right = [1-int(x) for x in bin(b)[2:][::-1]] if b else []
    l,n = len(left),len(rows)
    old = l+n+len(right)
    out = [0]*(old+1)
    for i in range(old):
        for j in range(i+1,old):
            forward = not (l <= i < j < l+n) or bool(rows[i-l] >> (j-l) & 1)
            u,v = (i,j) if forward else (j,i)
            out[u] |= 1 << v
    for v,bit in enumerate(left+[(cut >> j) & 1 for j in range(n)]+right):
        u,w = (old,v) if bit else (v,old)
        out[u] |= 1 << w
    return out


def physical_pair():
    rows = (64,121,35,37,13,81,28)
    signatures = [positive_signature(rows,c) for c in (12,20)]
    p,q = [s[0] for s in signatures]
    require(p[1:] == q[1:] and q[0]-p[0] == 2,signatures)
    for name in ('I','O'):
        require(signatures[0][2][name] == signatures[1][2][name],('polynomial',name))
    delta = [y-x for x,y in zip(signatures[0][2]['J'],signatures[1][2]['J'])]
    require(delta == [0,0,-2,-6,1,3,2,0],('central polynomial',delta))
    tests = 0
    for a,b in itertools.product(range(4),repeat=2):
        counts = []
        for cut,(sig,_,polynomials) in zip((12,20),signatures):
            h,d,l,i,o,s,e = sig
            expected = h+d*(g(a)+g(b)-2)+l*a*b+i*a+o*b+s*(a&1)+e*(b&1)
            actual = hp(literal(rows,cut,a,b))[0][-1]
            require(actual == expected,('physical pair',cut,a,b,actual,expected))
            counts.append(actual)
            tests += 1
        require(counts[1]-counts[0] == 2,('gap',a,b,counts))
    # The central matrix E13 has square zero and M E13=E13 M=d E13.
    # Test over integers independently of the binary-field orbit routine.
    def mm(a,b):
        return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    c = [[0,0,1],[0,0,0],[0,0,0]]
    for d,l,o,i,j in itertools.product(range(3),repeat=5):
        m = [[d,o,j],[0,l,i],[0,0,d]]
        target = [[d*x for x in row] for row in c]
        require(mm(m,c) == target == mm(c,m),'central propagation')
    require(mm(c,c) == [[0]*3 for _ in range(3)],'nilpotent')
    return {'rows':rows,'cuts':[12,20],'signatures':signatures,
            'actual_graph_checks':tests,'algebra_integer_checks':243}


def main():
    out = {'scope':'Complete stated F2 matrix universe; physical inherited B95 switch; all-field classification is proved in FIVE_MODULES.md.',
           'module_census':module_census(),'physical_pair':physical_pair()}
    path = Path(__file__).with_name('module_checks.json')
    path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'valid_pairs':out['module_census']['total_valid_pairs'],
                      'graph_checks':out['physical_pair']['actual_graph_checks'],
                      'status':'PASS'}))


if __name__ == '__main__':
    main()
