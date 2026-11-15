"""Independent exact depth-first traversal for the literal HYP-9732 witnesses."""
from collections import Counter


def paths(rows):
    n = len(rows)
    full = (1 << n) - 1
    def extend(word, seen):
        if seen == full:
            yield word
            return
        choices = rows[word[-1]] & ~seen
        while choices:
            bit = choices & -choices
            choices -= bit
            yield from extend(word + (bit.bit_length() - 1,), seen | bit)
    for v in range(n):
        yield from extend((v,), 1 << v)


def need(value, detail):
    if not value:
        raise ValueError(detail)


def audit(rows, x, expected):
    n = len(rows)
    for i in range(n):
        need(not (rows[i] >> i & 1), ('loop', i))
        for j in range(i):
            need((rows[i] >> j & 1) + (rows[j] >> i & 1) == 1, ('not tournament', i, j))
    q = lambda v: (x >> v) & 1
    extended = tuple(row | ((1 - q(v)) << n) for v, row in enumerate(rows)) + (x,)
    histogram = Counter()
    witness = None
    for path in paths(rows):
        k = sum(q(u)==1 and q(v)==0 for u,v in zip(path,path[1:]))
        if k>=2 and witness is None:
            witness = path
        histogram[k] += 1
    h = sum(histogram.values())
    hplus = reverse = 0
    for path in paths(extended):
        hplus += 1
        j = path.index(n)
        if 0<j<n:
            u,v=path[j-1],path[j+1]
            reverse += not (rows[u] >> v & 1)
    e = sum((k-1)*v for k,v in histogram.items() if k>=2)
    z = histogram[0]
    data = (h,hplus,hplus-2*h,e,z,reverse)
    need(data == expected,(data,expected))
    need(e + reverse - z == hplus-2*h,data)
    print('order',n,'cut',x,'rows',rows)
    print(' (H,Hplus,delta,E,Z,R) =',data,'crossing histogram',dict(sorted(histogram.items())))
    print(' double-crossing path',witness)


if __name__ == '__main__':
    audit((4030,4024,4090,1008,4064,4032,3971,3840,3584,3072,2056,8),4008,
          (2137,4269,-5,299,397,93))
    audit((1982,1976,2042,432,2016,1984,1931,1792,1536,1032,8),1960,
          (815,1629,-1,156,194,37))


def shelf_rows(r, s):
    n = 7 + r + s
    back = {(2, 1), (6, 0), (6, 1), (6, 3)}
    g = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            if (j, i) in back or (i == 3 and j >= 7+r):
                g[j] |= 1 << i
            else:
                g[i] |= 1 << j
    x = (1 << 3) + (1 << 5) + sum(1 << i for i in range(7,n))
    t = [a | (0 if x >> i & 1 else 1 << n) for i,a in enumerate(g)] + [x]
    return g, t


def shelf_stencil(rows, r, s):
    def color(v):
        if 7 <= v < 7+r:
            return 1
        if 7+r <= v < 7+r+s:
            return 2
        return 0
    results = Counter()
    for p in paths(rows):
        runs = [0]*3
        last = 0
        for v in p:
            c = color(v)
            if c != last:
                runs[c] += 1
            last = c
        results[runs[1], runs[2]] += 1
    need(all(1 <= a <= 2 and 1 <= b <= 2 for a,b in results), results)
    out = {}
    for (a,b), count in results.items():
        d = 2**((a==2)+(b==2))
        need(count % d == 0, (a,b,count))
        out[a,b] = count//d
    return out


def audit_shelves():
    stencils = [shelf_stencil(rows,2,2) for rows in shelf_rows(2,2)]
    print(' two-slot skeleton counts', stencils)
    expected = {(1,1):(173,361), (1,2):(407,831), (2,1):(309,627),
                (2,2):(815,1629), (1,3):(875,1771), (3,1):(581,1159),
                (3,2):(1631,3225), (2,3):(1827,3633), (3,3):(3731,7357)}
    for (r,s), result in expected.items():
        u,v = 2**r,2**s
        g,t=shelf_rows(r,s)
        n=len(g)
        x=t[-1]
        runs = Counter(sum((x >> a & 1) and not (x >> b & 1)
                           for a,b in zip(p,p[1:])) for p in paths(g))
        rev = total = 0
        for p in paths(t):
            total+=1
            j=p.index(n)
            if 0<j<n:
                rev += not (g[p[j-1]] >> p[j+1] & 1)
        direct = sum(runs.values()),total
        need(max(runs) == 2 and
             (runs[2],runs[0],rev) == (12*u*v-12*u+12,
                                        17*u*v-17*u-8*v+22,
                                        2*u*v-2*u-v+17), (r,s,runs,rev))
        fits = (68*u*v-68*u-19*v+75,133*u*v-133*u-31*v+157)
        skeleton = tuple(sum(c*(u-2)**(a-1)*(v-2)**(b-1)
                             for (a,b),c in m.items()) for m in stencils)
        need(direct == result == fits == skeleton,(r,s,direct,result,fits,skeleton))
        print(' shelf', (r,s), direct)


if __name__ == '__main__':
    audit_shelves()
