"""aud_forced10.py - classify the forced arcs of the order-10 2-arc-strong classes printed by aud_census (exh 10).

For every cut vertex v of each class, T - v = D_1 => ... => D_m; a forced arc (on every Hamiltonian cycle) is typed
D1 (inside D_1), Dm (inside D_m), junction (D_i -> D_(i+1)), middle (inside D_i, 1 < i < m), skip, or at-v.
The note (section 1) claims only D1, Dm and junction types occur.
"""
import re
import sys

rows_re = re.compile(r"forced-arc class n=(\d+) hc=(\d+) cutvertex=\d forced arcs:((?: \d+->\d+)+) \| scores[^|]*\| rows((?: [01]+)+)")


def strong(n, out, V):
    if V == 0:
        return True
    s = (V & -V).bit_length() - 1
    for direction in (0, 1):
        seen, fr = 1 << s, [s]
        while fr:
            x = fr.pop()
            for y in range(n):
                if (V >> y) & 1 and not (seen >> y) & 1:
                    arc = (out[x] >> y) & 1 if direction == 0 else (out[y] >> x) & 1
                    if arc:
                        seen |= 1 << y
                        fr.append(y)
        if seen != V:
            return False
    return True


def components(n, out, V):
    reach = {}
    for x in range(n):
        if (V >> x) & 1:
            seen, fr = 1 << x, [x]
            while fr:
                y = fr.pop()
                for z in range(n):
                    if (V >> z) & 1 and (out[y] >> z) & 1 and not (seen >> z) & 1:
                        seen |= 1 << z
                        fr.append(z)
            reach[x] = seen
    comps, left = [], V
    while left:
        c = 0
        for x in range(n):
            if (left >> x) & 1 and reach[x] & left == left:
                c |= 1 << x
        comps.append(c)
        left &= ~c
    return comps


tally = {}
lines = open(sys.argv[1]).read().splitlines()
ncls = 0
for line in lines:
    m = rows_re.search(line)
    if not m or int(m.group(1)) != 10:
        continue
    ncls += 1
    n = 10
    arcs = [tuple(map(int, a.split("->"))) for a in m.group(3).split()]
    rows = m.group(4).split()
    out = [sum(1 << w for w in range(n) if r[w] == "1") for r in rows]
    full = (1 << n) - 1
    cuts = [v for v in range(n) if not strong(n, out, full & ~(1 << v))]
    for (a, b) in arcs:
        types = set()
        for v in cuts:
            comps = components(n, out, full & ~(1 << v))
            ia = next((i for i, c in enumerate(comps) if (c >> a) & 1), None)
            ib = next((i for i, c in enumerate(comps) if (c >> b) & 1), None)
            if a == v or b == v:
                t = "at-v"
            elif ia == ib == 0:
                t = "D1"
            elif ia == ib == len(comps) - 1:
                t = "Dm"
            elif ia == ib:
                t = "middle"
            elif ib == ia + 1:
                t = "junction"
            else:
                t = "skip"
            types.add(t)
        key = "/".join(sorted(types))
        tally[key] = tally.get(key, 0) + 1
        if "junction" not in types:
            print("  non-junction forced arc %d->%d, hc=%s, cut vertices %s, rows %s" % (a, b, m.group(2), cuts, " ".join(rows)))
            for v in cuts:
                comps = components(n, out, full & ~(1 << v))
                print("     v=%d: components %s, N+(v) in D_1: %s" % (
                    v, [[x for x in range(n) if (c >> x) & 1] for c in comps],
                    [x for x in range(n) if (comps[0] >> x) & 1 and (out[v] >> x) & 1]))
print("order-10 forced-arc classes parsed: %d; forced arcs by type (over all cut vertices): %s" % (ncls, tally))
