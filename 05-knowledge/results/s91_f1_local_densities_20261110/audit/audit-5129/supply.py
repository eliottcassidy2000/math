# audit-5129: the product-law supply of THM-5129 (A2), written from the statement.
#
#   Slopes(D) = Slopes_0(D) * {2^j},
#   Slopes_0(D) = { d1 u w : D = d1 h d3, d1 attained, u in Root(h), w in Suf(d3) }.
#
# Every set is stored as a dict value -> least order of a realisation, so that the same tables give the least
# order of a rooted tournament with a given pair (D, U).


def build(strong, Dmax):
    """strong: iterable of (order, H, Pi_2, list of U_C(y)) over strong tournaments of order >= 3 (a list that
    contains every isomorphism class with H <= Dmax; repetitions are harmless)."""
    Root = {1: {1: 1}}          # the root is a trivial component: h = 1, u = 1, order 1
    PI = {}
    minord = {}
    for (o, h, p2, Ul) in strong:
        if h > Dmax:
            continue
        if minord.get(h, 99) > o:
            minord[h] = o
        d = Root.setdefault(h, {})
        for u in Ul:
            if d.get(u, 99) > o:
                d[u] = o
        d2 = PI.setdefault(h, {})
        if d2.get(p2, 99) > o:
            d2[p2] = o
    # attained values (H of some tournament = product of H of strong components), with least order; the empty
    # prefix has H = 1 and order 0
    Aord = {1: 0}
    for d in range(3, Dmax + 1, 2):
        best = None
        for h in minord:
            if d % h == 0 and (d // h) in Aord:
                c = minord[h] + Aord[d // h]
                if best is None or c < best:
                    best = c
        if best is not None:
            Aord[d] = best
    # Suf(d): products of Pi_2 over ordered factorisations of d into strong tournaments of order >= 3
    Suf = {1: {1: 0}}
    for d in range(3, Dmax + 1, 2):
        s = {}
        for h in PI:
            if d % h == 0 and (d // h) in Suf:
                for p2, o in PI[h].items():
                    for w, ow in Suf[d // h].items():
                        val = p2 * w
                        if s.get(val, 99) > o + ow:
                            s[val] = o + ow
        if s:
            Suf[d] = s
    S0 = {}
    for D in range(1, Dmax + 1, 2):
        s = {}
        for d1 in Aord:
            if D % d1:
                continue
            rest = D // d1
            for h in Root:
                if rest % h:
                    continue
                d3 = rest // h
                if d3 not in Suf:
                    continue
                for u, ou in Root[h].items():
                    for w, ow in Suf[d3].items():
                        val = d1 * u * w
                        o = Aord[d1] + ou + ow
                        if s.get(val, 99) > o:
                            s[val] = o
        if s:
            S0[D] = s
    return Root, PI, Aord, Suf, S0, minord


def odd_part(x):
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return x, v


def delta(S0D):
    """odd part -> least exponent j."""
    out = {}
    for s0 in S0D:
        d, v = odd_part(s0)
        if out.get(d, 99) > v:
            out[d] = v
    return out


def minorder(S0D, U):
    """Least order of a rooted tournament with offset D and slope U (None if U is not a slope)."""
    best = None
    for s0, o in S0D.items():
        if U % s0 == 0:
            q = U // s0
            if q & (q - 1) == 0:
                j = q.bit_length() - 1
                if best is None or o + j < best:
                    best = o + j
    return best


def is_slope(DL, U):
    d, v = odd_part(U)
    return d in DL and v >= DL[d]
