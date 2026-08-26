#!/usr/bin/env python3
"""Exact verifier for the rooted double-clone endpoint tensor.

This intentionally does not import the earlier pair probe.  It compares the
rooted tensor with a literal child subset DP and separately checks the closed
layer-moment and Johnson-gcd formulas.
"""

from itertools import combinations
from math import comb, gcd


def adj_from_label(n, z):
    out = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (z >> k) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            k += 1
    return out


def expand_pair(out, root):
    n = len(out)
    child = [0] * (n + 1)
    for i in range(n):
        child[i] = out[i]
    b = n
    for u in range(n):
        if u == root:
            continue
        if (out[root] >> u) & 1:
            child[b] |= 1 << u
        else:
            child[u] |= 1 << b
    child[root] |= 1 << b
    return child


def endpoint_dp(out):
    n = len(out)
    sz = 1 << n
    incoming = [0] * n
    for u in range(n):
        for v in range(n):
            if (out[u] >> v) & 1:
                incoming[v] |= 1 << u
    end = [[0] * n for _ in range(sz)]
    start = [[0] * n for _ in range(sz)]
    for v in range(n):
        end[1 << v][v] = start[1 << v][v] = 1
    for mask in range(1, sz):
        q = mask
        while q:
            v = (q & -q).bit_length() - 1
            q &= q - 1
            rest = mask ^ (1 << v)
            p = rest & incoming[v]
            while p:
                u = (p & -p).bit_length() - 1
                p &= p - 1
                end[mask][v] += end[rest][u]
            p = rest & out[v]
            while p:
                u = (p & -p).bit_length() - 1
                p &= p - 1
                start[mask][v] += start[rest][u]
    return end, start


def literal_caps(out):
    n = len(out)
    full = (1 << n) - 1
    end, start = endpoint_dp(out)

    def before(mask, x):
        if not mask:
            return 1
        return sum(end[mask][p] for p in range(n)
                   if (mask >> p) & 1 and (out[p] >> x) & 1)

    def after(mask, x):
        if not mask:
            return 1
        return sum(start[mask][p] for p in range(n)
                   if (mask >> p) & 1 and (out[x] >> p) & 1)

    cap = [[0] * n for _ in range(n)]
    for x in range(n):
        for y in range(x + 1, n):
            if not ((out[x] >> y) & 1):
                x0, y0 = y, x
            else:
                x0, y0 = x, y
            rem = full ^ (1 << x0) ^ (1 << y0)
            left = rem
            c = 0
            while True:
                right = rem ^ left
                c += before(left, x0) * after(right, y0)
                c += before(left, y0) * after(right, x0)
                if not left:
                    break
                left = (left - 1) & rem
            cap[x][y] = cap[y][x] = c
    h = sum(end[full])
    return h, cap


def rooted_tensor(out, root):
    """Return the virtual child endpoint oracle from the minimal stored state."""
    n = len(out)
    b = n
    child = expand_pair(out, root)
    end, start = endpoint_dp(out)
    u_mask = ((1 << n) - 1) ^ (1 << root)
    dend = [[0] * n for _ in range(1 << n)]
    dstart = [[0] * n for _ in range(1 << n)]

    def ea(mask):
        whole = mask | (1 << root)
        return sum(end[whole][p] for p in range(n)
                   if (mask >> p) & 1 and (out[p] >> root) & 1)

    def eb(mask):
        return ea(mask) + end[mask | (1 << root)][root]

    def sb(mask):
        whole = mask | (1 << root)
        return sum(start[whole][p] for p in range(n)
                   if (mask >> p) & 1 and (out[root] >> p) & 1)

    def sa(mask):
        return sb(mask) + start[mask | (1 << root)][root]

    # Only masks contained in U are meaningful.  Increasing numeric mask is
    # enough because every recurrence removes one bit.
    for mask in range(1, 1 << n):
        if mask & ~u_mask:
            continue
        q = mask
        while q:
            u = (q & -q).bit_length() - 1
            q &= q - 1
            rest = mask ^ (1 << u)
            p = rest
            while p:
                v = (p & -p).bit_length() - 1
                p &= p - 1
                if (out[v] >> u) & 1:
                    dend[mask][u] += dend[rest][v]
                if (out[u] >> v) & 1:
                    dstart[mask][u] += dstart[rest][v]
            if (out[root] >> u) & 1:
                dend[mask][u] += ea(rest) + eb(rest)
            else:
                dstart[mask][u] += sa(rest) + sb(rest)

    def split(mask):
        ca = bool(mask & (1 << root))
        cb = bool(mask & (1 << b))
        a_mask = mask & u_mask
        return a_mask, ca, cb

    def vend(mask, x):
        a_mask, ca, cb = split(mask)
        if not ca and not cb:
            return end[a_mask][x]
        if ca ^ cb:
            qmask = a_mask | (1 << root)
            return end[qmask][root if x in (root, b) else x]
        if x == root:
            return ea(a_mask)
        if x == b:
            return eb(a_mask)
        return dend[a_mask][x]

    def vstart(mask, x):
        a_mask, ca, cb = split(mask)
        if not ca and not cb:
            return start[a_mask][x]
        if ca ^ cb:
            qmask = a_mask | (1 << root)
            return start[qmask][root if x in (root, b) else x]
        if x == root:
            return sa(a_mask)
        if x == b:
            return sb(a_mask)
        return dstart[a_mask][x]

    return child, vend, vstart


def tensor_caps(out, root):
    child, vend, vstart = rooted_tensor(out, root)
    n = len(child)
    full = (1 << n) - 1

    def before(mask, x):
        if not mask:
            return 1
        return sum(vend(mask, p) for p in range(n)
                   if (mask >> p) & 1 and (child[p] >> x) & 1)

    def after(mask, x):
        if not mask:
            return 1
        return sum(vstart(mask, p) for p in range(n)
                   if (mask >> p) & 1 and (child[x] >> p) & 1)

    cap = [[0] * n for _ in range(n)]
    for x in range(n):
        for y in range(x + 1, n):
            x0, y0 = (x, y) if (child[x] >> y) & 1 else (y, x)
            rem = full ^ (1 << x0) ^ (1 << y0)
            left = rem
            c = 0
            while True:
                right = rem ^ left
                c += before(left, x0) * after(right, y0)
                c += before(left, y0) * after(right, x0)
                if not left:
                    break
                left = (left - 1) & rem
            cap[x][y] = cap[y][x] = c
    h = sum(vend(full, x) for x in range(n))
    return child, h, cap


def values_from_caps(out, h, cap):
    n = len(out)
    vals = []
    for mask in range(1 << n):
        f = h
        for u in range(n):
            if (mask >> u) & 1:
                for v in range(n):
                    if not ((mask >> v) & 1) and (out[u] >> v) & 1:
                        f += cap[u][v]
        vals.append(f)
    return vals


def binom(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def layer_closed(out, h, cap, m):
    """Return (sum, sumsq, anchor, gcd) without enumerating the layer."""
    n = len(out)
    edges = [(u, v, cap[u][v]) for u in range(n) for v in range(n)
             if (out[u] >> v) & 1]
    csum = sum(c for _, _, c in edges)
    qsum = sum(c * c for _, _, c in edges)
    outstar = instar = disjoint = 0
    for ei, (u, v, c) in enumerate(edges):
        for x, y, d in edges[ei + 1:]:
            if u == x:
                outstar += c * d
            elif v == y:
                instar += c * d
            elif len({u, v, x, y}) == 4:
                disjoint += c * d
            # The remaining adjacent case is a directed two-edge path and
            # cannot be simultaneously selected by a directed cut.
    count = binom(n, m)
    p1 = binom(n - 2, m - 1)
    total = count * h + p1 * csum
    sq = (count * h * h + p1 * (2 * h * csum + qsum)
          + 2 * (binom(n - 3, m - 1) * outstar
                 + binom(n - 3, m - 2) * instar
                 + binom(n - 4, m - 2) * disjoint))

    outcap = [sum(cap[i][j] for j in range(n) if (out[i] >> j) & 1)
              for i in range(n)]
    layer_gcd = 0
    k = m - 1
    for i in range(n):
        for j in range(i + 1, n):
            rest = [x for x in range(n) if x not in (i, j)]
            bvals = [cap[j][x] - cap[i][x] for x in rest]
            base = outcap[j] - outcap[i]
            if k == 0:
                pair_gcd = abs(base)
            elif k == len(rest):
                pair_gcd = abs(base - sum(bvals))
            else:
                r0 = bvals[:k]
                pair_gcd = abs(base - sum(r0))
                pivot = bvals[0]
                for x in bvals[1:]:
                    pair_gcd = gcd(pair_gcd, x - pivot)
            layer_gcd = gcd(layer_gcd, pair_gcd)

    anchor_mask = (1 << m) - 1
    anchor = h
    for u, v, c in edges:
        if (anchor_mask >> u) & 1 and not ((anchor_mask >> v) & 1):
            anchor += c
    return total, sq, anchor, layer_gcd


def direct_and_clone_decomposed_packet(out, cap, a, b):
    """Check the symmetric/antisymmetric clone packet decomposition."""
    n = len(out)
    degree = [sum(cap[i]) for i in range(n)]
    signed = [sum(cap[i][j] if (out[i] >> j) & 1 else -cap[i][j]
                  for j in range(n) if j != i) for i in range(n)]
    direct_c = sum(x * y for x, y in zip(degree, signed))
    direct_d = 0
    edges = [(i, j, cap[i][j]) for i in range(n) for j in range(i + 1, n)]
    for ei, (i, j, c) in enumerate(edges):
        for u, v, d in edges[ei + 1:]:
            if len({i, j, u, v}) == 4:
                direct_d += c * d

    vertices = [u for u in range(n) if u not in (a, b)]
    kappa = cap[a][b]
    e = {u: sum(cap[u][v] for v in vertices if v != u) for u in vertices}
    g = {
        u: sum(cap[u][v] if (out[u] >> v) & 1 else -cap[u][v]
               for v in vertices if v != u)
        for u in vertices
    }
    x = {u: cap[a][u] + cap[b][u] for u in vertices}
    y = {u: cap[a][u] - cap[b][u] for u in vertices}
    sign = {u: 1 if (out[a] >> u) & 1 else -1 for u in vertices}
    xsum = sum(x.values())
    ysum = sum(y.values())
    sx = sum(sign[u] * x[u] for u in vertices)
    sy = sum(sign[u] * y[u] for u in vertices)
    c2 = (2 * sum((g[u] - sign[u] * x[u]) * (e[u] + x[u])
                  for u in vertices)
          + sx * xsum + sy * ysum + 2 * kappa * (sx + ysum))
    assert c2 % 2 == 0
    clone_c = c2 // 2

    core_edges = [cap[u][v] for pos, u in enumerate(vertices)
                  for v in vertices[pos + 1:]]
    core_sum = sum(core_edges)
    core_sq = sum(z * z for z in core_edges)
    d4 = (2 * ((core_sum + xsum) ** 2 + 2 * kappa * core_sum
               + core_sq - sum((e[u] + x[u]) ** 2 for u in vertices))
          - xsum * xsum - ysum * ysum
          + sum(x[u] * x[u] + y[u] * y[u] for u in vertices))
    assert d4 % 4 == 0
    clone_d = d4 // 4
    assert (clone_c, clone_d) == (direct_c, direct_d)
    return direct_c, direct_d


def check_one(n, z, root):
    q = adj_from_label(n, z)
    child = expand_pair(q, root)
    h0, c0 = literal_caps(child)
    child1, h1, c1 = tensor_caps(q, root)
    assert child == child1
    assert (h0, c0) == (h1, c1), ('tensor', n, z, root, h0, h1)
    assert c1[root][n] == 2 * literal_caps(q)[0]
    direct_and_clone_decomposed_packet(child, c1, root, n)
    vals = values_from_caps(child, h1, c1)
    for m in range(1, n + 1):
        layer = [vals[s] for s in range(1 << (n + 1)) if s.bit_count() == m]
        got = layer_closed(child, h1, c1, m)
        want = (sum(layer), sum(x * x for x in layer),
                vals[(1 << m) - 1],
                gcd(*(x - layer[0] for x in layer)))
        assert got == want, ('layer', n, z, root, m, got, want)


def main():
    totals = []
    for n in range(2, 6):
        count = 0
        for z in range(1 << (n * (n - 1) // 2)):
            for root in range(n):
                check_one(n, z, root)
                count += 1
        totals.append((n, count))
        print('complete', n, count, 'PASS')

    # Hostile where the square-only Hamilton formula misses one repaired word.
    c3 = [1 << 1, 1 << 2, 1 << 0]
    child, h, cap = tensor_caps(c3, 0)
    print('C3_pair', h, 'internal_cap', cap[0][3], 'caps', cap)

    # Minimal strong hostile found in the strong isomorphism-class stream:
    # pair expansion does not force the interior Johnson lattice to be 2.
    q7 = adj_from_label(7, 1572047)
    t8, h8, c8 = tensor_caps(q7, 2)
    q7_bits = ''.join('1' if (1572047 >> k) & 1 else '0' for k in range(21))
    t8_bits = ''.join('1' if (t8[i] >> j) & 1 else '0'
                      for i in range(8) for j in range(i + 1, 8))
    gcds = tuple(layer_closed(t8, h8, c8, m)[3] for m in range(1, 8))
    assert q7_bits == '111100110011111111101'
    assert t8_bits == '1111001110011111111110100100'
    assert h8 == 387 and gcds == (6,) * 7
    print('q7_gcd6_hostile', q7_bits, 2, t8_bits, h8, gcds)
    print('q10_parent_endpoint_entries', 2 * 10 * (1 << 9))
    print('q10_new_double_clone_entries', 2 * 9 * (1 << 8))
    print('q10_cut_tensor_coordinates', 1 + 54, 'plus_inherited_internal_cap')
    print('TOTAL', totals, 'PASS')


if __name__ == '__main__':
    main()
