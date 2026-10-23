#!/usr/bin/env python3
"""Exact controls for restricted observer coupling and composite-base clocks.

No imports from producer packets; no floating point and no Python asserts.
Run from any directory. The printed receipt is deterministic under python -O.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import comb, factorial, prod
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def edges(n):
    return list(combinations(range(n), 2))


def masks(n, group):
    ee = edges(n)
    return [sum((p[a] > p[b]) << i for i, (a, b) in enumerate(ee)) for p in group]


def fwht(values):
    a = list(values)
    step = 1
    while step < len(a):
        for base in range(0, len(a), 2 * step):
            for j in range(base, base + step):
                x, y = a[j], a[j + step]
                a[j], a[j + step] = x + y, x - y
        step *= 2
    return a


def signed_sums(n, group):
    f = [0] * (1 << comb(n, 2))
    for x in masks(n, group):
        f[x] += 1
    return fwht(f)


def compositions(n):
    if n == 0:
        yield ()
    else:
        for k in range(1, n + 1):
            for rest in compositions(n - k):
                yield (k,) + rest


def young_group(parts):
    blocks, start = [], 0
    for k in parts:
        blocks.append(tuple(range(start, start + k)))
        start += k
    return [tuple(v for block in pp for v in block)
            for pp in product(*(permutations(b) for b in blocks))]


def image_mask(x, p, ee, lookup):
    out = 0
    for i, (a, b) in enumerate(ee):
        if (x >> i) & 1:
            out |= 1 << lookup[tuple(sorted((p[a], p[b])))]
    return out


def partitions(n, group, affine=False, paired=False):
    """Actual word orbits, independently of the Fourier calculation."""
    ee = edges(n)
    lookup = {e: i for i, e in enumerate(ee)}
    size = 1 << len(ee)
    transforms = []
    for p in group:
        # Relabel the reference tournament, so inversions live at IMAGE edges.
        eps = sum(1 << lookup[tuple(sorted((p[a], p[b])))]
                  for a, b in ee if p[a] > p[b]) if affine else 0
        transforms.append((p, eps))
    labels = [-1] * size
    orbits = []
    for x in range(size):
        if labels[x] >= 0:
            continue
        orbit = {image_mask(x, p, ee, lookup) ^ eps for p, eps in transforms}
        if paired:
            orbit |= {y ^ (size - 1) for y in list(orbit)}
        idx = len(orbits)
        for y in orbit:
            require(labels[y] in (-1, idx), "orbits overlap")
            labels[y] = idx
        orbits.append(sorted(orbit))
    return labels, orbits


def incidence(labels_g, labels_t, carrier):
    W = Counter((labels_g[x], labels_t[x]) for x in carrier)
    row, col = Counter(), Counter()
    for (g, t), count in W.items():
        row[g] += count
        col[t] += count
    return W, row, col


def apply_incidence(W, row, col, vector):
    middle = {t: Q(0) for t in col}
    for (g, t), count in W.items():
        middle[t] += Q(count, col[t]) * vector[g]
    result = {g: Q(0) for g in row}
    for (g, t), count in W.items():
        result[g] += Q(count, row[g]) * middle[t]
    return result


def audit_matrix(n, group, paired):
    gg, go = partitions(n, group, paired=paired)
    tt, _ = partitions(n, group, affine=True, paired=paired)
    W, row, col = incidence(gg, tt, range(len(gg)))
    sums = signed_sums(n, group)
    # Every H-orbit of masks gives an unsigned Walsh orbit vector.
    _, mask_orbits = partitions(n, group)
    count = 0
    for orbit in mask_orbits:
        F = orbit[0]
        if paired and F.bit_count() % 2:
            continue
        vector = {}
        for g, words in enumerate(go):
            x = words[0]
            value = Q(sum((-1) ** ((x & f).bit_count()) for f in orbit), len(orbit))
            for y in words:
                require(value == Q(sum((-1) ** ((y & f).bit_count()) for f in orbit), len(orbit)),
                        "Walsh vector not graph invariant")
            vector[g] = value
        lam = Q(sums[F], len(group)) ** 2
        require(apply_incidence(W, row, col, vector) ==
                {g: lam * v for g, v in vector.items()}, "subgroup spectral equation")
        count += 1
    require(count == len(go), "Walsh basis dimension")
    return count


def coupling_checks():
    full = {n: signed_sums(n, list(permutations(range(n)))) for n in range(1, 7)}
    product_checks = 0
    block_cases = 0
    for n in range(1, 7):
        ee = edges(n)
        for parts in compositions(n):
            group = young_group(parts)
            values = signed_sums(n, group)
            local_maps, start = [], 0
            for k in parts:
                local_maps.append((k, [ee.index((start + a, start + b)) for a, b in edges(k)]))
                start += k
            for f, val in enumerate(values):
                expected = 1
                for k, indices in local_maps:
                    ff = sum(((f >> j) & 1) << i for i, j in enumerate(indices))
                    expected *= full[k][ff]
                require(val == expected, "contiguous-block factorization")
                product_checks += 1
            residual = [Q(v, len(group)) ** 2 for v in values if abs(v) < len(group)]
            if max(parts) >= 3:
                require(max(residual) == Q(1, 9), "block relative gap")
                # Pairing keeps even masks; the P3 extremizer must survive.
                require(max(Q(v, len(group)) ** 2 for f, v in enumerate(values)
                            if f.bit_count() % 2 == 0 and abs(v) < len(group)) == Q(1, 9),
                        "paired block gap")
            elif max(parts) == 2:
                require(max(residual) == 0, "size-two blocks")
            else:
                require(not residual, "singleton block space")
            block_cases += 1
    fixed_port_rows = []
    for n in range(3, 7):
        ee = edges(n)
        for r in range(n):
            group = [p for p in permutations(range(n)) if p[r] == r]
            values = signed_sums(n, group)
            for k in range(n):
                # A particular star; signs depend on its reference positions.
                leaves = [v for v in range(n) if v != r][:k]
                f = sum(1 << ee.index(tuple(sorted((r, v)))) for v in leaves)
                total = sum((-1) ** j * comb(r, j) * comb(n - 1 - r, k - j)
                            for j in range(k + 1) if j <= r and 0 <= k - j <= n - 1 - r)
                require(abs(Q(values[f], len(group))) == abs(Q(total, comb(n - 1, k))),
                        "fixed-port Krawtchouk formula")
            residual = [Q(v, len(group)) ** 2 for v in values if abs(v) < len(group)]
            even_residual = [Q(v, len(group)) ** 2 for f, v in enumerate(values)
                             if f.bit_count() % 2 == 0 and abs(v) < len(group)]
            fixed_port_rows.append({"n": n, "port_1_based": r + 1,
                                    "largest_residual_eigenvalue": str(max(residual, default=Q(0))),
                                    "paired_largest_residual_eigenvalue": str(max(even_residual, default=Q(0)))})
    matrix_modes = 0
    for n in range(2, 5):
        groups = [young_group(c) for c in compositions(n)]
        groups += [[p for p in permutations(range(n)) if p[r] == r] for r in range(n)]
        # An arbitrary non-Young cyclic subgroup is a useful hostile control.
        groups.append([tuple((i + k) % n for i in range(n)) for k in range(n)])
        for group in groups:
            for paired in (False, True):
                matrix_modes += audit_matrix(n, group, paired)
    return {"contiguous_block_cases": block_cases, "mask_product_checks": product_checks,
            "independent_incidence_eigenvectors": matrix_modes, "fixed_port_rows": fixed_port_rows}


def component_sizes(W, row, col):
    adj = defaultdict(set)
    for g, t in W:
        adj[(0, g)].add((1, t))
        adj[(1, t)].add((0, g))
    unseen = set(adj)
    result = []
    while unseen:
        todo = [next(iter(unseen))]
        seen = set(todo)
        while todo:
            for v in adj[todo.pop()]:
                if v not in seen:
                    todo.append(v)
                    seen.add(v)
        unseen -= seen
        gs = [i for side, i in seen if side == 0]
        ts = [i for side, i in seen if side == 1]
        result.append((len(gs), len(ts), sum(row[g] for g in gs)))
    return sorted(result)


def path_checks():
    receipt = []
    for n in range(3, 7):
        ee = edges(n)
        group = list(permutations(range(n)))
        pathmask = sum(1 << ee.index((i, i + 1)) for i in range(n - 1))
        carrier = [x for x in range(1 << len(ee)) if not x & pathmask]
        for paired in (False, True):
            gg, _ = partitions(n, group, paired=paired)
            tt, _ = partitions(n, group, affine=True, paired=paired)
            W, row, col = incidence(gg, tt, carrier)
            require(row[gg[0]] == col[tt[0]] == W[gg[0], tt[0]] == 1,
                    "fixed-path absorbing transitive state")
            receipt.append({"n": n, "paired": paired, "words": len(carrier),
                            "component_sizes_graph_tournament_words": component_sizes(W, row, col)})
            if n >= 4:
                # Remove the isolated transitive state. The remaining carrier is
                # connected in these finite cases; do not assume it for all n.
                positive = [x for x in carrier if x]
                ww, rr, cc = incidence(gg, tt, positive)
                require(len(component_sizes(ww, rr, cc)) == 1, "nontransitive finite component")
                one_edge = next(x for x in positive if x.bit_count() == 1)
                g = gg[one_edge]
                mass = comb(n - 1, 2)
                require(rr[g] == mass, "one-edge row mass")
                diag = sum(Q(count * count, rr[g] * cc[t]) for (j, t), count in ww.items() if j == g)
                predicted = Q(n - 2, mass) + sum(Q(n - ell, mass * (2**(ell-1) + 1))
                                                for ell in range(3, n))
                require(diag == predicted, "one-edge return formula")
                pi = Q(mass, len(positive))
                rayleigh = (diag - pi) / (1 - pi)
                require(rayleigh > Q(1, 9), "even after transitive deletion the 8/9 gap fails")
                if n == 4:
                    expected = [[Q(11,15), Q(2,15), Q(1,15), Q(1,15)]] + [[Q(1,5), Q(2,5), Q(1,5), Q(1,5)]]*3
                    keys = sorted(rr)
                    actual = [[sum(Q(ww.get((j,t),0)*ww.get((k,t),0),rr[j]*cc[t]) for t in cc)
                               for k in keys] for j in keys]
                    require(actual == expected, "exact four-row path component")
                    f = {j: Q(j == g) - pi for j in rr}
                    require(apply_incidence(ww, rr, cc, f) == {j: Q(8,15)*v for j,v in f.items()},
                            "four-row extremal eigenvector")
                receipt[-1]['nontransitive_one_edge_return'] = str(diag)
                receipt[-1]['nontransitive_rayleigh_eigenvalue_lower_bound'] = str(rayleigh)
        # Converse transported back to the path chart is REFLECTION, not a free-edge flip.
        rev = tuple(reversed(range(n)))
        lookup = {e: i for i, e in enumerate(ee)}
        for x in carrier:
            reflected = image_mask(x, rev, ee, lookup)
            require(reflected in carrier, "path-chart converse")
    return receipt


def single_defect_paths():
    # A separate subset DP, using actual arcs, checks the inherited two-apex
    # count. The new return-probability formula uses this count as a dependency.
    result = []
    for k in range(3, 13):
        out = [(sum(1 << j for j in range(i + 1, k))) for i in range(k)]
        out[0] ^= 1 << (k-1)
        out[k-1] |= 1
        dp = [[0]*k for _ in range(1 << k)]
        for i in range(k):
            dp[1 << i][i] = 1
        for mask in range(1, 1 << k):
            for end in range(k):
                if dp[mask][end]:
                    nxt = out[end] & ~mask
                    while nxt:
                        bit = nxt & -nxt
                        nxt -= bit
                        dp[mask|bit][bit.bit_length()-1] += dp[mask][end]
        H = sum(dp[-1])
        require(H == 1 + 2**(k-2), "actual one-long-backedge path count")
        result.append([k, H])
    return result


def base_divisibility_checks():
    primes = (2, 3, 5, 7, 11, 13)
    receipt = []
    for exps in ((5, 7, 3, 4, 1, 2), (4, 7, 3, 4, 1, 2)):
        N = prod(p ** e for p, e in zip(primes, exps))
        hist = Counter()
        for powers in product(*(range(e + 1) for e in exps)):
            if not any(powers):
                continue
            a = prod(p ** r for p, r in zip(primes, powers))
            predicted = min(e // r for e, r in zip(exps, powers) if r)
            # Independent literal repeated-division oracle.
            b, actual = N, 0
            while b % a == 0:
                b //= a
                actual += 1
            require(predicted == actual, "composite-base valuation")
            hist[actual] += 1
        for k in range(1, max(exps) + 1):
            tau = prod(e // k + 1 for e in exps)
            nxt = prod(e // (k + 1) + 1 for e in exps)
            require(hist[k] == tau - nxt, "exact valuation level counts")
        def D(a):
            if a <= 1:
                return a
            b, value = N, 0
            while b % a == 0:
                b //= a
                value += 1
            return value
        core = {a: D(a) for a in range(8)}
        cyc = [4, 2, 5, 3, 7] if exps[0] == 5 else [2, 4]
        require(all(D(a) == b for a, b in zip(cyc, cyc[1:] + cyc[:1])), "clock cycle")
        require(D(13) == 2, "thirteen entrance")
        receipt.append({"N": N, "exponents": exps, "core": core,
                        "cycle": cyc, "levels": dict(sorted(hist.items())),
                        "nontrivial_cycle_basin_divisors": sum(v for k, v in hist.items() if k >= 2)})
    return receipt


def projective_checks():
    receipt = []
    for q in (2, 3, 5):
        def normalize(v):
            a = next(x for x in v if x)
            return tuple(x * pow(a, -1, q) % q for x in v)
        points = sorted({normalize(v) for v in product(range(q), repeat=3) if any(v)})
        M = [[int(sum(a * b for a, b in zip(p, line)) % q == 0) for line in points]
             for p in points]
        v = q * q + q + 1
        require(len(points) == v and all(sum(row) == q + 1 for row in M), "projective plane sizes")
        for i in range(v):
            for j in range(v):
                require(sum(a * b for a, b in zip(M[i], M[j])) == 1 + q * (i == j),
                        "incidence Gram")
        receipt.append({"q": q, "points": v, "line_size": q + 1,
                        "gap": str(Q(v, (q + 1) ** 2))})
    return receipt


if __name__ == "__main__":
    print(json.dumps({"status": "FINITE-EXACT PASS", "coupling": coupling_checks(),
                      "path_restrictions": path_checks(), "single_defect_paths": single_defect_paths(),
                      "arithmetic_clocks": base_divisibility_checks(),
                      "projective_incidence": projective_checks()}, indent=2, sort_keys=True))
