"""Exact R1-R3 controls and full bidegree (2,2) reconstruction; no packages."""
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path


def require(ok, context):
    if not ok:
        raise RuntimeError(context)


def tournament(n, mask):
    out = [0] * n
    for k, (u, v) in enumerate(combinations(range(n), 2)):
        if mask >> k & 1:
            out[u] |= 1 << v
        else:
            out[v] |= 1 << u
    return out


def reverse(out):
    full = (1 << len(out)) - 1
    return [full ^ r ^ (1 << u) for u, r in enumerate(out)]


def endpoints(out):
    n = len(out)
    dp = [[0] * n for _ in range(1 << n)]
    for u in range(n):
        dp[1 << u][u] = 1
    for s in range(1, 1 << n):
        for u in range(n):
            value = dp[s][u]
            if not value:
                continue
            targets = out[u] & ~s
            while targets:
                bit = targets & -targets
                targets -= bit
                dp[s | bit][bit.bit_length() - 1] += value
    return dp


def h(out):
    return sum(endpoints(out)[-1]) if out else 1


def permutation_h(out):
    return sum(all(out[u] >> v & 1 for u, v in zip(p, p[1:]))
               for p in permutations(range(len(out))))


def extend(out, masks):
    """Each mask bit is OLD -> PORT; ports are in their transitive order."""
    n = len(out)
    rows = out[:] + [0] * len(masks)
    for i, mask in enumerate(masks):
        port = n + i
        for v in range(n):
            if mask >> v & 1:
                rows[v] |= 1 << port
            else:
                rows[port] |= 1 << v
        for j in range(i + 1, len(masks)):
            rows[port] |= 1 << (n + j)
    return rows


def induce(out, labels):
    return [sum(1 << j for j, v in enumerate(labels) if out[u] >> v & 1)
            for u in labels]


def submasks(mask):
    sub = mask
    while True:
        yield sub
        if not sub:
            return
        sub = (sub - 1) & mask


def insert_zero(mask, v):
    low = mask & ((1 << v) - 1)
    return low | ((mask ^ low) << 1)


def cover_data(out):
    n = len(out)
    e = endpoints(out)
    b = endpoints(reverse(out))
    full = (1 << n) - 1
    c = [[0] * n for _ in out]
    for s in range(1, full):
        for u in range(n):
            if e[s][u]:
                for v in range(n):
                    c[u][v] += e[s][u] * b[full ^ s][v]
    return sum(e[-1]), b[-1], e[-1], c


def recover_cover(table, n):
    base = table[0][0]
    one = table[0]
    linear = [one[1 << u] - base for u in range(n)]
    potential = [table[1 << u][0] - one[1 << u] for u in range(n)]
    c = [[0] * n for _ in range(n)]
    for u, v in combinations(range(n), 2):
        symmetric = base + linear[u] + linear[v] - one[(1 << u) | (1 << v)]
        contrast = table[1 << u][1 << v] - table[1 << v][1 << u]
        current = potential[u] - potential[v] - contrast
        require((symmetric + current) % 2 == 0, "ordered cover parity")
        c[u][v] = (symmetric + current) // 2
        c[v][u] = (symmetric - current) // 2
    require(all((potential[u] - linear[u]) % 2 == 0 for u in range(n)), "start parity")
    b = [(potential[u] - linear[u]) // 2 for u in range(n)]
    e = [linear[u] + b[u] - sum(c[u]) for u in range(n)]
    return base, b, e, c


def coefficient_table(table, n):
    small = [s for s in range(1 << n) if s.bit_count() <= 2]
    coeff = {}
    for a in small:
        for b in small:
            value = sum((-1) ** (a.bit_count() + b.bit_count() - i.bit_count() - j.bit_count())
                        * table[i][j] for i in submasks(a) for j in submasks(b))
            if value:
                coeff[a, b] = value
    return coeff


def evaluate(coeff, a, b):
    return sum(value for (i, j), value in coeff.items() if i & a == i and j & b == j)


def basic_universe():
    counts = {"cores": 0, "full_two_port_entries": 0,
              "bidegree_reconstructions": 0, "cover_tomographies": 0,
              "deletion_contractions": 0, "permutation_counts": 0,
              "new_old_weight_identities": 0}
    for n in range(1, 5):
        for mask in range(1 << (n * (n - 1) // 2)):
            out = tournament(n, mask)
            width = 1 << n
            table = [[h(extend(out, [a, b])) for b in range(width)] for a in range(width)]
            counts["cores"] += 1
            counts["full_two_port_entries"] += width * width
            coeff = coefficient_table(table, n)
            for a in range(width):
                for b in range(width):
                    require(evaluate(coeff, a, b) == table[a][b], ("bidegree", n, mask, a, b))
                    counts["bidegree_reconstructions"] += 1
                    if n <= 3:
                        require(permutation_h(extend(out, [a, b])) == table[a][b], "independent H")
                        counts["permutation_counts"] += 1
            direct_data = cover_data(out)
            require(recover_cover(table, n) == direct_data, ("full ordered-cover data", n, mask))
            counts["cover_tomographies"] += 1
            _, starts, ends, c = direct_data
            for a in range(width):
                extended_c = cover_data(extend(out, [a]))[3]
                for u in range(n):
                    expected = starts[u] + ends[u] + sum(c[v][u] for v in range(n))
                    expected += sum((1 - (a >> v & 1)) * (c[u][v] - c[v][u]) for v in range(n))
                    require(extended_c[n][u] + extended_c[u][n] == expected, "new-old weight identity")
                    counts["new_old_weight_identities"] += 1
            for v in range(n):
                deleted = induce(out, [u for u in range(n) if u != v])
                for a in range(1 << (n - 1)):
                    common = insert_zero(a, v)
                    total = 0
                    for alpha in (-1, 1):
                        for beta in (-1, 1):
                            x = common | ((alpha == -1) << v)
                            y = common | ((beta == 1) << v)
                            total += alpha * beta * table[x][y]
                    require(total == 2 * h(extend(deleted, [a])), ("deletion contraction", n, mask, v, a))
                    counts["deletion_contractions"] += 1
    return counts


def three_port_universe():
    contractions = permutations_checked = 0
    for n in range(1, 4):
        for mask in range(1 << (n * (n - 1) // 2)):
            out = tournament(n, mask)
            for v in range(n):
                deleted = induce(out, [u for u in range(n) if u != v])
                for common_small in range(1 << (n - 1)):
                    common = insert_zero(common_small, v)
                    for other in range(1 << n):
                        other_small = (other & ((1 << v) - 1)) | ((other >> (v + 1)) << v)
                        for position in (0, 1):
                            total = 0
                            for alpha in (-1, 1):
                                for beta in (-1, 1):
                                    pair = [common | ((alpha == -1) << v), common | ((beta == 1) << v)]
                                    masks = pair + [other] if position == 0 else [other] + pair
                                    rows = extend(out, masks)
                                    value = h(rows)
                                    if n <= 2:
                                        require(value == permutation_h(rows), "independent three-port H")
                                        permutations_checked += 1
                                    total += alpha * beta * value
                            target = [common_small, other_small] if position == 0 else [other_small, common_small]
                            require(total == 2 * h(extend(deleted, target)), "three-port contraction")
                            contractions += 1
    return {"cores": "every labelled n=1,2,3", "adjacent_pair_positions": [0, 1],
            "contractions": contractions, "independent_permutation_counts": permutations_checked}


def reconstructed_deleted_response(out, deleted, global_masks):
    """Reverse R1 repeatedly; global masks retain original old labels."""
    n = len(out)
    if not deleted:
        return h(extend(out, global_masks)), 1
    v = deleted[-1]
    rest = deleted[:-1]
    total = queries = 0
    for alpha in (-1, 1):
        for beta in (-1, 1):
            first = global_masks[0] & ~(1 << v)
            x = first | ((alpha == -1) << v)
            y = first | ((beta == 1) << v)
            remaining = [mask & ~(1 << v) for mask in global_masks[1:]]
            value, calls = reconstructed_deleted_response(out, rest, [x, y] + remaining)
            total += alpha * beta * value
            queries += calls
    require(total % 2 == 0, "iterated contraction divisibility")
    return total // 2, queries


def repeated_controls():
    trials = queries = 0
    for mask in (0, 1, 12, 31, 63):
        out = tournament(4, mask)
        for deleted in combinations(range(4), 2):
            labels = [u for u in range(4) if u not in deleted]
            core = induce(out, labels)
            for small in range(4):
                global_mask = sum(1 << u for i, u in enumerate(labels) if small >> i & 1)
                got, calls = reconstructed_deleted_response(out, deleted, [global_mask])
                require(calls == 16, "query count per two deletions")
                require(got == h(extend(core, [small])), "iterated response")
                trials += 1
                queries += calls
    return {"n": 4, "core_masks": [0, 1, 12, 31, 63],
            "deleted_size": 2, "trials": trials, "actual_three_port_queries": queries}


def odd_converse_controls():
    summary = []
    for n in range(2, 6):
        symmetric = 0
        for mask in range(1 << (n * (n - 1) // 2)):
            out = tournament(n, mask)
            base, _, _, c = cover_data(out)
            require(sum(c[u][v] for u in range(n) for v in range(n) if out[u] >> v & 1)
                    == (n - 1) * base, "Hamilton incidence identity")
            allh = [sum(row) for row in endpoints(out)]
            allh[0] = 1
            full = (1 << n) - 1
            qtwo = sum(allh[s] * allh[full ^ s] for s in range(1 << n))
            require(qtwo == 2 * base + sum(map(sum, c)), "two-color cover identity")
            require(qtwo % 4 == 0, "complement pairing mod4")
            is_symmetric = all(c[u][v] == c[v][u] for u in range(n) for v in range(n))
            symmetric += is_symmetric
            if n % 2:
                require(not is_symmetric, "odd-order symmetric C")
        summary.append({"n": n, "cores": 1 << (n * (n - 1) // 2), "symmetric_C": symmetric})
    return summary


def hostile():
    out = [6, 12, 24, 17, 3]
    other = reverse(out)
    full = 31
    require(all(h(extend(out, [a])) == h(extend(other, [a])) for a in range(32)), "C5 one-port collision")
    values = [h(extend(t, [full ^ 1, full ^ 2])) for t in (out, other)]
    require(values == [89, 91], "C5 two-port split")
    require(all(row.bit_count() == 2 for row in out), "C5 degree obstruction")
    for v in range(5):
        card = induce(out, [u for u in range(5) if u != v])
        reach = [row | (1 << u) for u, row in enumerate(card)]
        for k in range(4):
            for u in range(4):
                if reach[u] >> k & 1:
                    reach[u] |= reach[k]
        require(all(row == 15 for row in reach), "C5 deletion not strong")
    return {"core": out, "converse": other, "old_to_port_masks": [30, 29], "counts": values}


def main():
    source = Path(__file__)
    digest = sha256(source.read_bytes()).hexdigest()
    result = {"status": "PASS", "source_sha256": digest,
              "basic": basic_universe(), "three_port": three_port_universe(),
              "repeated_contraction": repeated_controls(),
              "odd_converse": odd_converse_controls(), "cycle_five_hostile": hostile()}
    require(sha256(source.read_bytes()).hexdigest() == digest, "source changed during execution")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
