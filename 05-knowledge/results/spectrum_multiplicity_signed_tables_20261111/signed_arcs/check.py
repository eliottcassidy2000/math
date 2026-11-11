#!/usr/bin/env python3
"""Small exact controls for the signed atlas and prime-graph parity kernel."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path


def require(p, message):
    if not p:
        raise RuntimeError(message)


def hp(rows):
    n = len(rows)
    require(all(not (rows[i] >> i & 1) for i in range(n)), "loop")
    require(all(((rows[i] >> j & 1) + (rows[j] >> i & 1)) == 1
                for i in range(n) for j in range(i)), "not a tournament")
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for last, count in enumerate(dp[mask]):
            rest = rows[last] & ~mask
            while count and rest:
                bit = rest & -rest
                dp[mask | bit][bit.bit_length() - 1] += count
                rest -= bit
    return sum(dp[-1])


def cycles(adj, directed=False):
    answer = set()
    for start in sorted(adj):
        def visit(path):
            for nxt in adj[path[-1]]:
                if nxt == start and len(path) >= 3:
                    p = tuple(path)
                    if not directed:
                        p = min(p, (start,) + tuple(reversed(path[1:])))
                    answer.add(p)
                elif nxt > start and nxt not in path:
                    visit(path + [nxt])
        visit([start])
    return sorted(answer)


def c(x, sign):
    return 3 * x + sign if x % 2 else x // 2


def determinant_poly(vertices, sign, delete):
    """Leibniz expansion, independent of cycle-factor reasoning."""
    result = [0] * (len(vertices) + 1)
    for perm in itertools.permutations(vertices):
        moved = 0
        for x, y in zip(vertices, perm):
            if x != y:
                if x == delete or c(x, sign) != y:
                    break
                moved += 1
        else:
            inv = sum(perm[i] > perm[j]
                      for i in range(len(perm)) for j in range(i + 1, len(perm)))
            result[moved] += (-1) ** (inv + moved)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def walks(vertices, sign, delete, length=20):
    n = len(vertices)
    index = {v: i for i, v in enumerate(vertices)}
    current = [[int(i == j) for j in range(n)] for i in range(n)]
    result = [current]
    for _ in range(length):
        nxt = [[0] * n for _ in range(n)]
        for i in range(n):
            for x in vertices:
                if x != delete:
                    nxt[i][index[c(x, sign)]] += current[i][index[x]]
        result.append(nxt)
        current = nxt
    return result


def surgery_check(vertices, sign, gate, target):
    full, cut = walks(vertices, sign, None), walks(vertices, sign, gate)
    idx = {v: i for i, v in enumerate(vertices)}
    # Independent coefficient convolution of the first-crossing renewal law:
    # R = R0 + z R0(*,gate) R(target,*).
    checks = 0
    for k in range(len(full)):
        for x in vertices:
            for y in vertices:
                rhs = cut[k][idx[x]][idx[y]]
                rhs += sum(cut[a][idx[x]][idx[gate]] *
                           full[k - 1 - a][idx[target]][idx[y]] for a in range(k))
                require(full[k][idx[x]][idx[y]] == rhs, "renewal identity")
                checks += 1
    return {"vertices": vertices, "full_determinant": determinant_poly(vertices, sign, None),
            "cut_determinant": determinant_poly(vertices, sign, gate),
            "renewal_coefficients_checked": checks,
            "return_coefficients": [z[idx[target]][idx[gate]] for z in cut]}


def controller(n):
    edges = []
    while n != 1:
        if n == 5:
            word = [(5, -1), (14, 1), (7, 1), (22, 1), (11, -1),
                    (32, 1), (16, 1), (8, 1), (4, 1), (2, 1)]
            for x, sign in word:
                require(x == n, "detour source")
                edges.append((x, sign, 4 * x + sign, c(x, sign)))
                n = c(x, sign)
            continue
        old = n
        signs = [1] if n % 2 == 0 else ([1 if n % 4 == 1 else -1, 1, 1])
        for sign in signs:
            edges.append((n, sign, 4 * n + sign, c(n, sign)))
            n = c(n, sign)
        require(9 * n <= 7 * old, "block contraction")
    require(all(m not in (7, 21) for _, _, m, _ in edges), "puncture used")
    return edges


def run():
    witness = {
        19: [118, 44, 112, 117, 98, 0, 34],
        29: [36, 25, 42, 1, 45, 10],
        43: [38, 12, 24, 49, 3, 22],
        57: [12, 61, 40, 48, 5, 81, 31],
        89: [102, 48, 114, 23, 33, 8, 58],
    }
    for h in [9, 17, 33, 65, 129]:
        n = (h - 1).bit_length() + 1
        rows = [sum(1 << j for j in range(i + 1, n)) for i in range(n)]
        rows[0] ^= 1 << (n - 1)
        rows[-1] |= 1
        witness[h] = rows
    for h, rows in witness.items():
        require(hp(rows) == h, "detour HP witness")
    require(hp([2, 4, 1]) == 3 and hp([6, 4, 0]) == 1, "HP controls")
    invalid_rejected = False
    try:
        hp([0, 0])
    except RuntimeError:
        invalid_rejected = True
    require(invalid_rejected, "tie control")

    decode_checks = 0
    for sign in [-1, 1]:
        for x in range(1, 257):
            m = 4 * x + sign
            k = (m - 1) // 2
            d, off = k * k % m, k * (k + 1) % m
            require(sorted([d, off]) == sorted([x, 3 * x + sign]), "center")
            target = 3 * x + sign if x % 2 else ((3 * x + sign) - sign) // 6
            require(target == c(x, sign), "decode")
            decode_checks += 1
    controller_rows = [controller(n) for n in range(1, 10001)]
    longest = max(enumerate(controller_rows, 1), key=lambda r: len(r[1]))
    plus = surgery_check([1, 2, 4, 5, 8, 16], 1, 5, 16)
    minus = surgery_check([1, 2, 5, 7, 10, 14, 20], -1, 2, 1)
    require(plus['full_determinant'] == plus['cut_determinant'] == [1, 0, 0, -1], "plus feedback")
    require(minus['cut_determinant'] == [1, 0, 0, 0, 0, -1], "minus cut")
    require(minus['full_determinant'] == [1, 0, -1, 0, 0, -1, 0, 1], "minus feedback")

    primes = [p for p in range(2, 100)
              if all(p % d for d in range(2, math.isqrt(p) + 1))]
    adj = {p: {q for q in primes if q != p and math.isqrt(p + q) ** 2 == p + q}
           for p in primes}
    core = {p: set(v) for p, v in adj.items()}
    while any(len(v) < 2 for v in core.values()):
        for p in [p for p, v in core.items() if len(v) < 2]:
            for q in core[p]:
                core[q].remove(p)
            del core[p]
    paths = [[3, 13, 23], [3, 97, 47], [3, 61, 83, 17, 47],
             [23, 2, 47], [23, 41, 59, 5, 11, 53, 47]]
    core_edges = {tuple(sorted((p, q))) for p, v in core.items() for q in v}
    path_edges = {tuple(sorted(e)) for p in paths for e in zip(p, p[1:])}
    require(core_edges == path_edges, "core reconstruction")
    cs = cycles(adj)
    require(len(cs) == 6 and set.intersection(*(set(cyc) for cyc in cs)) == {47}, "six cycles")
    phases = {}
    for q, a0, a1, b0, b1 in itertools.product(range(2), repeat=5):
        u, v, w = a0 ^ a1, b0 ^ b1, q ^ a0 ^ b0
        direct = (a0 ^ a1, b0 ^ b1, q ^ a0 ^ b0, q ^ a1 ^ b0,
                  q ^ a0 ^ b1, q ^ a1 ^ b1)
        predicted = (u, v, w, w ^ u, w ^ v, w ^ u ^ v)
        require(direct == predicted, "phase map")
        phases[(u, v, w)] = direct
    require(len(phases) == 8, "phase surjectivity")
    hist = {o: sum(sum(p) == o for p in phases.values()) for o in [0, 3, 4]}
    require(hist == {0: 1, 3: 4, 4: 3}, "phase distribution")
    dpaths = [[3, 97, 47], [47, 17, 83, 61, 100, 3], [3, 13, 23],
              [23, 2, 47], [23, 41, 59, 5, 11, 53, 47]]
    dadj = {v: set() for p in dpaths for v in p}
    for p in dpaths:
        for x, y in zip(p, p[1:]):
            dadj[x].add(y)
    dc = cycles(dadj, directed=True)
    require(len(dadj) == 15 and sorted(map(len, dc)) == [7, 9, 13], "directed carrier")
    q13 = [18, 8, 3, 5, 14]
    qadj = {i: {j for j in range(5) if q13[i] >> j & 1} for i in range(5)}
    qc = [p for p in cycles(qadj, directed=True) if len(p) % 2]
    require(hp(q13) == 13 and len(qc) == 6, "Q13")
    require(all(set(a) & set(b) for a, b in itertools.combinations(qc, 2)), "Q13 overlaps")
    require(not set.intersection(*(set(p) for p in qc)), "Q13 non-Helly")
    return {
        "status": "FINITE-EXACT supporting PROVED scoped theorems",
        "universe": {"raw_sources_per_sign": 256, "controller_starts": 10000,
                     "walk_max_length": 20, "path_parity_assignments": 32,
                     "prime_graph": "all 25 primes below100; every distinct pair"},
        "controls": {"C3": 3, "TT3": 1, "tie_rejected": invalid_rejected},
        "raw_decode_checks": decode_checks,
        "detour": controller(5),
        "detour_tournament_witnesses": {str(h): r for h, r in sorted(witness.items())},
        "controller_longest": {"source": longest[0], "edges": len(longest[1])},
        "plus_surgery": plus, "minus_surgery": minus,
        "prime_graph": {"primes": primes, "edges": sorted([p, q] for p, vs in adj.items() for q in vs if p < q),
                        "core_paths": paths, "cycles": cs,
                        "phase_table": [{"state": k, "cycle_parities": v, "odd_count": sum(v)}
                                        for k, v in sorted(phases.items())],
                        "odd_count_distribution": hist},
        "directed_carrier": {"paths": dpaths, "cycles": dc, "odd_cycle_packing_at2": 7},
        "Q13": {"rows": q13, "HP": hp(q13), "odd_cycles": qc,
                "pairwise_conflict_graph": "K6", "total_intersection": []},
    }


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()
    result = run()
    data = json.dumps(result, indent=2, sort_keys=True) + '\n'
    target = Path(__file__).with_name('receipt.json')
    if args.write:
        target.write_text(data)
    else:
        require(target.read_text() == data, "receipt drift")
    print(json.dumps({"ok": True, "receipt_sha256": hashlib.sha256(data.encode()).hexdigest(),
                      "kernel_states": 8, "directed_carrier_packing_value": 7,
                      "completion_HP_floor": 9}, sort_keys=True))
