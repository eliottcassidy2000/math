#!/usr/bin/env python3
"""Fixed-horizon tensor interpolation; common finite arithmetic host."""
import hashlib
import itertools
import json
from verify import birth, compile_source, decode, h_dp, require, square, square_bank


def weights(mask):
    vertices = [i for i in range(mask.bit_length()) if mask >> i & 1]
    t = len(vertices)
    row = {0: (t - 1) * (t - 2) // 2}
    row.update({1 << i: 2 - t for i in vertices})
    row.update({(1 << i) | (1 << j): 1 for i, j in itertools.combinations(vertices, 2)})
    return [(m, c) for m, c in row.items() if c]


def replay(source, word):
    rows = source
    for mask in word:
        rows = birth(rows, mask)
    return rows


def checksum(values):
    return hashlib.sha256(json.dumps(values, separators=(",", ":")).encode()).hexdigest()


def tensor_test(source, horizon):
    domains = [range(1 << (len(source) + j)) for j in range(horizon)]
    sparse_domains = [[s for s in d if s.bit_count() <= 2] for d in domains]
    bank = {w: h_dp(replay(source, w)) for w in itertools.product(*sparse_domains)}
    values = []
    for word in itertools.product(*domains):
        total = 0
        for choices in itertools.product(*(weights(s) for s in word)):
            masks, cs = zip(*choices)
            coeff = 1
            for c in cs:
                coeff *= c
            total += coeff * bank[masks]
        direct = h_dp(replay(source, word))
        require(total == direct, ("tensor failure", source, word, total, direct))
        values.append(direct)
    return {"source": source, "horizon": horizon, "sparse_programs": len(bank),
            "all_programs": len(values), "counts_sha256": checksum(values)}


def common_host(source, horizon):
    pool, unused_clock = compile_source(source)
    sizes, charts, scales = [len(pool)], [], []
    for j in range(horizon):
        pool, candidates, scale = square_bank(pool)
        chart = {}
        for mask, label in candidates.items():
            chart[mask] = len(pool)
            pool.append(label)
        charts.append(chart)
        scales.append(scale)
        sizes.append(len(pool))
    witnesses = []
    domains = [[s for s in range(1 << (len(source) + j)) if s.bit_count() <= 2]
               for j in range(horizon)]
    for word in itertools.product(*domains):
        selected = list(range(len(source)))
        for j, mask in enumerate(word):
            global_mask = sum(1 << selected[i] for i in range(len(selected)) if mask >> i & 1)
            selected.append(charts[j][global_mask])
        labels = [pool[i] for i in selected]
        rows = replay(source, word)
        require(decode(labels) == rows, ("host incidence failure", word))
        witnesses.append({"word": word, "selected_pool_indices": selected, "H": h_dp(rows)})
    require(sizes == [3, 10, 66], "unexpected two-stage host sizes")
    return {"source": source, "sizes": sizes, "scales": [str(s) for s in scales],
            "labels": [str(v) for v in pool], "witnesses": witnesses}


def main():
    tests = [tensor_test((2, 4, 1), 2), tensor_test((4, 1, 2), 2), tensor_test((2, 4, 1), 3)]
    hostile = [h_dp(replay((2, 0), (s, 1))) for s in range(4)]
    require(hostile == [3, 1, 5, 1], "terminal privilege hostile")
    singleton = [h_dp(replay((0,), (s, 1))) for s in range(2)]
    require(singleton == [3, 1], "one-coordinate hostile")
    print(json.dumps({"status": "FINITE-EXACT; see HORIZON.md for all-size proof",
                      "tensor_tests": tests, "common_host": common_host((2, 4, 1), 2),
                      "nonterminal_pair_hostile": hostile, "nonterminal_single_hostile": singleton},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
