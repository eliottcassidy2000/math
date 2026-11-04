#!/usr/bin/env python3
"""Standalone exact controls for PROOF.md. Uses only the standard library."""
from collections import Counter, defaultdict
from itertools import permutations, product
from math import comb, factorial
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def check_tournament(rows):
    n = len(rows)
    for i in range(n):
        require(rows[i] >> i & 1 == 0, "loop")
        require(rows[i] >> n == 0, "out-of-range arc")
        for j in range(i):
            require(((rows[i] >> j) ^ (rows[j] >> i)) & 1 == 1,
                    "not a tournament")


def canonical(rows):
    cells = defaultdict(list)
    for v, row in enumerate(rows):
        cells[row.bit_count()].append(v)
    blocks = [list(permutations(cells[d])) for d in sorted(cells)]
    best = None
    best_rows = None
    for pblocks in product(*blocks):
        p = sum(pblocks, ())
        candidate = tuple(sum(((rows[p[i]] >> p[j]) & 1) << j
                              for j in range(len(rows)))
                          for i in range(len(rows)))
        if best is None or candidate < best:
            best, best_rows = candidate, candidate
    return best_rows


def frame_census(n):
    pairs = [(i, j) for i in range(n) for j in range(i+2, n)]
    seen = set()
    for mask in range(1 << len(pairs)):
        rows = [0] * n
        for i in range(n-1):
            rows[i] |= 1 << (i+1)
        for k, (i, j) in enumerate(pairs):
            if mask >> k & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
        seen.add(canonical(rows))
    return sorted(seen)


def strong(rows):
    n = len(rows)
    def closure(reverse):
        seen = 1
        changed = True
        while changed:
            old = seen
            for i in range(n):
                if seen >> i & 1:
                    if not reverse:
                        seen |= rows[i]
                    else:
                        seen |= sum(((rows[j] >> i) & 1) << j for j in range(n))
            changed = seen != old
        return seen == (1 << n)-1
    return closure(False) and closure(True)


def literal_count(rows, vertices=None):
    if vertices is None:
        vertices = tuple(range(len(rows)))
    return sum(all(rows[a] >> b & 1 for a, b in zip(p, p[1:]))
               for p in permutations(vertices))


def subset_counts(rows):
    """Counts all induced subsets by a last-vertex dynamic program."""
    n = len(rows)
    limit = 1 << n
    dp = [[0] * n for _ in range(limit)]
    totals = [0] * limit
    totals[0] = 1
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, limit):
        current = dp[mask]
        totals[mask] = sum(current)
        unused = (limit-1) ^ mask
        for v, count in enumerate(current):
            if not count:
                continue
            successors = rows[v] & unused
            while successors:
                bit = successors & -successors
                successors ^= bit
                w = bit.bit_length()-1
                dp[mask | bit][w] += count
    return totals


def ordinal_sum(parts):
    n = sum(map(len, parts))
    rows = []
    offset = 0
    for part in parts:
        later = ((1 << n)-1) ^ ((1 << (offset+len(part)))-1)
        rows.extend((row << offset) | later for row in part)
        offset += len(part)
    return tuple(rows)


def automorphism_count(rows):
    return sum(all(((rows[i] >> j) & 1) == ((rows[p[i]] >> p[j]) & 1)
                   for i in range(len(rows)) for j in range(len(rows)))
               for p in permutations(range(len(rows))))


def main():
    classes = {}
    profiles = {}
    literal_checks = 0
    expected_classes = [1, 1, 2, 4, 12, 56]
    for n, expected in enumerate(expected_classes, 1):
        reps = frame_census(n)
        require(len(reps) == expected, f"class count at n={n}")
        classes[n] = reps
        profiles[n] = []
        mass = 0
        for rows in reps:
            check_tournament(rows)
            counts = subset_counts(rows)
            full = (1 << n)-1
            h = literal_count(rows)
            require(h == counts[full], "literal/DP mismatch")
            literal_checks += 1
            deletes = []
            for v in range(n):
                hv = literal_count(rows, tuple(w for w in range(n) if w != v))
                literal_checks += 1
                require(hv == counts[full ^ (1 << v)], "deletion mismatch")
                deletes.append(hv)
            aut = automorphism_count(rows)
            mass += factorial(n) // aut
            if strong(rows) and n >= 3:
                profiles[n].append({"rows": rows, "h": h,
                                    "deletion_max": max(deletes), "aut": aut})
        require(mass == 1 << (n*(n-1)//2), "labelled orbit-mass mismatch")

    require(len(profiles[3]) == 1, "C3 classification")
    require(len(profiles[4]) == 1, "Q classification")
    c3, q = profiles[3][0], profiles[4][0]
    require((c3['h'], c3['deletion_max'], c3['aut']) == (3, 1, 3), "C3 data")
    require((q['h'], q['deletion_max'], q['aut']) == (5, 3, 1), "Q data")
    require(min(p['h'] for p in profiles[5]) == 9, "five-vertex floor")
    require(canonical(tuple(sum((1-((q['rows'][i] >> j) & 1)) << j
                                for j in range(4) if i != j)
                            for i in range(4))) == q['rows'], "Q converse")

    # Exact floor/cost inequalities for the closed-form Busch sequence.
    floor_checks = 0
    for s in range(3, 251):
        if s % 3 == 0:
            f = 3 * 5 ** (s//3-1)
        elif s % 3 == 1:
            f = 5 ** (s//3)
        else:
            f = 9 * 5 ** (s//3-1)
        require(f**3 * 81**(s-3) >= 3**s * 125**(s-3), "cost inequality")
        if s >= 5:
            require(f**3 >= 3**(s+1), "large-component inequality")
        floor_checks += 1
    require(125**3 > 3*81**3, "top-layer excess cutoff")

    # Theorem A and congruence reduce the entire top layer to these three
    # exceptional size patterns. Exhaust every actual strong class in them.
    exceptional_words = [(p,) for p in profiles[6]]
    exceptional_words += [(p, q) for p in profiles[5]]
    exceptional_words += [(q, p) for p in profiles[5]]
    exceptional_words += [(q, q, q)]
    pattern_checks = 0
    transitions = []
    for r in range(1, 41):
        survivors = []
        for word in exceptional_words:
            s = sum(len(p['rows']) for p in word)
            require(s % 3 == 0, "nonintegral template")
            a = r-s//3
            if a < 0:
                continue
            h = 3**a
            for p in word:
                h *= p['h']
            deletion_max = max(h // p['h'] * p['deletion_max'] for p in word)
            if a:
                deletion_max = max(deletion_max, h//3)
            low = max(3**(r-1), deletion_max)
            high = min(3**r-1, h-1)
            pattern_checks += 1
            if low <= high:
                survivors.append((tuple(len(p['rows']) for p in word), low, high))
        expected = [] if r < 4 else [((4, 4, 4), 25*3**(r-3), 3**r-1)]
        require(survivors == expected, f"unexpected top template at r={r}")
        if r <= 10:
            transitions.append({"r": r, "top_order": 3*r,
                                "extra_onset": None if r < 4 else 25*3**(r-3),
                                "additional_classes": 0 if r < 4 else comb(r-1, 3)})

    actual_controls = []
    named_words = [('C3^4', (c3,)*4), ('Q^3', (q,)*3)]
    for gap in range(4):
        word = [q, q, q]
        word.insert(gap, c3)
        named_words.append((f'Q^3 with triangle gap {gap}', tuple(word)))
    for name, word in named_words:
        rows = ordinal_sum([p['rows'] for p in word])
        check_tournament(rows)
        counts = subset_counts(rows)
        full = (1 << len(rows))-1
        h = counts[full]
        dmax = max(counts[full ^ (1 << v)] for v in range(len(rows)))
        predicted = 1
        for p in word:
            predicted *= p['h']
        predicted_deletion = max(predicted // p['h'] * p['deletion_max'] for p in word)
        require((h, dmax) == (predicted, predicted_deletion), "actual sum mismatch")
        actual_controls.append({"name": name, "n": len(rows), "h": h,
                                "maximum_deletion": dmax, "rows": rows})

    # Hostile controls: Q^3 starts at75, not74; source padding preserves H
    # and therefore destroys induced-minimality.
    rows = ordinal_sum([q['rows']]*3)
    counts = subset_counts(rows)
    require(max(counts[((1 << 12)-1) ^ (1 << v)] for v in range(12)) == 75,
            "threshold hostile")
    padded = ordinal_sum([(0,), rows])
    pc = subset_counts(padded)
    require(pc[-1] == pc[((1 << 13)-1) ^ 1] == 125, "padding hostile")

    output = {
        "status": "ALL CHECKS PASSED",
        "class_counts_1_to_6": expected_classes,
        "strong_class_counts_3_to_6": {str(n): len(profiles[n]) for n in range(3, 7)},
        "literal_path_and_deletion_checks": literal_checks,
        "floor_inequality_controls": floor_checks,
        "top_template_checks": pattern_checks,
        "threshold_transitions": transitions,
        "actual_adjacency_controls": actual_controls,
        "hostile_controls": ["Q^3 not minimal at B74", "singleton padding not minimal"],
        "scope": "Finite controls for the proved all-order top-layer and near-top-template theorems."
    }
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
