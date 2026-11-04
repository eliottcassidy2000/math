#!/usr/bin/env python3
"""All labelled seven-vertex possible unit3-cores, with explicit necessary filters."""
from itertools import combinations, permutations
import hashlib
import json
from pathlib import Path

PAIRS = list(combinations(range(7), 2))
INDEX = {e: i for i, e in enumerate(PAIRS)}
PERMS = list(permutations(range(7)))


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def adjacency(mask):
    A = [0] * 7
    for j, (u, v) in enumerate(PAIRS):
        if mask >> j & 1:
            A[u] |= 1 << v
            A[v] |= 1 << u
    return A


def passes(A):
    if min(a.bit_count() for a in A) < 3:
        return False
    if any((A[u] & A[v]).bit_count() > 2 for u, v in PAIRS):
        return False
    for v in range(7):
        left = A[v]
        while left:
            start = left & -left
            seen = start
            wave = start
            while wave:
                bit = wave & -wave
                wave -= bit
                u = bit.bit_length() - 1
                new = A[u] & A[v] & ~seen
                wave |= new
                seen |= new
            # Already maxdegree2 by the pairwise common-neighbour test.
            if all((A[u] & seen).bit_count() == 2 for u in range(7) if seen >> u & 1):
                if seen.bit_count() != 6:
                    return False
            left &= ~seen
    return True


def orbit(mask):
    E = [e for j, e in enumerate(PAIRS) if mask >> j & 1]
    result = set()
    for p in PERMS:
        result.add(sum(1 << INDEX[tuple(sorted((p[u], p[v])))] for u, v in E))
    return result


def count_colourings(A, q):
    colours = [-1] * 7
    order = sorted(range(7), key=lambda v: -A[v].bit_count())
    def rec(i):
        if i == 7:
            return 1
        u = order[i]
        bad = {colours[v] for v in range(7) if A[u] >> v & 1 and colours[v] >= 0}
        answer = 0
        for c in range(q):
            if c not in bad:
                colours[u] = c
                answer += rec(i + 1)
        colours[u] = -1
        return answer
    return rec(0)


def main():
    source = Path(__file__).read_bytes()
    A = [0] * 7
    degrees = [0] * 7
    old = 0
    labels = {}
    canonical = {}
    survived = []
    counts = dict(all_graphs=1 << 21, minimum_degree3=0, unit_local_survivors=0)
    for i in range(1 << 21):
        mask = i ^ (i >> 1)
        if i:
            bit = mask ^ old
            u, v = PAIRS[bit.bit_length() - 1]
            change = 1 if mask & bit else -1
            A[u] ^= 1 << v
            A[v] ^= 1 << u
            degrees[u] += change
            degrees[v] += change
        old = mask
        if min(degrees) < 3:
            continue
        counts["minimum_degree3"] += 1
        if not passes(A):
            continue
        counts["unit_local_survivors"] += 1
        survived.append(mask)
        if mask not in labels:
            O = orbit(mask)
            representative = min(O)
            labels.update({m: representative for m in O})
            B = adjacency(representative)
            canonical[representative] = dict(mask=representative,
                    edges=[e for j, e in enumerate(PAIRS) if representative >> j & 1],
                    degrees=sorted(a.bit_count() for a in B), labelled_orbit=len(O),
                    colourings3=count_colourings(B, 3), colourings4=count_colourings(B, 4))
    need(sum(c["labelled_orbit"] for c in canonical.values()) == len(survived) == len(labels), "orbit coverage")
    need(set(labels) == set(survived), "complete labelled coverage")
    counts["unlabelled_survivors"] = len(canonical)
    counts["survivor_mask_sha256"] = hashlib.sha256(
        b"".join(x.to_bytes(3, "little") for x in sorted(survived))).hexdigest()
    receipt = dict(universe=counts, representatives=[canonical[m] for m in sorted(canonical)],
                   source_sha256=hashlib.sha256(source).hexdigest())
    need(Path(__file__).read_bytes() == source, "source changed")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
