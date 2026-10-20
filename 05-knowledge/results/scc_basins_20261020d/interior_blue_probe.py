"""Exact affine tests for blue intermediate moves in the interior R4 ledger."""
import json
from itertools import combinations
from pathlib import Path
from moves import join, edge_set, need

B = (6, 4, 8, 3)


def problem(m, k, step):
    n = m + 4 + k
    aa = tuple(range(m))
    bb = tuple(range(m, m+4))
    ww = tuple(range(m+4, n))
    words = (aa+(bb[2], bb[3], bb[0], bb[1])+ww,
             ww+(bb[3], bb[0])+aa+(bb[2], bb[1]),
             ww+(bb[3], bb[1], bb[0])+aa+(bb[2],),
             aa+(bb[2], bb[0], bb[3], bb[1])+ww)
    transitive = lambda h: tuple(sum(1 << j for j in range(i+1, h)) for i in range(h))
    base = join([transitive(m), B, transitive(k)])
    free = [(u, v) for block in (aa, ww) for u, v in combinations(block, 2)
            if v != u+1]
    index = {e: i for i, e in enumerate(free)}
    constant = 1 << len(free)
    expressions = {(u, v): ((1 << index[u, v]) if (u, v) in index else 0)
                   ^ (constant if base[u] >> v & 1 else 0)
                   for u, v in combinations(range(n), 2)}
    for path in words[:step]:
        kept = edge_set(path)
        expressions = {e: value ^ (0 if e in kept else constant)
                       for e, value in expressions.items()}
    def arc(u, v):
        return expressions[min(u, v), max(u, v)] ^ (constant if u > v else 0)
    path = words[step]
    equations = [arc(u, v) ^ arc(path[n-1-j], path[n-1-i])
                 for i, u in enumerate(path) for j, v in enumerate(path) if i < j]
    pivots = {}
    for row in equations:
        while row and row != constant:
            pivot = (row & -row).bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
        if row == constant:
            return None
    bits = 0
    for pivot in sorted(pivots, reverse=True):
        row = pivots[pivot]
        if ((row & bits).bit_count() & 1) ^ bool(row & constant):
            bits |= 1 << pivot
    out = list(base)
    for i, (u, v) in enumerate(free):
        if bits >> i & 1:
            out[u] ^= 1 << v
            out[v] ^= 1 << u
    from moves import blue, complement
    actual = tuple(out)
    for p in words[:step]:
        actual = complement(actual, p)
    need(blue(actual, words[step]), "affine witness is not blue")
    return dict(m=m, k=k, step=step+1, free=len(free)-len(pivots),
                source=out, paths=words)


def main():
    rows = []
    for m in range(1, 7):
        for k in range(1, 15):
            for step in (1, 2):
                result = problem(m, k, step)
                if result is not None:
                    rows.append(result)
    Path(__file__).with_suffix(".json").write_text(json.dumps(rows, indent=2)+"\n")
    print(json.dumps([{a: r[a] for a in ("m", "k", "step", "free")} for r in rows],
                     indent=2))


if __name__ == "__main__":
    main()
