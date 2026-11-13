"""Independent literal cycle and triangular-join checks; no assert dependency."""
from itertools import product
from pathlib import Path
import json
import random


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def maps(m):
    return sorted({tuple((s*x+d) % m for x in range(m))
                   for s in (-1, 1) for d in range(m)})


def normal(m):
    ps = maps(m)
    out = {frozenset()}
    out.update(frozenset((x, p[x], q[x]) for x in range(m))
               for p, q in product(ps, repeat=2))
    if m % 2:
        out.update(frozenset([a]) for a in product(range(m), repeat=3))
    elif m > 2:
        out.update(frozenset([a, tuple((x+m//2) % m for x in a)])
                   for a in product(range(m), repeat=3))
    return out


def literal_triangles(m):
    out = set()
    ps = maps(m)
    for p, q, r in product(ps, repeat=3):
        out.add(frozenset((x, p[x], q[p[x]]) for x in range(m)
                          if r[q[p[x]]] == x))
    return out


def join(a, b, c):
    by_x = {}
    for x, bv, y in b:
        by_x.setdefault(x, []).append((bv, y))
    by_zy = {}
    for z, y, cv in c:
        by_zy.setdefault((z, y), []).append(cv)
    result = {}
    for av, x, z in a:
        for bv, y in by_x.get(x, ()):
            for cv in by_zy.get((z, y), ()):
                key = av, bv, cv
                result[key] = result.get(key, 0)+1
    need(all(v == 1 for v in result.values()), "multiple interior completions")
    return frozenset(result)


def main():
    report = {"status": "FINITE-EXACT", "literal_triangles": {}, "joins": {}}
    for m in range(2, 12):
        states = normal(m)
        expected = 1+2**2 if m == 2 else 1+(m**3 if m % 2 else m**3//2)+(2*m)**2
        need(len(states) == expected, (m, "count"))
        need(literal_triangles(m) == states, (m, "literal realization"))
        report["literal_triangles"][m] = len(states)
    for m in (2, 4):
        states = normal(m)
        for triple in product(states, repeat=3):
            need(join(*triple) in states, (m, "closure"))
        report["joins"][m] = {"exhaustive": True, "triples": len(states)**3}
    rng = random.Random(20261112)
    for m in (3, 5, 7, 10):
        states = normal(m)
        ordered = sorted(states, key=lambda x: tuple(sorted(x)))
        for _ in range(5000):
            need(join(*(rng.choice(ordered) for _ in range(3))) in states, "sample closure")
        report["joins"][m] = {"exhaustive": False, "triples": 5000}
    need({(x, x, x) for x in range(4)} != {(x, x, x) for x in (0, 2)}, "hostile")
    need(not [x for x in range(4) if (2*x-1) % 4 == 0], "even reflection hostile")
    report["large_modulus_formula_only"] = {m: 1+m**3//2+4*m*m for m in (30, 42)}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, sort_keys=True))


if __name__ == '__main__':
    main()
