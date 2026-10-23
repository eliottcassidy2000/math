"""Independent literal-tuple audit of the shared-clock compiler."""
from itertools import product
from math import isqrt
from pathlib import Path
from random import Random
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("clock_compiler", HERE / "check_shared_clocks.py")
source = importlib.util.module_from_spec(spec)
spec.loader.exec_module(source)
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def literal_paths(sums, N):
    # Deliberately enumerate tuples, not source affine constants or reflection recursion.
    return [p for p in product(range(1, N+1), repeat=len(sums)+1)
            if len(set(p)) == len(p)
            and all(p[i]+p[i+1] == s for i, s in enumerate(sums))]


def decoded(record):
    domains = []
    for c in record["components"]:
        candidates = [] if c["empty"] else [r for r in range(c["L"], c["U"]+1)
                       if r not in c["holes"] and (c["pin"] is None or r == c["pin"])]
        need(len(candidates) == c["count"], "component count matches literal interval")
        need(c["empty"] == (not candidates), "component empty flag exact")
        domains.append(candidates)
    answer = set()
    for roots in product(*domains):
        if any(f["left_sign"]*roots[f["left_component"]]
               - f["right_sign"]*roots[f["right_component"]] == f["rhs"]
               for f in record["cross_forbidden"]):
            continue
        answer.add(tuple(s*roots[c]+b for s, c, b in
                         zip(record["node_sign"], record["node_component"], record["node_offset"])))
    return answer


def main():
    rng = Random(42319)
    pools = {}
    for N in range(4, 9):
        for k in range(3):
            pools[N, k] = [p for p in product(range(1, N+1), repeat=k+1)
                           if len(set(p)) == len(p)
                           and all(isqrt(a+b)**2 == a+b for a, b in zip(p, p[1:]))]
    positive = with_loops = with_marks = with_distinct = cross = 0
    trials = 1000
    for test in range(trials):
        paths = []
        seeds = []
        for _ in range(rng.randrange(1, 5)):
            N = rng.randrange(4, 9)
            possible = [k for k in range(3) if pools[N, k]]
            seed = rng.choice(pools[N, rng.choice(possible)])
            paths.append(dict(sums=[a+b for a, b in zip(seed, seed[1:])], N=N))
            seeds.append(seed)
        owners = [(a, i) for a, p in enumerate(paths) for i in range(len(p["sums"])+1)]
        ports = [] if test % 10 == 0 else [tuple(rng.choices(owners, k=2))
                                          for _ in range(rng.randrange(5))]
        distinct = [] if test % 10 == 0 else [tuple(rng.choices(owners, k=2))
                                             for _ in range(rng.randrange(4))]
        marks = []
        if test % 10:
            for _ in range(rng.randrange(3)):
                a, i = rng.choice(owners)
                value = seeds[a][i] if rng.randrange(4) else rng.choice((0, paths[a]["N"]+1))
                marks.append(((a, i), value))
        local = [literal_paths(p["sums"], p["N"]) for p in paths]
        wanted = set()
        for ps in product(*local):
            if not all(ps[a][i] == ps[b][j] for (a, i), (b, j) in ports):
                continue
            if not all(ps[a][i] != ps[b][j] for (a, i), (b, j) in distinct):
                continue
            if not all(ps[a][i] == value for (a, i), value in marks):
                continue
            wanted.add(tuple(p[0] for p in ps))
        record = source.compile_network(paths, ports, marks, distinct)
        need(decoded(record) == wanted, "literal tuple network comparison")
        positive += bool(wanted)
        with_loops += any(a == b for (a, i), (b, j) in ports)
        with_marks += bool(marks)
        with_distinct += bool(distinct)
        cross += bool(record["cross_forbidden"])
    local_empty = source.local_domain([4], 2)
    need(local_empty["empty"] and not literal_paths([4], 2), "half-diagonal singleton domain is empty")
    out = dict(status="PASS", seed=42319, networks=trials, positive=positive,
               with_loop_equalities=with_loops, with_absolute_marks=with_marks,
               with_distinctness=with_distinct, with_cross_component_factors=cross,
               scope="1..4 actual square-clock paths, length0..2, hosts4..8; literal tuple enumeration",
               checks=CHECKS)
    text = json.dumps(out, indent=2, sort_keys=True)+"\n"
    target = HERE / "shared_clocks_independent.json"
    if "--write" in sys.argv:
        target.write_text(text, encoding="utf-8")
    else:
        need(json.loads(target.read_text(encoding="utf-8")) == json.loads(text), "frozen independent receipt")
    print(json.dumps(dict(status="PASS", checks=CHECKS, networks=trials, positive=positive)))


if __name__ == "__main__":
    main()
