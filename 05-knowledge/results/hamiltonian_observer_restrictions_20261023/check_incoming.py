"""Exact carrier and representative-cut controls for the late integration."""
from itertools import combinations
from pathlib import Path
import json
import sys

checks = 0


def need(ok, why):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(why)


def graph_on(states):
    return {x: {y for y in states if sum(a != b for a, b in zip(x, y)) == 1}
            for x in states}


def connected(adj):
    reached = {next(iter(adj))}
    todo = list(reached)
    while todo:
        x = todo.pop()
        for y in adj[x] - reached:
            reached.add(y)
            todo.append(y)
    return len(reached) == len(adj)


def degrees(edges, n):
    return sorted(sum(v in e for e in edges) for v in range(n))


def main():
    reports = []
    for p in range(3, 41):
        s = {(a, b) for a in range(p) for b in range(p) if (b-a) % p in (0, 1)}
        adj = graph_on(s)
        need(len(s) == 2*p and connected(adj), "connected phase carrier")
        need(all(len(x) == 2 for x in adj.values()), "cycle degrees")
        ordered = [(a, b) for a in range(p) for b in (a, (a+1) % p)]
        need(set(ordered) == s, "all actual phase labels retained")
        need(all(ordered[(i+1) % len(ordered)] in adj[x]
                 for i, x in enumerate(ordered)), "displayed cycle")
        exact = {x for x in s if x[0]-x[1] in (0, -1)}
        need(s-exact == {(p-1, 0)}, "only wrap representative removed")
        cut = graph_on(exact)
        need(connected(cut), "cut remains connected")
        need(sorted(len(x) for x in cut.values()) == [1, 1]+[2]*(2*p-3), "cut is a path")
        for modulus in (p+1, 2*p, 3*p+2):
            modular = {x for x in s if (x[0]-x[1]) % modulus in (0, modulus-1)}
            need(modular == exact, "bounded-range modular guard")
        need({x for x in s if (x[0]-x[1]) % p in (0, p-1)} == s,
             "threshold hostile: modulus p does not cut")
        for modulus in range(2, 47):
            modular = {x for x in s if (x[0]-x[1]) % modulus in (0, modulus-1)}
            wrap = p % modulus == 0 or (p-1) % modulus == 0
            need(modular == (s if wrap else exact), "exact all-modulus alias criterion")
        # Move-coordinate labels survive the diagonal-phase quotient.
        for a in range(p):
            d, u, following = (a,a), (a,(a+1)%p), ((a+1)%p,(a+1)%p)
            need(d[0] == u[0] and d[1] != u[1] and u in adj[d], "first edge changes b")
            need(u[1] == following[1] and u[0] != following[0] and following in adj[u],
                 "second edge changes a")
        need([k for k in range(1,p+1) if k % p == 0] == [p], "two-edge circuit return clock")
        reports.append({"p": p, "cycle_vertices": len(s), "cut_vertices": len(exact)})
    three = {(a,b) for a in range(3) for b in range(3)
             if (b-a) % 3 in (0,1) and (a-b) % 8 in (0,7)}
    need(len(three) == 5 and (2,0) not in three, "incoming exact p3 mod8 control")
    # Explicit isomorphisms for equal shapes; degree multisets reject unequal ones.
    ga, gb, gc = [{(0,1),(0,2)}, {(0,1),(1,2)}, {(0,1),(2,3)}]
    expand = lambda g: {(i, j+1) for i, j in g}
    for n in range(4, 12):
        need(degrees(ga,n-1) == degrees(gb,n-1), "same compressed path shape")
        swap01 = lambda v: 1-v if v < 2 else v
        need({tuple(sorted((swap01(a),swap01(b)))) for a,b in ga} == gb,
             "actual compressed graph isomorphism")
        need(degrees(expand(ga),n) != degrees(expand(gb),n), "full graph distinguishes")
        if n >= 5:
            need(degrees(expand(ga),n) == degrees(expand(gc),n), "same full path shape")
            perm = [2,1,0,4,3]+list(range(5,n))
            need({tuple(sorted((perm[a],perm[b]))) for a,b in expand(ga)} == expand(gc),
                 "actual full graph isomorphism")
            need(degrees(ga,n-1) != degrees(gc,n-1), "compressed graph distinguishes")
    result = {"status":"FINITE-EXACT", "checks":checks, "phase_cases":reports,
              "p3_mod8_survivors":[list(x) for x in sorted(three)],
              "scope":"quotient carriers and representative cut only"}
    path = Path(__file__).with_name("incoming_receipt.json")
    if "--write" in sys.argv:
        path.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    else:
        need(json.loads(path.read_text(encoding="utf-8")) == result, "frozen receipt")
    print(f"PASS incoming carrier and phase-cut checks: {checks}")


if __name__ == "__main__":
    main()
