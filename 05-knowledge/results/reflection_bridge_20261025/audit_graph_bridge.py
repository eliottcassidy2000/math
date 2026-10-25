"""Independent small native audit of the exponent-three threshold map."""
from itertools import product, combinations
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


D = [x for x in product(range(4), range(2), range(2))
     if x not in [(0, 0, 0), (3, 1, 1)]]
E = list(combinations(range(4), 2))


def phi(x):
    a, b, c = x
    predicates = [a+2*b >= 3, b == 1, a >= 2,
                  a+b+c >= 3, c == 1, a+2*c >= 3]
    return {e for e, yes in zip(E, predicates) if yes}


def arc(x, i, j):
    if i > j:
        return not arc(x, j, i)
    return j == i+1 or (i, j-1) not in phi(x)


def reach(x, start):
    seen, todo = {start}, [start]
    while todo:
        i = todo.pop()
        for j in range(5):
            if j != i and arc(x, i, j) and j not in seen:
                seen.add(j)
                todo.append(j)
    return seen


def main():
    for x in D:
        a, b, c = x
        g = phi(x)
        cx, qx = (3-a, 1-b, 1-c), (a, c, b)
        require(phi(cx) == set(E)-g, "complement")
        require(phi(qx) == {tuple(sorted((3-i, 3-j))) for i, j in g}, "reflection")
        require(arc(x, 4, 0) == (a >= 2), "marked closure")
        strong = all(len(reach(x, i)) == 5 for i in range(5))
        require(strong == (a >= 2 or b == c == 1), "strong")
        blue = all(arc(x, i, j) == arc(x, 4-j, 4-i)
                   for i, j in combinations(range(5), 2))
        require(blue == (b == c), "marked blue")
        require(int((0, 2) in g)-int((1, 3) in g) == b-c, "g coordinate")
        cycle = any((arc(x, i, j) and arc(x, j, ell) and arc(x, ell, i))
                    or (arc(x, j, i) and arc(x, ell, j) and arc(x, i, ell))
                    for i, j, ell in combinations(range(5), 3))
        require((not cycle) == (x == (1, 0, 0)), "transitive")
        for y in D:
            if all(i <= j for i, j in zip(x, y)):
                require(g <= phi(y), "monotone")
    require(len({frozenset(phi(x)) for x in D}) == 14, "injection")
    require(phi((1, 0, 0)) <= phi((0, 1, 0)), "extra-comparison hostile")
    print(json.dumps(dict(status="FINITE-EXACT", sources=14, ordered_source_pairs=196,
                          literal_tournament_triples=140, tournament_predicates=5,
                          verdict="pass"), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
