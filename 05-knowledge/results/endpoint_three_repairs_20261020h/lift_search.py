"""Small symbolic block-word probe, not a tournament census."""
from itertools import combinations, permutations
import json
from pathlib import Path


def edges(word):
    return {tuple(sorted(pair)) for pair in zip(word, word[1:])}


def val(fixed, u, v):
    pair = tuple(sorted((u, v)))
    return None if pair not in fixed else (fixed[pair] if u < v else 1 - fixed[pair])


def move(fixed, word):
    retained = edges(word)
    return {pair: value if pair in retained else 1 - value for pair, value in fixed.items()}


def legal(fixed, word):
    return all(val(fixed, u, v) == 1 for u, v in zip(word, word[1:]))


def black(fixed, word, standard_anti=False):
    n = len(word); rho = dict(zip(word, word[::-1]))
    for u, v in combinations(range(n), 2):
        a, b = val(fixed, u, v), val(fixed, rho[v], rho[u])
        if a is not None and b is not None and a != b:
            return dict(type="fixed_arc", pair=[u, v, rho[v], rho[u]])
    lower = [sum(val(fixed, u, v) == 1 for v in range(n) if u != v) for u in (word[0], word[-1])]
    upper = [n - 1 - sum(val(fixed, u, v) == 0 for v in range(n) if u != v) for u in (word[0], word[-1])]
    if sum(lower) > n - 1 or sum(upper) < n - 1:
        return dict(type="degree_bound", lower=lower, upper=upper)
    if standard_anti:
        g = {v: rho[n - 1 - v] for v in range(n)}; unseen = set(range(n))
        while unseen:
            u = min(unseen); orbit = []; v = u
            while v not in orbit:
                orbit.append(v); unseen.remove(v); v = g[v]
            if len(orbit) % 2 == 0:
                return dict(type="even_order", orbit=orbit)
            for a in (0, n - 1):
                if a in orbit:
                    out = sum(val(fixed, a, v) == 1 for v in orbit if a != v)
                    if 2 * out != len(orbit) - 1:
                        return dict(type="orbit_capacity", orbit=orbit, vertex=a, out=out)
    return None


def xor_words(target_edges, first, second, n):
    graph = {v: [] for v in range(n)}
    for u, v in target_edges ^ edges(first) ^ edges(second):
        graph[u].append(v); graph[v].append(u)
    ends = [v for v in range(n) if len(graph[v]) == 1]
    if len(ends) != 2 or any(len(graph[v]) not in (1, 2) for v in range(n)):
        return []
    word = [ends[0]]; prev = None
    while True:
        choices = [v for v in graph[word[-1]] if v != prev]
        if not choices:
            break
        prev, nxt = word[-1], choices[0]; word.append(nxt)
        if len(word) > n:
            return []
    return [tuple(word), tuple(word[::-1])] if len(word) == n else []


def search(n, x_to_e, b_to_e, y_to_x):
    h, c, z, e = n - 3, n - 2, n - 1, n - 4
    blocks = [(0,), (1,), (2,), tuple(range(3, n - 3)), (h,), (c,), (z,)]
    fixed = {(i, i + 1): 1 for i in range(n - 1)}
    fixed.update({(0, j): int(j not in (h, c, z)) for j in range(1, n)})
    fixed.update({(1, c): 0, (2, c): 0, (2, e): x_to_e, (1, e): b_to_e, (2, h): 1 - y_to_x})
    for (u, v), bit in list(fixed.items()):
        pair = (n - 1 - v, n - 1 - u)
        if pair in fixed and fixed[pair] != bit:
            return dict(incompatible=True)
        fixed[pair] = bit
    words = [tuple(v for i in perm for v in blocks[i]) for perm in permutations(range(7))]
    P = tuple(range(n)); B = move(fixed, P); EP = edges(P)
    firsts = [(w, black(fixed, w, True)) for w in words if legal(fixed, w)]
    firsts = [(w, witness) for w, witness in firsts if witness]
    for Q1, witness1 in firsts:
        S = move(fixed, Q1)
        for Q2 in words:
            if not legal(S, Q2):
                continue
            for Q3 in xor_words(EP, Q1, Q2, n):
                if not legal(B, Q3):
                    continue
                witness2 = black(S, Q2); witness3 = black(B, Q3, True)
                if witness2 and witness3:
                    return dict(words=[Q1, Q2, Q3], witnesses=[witness1, witness2, witness3], first_words=len(firsts))
    return dict(words=None, first_words=len(firsts), candidate_words=len(words))


def main():
    result = []
    for x_to_e, b_to_e, y_to_x in [(x, b, y) for x in (0, 1) for b in (0, 1) for y in (0, 1)]:
        item = dict(n=12, x_to_e=x_to_e, b_to_e=b_to_e, y_to_x=y_to_x,
                    result=search(12, x_to_e, b_to_e, y_to_x))
        result.append(item); print(json.dumps(item))
    Path(__file__).with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
