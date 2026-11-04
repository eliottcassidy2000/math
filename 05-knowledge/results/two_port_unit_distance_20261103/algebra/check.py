"""Permutation-based exact controls for the two-port algebraic obstructions."""
from argparse import ArgumentParser
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import json

SOURCE = Path(__file__)
SOURCE_HASH = sha256(SOURCE.read_bytes()).hexdigest()
CHECKS = 0


def require(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(why)


def make(n, code):
    out = [0] * n
    for bit, (u, v) in enumerate(combinations(range(n), 2)):
        if code >> bit & 1:
            out[u] |= 1 << v
        else:
            out[v] |= 1 << u
    return tuple(out)


def data(out):
    n = len(out)
    H, start, end = 0, [0] * n, [0] * n
    C = [[0] * n for _ in range(n)]
    for word in permutations(range(n)):
        legal = [bool(out[u] >> v & 1) for u, v in zip(word, word[1:])]
        bad = sum(not x for x in legal)
        if bad == 0:
            H += 1
            start[word[0]] += 1
            end[word[-1]] += 1
        for i, (u, v) in enumerate(zip(word, word[1:])):
            if bad == (not legal[i]):
                C[u][v] += 1
    return H, tuple(start), tuple(end), tuple(map(tuple, C))


def subset_h(out):
    n = len(out)
    dp = [[0] * n for _ in range(1 << n)]
    result = [0] * (1 << n)
    result[0] = 1
    for mask in range(1, 1 << n):
        for u in range(n):
            if mask >> u & 1:
                rest = mask ^ (1 << u)
                dp[mask][u] = sum(dp[rest][v] for v in range(n) if out[v] >> u & 1) if rest else 1
        result[mask] = sum(dp[mask])
    return result


def triangle_word_counts(out, triple):
    edges = {frozenset(pair) for pair in combinations(triple, 2)}
    coeff = [0, 0]
    for word in permutations(range(len(out))):
        signs = []
        legal = True
        for u, v in zip(word, word[1:]):
            if frozenset((u, v)) in edges:
                signs.append(bool(out[u] >> v & 1))
            elif not out[u] >> v & 1:
                legal = False
                break
        if legal and len(signs) == 2:
            require(signs[0] == signs[1], "cyclic triangle's two edges are coherent")
            coeff[0 if signs[0] else 1] += 1
    return coeff


def at_half(out):
    # Every weighted orientation equals1/2 at the global converse midpoint.
    n = len(out)
    from fractions import Fraction
    from math import factorial
    return Fraction(factorial(n), 2 ** (n - 1))


def global_converse_polynomial(out):
    """Power coefficients of P((1-t)A+tA^t), by literal words."""
    from math import comb
    n = len(out)
    coeff = [0] * n
    for word in permutations(range(n)):
        forward = sum(bool(out[u] >> v & 1) for u, v in zip(word, word[1:]))
        backward = n - 1 - forward
        for j in range(forward + 1):
            coeff[backward + j] += (-1) ** j * comb(forward, j)
    return coeff


def run():
    census = []
    total_triangle_tests = 0
    for n in range(1, 6):
        edges = list(combinations(range(n), 2))
        codes = range(1 << len(edges))
        tables = [data(make(n, code)) for code in codes]
        symmetric = 0
        odd_separations = 0
        triangle_tests = 0
        for code in codes:
            out = make(n, code)
            H, S, E, C = tables[code]
            h = subset_h(out)
            require(H == h[-1] and H % 2 == 1, "independent path count and parity")
            require(sum(C[u][v] for u, v in edges if out[u] >> v & 1)
                    + sum(C[v][u] for u, v in edges if out[v] >> u & 1)
                    == (n - 1) * H, "Euler identity")
            cover2 = sum(map(sum, C))
            Q2 = sum(h[mask] * h[(1 << n) - 1 - mask] for mask in range(1 << n))
            require(Q2 == 2 * H + cover2, "two-colour/cover dictionary")
            if n >= 2:
                require(Q2 % 4 == 0, "complement-pair divisibility")
            symmetric += all(C[u][v] == C[v][u] for u, v in edges)
            if n >= 3 and n % 2:
                require(any(C[u][v] != C[v][u] for u, v in edges), "odd converse separated")
                odd_separations += 1
            reversed_C = tables[((1 << len(edges)) - 1) ^ code][3]
            require(all(reversed_C[u][v] == C[v][u] for u in range(n) for v in range(n)),
                    "converse transposes C")
            for triple in combinations(range(n), 3):
                if any(sum(out[u] >> v & 1 for v in triple if v != u) != 1 for u in triple):
                    continue
                bits = [i for i, edge in enumerate(edges) if set(edge) <= set(triple)]
                flipped = code ^ sum(1 << i for i in bits)
                C1 = tables[flipped][3]
                oriented = [(u, v) if out[u] >> v & 1 else (v, u)
                            for u, v in combinations(triple, 2)]
                delta = sum((C1[v][u] - C1[u][v]) - (C[v][u] - C[u][v])
                            for u, v in oriented)
                outside = ((1 << n) - 1) ^ sum(1 << u for u in triple)
                A, B = triangle_word_counts(out, triple)
                require(A >= 3 * h[outside] and B >= 3 * h[outside], "three block insertions")
                require(delta == 2 * (A + B) >= 12 * h[outside], "strict directional gradient change")
                require(C != C1, "triangle flip separates gradient")
                triangle_tests += 1
        total_triangle_tests += triangle_tests
        census.append({"n": n, "labelled": len(tables), "symmetric_C": symmetric,
                       "odd_converse_separations": odd_separations,
                       "triangle_flips": triangle_tests})
    transitive3 = (6, 4, 0)
    cycle3 = (2, 4, 1)
    cycle_sink = (10, 12, 9, 0)
    require(global_converse_polynomial(transitive3) == [1, 2, -2], "nonconvex hostile")
    require(global_converse_polynomial(cycle3) == [3, -6, 6], "nonconcave hostile")
    require(global_converse_polynomial(cycle_sink) == [3, 0, 0, 0], "flat converse direction")
    require(data(cycle_sink)[3] != tuple(zip(*data(cycle_sink)[3])), "flat scalar direction still has different C")
    require(sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_HASH, "source stable")
    return {"status": "FINITE-EXACT PASS", "source_sha256": SOURCE_HASH,
            "checks": CHECKS, "census": census, "directed_triangle_flip_tests": total_triangle_tests,
            "hostile_converse_polynomials": [[1, 2, -2], [3, -6, 6], [3, 0, 0, 0]],
            "scope": "All1099 labelled tournaments n=1..5, every cyclic triple, complete permutation-derived H/S/E/C. Universal statements have separate proofs."}


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    text = json.dumps(run(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8", newline="\n")
    else:
        print(text, end="")
