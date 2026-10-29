"""Recovered successor-rook recurrence with explicit directed-cycle removal."""
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial
from pathlib import Path
import json
import renewal

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def rook_polynomial(D, vertices):
    nodes = tuple(i for i in range(len(D)) if vertices >> i & 1)
    @lru_cache(None)
    def recur(row, columns):
        if row == len(nodes):
            return (1,)
        out = list(recur(row + 1, columns)) + [0]
        choices = D[nodes[row]] & columns
        while choices:
            bit = choices & -choices
            choices -= bit
            for i, c in enumerate(recur(row + 1, columns ^ bit)):
                out[i + 1] += c
        return tuple(out)
    return recur(0, vertices)


def cycle_free_rooks(D):
    n = len(D)
    cycles = []
    for size in range(3, n + 1):
        for nodes in combinations(range(n), size):
            for tail in permutations(nodes[1:]):
                word = (nodes[0],) + tail
                if all(D[word[i]] >> word[(i + 1) % size] & 1 for i in range(size)):
                    cycles.append(sum(1 << x for x in nodes))
    full = (1 << n) - 1
    out = [0] * (n + 1)
    def packs(start, used, sign):
        size = used.bit_count()
        for i, c in enumerate(rook_polynomial(D, full ^ used)):
            out[size + i] += sign * c
        for j in range(start, len(cycles)):
            if not cycles[j] & used:
                packs(j + 1, used | cycles[j], -sign)
    packs(0, 0, 1)
    return out


def literal_forests(D):
    arcs = [(i, j) for i in range(len(D)) for j in range(len(D)) if D[i] >> j & 1]
    out = [0] * (len(D) + 1)
    for mask in range(1 << len(arcs)):
        succ, incoming = {}, set()
        legal = True
        for e, (i, j) in enumerate(arcs):
            if mask >> e & 1:
                if i in succ or j in incoming:
                    legal = False
                    break
                succ[i] = j
                incoming.add(j)
        if not legal:
            continue
        for i in succ:
            seen = set()
            while i in succ:
                if i in seen:
                    legal = False
                    break
                seen.add(i)
                i = succ[i]
            if not legal:
                break
        if legal:
            out[mask.bit_count()] += 1
    return out


def main(write=True):
    global CHECKS
    CHECKS = 0
    counts = []
    for n in range(3, 6):
        totals = []
        for x in range(1 << ((n - 1) * (n - 2) // 2)):
            T = renewal.presentation(n, x)
            partner = renewal.ports.core.complement(T, tuple(range(n)))
            forbidden = renewal.ports.core.converse(partner)
            coefficients = cycle_free_rooks(forbidden)
            need(coefficients == literal_forests(forbidden), 'cycle-pack correction versus literal acyclic arc subsets')
            H = sum((-1) ** e * c * factorial(n - e) for e, c in enumerate(coefficients))
            need(H == renewal.oracle.hamilton(partner), 'rook recurrence counts actual native partner paths')
            totals.append(H)
        counts.append({'vertices': n, 'presentations': len(totals), 'sum_partner_H': sum(totals)})
    triangle = (2, 4, 1)
    R = rook_polynomial(triangle, 7)
    A = cycle_free_rooks(triangle)
    naive = sum((-1) ** e * c * factorial(3 - e) for e, c in enumerate(R))
    exact = sum((-1) ** e * c * factorial(3 - e) for e, c in enumerate(A))
    need((naive, exact) == (2, 3), 'uncorrected rook count fails already on a directed triangle')
    out = {'status': 'FINITE-EXACT n=3..5 fixed-path presentations, recovered classical inclusion-exclusion',
           'checks': CHECKS, 'rows': counts, 'triangle_rooks': R, 'triangle_cycle_free_rooks': A,
           'triangle_naive_H': naive, 'triangle_actual_H': exact}
    if write:
        (HERE / 'rooks.json').write_text(json.dumps(out, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'checks': CHECKS, 'rows': counts, 'triangle': [naive, exact]}))
    return out


if __name__ == '__main__':
    main()
