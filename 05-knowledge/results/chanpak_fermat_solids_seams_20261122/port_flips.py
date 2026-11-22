"""Exact changes under one reversed interport or intraport arc in the strong lower host."""
from collections import Counter, defaultdict
from itertools import combinations, product
from fractions import Fraction
from math import comb
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'chanpak_support_birth_grammar_20261122'))
from check_perfect_port import LOWER, LPORTS, H, inflate


def bridge_tensor(source, destination):
    ordinary = sum(1 << i for i in range(8) if i not in LPORTS)
    dp = {(0, (0, 0, 0), -1, False): 1}
    table = Counter()
    while dp:
        next_layer = defaultdict(int)
        for (used, counts, last, crossed), ways in dp.items():
            if crossed and used == ordinary and all(counts):
                table[counts] += ways
            for v in range(8):
                if v == last:
                    continue
                normal = last == -1 or bool(LOWER[last] >> v & 1)
                special = not crossed and (last, v) == (source, destination)
                if not (normal or special):
                    continue
                mask = used
                new_counts = counts
                if v in LPORTS:
                    i = LPORTS.index(v)
                    if counts[i] == 3:
                        continue
                    new_counts = counts[:i] + (counts[i] + 1,) + counts[i+1:]
                elif used >> v & 1:
                    continue
                else:
                    mask |= 1 << v
                if normal:
                    next_layer[mask, new_counts, v, crossed] += ways
                if special:
                    next_layer[mask, new_counts, v, True] += ways
        dp = next_layer
    return table


def end_weight(width, index, runs):
    return sum((-1) ** (runs-1-k) * comb(runs-1, k) * (k+1) ** index * k ** (width-1-index)
               for k in range(runs))


def onto(width, runs):
    return sum((-1) ** (runs-k) * comb(runs, k) * k ** width for k in range(runs+1))


def bridge_grid(widths, a, b, tensor):
    c = 3-a-b
    weights = Counter()
    for runs, ways in tensor.items():
        weights[runs[a], runs[b]] += ways * onto(widths[c], runs[c])
    left = [[end_weight(widths[a], i, r) for r in range(1, 4)] for i in range(widths[a])]
    right = [[end_weight(widths[b], widths[b]-1-j, r) for r in range(1, 4)] for j in range(widths[b])]
    return [[sum(left[i][r-1] * right[j][s-1] * coef for (r,s), coef in weights.items())
             for j in range(widths[b])] for i in range(widths[a])]


def cross_changes(widths, a, b, banks):
    if not (LOWER[LPORTS[a]] >> LPORTS[b] & 1):
        a,b = b,a
    old = bridge_grid(widths, a, b, banks[a,b])
    new = bridge_grid(widths, b, a, banks[b,a])
    return {(LPORTS[a],i,LPORTS[b],j):new[j][i]-old[i][j] for i,j in product(range(widths[a]),range(widths[b]))}


def check_corner_polynomials(banks):
    for x, y in product(range(2, 7), repeat=2):
        twos = (8*x*x, x, x)
        threes = (27*y*y, y, y)
        def full(i, r):
            return (0, 1, twos[i]-2, threes[i]-3*twos[i]+3)[r]
        def endpoint(i, r, new):
            if new:
                return (0, 0, 1, Fraction(twos[i], 2)-2)[r]
            return Fraction(full(i, r), r)
        predicted_twice = (24*x**3-4*x**2+9*x-4,
                           80*x**4-8*x**3+28*x**2-x-12 + y*(72*x**3-80*x**2+18*x-20),
                           64*x**4+8*x**2)
        for (a, b), twice in zip(combinations(range(3), 2), predicted_twice):
            c = 3-a-b
            def corner(table, new):
                return sum(ways*full(c, r[c])*endpoint(a, r[a], new)*endpoint(b, r[b], new)
                           for r, ways in table.items())
            assert 2*(corner(banks[b,a], True)-corner(banks[a,b], False)) == twice


def main():
    banks = {(i,j):bridge_tensor(LPORTS[i], LPORTS[j]) for i,j in product(range(3),repeat=2) if i!=j}
    check_corner_polynomials(banks)
    for p in (5,7,61):
        widths=(p,(p-3)//2,(p-3)//2)
        for a,b in combinations(range(3),2):
            deltas = cross_changes(widths, a, b, banks)
            if p==5:
                assert min(deltas.values()) == {(0,1):94, (0,2):860, (1,2):426}[a,b]
                rows = inflate(LOWER, LPORTS, widths)
                offsets = [sum(dict(zip(LPORTS,widths)).get(i,1) for i in range(v)) for v in range(8)]
                h = H(rows)
                for (u,i,v,j), delta in deltas.items():
                    x,y = offsets[u]+i,offsets[v]+j
                    other = list(rows);other[x] ^= 1<<y;other[y] ^= 1<<x
                    assert H(other)-h == delta
            if p == 61:
                V = list(deltas.values())
                closest = min(deltas, key=lambda key: abs(deltas[key]-2))
                print(LPORTS[a],LPORTS[b],len(V),'negative/zero/positive',sum(v<0 for v in V),V.count(0),sum(v>0 for v in V),
                      'closest',closest,deltas[closest], 'unit',V.count(2),flush=True)


if __name__=='__main__':
    main()
