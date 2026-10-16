#!/usr/bin/env python3
"""Source-preserving packets: exact residuals, reset controllers, native clocks.

Run from the repository root. No eventual ROOT labels or empirical signs enter
controller selection. All acceptance decisions use integers or Fractions.
"""
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path
import numpy as np

from collatz_odd_ray_overlap_20261015d import rows, strip, w
from collatz_gram_slack_20261015b import T, close_graph, need

ROOT = Path(__file__).resolve().parents[2]


def image(row, action):
    out = defaultdict(F)
    for n, v in row.items():
        e, s = action[n]
        out[e] += s*v
    return {e: v for e, v in out.items() if v}


def packet(row, action):
    """The minimum conformal, rowwise 7/10 residual, in source coordinates."""
    fibres = defaultdict(lambda: [[], []])
    for n, v in sorted(row.items()):
        if v:
            e, s = action[n]
            fibres[e][int(s*v < 0)].append([n, abs(F(v))])
    survivors = {}
    for pos, neg in fibres.values():
        i = j = 0
        while i < len(pos) and j < len(neg):
            amount = min(pos[i][1], neg[j][1])
            pos[i][1] -= amount
            neg[j][1] -= amount
            i += pos[i][1] == 0
            j += neg[j][1] == 0
        for n, amount in pos+neg:
            if amount:
                survivors[n] = amount
    u = sum(abs(v) for v in row.values())
    v = sum(abs(v) for v in image(row, action).values())
    r = max(F(0), F(10*v-7*u, 3))
    left, residual = r, {}
    for n, amount in sorted(survivors.items()):
        take = min(left, amount)
        if take:
            residual[n] = take if row[n] > 0 else -take
            left -= take
    need(left == 0, 'survivor capacity')
    core = {n: F(v)-residual.get(n, 0) for n, v in row.items()}
    core_image, res_image = image(core, action), image(residual, action)
    need(sum(abs(v) for v in residual.values()) == r, 'exact residual norm')
    need(sum(abs(v) for v in res_image.values()) == r, 'residual never cancels')
    need(sum(abs(v) for v in core.values()) == u-r, 'conformal source split')
    need(sum(abs(v) for v in core_image.values()) == v-r, 'all old cancellation retained')
    need(10*(v-r) <= 7*(u-r), 'literal packet contraction')
    need(all(core.get(n, 0)+residual.get(n, 0) == value for n, value in row.items()),
         'each original source weight preserved')
    return core, residual, r


def annular_rows(k, cache, materialize=False):
    """O(X + number_of_bins * number_of_terminals) exact event sweep.

    This is checked against the independently implemented c_I(n) formula.
    Matrix entries and the sweep are int64; squared sums are Python integers.
    """
    X, b = 1 << k, 1 << (k//2)
    need(32*b*X < 2**63, 'int64 absolute row/sweep bound')
    cols = sorted({cache[n][0] for n in range(1, 8*X+1)})
    index = {e: i for i, e in enumerate(cols)}
    target = np.array([0]+[index[cache[n][0]] for n in range(1, 8*X+1)], dtype=np.int64)
    value = np.array([0]+[w(n)*cache[n][1] for n in range(1, 8*X+1)], dtype=np.int64)
    delta = np.zeros(len(cols), dtype=np.int64)
    np.add.at(delta, target[X+1:4*X+1], value[X+1:4*X+1])
    norms, matrix = [], []
    for a in range(X, 2*X, b):
        row = b*delta
        nn = np.arange(a, a+b, dtype=np.int64)
        remaining = a+b-1-nn
        for j in (1, 2, 3, 4):
            ids = 4*nn+j
            np.add.at(row, target[ids], value[ids]*remaining)
            np.add.at(delta, target[ids], value[ids])
        ids = nn+1
        np.add.at(row, target[ids], -value[ids]*remaining)
        np.add.at(delta, target[ids], -value[ids])
        norms.append(int(np.sum(np.abs(row))))
        if materialize:
            matrix.append(row)
    return norms, cols, matrix


def extend(edges, terminals, reset):
    """One odd block, or through the first valuation >=2; then graph glue."""
    added, histogram = [], Counter()
    for source in sorted(set(terminals)-{1}):
        x, ones = source, 0
        while x != 1:
            y = T(x)
            if x not in edges:
                added.append((x, y))
            edges[x] = y
            a = 1
            while y % 2 == 0:
                z = y//2
                if y not in edges:
                    added.append((y, z))
                edges[y] = z
                y, a = z, a+1
            x = y
            if a >= 2 or not reset:
                histogram[ones] += 1
                break
            ones += 1
    need(all(x != 1 and y == T(x) for x, y in added), 'every new native edge authenticated')
    return len(added), dict(sorted(histogram.items()))


def ceil_sqrt_fraction(numerator, denominator):
    a = math.isqrt(numerator//denominator)
    return a if a*a*denominator == numerator else a+1


def run(k, reset, stop_below_seven=True, controls=False):
    X, b = 1 << k, 1 << (k//2)
    edges = {n: n//2 for n in range(2, 8*X+1, 2)}
    cache, cycles = close_graph(edges, range(1, 8*X+1))
    norms, cols, matrix = annular_rows(k, cache, controls)
    initial = sum(v*v for v in norms)
    need(initial <= 64*b*X**3, 'initial 8 X^(3/2) bound')
    bound = F(ceil_sqrt_fraction(initial, b))
    out = []
    for stage in range(1, k+9):
        old_cache = cache if controls else None
        old_cols, old_matrix = cols, matrix
        new_edges, histogram = extend(edges, cols, reset)
        cache, cycles = close_graph(edges, range(1, 8*X+1))
        need(not cycles, 'no witnessed nonROOT cycle in this finite universe')
        newer, cols, matrix = annular_rows(k, cache, controls)
        need(all(v <= u for u, v in zip(norms, newer)), 'row norm native monotonicity')
        rnum = [max(0, 10*v-7*u) for u, v in zip(norms, newer)]
        residual_square_num = sum(z*z for z in rnum)
        Qsquare_num = sum(z*z for z in newer)
        bound = F(7, 10)*bound+F(3, 10)*ceil_sqrt_fraction(residual_square_num, 9*b)
        if controls:
            direct, direct_cols = rows(k, cache)
            pos = {e: j for j, e in enumerate(direct_cols)}
            need(np.array_equal(np.array(matrix), direct[:, [pos[e] for e in cols]]),
                 'event sweep equals original integer c_I(n) rows')
            transport, _ = close_graph(edges, old_cols)
            action = {n: transport[n][:2] for n in old_cols}
            for old, new in zip(old_matrix, matrix):
                literal = {e: F(int(v)) for e, v in zip(old_cols, old) if v}
                pushed = image(literal, action)
                need(pushed == {e: F(int(v)) for e, v in zip(cols, new) if v},
                     'same source rows under compiled native identification')
                packet(literal, action)
            # An independent literal T replay checks exact, not just parity, clocks.
            for n in sorted({1, X+1, 2*X-1, 4*X+3, 8*X}):
                e, sign, clock = cache[n]
                x = n
                for _ in range(clock):
                    need(x != 1, 'strict ROOT stopping')
                    x = T(x)
                need((x, (-1)**clock) == (e, sign), 'literal native clock replay')
            del old_cache
        record = dict(stage=stage, residual_square=[residual_square_num, 9*b],
                      residual_over_X=math.sqrt(residual_square_num/(9*b))/X,
                      row_envelope=[Qsquare_num, b], Q_over_X=math.sqrt(Qsquare_num/b)/X,
                      pass_seven=residual_square_num < 441*b*X*X,
                      row_norms=newer, residual_numerators=rnum, terminals=len(cols),
                      new_edges=new_edges, edges=len(edges), reset_ones_histogram=histogram,
                      unresolved_sources=sum(cache[n][0] != 1 for n in range(1, 8*X+1)),
                      max_native_clock=max(cache[n][2] for n in range(1, 8*X+1)))
        out.append(record)
        norms = newer
        if stop_below_seven and Qsquare_num < 49*b*X*X:
            break
    # The remaining k+8 stages can be identities, C=0, R=B, with the same clock.
    padded = k+8-len(out)
    if padded:
        need(out[-1]['row_envelope'][0] < 49*b*X*X, 'identity residual is below 7X')
    rec = dict(k=k, controller='first_reset' if reset else 'one_odd',
               source_interval=[1, 8*X], bin_width=b, initial_envelope=[initial, b],
               stages=out, identity_padding=padded,
               discounted_upper_bound=[bound.numerator, bound.denominator],
               discounted_bound_below_7p5X=2*bound < 15*X,
               all_stage_residuals_below_7X=all(v['pass_seven'] for v in out))
    print(json.dumps(dict(k=k, controller=rec['controller'], stages=len(out),
                          max_R_over_X=max(v['residual_over_X'] for v in out),
                          Q_over_X=out[-1]['Q_over_X'], padded=padded)), flush=True)
    return rec


def controls():
    count = 0
    actions = [{j: (j % 2, 1) for j in range(4)},
               {j: (j % 2, (-1)**j) for j in range(4)},
               {j: (0, (-1)**(j//2)) for j in range(4)}]
    for values in itertools.product(range(-2, 3), repeat=4):
        for action in actions:
            packet(dict(enumerate(map(F, values))), action)
            count += 1
    # No collision: all mass must remain in the conformal residual.
    _, _, r = packet({0: F(2), 1: F(-3)}, {0: (0, 1), 1: (1, 1)})
    need(r == 5, 'injective hostile')
    # Inherited reset-ray recurrence, checked on actual paths and with clock.
    rays = []
    for p in range(4):
        j = p+1
        period = 2*3**p
        multiplier = 2**period
        intercept = (multiplier-1)*(3**j-2**j)//3**j
        need((multiplier-1)*(3**j-2**j) % 3**j == 0, 'integral reset-ray intercept')
        cases = 0
        for n in range(3, 4096, 2):
            if (n+1 & -(n+1)).bit_length()-2 != p:
                continue
            x, word = n, []
            for _ in range(j):
                x, a = strip(3*x+1)
                word.append(a)
            need(word[:-1] == [1]*p and word[-1] >= 2, 'actual first reset')
            lifted = multiplier*n+intercept
            y, newword = lifted, []
            for _ in range(j):
                y, a = strip(3*y+1)
                newword.append(a)
            need(y == x and newword == word[:-1]+[word[-1]+period],
                 'same endpoint, unchanged head, exact clock increment')
            cases += 1
        rays.append(dict(p=p, period=period, multiplier=multiplier,
                         intercept=intercept, checked_sources=cases))
    need(2**18-1 == 27*7*19*73 and (2**9+1)//27 == 19 and 3**3-2**3 == 19,
         '18/19 arithmetic, with factors and anchors distinguished')
    need(F(7,10)**2 < F(1,2) and 8*F(7,10)**8 < F(1,2),
         'all-k stage bound via (sqrt(2)*7/10)^k <= 1')
    return dict(packet_checks=count, reset_rays=rays)


def main():
    checks = controls()
    census = [run(k, True, controls=(k in (8, 9))) for k in range(8, 19)]
    hostile = run(18, False, stop_below_seven=False)
    need(all(r['all_stage_residuals_below_7X'] for r in census), 'declared reset census')
    need(not hostile['all_stage_residuals_below_7X'], 'one-odd schedule hostile')
    result = dict(status='PROVED packet reduction; FINITE-EXACT census; all-scale residual OPEN',
                  controls=checks, first_reset=census, one_odd_hostile=hostile)
    data = json.dumps(result, indent=2, sort_keys=True)+'\n'
    path = ROOT/'05-knowledge/results/collatz_survivor_packets_20261015e.json'
    path.write_text(data)
    print(json.dumps(dict(output=str(path), sha256=hashlib.sha256(data.encode()).hexdigest())))


if __name__ == '__main__':
    main()
