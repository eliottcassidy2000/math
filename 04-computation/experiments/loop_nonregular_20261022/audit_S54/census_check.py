"""Auditor: completeness of the almost-regular class lists ar{4,6,8,10}_classes.txt.
(a) labelled number of tournaments with a prescribed score vector, by an independent memoised recursion
    (vertex 0 chooses its out-set; the rest must realise the reduced vector), for regular orders 5, 7, 9, 11;
(b) brute-force canonical forms over score-respecting relabellings (canon_brute.exe -> canon_ar*.out):
    classes pairwise distinct, |Aut| equal to the file's value, mass sum n!/|Aut| equal to (a);
(c) self-converse counts.   Usage: python census_check.py > census_check.out"""
import os, sys, math
from functools import lru_cache
from itertools import combinations
HERE = os.path.dirname(os.path.abspath(__file__))
UP = os.path.join(HERE, '..')


@lru_cache(maxsize=None)
def count(scores):
    """number of labelled tournaments on len(scores) vertices with out-degree vector `scores` (sorted tuple)."""
    m = len(scores)
    if m == 0:
        return 1
    if any(s < 0 or s > m - 1 for s in scores) or sum(scores) != m * (m - 1) // 2:
        return 0
    s0, rest = scores[0], scores[1:]
    tot = 0
    # vertex 0 beats a subset X of size s0 of the others; each y in X keeps its score, each y not in X loses one
    # (it used one win on vertex 0).  Count subsets by multiplicity pattern over equal scores.
    vals = sorted(set(rest))
    mult = [rest.count(v) for v in vals]

    def rec(i, chosen, acc):
        nonlocal tot
        if i == len(vals):
            if chosen == s0:
                new = []
                for v, k, x in acc:
                    new += [v] * x + [v - 1] * (k - x)
                tot += math.prod(math.comb(k, x) for v, k, x in acc) * count(tuple(sorted(new)))
            return
        for x in range(0, mult[i] + 1):
            if chosen + x > s0:
                break
            rec(i + 1, chosen + x, acc + [(vals[i], mult[i], x)])
    rec(0, 0, [])
    return tot


def main():
    for q in (5, 7, 9, 11):
        print(f'labelled regular {q}-tournaments: {count(tuple([(q - 1) // 2] * q))}')
    for n in (4, 6, 8):
        a = count(tuple(sorted([n // 2 - 1] * (n // 2) + [n // 2] * (n // 2))))
        print(f'labelled almost-regular {n}-tournaments: {a} per fixed score vector x C({n},{n // 2}) = {a * math.comb(n, n // 2)}')
    for n in (4, 6, 8, 10):
        fn = os.path.join(UP, f'ar{n}_classes.txt')
        cf = os.path.join(HERE, f'canon_ar{n}.out')
        if not os.path.exists(cf):
            print(f'n={n}: canon file missing'); continue
        file_aut = [int(l.split()[2]) for l in open(fn)]
        mine = [l.split() for l in open(cf)]
        if len(mine) != len(file_aut):
            print(f'n={n}: canon file incomplete ({len(mine)} of {len(file_aut)})'); continue
        canon = [m[0] for m in mine]; aut = [int(m[1]) for m in mine]; sc = [int(m[2]) for m in mine]
        mass = sum(math.factorial(n) // a for a in aut)
        target = count(tuple([n // 2] * (n + 1)))
        print(f'n={n}: classes {len(canon)}, distinct canonical forms {len(set(canon))}, '
              f'|Aut| agrees with file {sum(a == b for a, b in zip(aut, file_aut))}/{len(aut)}, '
              f'mass {mass} vs labelled regular {n + 1}-tournaments {target}: {mass == target}; '
              f'self-converse classes {sum(sc)}; |Aut| distribution { {a: aut.count(a) for a in sorted(set(aut))} }')


if __name__ == '__main__':
    main()
