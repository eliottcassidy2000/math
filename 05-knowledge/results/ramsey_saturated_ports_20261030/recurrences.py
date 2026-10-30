"""Small complete universes checking exact labelled and rooted growth laws."""
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent


def require(ok, why):
    if not ok:
        raise RuntimeError(why)


def legal_words(C, targets):
    n = len(C)
    for f in product(range(3), repeat=n):
        if all(not any(all(C[u][v] == c for u, v in combinations(S, 2))
                       for S in combinations([u for u in range(n) if f[u] == c], targets[c]-1))
               for c in range(3)):
            yield f


def extend(C, f):
    return tuple(tuple(row)+(f[i],) for i, row in enumerate(C))+(tuple(f)+(-1,),)


def word(C, p):
    return tuple(C[p[i]][p[j]] for i, j in combinations(range(len(C)), 2))


def aut(C):
    original = word(C, tuple(range(len(C))))
    return [p for p in permutations(range(len(C))) if word(C, p) == original]


def reps(level):
    out = {}
    for C in level:
        key = min(word(C, p) for p in permutations(range(len(C))))
        out.setdefault(key, C)
    return list(out.values())


def vertex_orbits(C):
    group = aut(C)
    return len({tuple(sorted({p[v] for p in group})) for v in range(len(C))})


def check(k):
    targets = (3, 3, k)
    level, counts, small = [()], [1], {}
    for n in range(1, 6):
        level = [extend(C, f) for C in level for f in legal_words(C, targets)]
        counts.append(len(level))
        if n in (3, 4):
            small[n] = reps(level)
    rooted_parent = 0
    mass_parent = Fraction(0)
    for C in small[3]:
        group = aut(C)
        words = list(legal_words(C, targets))
        mass_parent += Fraction(len(words), len(group))
        rooted_parent += len({min(tuple(f[p[i]] for i in range(3)) for p in group) for f in words})
    rooted_child = sum(vertex_orbits(C) for C in small[4])
    mass_child = sum((Fraction(1, len(aut(C))) for C in small[4]), Fraction(0))
    require(rooted_parent == rooted_child, 'Unmarked parent attachments count rooted child orbits')
    require(mass_parent/4 == mass_child, 'Automorphism-weighted recurrence')
    require(mass_child*24 == counts[4], 'Labelled count from group mass')
    return {'targets': targets, 'labelled_counts_n0_to_n5': counts,
            'unlabelled_n3_n4': [len(small[3]), len(small[4])],
            'rooted_child_orbits_n4': rooted_child, 'child_group_mass_n4': str(mass_child),
            'parent_weighted_attachment_mass': str(mass_parent)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()
    result = {str(k): check(k) for k in [3, 5, 6]}
    require(result['5']['labelled_counts_n0_to_n5'][-1]+1 == result['6']['labelled_counts_n0_to_n5'][-1],
            'Only the all-blue K5 separates the two five-vertex universes')
    text = json.dumps(result, indent=2, sort_keys=True)+'\n'
    target = HERE/'recurrences.json'
    if args.write:
        target.write_text(text, encoding='utf-8', newline='\n')
    else:
        require(json.loads(target.read_text()) == json.loads(text), 'Growth checks reproduce')
    print(text, end='')


if __name__ == '__main__':
    main()
