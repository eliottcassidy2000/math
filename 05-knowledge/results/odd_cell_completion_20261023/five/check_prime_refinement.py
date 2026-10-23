"""All binary colourings of selected prime cyclic tournaments; exact controls."""
from itertools import product
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent


def refine(A, cells):
    while True:
        new = []
        for C in cells:
            buckets = {}
            for i in C:
                word = tuple(sum(A[i][j] for j in D) for D in cells)
                buckets.setdefault(word, []).append(i)
            new.extend(buckets.values())
        if len(new) == len(cells):
            return new
        cells = new


def main():
    outputs = []
    total = 0
    for p in (3, 5, 7, 11, 13):
        skew = [{i if sign else p-i for i, sign in enumerate(bits, 1)}
                for bits in product((0, 1), repeat=(p-1)//2)]
        trivial_multiplier = [S for S in skew if all(
            {a*s % p for s in S} != S for a in range(2, p))]
        # All trivial-multiplier Cayley types at p<=7, and the cyclic R_p
        # representative at11/13. This is not an exhaustive p13 tournament test.
        selected = trivial_multiplier if p <= 7 else [set(range(1, (p-1)//2+1))]
        cases = 0
        for S in selected:
            A = [[int((j-i) % p in S) for j in range(p)] for i in range(p)]
            for mask in range(1 << (p-1)):
                first = [0]+[i for i in range(1, p) if mask >> (i-1) & 1]
                second = [i for i in range(p) if i not in first]
                if not second:
                    continue
                if any(len(C) != 1 for C in refine(A, [first, second])):
                    raise ValueError(('proper nondiscrete prime partition', p, sorted(S), first))
                cases += 1
        outputs.append({'prime': p, 'all_skew_connection_sets': len(skew),
                        'trivial_multiplier_connection_sets': len(trivial_multiplier),
                        'tested_connection_sets': [sorted(S) for S in selected],
                        'all_two_colour_splits_for_tested_sets': cases})
        total += cases
    result = {'scope': 'all binary splits for declared cyclic prime tournaments; general statement proved separately',
              'cases': outputs, 'total_binary_split_controls': total}
    destination = HERE/'prime_refinement.json'
    if '--write' in sys.argv:
        destination.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    elif json.loads(destination.read_text()) != result:
        raise ValueError('receipt mismatch')
    print(json.dumps({'binary_split_controls': total}, sort_keys=True))


if __name__ == '__main__':
    main()
