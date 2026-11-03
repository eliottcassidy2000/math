"""Independent four-coordinate fresh-block audit; no repository imports.

All rooted labelled tournaments through order3, plus regularC5 and the actual
adaptive H37 order6 source; k=1,2,3 and every position of the new marked vertex.
"""
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import json


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def from_code(n, code):
    rows = [0]*n
    for bit, (i, j) in enumerate(combinations(range(n), 2)):
        if code >> bit & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return tuple(rows)


def birth(rows, v):
    n = len(rows)
    return tuple(r | ((1 << n) if i == v else 0) for i, r in enumerate(rows)) + (((1 << n)-1) ^ (1 << v),)


def profile(rows, v):
    n, full = len(rows), (1 << len(rows))-1

    @lru_cache(None)
    def ending(mask, last):
        rest = mask ^ (1 << last)
        if not rest:
            return 1
        return sum(ending(rest, u) for u in range(n) if rest >> u & 1 and rows[u] >> last & 1)

    @lru_cache(None)
    def H(mask):
        return sum(ending(mask, u) for u in range(n) if mask >> u & 1) if mask else 1

    @lru_cache(None)
    def starting(first, remaining):
        if not remaining:
            return 1
        return sum(starting(u, remaining ^ (1 << u)) for u in range(n)
                   if remaining >> u & 1 and rows[first] >> u & 1)

    F = sum(ending(mask, v)*H(full ^ mask) for mask in range(1, full+1) if mask >> v & 1)
    return H(full), starting(v, full ^ (1 << v)), H(full ^ (1 << v)), F


def formula(z, A, B):
    h, s, q, F = z
    a, d, k = h-s, F-q, A+B+1
    return (a+q+(2**k-1)*d,
            2**B*d if A else a+(2**B-1)*d,
            h if k == 1 else a+q+(2**(k-1)-1)*d,
            2**A*a+(2**(B+1)*3**A-2**A)*d+2**(A+1)*q+(s-q if B == 0 else 0))


def main():
    sources = [from_code(n, c) for n in range(1, 4) for c in range(1 << (n*(n-1)//2))]
    sources += [(10, 20, 33, 6, 13, 27),
                tuple(sum(1 << ((i+d) % 5) for d in (1, 2)) for i in range(5))]
    controls = []
    for T in sources:
        for v in range(len(T)):
            z = profile(T, v)
            for k in range(1, 4):
                U = T
                for _ in range(k):
                    U = birth(U, v)
                for B in range(k):
                    A, w = k-1-B, len(T)+B
                    predicted, actual = formula(z, A, B), profile(U, w)
                    need(actual == predicted, 'four literal profile coordinates')
                    if A == 1 and B == 0:
                        need(actual[0]-actual[2] == 2*actual[1], 'singular single-command image')
                    controls.append({'source_rows': T, 'old_root': v, 'source_profile': z,
                                     'newer_A': A, 'older_B': B, 'child_rows': U,
                                     'new_mark': w, 'actual_profile': actual})
    need(len(controls) == 240, 'declared complete small universe')
    result = {
        'status': 'PASS', 'checks': len(controls), 'maximum_child_order': 9,
        'scope': 'All roots of all labelled tournaments orders1..3, plus regularC5 and adaptive '
                 'H37 order6; k1..3 and every new mark. Direct subset/end recursion counts all four '
                 'profile coordinates. This finite check is separate from the all-order combinatorial proof.',
        'profile_convention': '(h,s,q,F): Hamilton paths; paths starting at root; H after deleting root; '
                              'ordered two-slot spanning covers with first nonempty path ending at root, '
                              'second path possibly empty. H(empty)=1.',
        'controls': controls,
    }
    path = Path(__file__).with_suffix('.json')
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'PASS', 'checks': len(controls), 'maximum_child_order': 9}))


if __name__ == '__main__':
    main()
