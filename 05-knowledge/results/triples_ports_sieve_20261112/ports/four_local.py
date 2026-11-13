"""Finite local lemmas with an explicit all-order reduction in FOUR_RETURNS.md."""
from collections import Counter
import json
import sys
import verify as v


def trace(B, word):
    states = []
    for k in word:
        B = v.mutate(B, k)
        states.append(B)
    return states


def cycle_block(B):
    n = len(B)
    cycle = B[0][1] == B[1][2] == B[2][0]
    sides = {B[i][j] for i in range(3) for j in range(3, n)}
    return cycle and len(sides) == 1


def main():
    rows = []
    for word, prefix in [((0, 1, 0, 2), 2), ((0, 1, 2, 1), 0)]:
        for n in (3, 4, 5):
            returns = 0
            for mask in range(1 << (n*(n-1)//2)):
                ss = trace(v.matrix(n, mask), word)
                if v.tournament(ss[-1]):
                    returns += 1
                    v.need(v.tournament(ss[prefix]), 'prefix-tournament local lemma')
            rows.append(dict(word=word, n=n, returns=returns, necessary_native_prefix=prefix+1))
    for n in (4, 5):
        valid = []
        for mask in range(1 << (n*(n-1)//2)):
            ss = trace(v.matrix(n, mask), (0, 1, 2, 3))
            if v.tournament(ss[-1]):
                valid.append(mask)
        v.need(valid == [0, (1 << (n*(n-1)//2))-1], 'four-distinct extreme total order')
        rows.append(dict(word=(0, 1, 2, 3), n=n, return_masks=valid))
    triple_returns = []
    for n in (3, 4, 5):
        valid = []
        for mask in range(1 << (n*(n-1)//2)):
            B = v.matrix(n, mask)
            F = trace(B, (0, 1, 2, 0))[-1]
            yes = v.tournament(F)
            if n >= 4:
                v.need(yes == cycle_block(B), 'extreme cyclic triple iff')
                if yes:
                    expected = tuple(tuple(-B[i][j] if i < 3 or j < 3 else B[i][j]
                                           for j in range(n)) for i in range(n))
                    v.need(F == expected, 'exact global block-transfer local matrix')
            elif yes:
                v.need(v.hamilton(B) == v.hamilton(F), 'three-vertex exceptional H')
                triple_returns.append(dict(initial=mask, final=v.maskof(F), H=v.hamilton(B)))
            if yes:
                valid.append(mask)
        rows.append(dict(word=(0, 1, 2, 0), n=n, return_masks=valid))
    out = dict(status='PASS: finite lemmas used in all-order proof',
               local_universe='All labelled tournaments through five for the displayed words',
               rows=rows, exceptional_triples=triple_returns)
    name = sys.argv[1] if len(sys.argv) > 1 else 'four_local_receipt.json'
    path = (v.HERE/name).resolve()
    v.need(path.parent == v.HERE, 'owned directory')
    path.write_text(json.dumps(out, sort_keys=True, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status='PASS', lemma_rows=len(rows), triples=triple_returns), sort_keys=True))


if __name__ == '__main__':
    main()
