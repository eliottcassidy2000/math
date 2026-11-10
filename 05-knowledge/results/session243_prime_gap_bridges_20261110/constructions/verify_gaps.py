#!/usr/bin/env python3
"""Exact two universal-even-gap compilers; imports only sibling stdlib verifier."""
from itertools import product
import json
import sys
from pathlib import Path
# Resolve only this committed sibling, including under isolated mode.
import importlib.util
spec = importlib.util.spec_from_file_location('physical_verify', Path(__file__).with_name('verify.py'))
pv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pv)


def act(rows, root, word):
    for letter in word:
        n = len(rows)
        rows = pv.join(rows, [0]) if letter == 'O' else pv.graft(rows, root, [0])
        if letter == 'J':
            root = n
    return rows, root


def state_act(z, word):
    D, E, U, d = z
    for letter in word:
        if letter == 'O':
            D, E, U, d = D, 0, 2*U, d
        elif letter == 'Z':
            D, E, U, d = D+U-E, d, 2*U-d, d
        else:
            D, E, U, d = D+U-E, D-E, D+U, D
    return D, E, U, d


def response_word(U):
    pv.need(U >= 2, 'response domain')
    letters = []
    while U > 2:
        if U % 2 == 0:
            letters.append('O')
            U //= 2
        else:
            letters.append('Z')
            U = (U+1)//2
    return ''.join(reversed(letters))


def delete_vertex(rows, deleted):
    kept = [i for i in range(len(rows)) if i != deleted]
    return [sum(1 << j for j,v in enumerate(kept) if rows[u] >> v & 1) for u in kept]


def cheap_pair(N):
    pv.need(N >= 1, 'gap domain')
    source_word = '' if N == 1 else response_word(N) + 'O'
    rows, root = act([2, 0], 0, source_word)
    lower, _ = act(rows, root, 'Z')
    upper, _ = act(rows, root, 'OZ')
    return source_word, rows, root, lower, upper


def endpoint_source(N, extra=0):
    # n = max(4, ceil(log2(N+1))+2) + extra; ceil(log2(N+1))=N.bit_length().
    n = max(4, N.bit_length()+2) + extra
    base = pv.endpoint(n-2)
    root = n
    selected = [i for i in range(1, n-1) if N >> (i-1) & 1]
    rows = [x | ((1 << n) if i in selected else 0) for i, x in enumerate(base)]
    rows.append(sum(1 << i for i in range(n) if i not in selected))
    return rows, root


def main():
    if '--fail' in sys.argv:
        pv.need(False, 'intentional failure gate')
    # Exact linear identities on a basis, independent of realizability.
    basis_rows = {}
    for word in ['JOJOZ', 'JOOOZ']:
        basis_rows[word] = [state_act(tuple(int(i==j) for i in range(4)), word)[0] for j in range(4)]
    pv.need(basis_rows == {'JOJOZ':[9,-3,9,0], 'JOOOZ':[9,-1,9,0]}, 'endpoint word row')
    # Verify transfer on actual small graph states, not just matrix arithmetic.
    transfer_tests = 0
    for n in range(1,5):
        for code in range(1 << (n*(n-1)//2)):
            rows = pv.labelled(n,code)
            for root in range(n):
                z = pv.state(rows,root)
                for letter in 'OZJ':
                    result, newroot = act(rows,root,letter)
                    pv.need(pv.state(result,newroot) == state_act(z,letter), 'literal one-letter law')
                    transfer_tests += 1

    cheap = []
    for N in range(1,65):
        word, rows, root, lower, upper = cheap_pair(N)
        D,E,U,d = pv.state(rows,root)
        pv.need(U==2*N and E==0 and d==1 and D<=N, 'cheap source')
        hlow,hhigh = pv.tables(lower)[0][-1],pv.tables(upper)[0][-1]
        pv.need((hlow,hhigh)==(D+2*N,D+4*N) and hhigh-hlow==2*N and hhigh<=5*N, 'cheap counts')
        pv.need(pv.strong(lower) and pv.strong(upper), 'cheap strongness')
        pv.need(delete_vertex(upper,len(upper)-2)==lower, 'literal induced deletion')
        cap=(N-1).bit_length()+4
        pv.need(max(len(lower),len(upper))<=cap, 'cheap size')
        cheap.append(dict(N=N,source_word=word,D=D,lower=hlow,upper=hhigh,orders=[len(lower),len(upper)]))
    # Pure integer recurrence box beyond the graph-DP regression range.
    for N in range(1,10001):
        word = '' if N==1 else response_word(N)+'O'
        D,E,U,d = state_act((1,0,2,1),word)
        pv.need((E,U,d)==(0,2*N,1) and D<=N, 'response inverse')
        expected_equal = N==1 or (N>=3 and (N-1)&(N-2)==0)
        pv.need((D==N)==expected_equal, 'sharp equality boundary')
        pv.need(len(word)+4<=(N-1).bit_length()+4, 'response length')

    endpoint = []
    for N, extra in [(1,0),(2,0),(3,0),(7,0),(13,0),(21,0),(1,1),(1,2),(7,1)]:
        rows, root = endpoint_source(N,extra)
        z = pv.state(rows,root)
        D,E,U,d=z
        pv.need(E==N and pv.strong(rows), 'endpoint supplier')
        lower,rlow=act(rows,root,'JOJOZ')
        upper,rhigh=act(rows,root,'JOOOZ')
        hlow,hhigh=pv.tables(lower)[0][-1],pv.tables(upper)[0][-1]
        pv.need((hlow,hhigh)==(9*D-3*N+9*U,9*D-N+9*U), 'endpoint word literal')
        pv.need(hhigh-hlow==2*N and len(lower)==len(upper) and pv.strong(lower) and pv.strong(upper), 'same-order gap')
        pv.need(hlow%3==0 and hlow>3, 'endpoint compiler composite obstruction')
        endpoint.append(dict(N=N,extra=extra,source_rows=rows,source_state=z,
            lower=hlow,upper=hhigh,common_order=len(lower)))
    # Short universal coefficient search records scope, not a global minimality claim.
    rows_by_word={}
    for length in range(5):
        for letters in product('OZJ',repeat=length):
            word=''.join(letters)
            row=tuple(state_act(tuple(int(i==j) for i in range(4)),word)[0] for j in range(4))
            rows_by_word.setdefault(row,word)
    no_short_pair=all((r[0],r[1]+2,r[2],r[3]) not in rows_by_word for r in rows_by_word)
    pv.need(no_short_pair,'short word hostile')
    print(json.dumps(dict(status='PASS',one_letter_literal_state_checks=transfer_tests,
        universal_word_rows=basis_rows,cheap_graph_checks=cheap,cheap_integer_inverse_targets=10000,
        endpoint_literal_pairs=endpoint,universal_row_pair_absent_for_both_word_lengths_at_most=4,
        scope='Even differences of attainable total counts, not prime pairs, consecutive attainable values, or arbitrary targets.'),indent=2,sort_keys=True))


if __name__=='__main__':
    main()
