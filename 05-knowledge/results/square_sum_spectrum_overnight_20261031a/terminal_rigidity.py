#!/usr/bin/env python3
"""Universal prefix identity, independently witnessed before 17-state reduction."""
from itertools import product
import json
from pathlib import Path
import sys
import sympy as sp

PREVIOUS = Path(__file__).resolve().parent.parent / 'square_sum_spectrum_overnight_20261030'
sys.path.insert(0, str(PREVIOUS))
from minimal_feedback import model
from audit_feedback import gaps, require, subset_dp, tournament

TERMS = [((),2), ((1,),-2), ((2,),-3), ((3,),3),
         ((1,0),-2), ((1,1),2), ((1,2),1), ((1,3),-1),
         ((2,2),1), ((2,3),-1), ((1,0,2),2), ((1,0,3),-2),
         ((1,2,2),-1), ((1,2,3),1)]


RESTRICTED = [
    ('03',24,[((),-33),((3,),22),((3,0),51),((3,3),12),
               ((3,0,0),-21),((3,0,3),-51),((3,3,3),2),
               ((3,0,0,3),15),((3,0,3,0),3)]),
    ('12',4,[((1,),4),((2,),-2),((1,1),-4),((1,2),9),
             ((2,1),-4),((2,2),2),((1,1,2),-7),((1,2,1),-2),
             ((1,2,2),-3),((2,1,1),4),((1,1,2,1),6),
             ((1,1,2,2),1),((1,2,1,1),-4)])]


def main():
    M, u, v = model()
    D = sp.zeros(1, 36)
    for word, coefficient in TERMS:
        row = u
        for a in word:
            row = row * M[a]
        D += coefficient * row
    require(sum(c for _, c in TERMS) == 0, 'nonzero coefficient mass')
    require((D*v)[0] == 2, 'wrong constant term')
    cols, words = [v], [()]
    i = 0
    while i < len(cols):
        for a in range(4):
            candidate = M[a] * cols[i]
            if sp.Matrix.hstack(*cols, candidate).rank() > len(cols):
                cols.append(candidate)
                words.append((a,)+words[i])
        i += 1
    C = sp.Matrix.hstack(*cols)
    for a in range(4):
        require((D*M[a]-D)*C == sp.zeros(1, len(cols)), 'invariance failed')
        require(C.row_join(M[a]*C).rank() == len(cols), 'observable closure failed')
    # An unrelated vertex-subset algorithm verifies the 14-term identity on
    # all suffixes through length two; the matrix closure proves all lengths.
    literal_controls = 0
    for n in range(3):
        for suffix in product(range(4), repeat=n):
            value = sum(c*subset_dp(tournament(p+suffix)) for p,c in TERMS)
            require(value == 2, 'independent prefix identity failed')
            literal_controls += 1
    restricted_receipts=[]
    for alphabet,constant,terms in RESTRICTED:
        row=sp.zeros(1,36)
        require(sum(c for _,c in terms)==0,'restricted coefficient sum')
        for word,coefficient in terms:
            require(set(word)<=set(map(int,alphabet)),'restricted prefix alphabet')
            r=u
            for a in word:r=r*M[a]
            row += coefficient*r
        require((row*v)[0]==constant,'restricted constant')
        for a in range(4):
            require((row*M[a]-row)*C==sp.zeros(1,len(cols)),
                    'restricted identity observable invariance')
        for n in range(3):
            for suffix in product(range(4),repeat=n):
                require(sum(c*subset_dp(tournament(p+suffix)) for p,c in terms)==constant,
                        'independent restricted prefix identity')
        restricted_receipts.append({'prefix_alphabet':alphabet,'constant':constant,
                                    'coefficient_sum':0,'terms':terms,
                                    'independent_suffix_controls':21,
                                    'suffix_alphabet':'0123'})
    short = [w for n in range(3) for w in product(range(4), repeat=n)]
    debt_controls = 0
    for A in short:
        for B in short:
            if not B:
                continue
            nontrivial = any(a != 3 for a in A+B)
            for m in range(1, 5):
                difference = gaps((3,)*m+B+A)-gaps((3,)*m+A)
                require(difference >= 2**m if nontrivial else difference == 0,
                        'odd-cycle insertion bound failed')
                debt_controls += 1
    print(json.dumps({'status':'passed', 'prefix_coefficients':TERMS,
                      'coefficient_sum':0, 'universal_constant':2,
                      'restricted_prefix_identities':restricted_receipts,
                      'observable_dimension':len(cols), 'observable_basis':words,
                      'independent_subset_dp_identity_controls':literal_controls,
                      'insertion_debt_controls':debt_controls,
                      'scope':'All-length identity via exact observable closure; quantitative insertion bound additionally uses the proved OCF formula.'},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
