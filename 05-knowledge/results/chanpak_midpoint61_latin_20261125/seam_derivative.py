"""Collect exact exponential coefficients for a contact with an interior chain rank."""
from collections import Counter
from itertools import product
from math import comb, prod
from seam_quarter import ROWS, two_slot_kernel

PORTS = (0, 4, 6, 9)


def ordinary_filling(k):
    # Positive-width chains annihilate the omitted base zero.
    return [(s, (-1)**(k-s) * comb(k, s)) for s in range(1, k+1)]


def pinned(k, first):
    # Interior rank leaves positive numbers of ranks on both sides.
    result = []
    for s in range(1,k):
        bases = (s,s+1) if first else (s+1,s)
        result.append((bases,(-1)**(k-1-s)*comb(k-1,s)))
    return result


def coefficients(kernel, mark, singleton):
    total=Counter()
    other=[p for p in PORTS if p!=mark]
    for left,right in [(mark,singleton),(singleton,mark)]:
        # Changes give count in the new order minus count in the old order.
        sign=-1 if ROWS[left]>>right&1 else 1
        for runs,words in kernel[left,right].items():
            count=dict(zip(PORTS,runs))
            terms=[ordinary_filling(count[p]) for p in other]
            for assignment in product(*terms):
                base=tuple(a for a,b in assignment)
                weight=prod(b for a,b in assignment)
                for tail,factor in pinned(count[mark],first=(mark==right)):
                    total[base+tail]+=sign*words*weight*factor
    return {key:value for key,value in sorted(total.items()) if value}


if __name__ == "__main__":
    K=two_slot_kernel(ROWS,PORTS)
    for mark in [6,9]:
        print(mark, coefficients(K,mark,8))
