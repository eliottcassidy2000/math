"""Exact parity-dependent affine obstructions in the complete small receiver census."""
from collections import Counter
from math import gcd

from check import graft_chart, strong_readers


def branches(target, chart):
    c, h, even, odd = chart
    return tuple(j for j, k in enumerate((even, odd))
                 if (target-k-j*c) % (2*gcd(c,h)) == 0)


def main():
    n = (1 << 60)*((1 << 61)-1)-1
    classes = Counter()
    uniformly_rejected = 0
    for reader in strong_readers(9):
        c, h, k0, k1 = graft_chart(reader)
        classes[branches(n, (c, h, k0, k1))] += 1
        if c % 3 == h % 3 == 0 and k0 % 3 != 0 and k1 % 3 != 0:
            uniformly_rejected += 1
    assert classes == {(): 214, (0,): 24, (0,1): 372}
    assert uniformly_rejected == 99
    print('610 receivers, parity branches and numbers:', classes)
    print('All odd exponents excluded by modulo three:', uniformly_rejected)
    print('Triangle chart:', graft_chart((2,4,1)))


if __name__ == '__main__':
    main()
