#!/usr/bin/env python3
"""A11 countdown ends (S56).  Round map on (a0, b, c) with a0 in {0,1} (paper's rules A.11, a-loop collapsed):
C(a0, b, c+2, d) -> C(a0+b+4, 3, c, d) -> C((a0+b+4) mod 2, 3 + 3*floor((a0+b+4)/2), c, d).
At c in {0,1}, d = 3: HALT iff a0 = 0 and (c = 0 and b odd, or c = 1 and b even).  Every other end case starts
a new countdown with c ~ b/2 (rules C(a,2k,0,3)->C(a,0,k+3,3), C(1,2k+1,0,3)->C(6,3,k,3), ...).
Usage: python a11_countdown_20261025.py seed cmin cmax trials
"""
import sys, random
from collections import Counter


def countdown(a0, b, c):
    while c >= 2:
        a = a0 + b + 4
        c -= 2
        a0, b = a % 2, 3 + 3 * (a // 2)
    return a0, b, c


if __name__ == '__main__':
    seed, cmin, cmax, trials = (int(x) for x in sys.argv[1:5])
    rng = random.Random(seed)
    for cpar in (0, 1):
        cnt = Counter()
        for _ in range(trials):
            a0, b, c = rng.randrange(2), rng.randrange(1, 1 << 64), 2 * rng.randrange(cmin, cmax) + cpar
            a0, b, c = countdown(a0, b, c)
            cnt[a0 == 0 and ((c == 0 and b % 2 == 1) or (c == 1 and b % 2 == 0))] += 1
        print(f'c parity {cpar}: halting fraction {cnt[True] / trials:.4f} over {trials} countdowns '
              f'(c/2 in [{cmin},{cmax}))', flush=True)
