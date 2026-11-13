"""Exact small-prefix cross-check; the exponential rate is proved separately."""
from collections import Counter
from fractions import Fraction

for K in range(2, 13):
    c = {0: 1}
    for j in range(1, K+1):
        d = Counter()
        for s, m in c.items():
            for t in (s, s+1):
                if 3**t >= 2**j:
                    d[t] += m
        c = d
    count = 0
    words = set()
    for r in range(1, 2**K+1):
        x, a, b, threshold, w = r, 1, 0, None, []
        for j in range(1, K+1):
            z = x & 1
            w.append(z)
            x = (3*x+1)//2 if z else x//2
            if z:
                a, b = 3*a, 3*b + 2**(j-1)
            if a < 2**j:
                h = Fraction(b, 2**j-a)
                threshold = h if threshold is None else min(h, threshold)
        words.add(tuple(w))
        count += threshold is None
        for n in (r, r+2**K, r+3*2**K):
            if n <= 1:
                continue
            x = n
            for j in range(K):
                x = (3*x+1)//2 if x & 1 else x//2
                if x < n:
                    break
            if (x < n) != (threshold is not None and n > threshold):
                raise ValueError((K,n,x,threshold))
            if x >= n and not (n % 4 == 3 and (3*n-1)//4<n):
                raise ValueError('invalid repair')
    if len(words)!=2**K or count!=sum(c.values()):
        raise ValueError('parity word census')
print('K=2..12: parity cylinders, exact descent thresholds, and repair controls passed')
