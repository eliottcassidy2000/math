#!/usr/bin/env python3
"""A11: countdown-end halting frequency measured directly on the paper's rule system (S56).
Random starts C(a0, b, c, d) with a0 in {0,1}, d in {3,6}, c in [2*cmin, 2*cmax); run the rule system (step) until
HALT or until the first restart of a new countdown (c jumps up); tally HALT vs restart, by c parity and d.
Usage: python a11_rule_check_20261025.py seed trials cmin cmax
"""
import sys, random
from collections import Counter
from deepmind_cryptids_20261025 import step, HALT

if __name__ == '__main__':
    seed, trials, cmin, cmax = (int(x) for x in sys.argv[1:5])
    rng = random.Random(seed)
    cnt = Counter()
    for _ in range(trials):
        a0 = rng.randrange(2); b = rng.randrange(1, 1 << 40); d = rng.choice((3, 6))
        c = 2 * rng.randrange(cmin, cmax) + rng.randrange(2)
        p = (a0, b, c, d)
        key = (c % 2, d)
        cmin_seen = c
        for _ in range(10 * c + 100):
            q = step('A11', p)
            if q == HALT:
                cnt[key + ('HALT',)] += 1
                break
            if q[2] > cmin_seen + 2:      # a new countdown started (c reset upward)
                cnt[key + ('restart',)] += 1
                break
            cmin_seen = min(cmin_seen, q[2])
            p = q
        else:
            cnt[key + ('undecided',)] += 1
    for k in sorted(cnt):
        print(k, cnt[k])
    for cp in (0, 1):
        for d in (3, 6):
            h, r = cnt[(cp, d, 'HALT')], cnt[(cp, d, 'restart')]
            if h + r:
                print(f'c parity {cp}, d = {d}: halting fraction {h / (h + r):.3f} ({h}/{h + r})')
