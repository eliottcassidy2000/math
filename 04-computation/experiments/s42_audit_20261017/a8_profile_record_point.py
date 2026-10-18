#!/usr/bin/env python3
"""S42 audit a8: scale profile E_Y(N)/N (brute force from the definition, a1.brute) at the E_100 record N = 48,251,574 and at
N = 133,970,036 (R^2/N = 7.48), to test the note's claim that every large climb happens below Y = 30-100."""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from a1_frontier_bruteforce import brute
out = {}
for N in (48251574, 133970036):
    row = {}
    for Y in (1, 3, 10, 30, 100, 300, 1000, 3000, 10000):
        t0 = time.time()
        c, E, R, rdir, bad = brute(N, Y)
        assert R == rdir and bad == 0
        row[Y] = E / N
        print('N = %d  Y = %5d: E_Y/N = %.4f  (R = %d) [%.0f s]' % (N, Y, E / N, R, time.time() - t0), flush=True)
    out[N] = row
json.dump({str(k): v for k, v in out.items()}, open(os.path.join(HERE, 'a8_profile_record_point.json'), 'w'), indent=1)
