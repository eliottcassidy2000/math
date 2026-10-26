#!/usr/bin/env python3
"""Deep certified-accelerated runs of the slow Cryptids A1, A2, A5, A10 from their blank-tape start (S56).

Stops on HALT, on a certified translated cycler (CAPTURE: proof of non-halting from that configuration, hence of
the machine), on an infinite segment cycle, or on the time limit.  Progress lines every 60 s.
Usage: python rs_deep_20261025.py NAME seconds
"""
import sys, os, time, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rs_accel_20261025 as R
from deepmind_cryptids_20261025 import step, START as DSTART, HALT

if __name__ == '__main__':
    nm, secs = sys.argv[1], float(sys.argv[2])
    eng = R.Engine(nm)
    p = DSTART[nm]
    pre = 0
    while not R.ANCHOR[nm](p):
        p = step(nm, p)
        pre += 1
    t0 = time.time()
    last = [t0]
    hist = []

    def stop(s, n):
        now = time.time()
        if now - last[0] > 60:
            last[0] = now
            line = f'[{now - t0:.0f}s] segments {n:.4e} state bits {max(s).bit_length()} stats {eng.stats}'
            print(line, flush=True)
            hist.append((round(now - t0), n, max(s).bit_length()))
        return now - t0 > secs
    st = eng.orbit(p, stop=stop)
    status, s, nseg, info = st
    print('RESULT', nm, status, 'segments', nseg, 'state bits', max(s).bit_length() if isinstance(s, tuple) else None,
          'stats', eng.stats, flush=True)
    if info is not None:
        print('INFO', str(info)[:3000], flush=True)
    if status in ('CAPTURE', 'HALT', 'INFINITE-SEGMENT-CYCLE'):
        print('STATE', s, flush=True)
    json.dump({'machine': nm, 'status': status, 'segments': nseg, 'state': [str(x) for x in s] if isinstance(s, tuple) else str(s),
               'info': str(info)[:20000], 'progress': hist, 'stats': eng.stats},
              open(os.path.join(HERE, f'rs_deep_{nm}.json'), 'w'), indent=1)
