"""AUDIT G1 (independent): extract S(n) and per-value class counts from klein's census T-lines.

Census line format (klein README): 'T h c3 strong count'.  S(n) = set of h over all T-lines of census_n.
Writes S_<n>.txt (one value per line) for n = 1..11 into this directory, and prints per-order totals,
first-creator class counts (h not in S(n-1)) and shares.  Shares no code with the session's spectra.py/fc_share.py.
"""
import gzip
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, '..', '..', 'path_spectrum_two_angles_metagraph_20261106', 'data'))


def read_T(n):
    fn = os.path.join(DATA, 'census_%d.txt' % n)
    if os.path.exists(fn):
        f = open(fn, 'rt')
    else:
        f = gzip.open(fn + '.gz', 'rt')
    per = {}
    header = None
    with f:
        for line in f:
            if line.startswith('# order'):
                header = line.strip()
            if not line.startswith('T '):
                continue
            t = line.split()
            assert len(t) == 5, line
            h, cnt = int(t[1]), int(t[4])
            per[h] = per.get(h, 0) + cnt
    return per, header


def main():
    S = {1: {1}, 2: {1}}
    per = {}
    for n in range(3, 12):
        per[n], hdr = read_T(n)
        S[n] = set(per[n])
        print('order %2d: %s' % (n, hdr))
    for n in range(1, 12):
        with open(os.path.join(HERE, 'S_%d.txt' % n), 'w') as g:
            for h in sorted(S[n]):
                g.write('%d\n' % h)
    print()
    print('n  classes  FC-classes  share  |S(n)|  maxH  new-values  S(n-1)<=S(n)')
    for n in range(3, 12):
        tot = sum(per[n].values())
        fc = sum(c for h, c in per[n].items() if h not in S[n - 1])
        newv = len(S[n] - S[n - 1])
        print('%2d %10d %10d %.4f %6d %6d %6d %s' % (n, tot, fc, fc / tot, len(S[n]), max(S[n]), newv,
                                                    S[n - 1] <= S[n]))
    # solid prefixes: largest x such that every odd h <= x except 7, 21 is in S(n)
    for n in (9, 10, 11):
        x = 1
        while True:
            h = x + 2
            if h in (7, 21) or h in S[n]:
                x = h
            else:
                break
        holes = [h for h in range(1, max(S[n]) + 1, 2) if h not in S[n] and h not in (7, 21)]
        print('n=%d: every odd h <= %d except 7,21 in S(n); %d odd non-values below max (least %s); evens in S: %d'
              % (n, x, len(holes), holes[0] if holes else None, sum(1 for h in S[n] if h % 2 == 0)))


if __name__ == '__main__':
    main()
