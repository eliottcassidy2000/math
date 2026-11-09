"""Exact path spectra S(n) for n <= 11 (klein-2026-11-06 census, FINITE-EXACT), as membership oracles.

S(n) = {H(T) : |T| = n}.  For n <= 8 the set is read from census_n.txt ("T h ..." lines).
For n = 9, 10, 11: S(n) = odd h <= max H(n), minus {7, 21}, minus the listed unattained values.
For n >= 12 only the inclusion S(11) <= S(n) is known; in_S returns None (unknown) outside S(11).
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
KLEIN = os.path.join(HERE, '..', '..', '..', '05-knowledge', 'results',
                     'path_spectrum_two_angles_metagraph_20261106', 'data')
MAXH = {1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189, 8: 661, 9: 3357, 10: 15745, 11: 95095}

_small = {}
_holes = {}


def _load_small(n):
    if n in _small:
        return _small[n]
    if n <= 2:
        s = {1}
    else:
        s = set()
        with open(os.path.join(KLEIN, 'census_%d.txt' % n)) as f:
            for line in f:
                if line.startswith('T '):
                    s.add(int(line.split()[1]))
    _small[n] = s
    return s


def _load_holes(n):
    if n in _holes:
        return _holes[n]
    s = set()
    with open(os.path.join(HERE, 'unattained_%d.txt' % n)) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                s.add(int(line))
    _holes[n] = s
    return s


def in_S(n, h):
    """True/False if decided by the exact spectra, None if unknown (n >= 12 and h outside S(11))."""
    if n <= 0:
        return False
    if n <= 8:
        return h in _load_small(n)
    if n <= 11:
        return (h % 2 == 1 and 1 <= h <= MAXH[n] and h not in (7, 21) and h not in _load_holes(n))
    if in_S(11, h):
        return True
    return None


if __name__ == '__main__':
    for n in range(1, 12):
        print(n, sum(1 for h in range(1, MAXH[n] + 1) if in_S(n, h)))
