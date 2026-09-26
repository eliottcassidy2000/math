#!/usr/bin/env python3
"""gilbreath_certificate_20260926.py -- the wall theorem (THM-4511, any number of size-4 defects) as a
finite certificate for Gilbreath's conjecture (session gilbreath-fermat-platonic-20260926, opus, 2026-09-26).

For the primes below P, compute row by row: r_2(P) = first row whose entries after the leading 1 are all in
{0,2}; r_4(P) = first row whose entries after the leading 1 are all in {0,2,4}. By THM-4511 (multi-defect form)
the finite triangle is decided at row r_4: the leading 1 survives all later rows iff, in row r_4, the first 4
(if any) is preceded by at least one 2. The rows between r_4 and r_2 need not be computed.
Also: the light-cone certificate per row -- t(r) = number of leading cells (after column 0) in {0,2,4}; if the
first 4 among them is preceded by a 2 (or there is none), the leading 1 is safe for the next t(r) rows.
Also: the exact necessary condition 'no row has its first entry >= 4 equal to 4 with only zeros before it'
(else the conjecture fails) is checked on every computed row.
Usage: python3 gilbreath_certificate_20260926.py [P ...]
"""
import sys
import numpy as np


def primes_below(P):
    s = np.ones(P, dtype=bool); s[:2] = False
    for i in range(2, int(P ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


def main():
    Ps = [int(x) for x in sys.argv[1:]] or [200000, 1000000, 10000000]
    for P in Ps:
        row = primes_below(P)
        n0 = len(row)
        r_4 = None; r_2 = None; cert = None
        zero_prefix_rows = 0
        tcert = []
        r = 0
        while r_2 is None:
            r += 1
            row = np.abs(np.diff(row))
            assert row[0] == 1
            tail = row[1:]
            big6 = np.nonzero(tail >= 6)[0]
            big4 = np.nonzero(tail >= 4)[0]
            if len(big4) and tail[big4[0]] == 4 and not tail[:big4[0]].any():
                zero_prefix_rows += 1  # would destroy the leading 1 (never happens if the conjecture holds)
            t = int(big6[0]) if len(big6) else len(tail)  # cells 1..t are in {0,2,4}
            if len(big4) and big4[0] < t:
                safe = bool((tail[:big4[0]] == 2).any())
            else:
                safe = True
            if r <= 12 or r % 10 == 0:
                tcert.append((r, t, safe))
            if r_4 is None and len(big6) == 0:
                r_4 = r
                cert = safe
            if len(big4) == 0:
                r_2 = r
        print("primes below %d (%d primes): r_4 = %d (first all-{0,2,4} row), certificate at r_4: %s; r_2 = %d (first all-{0,2} row); rows with a zero-only prefix before a first 4: %d" % (P, n0, r_4, 'first 4 preceded by a 2 -> leading 1 safe for ever' if cert else 'NOT CERTIFIED', r_2, zero_prefix_rows))
        print("  light-cone certificates (row, t = cells in {0,2,4} from column 1, safe for t rows):", tcert[:16])


if __name__ == '__main__':
    main()
