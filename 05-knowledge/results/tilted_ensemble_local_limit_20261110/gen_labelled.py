#!/usr/bin/env python3
"""Print every labelled tournament on m vertices in gentourng's format (upper triangle row by row, '1' = i->j, i<j)."""
import sys
m = int(sys.argv[1]); L = m*(m - 1)//2
w = sys.stdout.write
for c in range(1 << L):
    w((format(c, '0%db' % L) if L else '') + '\n')
