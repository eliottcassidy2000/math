#!/usr/bin/env python3
"""Independent check of a colouring certificate for G_q = Cay(F_(q^2), mu_(q+1)) written by tabucol_gq.c.
Usage: python3 verify_colouring.py FILE"""
import sys
lines = open(sys.argv[1]).read().split('\n')
hdr = lines[0]; q = int(hdr.split('q=')[1].split()[0]); r = int(hdr.split('r=')[1].split()[0])
col = list(map(int, lines[1].split()))
assert len(col) == q * q
assert pow(r, (q - 1) // 2, q) == q - 1, 'r must be a non-residue'
mu = [(a, b) for a in range(q) for b in range(q) if (a * a - r * b * b) % q == 1]
assert len(mu) == q + 1
bad = sum(1 for a in range(q) for b in range(q) for (x, y) in mu if col[a * q + b] == col[((a + x) % q) * q + (b + y) % q])
print('%s: q = %d, colours used = %d, |mu| = %d, monochromatic edges = %d -> %s' % (sys.argv[1], q, len(set(col)), len(mu), bad // 2, 'PROPER' if bad == 0 else 'NOT PROPER'))
