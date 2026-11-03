#!/usr/bin/env python3
"""Exact certificate: nine square entries, eight square line sums, not magic."""
from math import isqrt

roots = ((2655, 2976, 320), (2624, 2520, 1665), (1440, 895, 3624))
grid = tuple(tuple(x*x for x in row) for row in roots)
flat = sum(grid, ())
assert len(set(flat)) == 9 and min(flat) > 0
lines = list(grid) + list(zip(*grid)) + [tuple(grid[i][i] for i in range(3)), tuple(grid[i][2-i] for i in range(3))]
totals = tuple(sum(line) for line in lines)
assert totals == (4001**2,)*6 + (5151**2, 2920**2)
assert all(isqrt(s)**2 == s for s in totals)
assert len(set(totals)) == 3  # This is expressly not a magic square.

# Independently verify the quaternion derivation without floating-point arithmetic.
a,b,c,d = 48,32,23,12
N = a*a+b*b+c*c+d*d
M = ((a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)),
     (2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)),
     (2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d))
assert N == 4001
assert all(sum(M[i][k]*M[j][k] for k in range(3)) == (N*N if i==j else 0)
           for i in range(3) for j in range(3))
assert roots == tuple(tuple(abs(row[j]) for j in (0,2,1)) for row in M)
print('PASS: all nine entries are distinct positive integer squares.')
for row in grid: print(*row)
print('Row, column, and diagonal sum roots:', *(isqrt(s) for s in totals))
print('Not magic: the eight line sums are not equal.')
