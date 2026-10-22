# n = 25 complete class table (circ25_H.txt): additive energy E(S) extremes against the H-maximum.  mac-mini fanopsp.
import re
from math import gcd
n = 25; m = 12
H = {}
for line in open('circ25_H.txt'):
    mm = re.match(r'S=\{([\d,]+)\} H=(\d+)', line)
    if mm: H[tuple(int(x) for x in mm.group(1).split(','))] = int(mm.group(2))
def E(S):
    Ss = set(S); return sum(1 for a in S for b in S for c in S if (a + b - c) % n in Ss)
data = sorted((E(S), h, S) for S, h in H.items())
R = tuple(range(1, m + 1))
print(f"n=25: {len(data)} classes; E range [{data[0][0]}, {data[-1][0]}]")
print(f"  E-max: E={data[-1][0]} S={list(data[-1][2])} {'[R_25]' if data[-1][2] == R else ''} (next E={data[-2][0]}); H={data[-1][1]}")
print(f"  E-min: E={data[0][0]} S={list(data[0][2])} (next E={data[1][0]}); H={data[0][1]}")
hb = max(data, key=lambda d: d[1])
print(f"  H-max: E={hb[0]} S={list(hb[2])}; H-max is the E-max: {hb[2] == data[-1][2]}")
