# Margin of the Paley peak: H(P_q with the difference class {1, -1} reversed) / H(P_q), from the complete circulant tables
# (q <= 23) and P31_oneflip.out (q = 31).  All one-flip neighbours of P_q are isomorphic (residue multipliers act simply
# transitively on the (q-1)/2 difference classes).  mac-mini-2026-10-21-fanopsp.
import re, glob
from math import gcd
from fractions import Fraction
rows = {}
for f in glob.glob('../tournament_pseudoprimes_20261021/circ_*.out'):
    for line in open(f):
        m_ = re.match(r'n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+)', line)
        if m_: rows.setdefault(int(m_.group(1)), {})[frozenset(int(x) for x in m_.group(2).split(','))] = int(m_.group(4))
def canon(n, S):
    return min(tuple(sorted((a * s) % n for s in S)) for a in range(1, n) if gcd(a, n) == 1)
for q in (7, 11, 19, 23):
    H = {canon(q, S): h for S, h in rows[q].items()}
    QR = {x * x % q for x in range(1, q)}
    flip = (QR - {1}) | {q - 1}
    hP, hF = H[canon(q, QR)], H[canon(q, flip)]
    print(f"q={q}: H(P_q) = {hP}, H(P_q flipped) = {hF}, margin = {hF / hP:.6f}")
hP31 = 19876729650849384218015865
hF31 = int(open('P31_oneflip.out').read().split('H=')[1])
print(f"q=31: H(P_q) = {hP31}, H(P_q flipped) = {hF31}, margin = {hF31 / hP31:.6f}")
