"""
graph_carriers_line_classification_20261013.py -- S30 companion to graph_carriers_level_factorisation_20261013.py.
Where do the lines of THM-4629 sit in the Farey graph?  The m = 2 census of THM-4629 is the set of pairs (Q, P) with
|Q theta - P| < 0.07, Q <= 768, (Q, P) not both even (collatz_parity_census_m2_m3_20261013.out, header of its section
m = 2).  Each line is classified as: a convergent of theta = log2 3; a semiconvergent (intermediate fraction
(P_(n-1) + t P_n)/(Q_(n-1) + t Q_n), 1 <= t < a_(n+1)); a non-reduced multiple j (Q0, P0) of a primitive pair; or another
reduced fraction, in which case we record whether it is a Farey neighbour (|P Q' - P' Q| = 1) of some convergent or
semiconvergent with Q' <= 10^7.  Also: the characters (m = 2) of census lines consecutive in Q, the fan around 84/53,
and the record lines (minimal |Q theta - P| up to Q) against the convergents (Lagrange).
"""
import math
from mpmath import mp, log, floor

mp.dps = 120
theta = log(3) / log(2)
TH = float(theta)

# continued fraction and convergents up to Q <= 1e7
a = []; x = theta
for _ in range(80):
    fl = int(floor(x)); a.append(fl); x = 1 / (x - fl)
P = [1, a[0]]; Q = [0, 1]
for n in range(1, len(a)):
    P.append(a[n] * P[-1] + P[-2]); Q.append(a[n] * Q[-1] + Q[-2])
conv = [(Q[i + 1], P[i + 1]) for i in range(len(a))]          # (Q_n, P_n), n = 0, 1, ...
conv = [c for c in conv if c[0] <= 10 ** 7]
semi = set()
for n in range(1, len(conv) - 1):
    (q0, p0), (q1, p1) = conv[n - 1], conv[n]
    an1 = a[n + 1]
    for t in range(1, an1):
        semi.add((q0 + t * q1, p0 + t * p1))
semi = {s for s in semi if s[0] <= 10 ** 7}
convset = set(conv)
geo = sorted(convset | semi)
print(f"# convergents with Q <= 1e7: {len(convset)}; semiconvergents: {len(semi)}; geodesic vertices: {len(geo)}")
print(f"  convergents (Q, P): {conv[:12]}")

# the census
census = []
for Qv in range(1, 769):
    Pv = int(round(Qv * TH))
    if abs(Qv * TH - Pv) < 0.07 and not (Qv % 2 == 0 and Pv % 2 == 0):
        census.append((Qv, Pv))
print(f"\n# m = 2 census: |Q theta - P| < 0.07, Q <= 768, (Q, P) not both even: {len(census)} lines")

def neighbour_of_geodesic(q, p):
    return [(q2, p2) for (q2, p2) in geo if abs(p * q2 - p2 * q) == 1]

cls = {"convergent": [], "semiconvergent": [], "non-reduced": [], "other reduced": []}
nbr_yes, nbr_no = [], []
for (q, p) in census:
    g = math.gcd(q, p)
    if g > 1:
        cls["non-reduced"].append(((q, p), g, (q // g, p // g)))
    elif (q, p) in convset:
        cls["convergent"].append((q, p))
    elif (q, p) in semi:
        cls["semiconvergent"].append((q, p))
    else:
        cls["other reduced"].append((q, p))
        nb = neighbour_of_geodesic(q, p)
        (nbr_yes if nb else nbr_no).append(((q, p), nb[:3]))
print(f"  convergents: {len(cls['convergent'])} {cls['convergent']}")
print(f"  semiconvergents: {len(cls['semiconvergent'])} {cls['semiconvergent']}")
print(f"  non-reduced multiples: {len(cls['non-reduced'])} {[(c[0], f'{c[1]}x{c[2]}') for c in cls['non-reduced']]}")
print(f"  other reduced fractions: {len(cls['other reduced'])}; Farey neighbours of a geodesic vertex (Q' <= 1e7): {len(nbr_yes)}; of none: {len(nbr_no)}")
print(f"    neighbours of none: {[c[0] for c in nbr_no]}")
print(f"    examples of neighbours: {[(c[0], c[1][0]) for c in nbr_yes[:8]]}")

# fan around 84/53: (41 + 53 t, 65 + 84 t)
fan = [(41 + 53 * t, 65 + 84 * t) for t in range(1, 14)]
names = {(1, 0): "tau", (0, 1): "omega", (1, 1): "C"}
print(f"\n# fan (41 + 53t, 65 + 84t), t = 1..13: in census? {[((q, p) in set(census)) for (q, p) in fan]}")
print(f"  characters: {[names[(q % 2, p % 2)] for (q, p) in fan]} (alternate, as Farey neighbours must)")

# consecutive census lines sharing a character
chars = [names[(q % 2, p % 2)] for (q, p) in census]
same = [(census[i], census[i + 1], chars[i]) for i in range(len(census) - 1) if chars[i] == chars[i + 1]]
print(f"\n# census lines consecutive in Q with the same character: {len(same)} of {len(census) - 1}: {same}")

# record lines = convergents (Lagrange)
records = []; best = 1.0
for Qv in range(1, 20001):
    d = abs(Qv * TH - round(Qv * TH))
    if d < best:
        best = d; records.append((Qv, int(round(Qv * TH))))
print(f"\n# record lines (minimal |Q theta - P| up to Q, Q <= 20000): {records}")
print(f"  all convergents? {all(r in convset for r in records)}; all convergents with Q <= 20000 are records? {all(c in records for c in conv if c[0] <= 20000 and c[0] >= 1)}")
print("done")
