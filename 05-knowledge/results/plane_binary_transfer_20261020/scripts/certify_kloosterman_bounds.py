"""Exact rational enclosures: Machin pi and Taylor cos, no floating proof steps."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json

def atan_inv(n, terms=40):
    s = sum((F((-1) ** k, (2 * k + 1) * n ** (2 * k + 1)) for k in range(terms)), F(0))
    nxt = F((-1) ** terms, (2 * terms + 1) * n ** (2 * terms + 1))
    return (min(s, s + nxt), max(s, s + nxt))
a, b = atan_inv(5)
c, d = atan_inv(239)
PI = (16 * a - 4 * d, 16 * b - 4 * c)

def cos_frac(n, p, N=24):
    n = min(n % p, -n % p)
    lo, hi = (2 * n * PI[0] / p, 2 * n * PI[1] / p)
    L = H = F(0)
    for k in range(N + 1):
        a, b = (lo ** (2 * k) / factorial(2 * k), hi ** (2 * k) / factorial(2 * k))
        if k % 2:
            L -= b
            H -= a
        else:
            L += a
            H += b
    r = hi ** (2 * N + 2) / factorial(2 * N + 2)
    return (L - r, H + r)

def kl(p, c):
    pairs = [cos_frac(x + c * pow(x, -1, p), p) for x in range(1, p)]
    return (sum((a for a, b in pairs)), sum((b for a, b in pairs)))

def dec_interval(bounds, den=10 ** 12):
    a, b = bounds
    lo = a.numerator * den // a.denominator
    hi = -(-b.numerator * den // b.denominator)
    return [str(F(lo, den)), str(F(hi, den))]
out = {'method': 'Exact Fractions; pi=16 atan(1/5)-4 atan(1/239), 40 terms each; cos degree 48 and Lagrange remainder', 'p17_c2': dec_interval(kl(17, 2))}
if not kl(17, 2)[0] > 6:
    raise RuntimeError('Exact verification gate failed: scripts/certify_kloosterman_bounds.py:36')
bounds23 = {c: kl(23, c) for c in range(1, 23)}
if not all((b < 8 for a, b in bounds23.values())):
    raise RuntimeError('Exact verification gate failed: scripts/certify_kloosterman_bounds.py:38')
if not bounds23[11][0] > max((b for c, (a, b) in bounds23.items() if c != 11)):
    raise RuntimeError('Exact verification gate failed: scripts/certify_kloosterman_bounds.py:39')
K = F(796068719, 10 ** 8)
if not all((b < K for a, b in bounds23.values())):
    raise RuntimeError('Exact verification gate failed: scripts/certify_kloosterman_bounds.py:41')
out.update({'p23_all_c': {str(c): dec_interval(v) for c, v in bounds23.items()}, 'p23_maximizer': 11, 'p23_upper_envelope': str(K), 'p23_all_triangles_monochromatic_density_lower': str((1 - 3 * K / 24) / 4), 'p23_placement_count': 23 ** 2 * 24, 'p23_monochromatic_placements_at_least': 16})
if not (1 - 3 * K / 24) / 4 * (23 ** 2 * 24) > 15:
    raise RuntimeError('Exact verification gate failed: scripts/certify_kloosterman_bounds.py:44')
print(json.dumps(out, indent=2, sort_keys=True))
