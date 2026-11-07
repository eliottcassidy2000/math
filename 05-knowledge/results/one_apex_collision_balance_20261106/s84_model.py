"""Poisson-pair model predictions (S84) for one-apex counts on reduced words Y_n.
Limit law of g(v)/2^n for v uniform in Y_n is XY with X, Y independent uniform on [1/2,1).
Converse classes {v, v*} form a Poisson process of intensity mu(t) = f_Y(t)/4 per odd value,
f_Y(t) = 4 ln(4t) on [1/4,1/2], 4 ln(1/t) on [1/2,1].
Predictions per |Y|:
  values with exactly q classes: 2 * [ (ln 2)^(q+1)/(4 (q+1)!) + int_0^{ln 2} e^{-2u} u^q / q! du ]
  non-converse ordered pairs: 12 - 16 ln 2;  energy E^Y/|Y| -> 14 - 16 ln 2
  smoothed slope (E^Y(K)-E^Y(0))/(K |Y|) -> 12 - 16 ln 2 (K even)
For X_n (all words, layered by trailing ones, layers independent): energy -> 6 - 4 ln 2.
"""
import math
from scipy.integrate import quad

L2 = math.log(2)
print('energy Y: 14-16ln2 = %.6f   nonconverse per |Y| = 12-16ln2 = %.6f' % (14 - 16 * L2, 12 - 16 * L2))
print('energy X (independent layers): 6-4ln2 = %.6f' % (6 - 4 * L2))
tot = 0.0
for q in range(1, 9):
    a = L2 ** (q + 1) / (4 * math.factorial(q + 1))
    b, _ = quad(lambda u: math.exp(-2 * u) * u ** q / math.factorial(q), 0, L2)
    frac = 2 * (a + b)
    tot += frac * (2 * q) ** 2
    print('q=%d (multiplicity %2d): fraction of |Y| = %.6f' % (q, 2 * q, frac))
print('energy from histogram (check) = %.6f' % tot)
# fraction of odd values in [2^(n-2), 2^n] hit (per |Y|): 2*int (1 - e^{-mu})
a0, _ = quad(lambda t: 1 - math.exp(-math.log(4 * t)), 0.25, 0.5)
b0, _ = quad(lambda t: 1 - math.exp(math.log(t)), 0.5, 1.0)
print('values attained per |Y| (model) = %.6f' % (2 * (a0 + b0)))
