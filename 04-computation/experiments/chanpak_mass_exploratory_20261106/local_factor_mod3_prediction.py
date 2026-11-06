# local_factor_mod3_prediction.py (opus-2026-11-06-S82): predicted mean of r2 over odd N = r (mod 3), r = 0, 1, 2,
# for the order-<=7 bank, assuming g(v) mod 3 is asymptotically uniform and independent of (v mod 3, parity, rho)
# (Session 193 modular_sync Theorems 3-4). Triples with 3 | y_h and 3 | y_R give Phi = -y_e e(v) (mod 3): v odd (rho < 1/2)
# lands in class 0 and v even (rho >= 1/2) in class -y_e; all other triples spread evenly over the three classes.
import math, sys
fn = sys.argv[1] if len(sys.argv) > 1 else 'rooted_triple_bank_order7_20261104.txt'
Y = [tuple(map(int, l.split())) for l in open(fn)]
base = 0.0; extra = [0.0, 0.0, 0.0]
for a, b, c in Y:
    if a % 3 == 0 and c % 3 == 0:
        extra[0] += 3 * (2.0 / a) * math.log(1 + a / (2.0 * c))
        extra[(-b) % 3] += 3 * (2.0 / a) * math.log((c + a) / (c + a / 2.0))
    else:
        base += (2.0 / a) * math.log(1 + a / c)
pred = [base + extra[r] for r in range(3)]
print('predicted mean r2 by N mod 3 (classes 0, 1, 2):', ['%.4f' % x for x in pred])
print('average of the three classes: %.4f (the mass M_Y)' % (sum(pred) / 3))
