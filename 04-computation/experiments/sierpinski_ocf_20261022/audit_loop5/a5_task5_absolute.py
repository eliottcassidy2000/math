# Absolute z^10 coefficients of X_3 and X_4 for every DRT(23) class: own exact values (a5_task5_z10_exact.json) vs the
# engine's DRT-mode values (drt23_five_loop.json 'coeffs'): a check of the universal (class-independent) part that the
# identity does not test.  Output: a5_task5_absolute.out.  Run from audit_loop5/.
import json
from fractions import Fraction
own = json.load(open('a5_task5_z10_exact.json'))
eng = json.load(open('../drt23_five_loop.json'))
worst = {3: 0, 4: 0}
for k, lab in enumerate(eng['labels']):
    for e in (3, 4):
        x = Fraction(*own[lab][str(e)])
        d = abs(float(x) - eng[f'X{e}'][k]['coeffs'][10])
        worst[e] = max(worst[e], d / abs(float(x)))
x3 = Fraction(*own['P_23']['3']); x4 = Fraction(*own['P_23']['4']); x5 = Fraction(*own['P_23']['5'])
print(f"P_23 exact: X_3^(10) = {float(x3):+.15e}, X_4^(10) = {float(x4):+.15e}, X_5^(10) (all V <= 5 cores) = {float(x5):+.15e}")
print(f"max relative |own - engine| over the 14 tournaments, absolute z^10 coefficients: X_3 {worst[3]:.1e}, X_4 {worst[4]:.1e}")
for e in range(6):
    print(f"   P_23 X_{e}^(10) = {Fraction(*own['P_23'][str(e)])}")
