"""Verify quoted representations exactly (sympy, omega = exp(2 pi i/3))."""
import sympy as sp
w = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2
def val(reps):
    return sp.nsimplify(sp.expand(sum((a + b * w) ** 4 for a, b in reps)))
quoted = {
 '47 (note/HYP)': [(-20, -10), (-13, -17), (-8, -15), (-8, 6)],
 '47 (small box, audit)': [(4, 2), (2, -2), (4, 1), (3, -1)],
 '102 (note/HYP)': [(113, 99), (51, -41), (104, 31), (106, 13)],
 '102 (second, audit)': [(92, 41), (14, -99), (73, -31), (93, -13)],
 '186 (note/HYP)': [(37, -44), (16, 40), (54, 36), (77, 14)],
 '186 (audit)': [(18, -36), (24, 40), (81, 44), (63, -14)],
 '30 (HYP window)': [(2, -3), (1, -1), (3, 1), (2, 5)],
 '22 (box validation)': [(0, 1), (6, 7), (1, 1), (1, 7)],
 '38 (box validation)': [(4, -1), (2, -3), (5, 1), (5, 3)],
}
for k, v in quoted.items():
    print(f'{k:24s} -> {val(v)}')
# negative control: perturb one coordinate
print('perturbed 102 (113 -> 112) ->', val([(112, 99), (51, -41), (104, 31), (106, 13)]))
