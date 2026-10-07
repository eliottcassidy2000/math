from fractions import Fraction
import sympy as sp
a = [2, 40, 544, 6912, 87552, 1094144, 13534208, 165978112, 2022215680, 24520433664, 296325611520, 3572784594944]
def fit(seq, d, start=0):
    # a_n = sum_{i=1..d} c_i a_{n-i}, use equations from n=start+d .. start+2d-1, verify rest
    rows = []; rhs = []
    for n in range(start + d, start + 2 * d):
        rows.append([seq[n - i] for i in range(1, d + 1)]); rhs.append(seq[n])
    M = sp.Matrix(rows); b = sp.Matrix(rhs)
    if M.det() == 0: return None
    c = M.LUsolve(b)
    ok = all(sum(c[i - 1] * seq[n - i] for i in range(1, d + 1)) == seq[n] for n in range(start + d, len(seq)))
    return list(c), ok
for d in range(1, 6):
    for st in range(0, 3):
        r = fit(a, d, st)
        if r and r[1]:
            print("order", d, "start", st, "coeffs", r[0])
            x = sp.symbols('x'); print("  charpoly", sp.factor(x**d - sum(r[0][i]*x**(d-1-i) for i in range(d))))
