# S43 referee: the congruence/integrality test of section 2.4 cannot pin the cusp coefficients:
# every c in -22/5 + 11Z (11A) and (c1,c2) in (-46/11,-23/11) + 23Z^2 (23A) passes.
import sys
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from audit_k3_twining import H_series, phi2, eta_prod, shift, mul
N = 60
H, ie3 = H_series(N)
A = [x / 2 for x in H]
def traces(chi, T):
    t = mul(T, ie3, N)
    return [Fr(chi, 24) * H[i] - t[i] for i in range(N + 1)]   # = Tr(g|K_n) = 2 A_n(g)
f11 = shift(eta_prod({1: 2, 11: 2}, N), 1, N)
base11 = [Fr(11, 6) * x for x in phi2(11, N)]
for c in (Fr(-22, 5), Fr(-22, 5) + 11, Fr(-22, 5) - 11, Fr(-22, 5) + 110):
    T = [base11[i] + c * f11[i] for i in range(N + 1)]
    tr = traces(2, T)
    ok = all(tr[n].denominator == 1 and (tr[n] - 2 * A[n]) % 11 == 0 for n in range(1, N + 1))
    print('11A  c = %-8s integral & == dim mod 11 for n <= %d: %s' % (c, N, ok))
u = shift(eta_prod({1: 2, 23: 2}, N), 2, N)
e1 = shift(eta_prod({1: 1, 23: 1}, N), 1, N)
th = [Fr(0)] * (N + 1)
for x in range(-20, 21):
    for y in range(-10, 11):
        v = x * x + x * y + 6 * y * y
        if v <= N:
            th[v] += 1
w = mul(e1, th, N)
base23 = [Fr(23, 12) * x for x in phi2(23, N)]
for c1, c2 in ((Fr(-46, 11), Fr(-23, 11)), (Fr(-46, 11) + 23, Fr(-23, 11)), (Fr(-46, 11), Fr(-23, 11) - 23), (Fr(-46, 11) - 23, Fr(-23, 11) + 46)):
    T = [base23[i] + c1 * u[i] + c2 * w[i] for i in range(N + 1)]
    tr = traces(1, T)
    ok = all(tr[n].denominator == 1 and (tr[n] - 2 * A[n]) % 23 == 0 for n in range(1, N + 1))
    print('23A  (c1, c2) = (%s, %s) integral & == dim mod 23 for n <= %d: %s' % (c1, c2, N, ok))
