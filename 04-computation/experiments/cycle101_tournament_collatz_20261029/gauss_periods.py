"""THM-4849 (F): the binary values c_S = sum_{j in S} 2^j of the two Paley base points satisfy
c_QR + c_NQR = -1 and c_QR * c_NQR = (q+1)/4 modulo 2^q - 1, so both are roots of X^2 + X + (q+1)/4 there
(the Gauss periods evaluated at the q-th root of unity 2).  Also the trace relation T-trace = 2 U-trace + k on the
five (3,7) cycles."""
from fractions import Fraction
for q in (7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83):
    QR = sorted({x * x % q for x in range(1, q)}); N = sorted(set(range(1, q)) - set(QR))
    M = 2 ** q - 1
    cS = sum(2 ** j for j in QR); cN = sum(2 ** j for j in N)
    ok = (cS + cN) % M == M - 1 and (cS * cN - (q + 1) // 4) % M == 0 and all((c * c + c + (q + 1) // 4) % M == 0 for c in (cS, cN))
    print('q=%d  Gauss-period identity mod 2^q - 1: %s' % (q, ok))
def T(x):
    return x / 2 if x.numerator % 2 == 0 else (3 * x + 1) / 2
for start in (19, 23, 29, 31, 37):
    x0 = Fraction(start, 101); x = x0; pts = []
    while True:
        pts.append(x); x = T(x)
        if x == x0:
            break
    odd = [p for p in pts if p.numerator % 2]
    print('cycle of %d/101: T-length %d, T-trace %s, U-trace %s, 2*U-trace + k = %s' % (start, len(pts), sum(pts), sum(odd), 2 * sum(odd) + len(odd)))
