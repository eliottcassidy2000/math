from mpmath import mp, mpf, power
mp.dps = 80
import sys
LMAX = int(sys.argv[1])
bestN = (mpf(0), None); bestP=(mpf(0),None)
p3 = 1; k = 0
for L in range(1, LMAX + 1):
    twoL = 1 << L
    # smallest k with 3^k > 2^L (negative side)
    while p3 <= twoL:
        p3 *= 3; k += 1
    r = mpf(twoL) / mpf(p3)
    b = 1 / (3 * (1 - power(r, mpf(1) / k)))
    if b > bestN[0]: bestN = (b, (L, k))
    # largest k with 3^k < 2^L (positive side) is k-1 (if 3^(k-1) < 2^L, and 2^L <= 4^(k-1))
    kk = k - 1; q = p3 // 3
    if kk >= 1 and q < twoL and twoL <= (1 << (2*kk)):
        r = mpf(twoL)/mpf(q); b = 1/(3*(power(r, mpf(1)/kk)-1))
        if b > bestP[0]: bestP = (b, (L, kk))
print(LMAX, "neg:", mp.nstr(bestN[0], 15), bestN[1], " pos:", mp.nstr(bestP[0], 15), bestP[1])
