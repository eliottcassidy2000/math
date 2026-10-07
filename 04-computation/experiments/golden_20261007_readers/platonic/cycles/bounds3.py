# Fast version: bounds via log form with mpmath only at candidate maxima.
# For each L: kneg = smallest k with 3^k > 2^L, kpos = kneg-1.
# eps_neg = kneg*log2(3) - L > 0 ; bound_neg = 1/(3(1 - 2^(-eps_neg/kneg)))
# eps_pos = L - kpos*log2(3) > 0 ; bound_pos = 1/(3(2^(eps_pos/kpos) - 1))  (needs L <= 2 kpos)
from mpmath import mp, mpf, log, power, floor
import sys
mp.dps = 50
LMAX = int(sys.argv[1])
a = log(3)/log(2)
bestN=(0,None); bestP=(0,None)
for L in range(1, LMAX+1):
    kneg = int(floor(L/a)) + 1
    eN = kneg*a - L
    bN = 1/(3*(1 - power(2, -eN/kneg)))
    if bN > bestN[0]: bestN=(bN,(L,kneg))
    kpos = kneg - 1
    if kpos >= 1 and L <= 2*kpos:
        eP = L - kpos*a
        if eP > 0:
            bP = 1/(3*(power(2, eP/kpos) - 1))
            if bP > bestP[0]: bestP=(bP,(L,kpos))
print(LMAX, "neg:", mp.nstr(bestN[0],15), bestN[1], " pos:", mp.nstr(bestP[0],15), bestP[1])
