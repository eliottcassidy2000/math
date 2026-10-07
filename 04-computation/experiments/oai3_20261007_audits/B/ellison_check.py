# Ellison (1971), via Waldschmidt arXiv:0908.4031 p.6-7: |2^x - 3^y| > 2^x e^{-x/10} for x >= 12, x != 13,14,16,19,27, all y.
# Exact check: for each x, the minimising y is floor or ceil of x/log2(3); compare |2^x-3^y| e^{x/10} vs 2^x using
# exact integers and a rigorous rational enclosure of e^{x/10} (via mpmath interval arithmetic).
import mpmath
from mpmath import iv
iv.dps = 60; mpmath.mp.dps = 60
exc = []
L = mpmath.log(3)/mpmath.log(2)
p3 = {}
for x in range(12, 20001):
    y0 = int(mpmath.floor(x / L))
    two = 1 << x
    best = None
    for y in (y0 - 1, y0, y0 + 1, y0 + 2):
        if y < 0: continue
        d = abs(two - 3**y)
        if best is None or d < best: best = d
    # exception iff best <= 2^x e^{-x/10}  <=>  best * e^{x/10} <= 2^x
    E = iv.exp(iv.mpf(x)/10)
    lo = iv.mpf(best) * E
    if lo.b <= two: exc.append(x)            # certainly <=
    elif lo.a <= two: print("UNDECIDED at", x)  # enclosure straddles
print("exceptions in [12, 20000]:", exc)
# floating scan further, 20000 < x <= 10^6 (log-scale; margin printed)
mpmath.mp.dps = 40
worst = None
for x in range(20001, 1000001):
    y0 = mpmath.floor(x / L)
    for y in (y0, y0 + 1):
        dlt = y*L - x
        r = abs(1 - mpmath.power(2, dlt))
        ratio = r / mpmath.exp(-mpmath.mpf(x)/10)
        if worst is None or ratio < worst[0]: worst = (ratio, x, int(y))
print("min over 20000<x<=1e6 of |1-3^y/2^x| / e^{-x/10}:", mpmath.nstr(worst[0], 5), "at x =", worst[1])
# THM-4512 closure constants: for j >= 65, A = bitlen(3^j), N(w)/2^A < 2^{-j} e^{A/10} <= e^{0.1} e^{-(ln2 - log2(3)/10) j}
c = mpmath.log(2) - (mpmath.log(3)/mpmath.log(2))/10
print("exponent constant ln2 - log2(3)/10 =", mpmath.nstr(c, 8), " (THM writes 0.535)")
print("bound at j=65: e^{0.1} e^{-c*65} =", mpmath.nstr(mpmath.e**0.1 * mpmath.e**(-c*65), 6))
# exact worst case over 65 <= j <= 3000 of 2^{-j} e^{A/10}
mx = max((mpmath.mpf(2)**(-j) * mpmath.e**(mpmath.mpf((3**j).bit_length())/10), j) for j in range(65, 3001))
print("max_{65<=j<=3000} 2^-j e^(A/10) =", mpmath.nstr(mx[0], 6), "at j =", mx[1])
# does '1.11 e^{-0.535 j}' bound 1.10517 e^{-c j}? ratio grows like e^{(0.535-c) j}
for j in (65, 1000, 10**4):
    print(f"  j={j}: [2^-j e^(A/10)] / [1.11 e^(-0.535 j)] =", mpmath.nstr(mpmath.mpf(2)**(-j) * mpmath.e**(mpmath.mpf((3**j).bit_length())/10) / (mpmath.mpf('1.11')*mpmath.e**(-mpmath.mpf('0.535')*j)), 6))
