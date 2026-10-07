# Print the integer cycles of 3x+d (Terras form x/2 | (3x+d)/2) realised by Ellison shapes with small d.
from math import gcd
exec(open('nearcyc.py').read().split("for (L, k) in")[0])
def cyc_min(x0, d, L):
    xs = [x0]; x = x0
    for _ in range(L - 1):
        x = x // 2 if x % 2 == 0 else (3 * x + d) // 2
        xs.append(x)
    return min(xs, key=abs), xs
for (L, k, dd) in [(27, 17, 5), (27, 17, 71), (19, 12, 23)]:
    m = 2**L - 3**k
    mins = set()
    for w in lyndon_words_k(L, k):
        c = cw(w)
        if (abs(m) // gcd(c, abs(m))) == dd:
            x = c * dd // m  # integer point of the 3x+dd cycle
            mins.add(cyc_min(x, dd, L)[0])
    print("shape (%d,%d): cycles of 3x+%d with min elements %s" % (L, k, dd, sorted(mins)))
