"""Independent PARI check of the order-5 strict-Fricke counterexample (audit S53).

For y in Lambda (coordinates in the LLL Leech basis of audit_common), v = y/5. Build the lattice L = Lambda + Z v,
scale its Gram matrix by 25 (integral), enumerate ALL its vectors of norm <= 2 with PARI qfminim (exact for integral
forms), and sort them by coset k (x = lambda + k v mod Lambda, read off from the last coordinate mod 5).
Prints the norm profile of Lambda + k v for k = 1, 2 (k = 3, 4 are the negatives)."""
import sys
from fractions import Fraction
from cypari import pari
from audit_common import GLAM

EX = {
    'c-like': [-4, 5, -1, 0, 8, -3, -4, 4, 1, 4, 2, 12, 0, -1, -3, -2, -6, 7, -8, -5, 3, 2, 6, 3],
    'irrational': [-4, 5, -3, 1, 10, -3, -4, 6, 2, 3, 3, 14, 1, 0, -3, -2, -7, 9, -10, -6, 4, 2, 6, 4],
}
P = 5
for name, y in EX.items():
    G = [[Fraction(x) for x in row] for row in GLAM]
    yy = sum(y[i] * G[i][j] * y[j] for i in range(24) for j in range(24))
    # generators of L = Lambda + Z v: e_1..e_24 and v = y/5 -> 25 generators, rank 24
    # basis via HNF of the generator matrix in coordinates (e_i ; y/5), scaled by 5
    gensM = pari.matrix(24, 25, [ (5 if i == j else 0) if j < 24 else y[i] for i in range(24) for j in range(25)])
    H = pari.mathnf(gensM)                        # columns: basis of 5L in Lambda-coordinates
    B = [[int(H[i, j]) for i in range(24)] for j in range(24)]   # 5 * basis vectors, Lambda coords
    Gram = [[sum(B[a][i] * GLAM[i][j] * B[b][j] for i in range(24) for j in range(24)) for b in range(24)] for a in range(24)]
    # Gram of 5L = 25 * Gram(L); integral
    GP = pari.matrix(24, 24, [x for row in Gram for x in row])
    res = pari.qfminim(GP, 50, 10 ** 6, 0)        # vectors of L with norm <= 2  <=> 25 x^2 <= 50
    cnt, V = int(res[0]), res[2]
    prof = {}
    for j in range(cnt // 2):
        c = [int(V[i, j]) for i in range(24)]
        x5 = [sum(c[a] * B[a][i] for a in range(24)) for i in range(24)]   # 5x in Lambda coords
        nr = Fraction(sum(x5[i] * GLAM[i][k] * x5[k] for i in range(24) for k in range(24)), 25)
        # coset: 5x = 5 lambda + k y  =>  k determined mod 5 by x5 - k y in 5 Z^24
        k = next(k for k in range(5) if all((x5[i] - k * y[i]) % 5 == 0 for i in range(24)))
        kk = min(k, (5 - k) % 5)                   # +-x pairs: x in coset k, -x in coset -k
        prof.setdefault(kk, {}).setdefault(nr, 0)
        prof[kk][nr] += 2 if k == 0 else 1
        if k != 0:
            prof[kk][nr] += 1                      # the partner -x lies in coset -k, same norm profile class
    print(f'{name}: y^2 = {yy}; vectors of Lambda + Zv with norm <= 2: {cnt}')
    for kk in (1, 2):
        d = prof.get(kk, {})
        # each class k and -k counted together above: divide by 2 to get the profile of the single coset k
        print(f'   coset {kk} v: ' + ', '.join(f'{nr}: {c // 2}' for nr, c in sorted(d.items())))
