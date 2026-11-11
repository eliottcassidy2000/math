# Test the tower: sum_{j=1..k} tr(A^(2j+1)) = C(|A^k 1|, 2) + floor(n/4) (mod 2) for tournament matrices over F_2?
import random
from math import comb
rng = random.Random(11)
def rand_rows(n):
    row = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5: row[i] |= 1 << j
            else: row[j] |= 1 << i
    return row
def mul(row, M, n):  # rows of A*M over F_2
    out = []
    for a in range(n):
        acc = 0; r = row[a]
        while r:
            b = (r & -r).bit_length() - 1; acc ^= M[b]; r &= r - 1
        out.append(acc)
    return out
def matvec(row, v, n):
    w = 0
    for a in range(n):
        if bin(row[a] & v).count('1') & 1: w |= 1 << a
    return w
res = {}
for n in range(3, 15):
    for trial in range(400):
        row = rand_rows(n)
        P = row[:]            # A^1
        vec = (1 << n) - 1    # A^0 1
        trsum = 0
        for k in range(1, 7):
            # P = A^(2k-1) -> A^(2k+1)
            P = mul(row, mul(row, P, n), n)
            tr = sum((P[a] >> a) & 1 for a in range(n)) & 1
            trsum ^= tr
            vec = matvec(row, vec, n)   # A^k 1
            w = bin(vec).count('1')
            ok = ((comb(w, 2) + n // 4 + trsum) % 2 == 0)
            res.setdefault(k, [0, 0])[0 if ok else 1] += 1
for k in sorted(res): print("k =", k, "ok / fail:", res[k])
