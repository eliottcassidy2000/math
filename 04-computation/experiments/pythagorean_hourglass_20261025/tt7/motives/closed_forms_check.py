# Check the explicit closed forms of the four triple products (sums over s != t, s, t not in {0, 1}; phi the quadratic
# character, L(l) = sum_w phi(w(w-1)(w-l))):
#   Y02 = - sum phi(s) phi(t) phi(t-1) L(s)^2 L(t/s)
#   Y05 = - sum phi(s-1) phi(t) phi(t-1) L(s)^2 L(t/s)
#   Y09 =   sum phi(s) phi(s-1) phi(t) L(s) L(t) L(t/s)
#   Y18 =   sum phi(s-1) phi(t-1) L(s) L(t) L(t/s)
import numpy as np
from sympy import primerange
def closed(q):
    chi = np.full(q, -1, dtype=np.int64); chi[0] = 0; chi[[(x * x) % q for x in range(1, q)]] = 1
    w = np.arange(q); base = chi[(w * (w - 1)) % q]
    L = np.array([int((base * chi[(w - l) % q]).sum()) for l in range(q)], dtype=np.int64)
    inv = np.zeros(q, dtype=np.int64); inv[1:] = [pow(int(x), q - 2, q) for x in range(1, q)]
    t = np.arange(2, q)
    Y = np.zeros(4, dtype=object)
    for s in range(2, q):
        Lr = L[(t * inv[s]) % q]; ok = (t != s)
        common = Lr * ok
        Y[0] += -int((chi[s] * chi[t] * chi[(t - 1) % q] * L[s] ** 2 * common).sum())
        Y[1] += -int((chi[(s - 1) % q] * chi[t] * chi[(t - 1) % q] * L[s] ** 2 * common).sum())
        Y[2] += int((chi[s] * chi[(s - 1) % q] * chi[t] * L[s] * L[t] * common).sum())
        Y[3] += int((chi[(s - 1) % q] * chi[(t - 1) % q] * L[s] * L[t] * common).sum())
    return list(Y)
rows = {int(l.split()[0]): list(map(int, l.split()[1:])) for l in open('cheapY_7_999.txt')}
ok = all(closed(q) == [rows[q][0], rows[q][2], rows[q][4], rows[q][5]] for q in primerange(7, 400) if q % 4 == 3)
print("closed forms reproduce the cone-reduced Y02, Y05, Y09, Y18 at every prime q = 3 mod 4 below 400:", ok)
