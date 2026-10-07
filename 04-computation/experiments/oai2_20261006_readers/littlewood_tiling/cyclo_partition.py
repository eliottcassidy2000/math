"""Explicit replacement for the random partition in the stacking lemma of the 3D tiling preprint:
Z_q = (C_0 u {0}) u C_1 u ... u C_{s-1} (cyclotomic classes of order s, q = 1 mod s prime).
Check E - E = Z_q for every class.  s = 2: the Paley partition ({0} u QR, NQR)."""
import sympy, numpy as np
def classes(q, s):
    g = sympy.primitive_root(q)
    ind = {}
    x = 1
    for e in range(q - 1):
        ind[x] = e; x = x * g % q
    C = [[] for _ in range(s)]
    for x in range(1, q):
        C[ind[x] % s].append(x)
    C[0].append(0)
    return C
def full_diff(E, q):
    E = np.array(E)
    D = (E[:, None] - E[None, :]) % q
    return len(np.unique(D)) == q
for s in range(2, 9):
    fails = []; ok = 0
    for q in sympy.primerange(3, 6000):
        if (q - 1) % s: continue
        C = classes(q, s)
        if all(full_diff(E, q) for E in C): ok += 1
        else: fails.append(q)
    print(f"s={s}: primes q=1 mod s below 6000: {ok} OK; failures: {fails[-8:]} (largest {max(fails) if fails else None}); s^4 = {s**4}")
