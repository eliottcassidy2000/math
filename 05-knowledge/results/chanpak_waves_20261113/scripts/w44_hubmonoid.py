"""Wave 44: the hub-seam monoid.  State of a tournament T with hub A:
   H = H(T), U = sum over nonempty paths Q of T-A ending at an in-neighbour of A of H(T-A-Q),
   Psi = H(T-A) + sum over those Q of (-1)^|Q| H(T-A-Q),  K = H(T-A).
Operations (all exact, proved in the README):
   S = triangle over T seamed at the hub (A beats x only):  H'=3H+3U-Psi+3K, U'=12U+6K, Psi'=K,  K'=3K
   P = plain triangle over T (triangle beats A):              H'=3H,          U'=12U+9K, Psi'=0,  K'=3K
   E = extend: new vertex b, A -> b, b beats T-A:              H'=H+U-Psi+K,   U'=2U,     Psi'=K,  K'=K
A reusable bounded digit = two words with identical matrices except the (H, K) entry: the difference is a constant
times K, and every later operation sees identical states.  Search all words up to a length."""
import itertools, sys
import numpy as np

S = np.array([[3, 3, -1, 3], [0, 12, 0, 6], [0, 0, 0, 1], [0, 0, 0, 3]], dtype=object)
P = np.array([[3, 0, 0, 0], [0, 12, 0, 9], [0, 0, 0, 0], [0, 0, 0, 3]], dtype=object)
E = np.array([[1, 1, -1, 1], [0, 2, 0, 0], [0, 0, 0, 1], [0, 0, 0, 1]], dtype=object)
ops = {'S': S, 'P': P, 'E': E}

def mat(word):          # word applied left to right
    M = np.identity(4, dtype=object)
    for ch in word:
        M = ops[ch].dot(M)
    return M

def key(M):             # everything except the (H, K) entry
    return tuple(M[i, j] for i in range(4) for j in range(4) if not (i == 0 and j == 3))

L = int(sys.argv[1]) if len(sys.argv) > 1 else 6
found = 0
for n in range(1, L + 1):
    groups = {}
    for w in itertools.product('SPE', repeat=n):
        w = ''.join(w)
        M = mat(w)
        groups.setdefault(key(M), []).append((M[0, 3], w))
    for k, lst in groups.items():
        vals = sorted(set(v for v, _ in lst))
        if len(vals) > 1:
            found += 1
            if found <= 12:
                print(f"length {n}: H-scale {k[0]}, deposits (times K): {vals}  words: {[w for _, w in sorted(lst)][:6]}")
print(f"words up to length {L}: {found} classes with two or more distinct bounded deposits")
