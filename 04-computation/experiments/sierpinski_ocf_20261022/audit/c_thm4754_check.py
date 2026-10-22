# Audit of THM-4754 statements 1-6 (auditor's own brute force).
import random
from math import comb
rng = random.Random(4754)

def sign3(z):
    """3-adic sign (unit part mod 3) of a nonzero integer or residue representative."""
    assert z != 0
    while z % 3 == 0: z //= 3
    return z % 3

def sgn_mod(z, m):
    z %= m
    assert z != 0
    while z % 3 == 0: z //= 3
    return z % 3

def S_adj(k):
    m = 3 ** k
    return [[1 if (i != j and sgn_mod(j - i, m) == 1) else 0 for j in range(m)] for i in range(m)]

def lex_power(k):
    """C3^[k] on F_3^k encoded as integers: vertex x <-> digits; top level = lowest digit.
       x -> y iff at the lowest digit position where they differ, digit_y = digit_x + 1 mod 3."""
    m = 3 ** k
    A = [[0] * m for _ in range(m)]
    for x in range(m):
        for y in range(m):
            if x == y: continue
            a, b = x, y
            while a % 3 == b % 3: a //= 3; b //= 3
            A[x][y] = 1 if (b % 3 - a % 3) % 3 == 1 else 0
    return A

# (1) tournament + equality with the lexicographic power under the base-3 digit map (identity on integers)
ok1 = True
for k in range(1, 6):
    A = S_adj(k); m = 3 ** k
    if any(A[i][j] + A[j][i] != 1 for i in range(m) for j in range(m) if i != j): ok1 = False
    if A != lex_power(k): ok1 = False
print("(1) S_k = Cay(Z/3^k, sign 1) is a tournament and equals C3^[k] (digit map), k<=5:", ok1)
# lexicographic power is C3[S_{k-1}] with top level = residue mod 3: covered by the digit rule.
# (2) self-similarity
ok2 = True
for k in range(2, 7):
    m = 3 ** k; mp_ = 3 ** (k - 1)
    for _ in range(3000):
        x, y = rng.randrange(mp_), rng.randrange(mp_)
        b = rng.randrange(3)
        if x != y:
            if (sgn_mod(y - x, mp_) == 1) != (sgn_mod((3 * y + b) - (3 * x + b), m) == 1): ok2 = False
        u, v = rng.randrange(m), rng.randrange(m)
        if u % 3 != v % 3 and (v - u) % 3 == 1 and sgn_mod(v - u, m) != 1: ok2 = False
print("(2) x->3x+b embeds S_(k-1) on class b; class b -> class b+1 (k<=6, random):", ok2)
ok2b = all(((-pow(2, -1, 3 ** k)) % 3 ** k) == (3 ** k - 1) // 2 for k in range(1, 30)) and all(int('1' * k, 3) == (3 ** k - 1) // 2 for k in range(1, 30))
it = all((lambda x, k: (lambda f: f(f, x, k))(lambda f, x, k: x if k == 0 else f(f, 3 * x + 1, k - 1)))(x, k) == 3 ** k * x + (3 ** k - 1) // 2 for x in range(-5, 6) for k in range(0, 8))
print("    -1/2 = ...111 (base 3) mod 3^k and (3x+1)^k(x) = 3^k x + (3^k-1)/2:", ok2b and it)
# (3) multipliers
ok3 = True
for k in range(1, 6):
    m = 3 ** k; A = S_adj(k)
    for lam in range(1, m):
        if lam % 3 == 0: continue
        for i in range(m):
            for j in range(m):
                if i == j: continue
                a = A[i][j]; b = A[(lam * i) % m][(lam * j) % m]
                if lam % 3 == 1 and a != b: ok3 = False
                if lam % 3 == 2 and a == b: ok3 = False
print("(3) unit lambda = 1 mod 3 automorphism, = 2 mod 3 reverses every arc (all units, k<=5):", ok3, "; 1/2 mod 3 =", pow(2, -1, 3))
# (4) Collatz arc-flip law
def C(x): return 3 * x + 1 if x % 2 else x // 2
def T(x): return (3 * x + 1) // 2 if x % 2 else x // 2
bad4 = badT = 0; tested = 0; fail_wo = 0
for _ in range(20000):
    t = rng.randint(1, 40)
    x = rng.randint(1, 10 ** 12)
    y = x + rng.randint(1, 10 ** 6) * 2 ** t
    if rng.random() < 0.5: y = x - rng.randint(1, 10 ** 6) * 2 ** t
    if y == x: continue
    a, b, e = x, y, 0
    for _s in range(t):
        if a % 2 != b % 2: bad4 += 1; break
        if a % 2 == 0: e += 1
        a, b = C(a), C(b)
    tested += 1
    if sign3(b - a) != (sign3(y - x) if e % 2 == 0 else 3 - sign3(y - x)): bad4 += 1
    a, b = x, y
    for _s in range(t): a, b = T(a), T(b)
    if sign3(b - a) != (sign3(y - x) if t % 2 == 0 else 3 - sign3(y - x)): badT += 1
# sanity: without the 2-adic condition the law can fail
for _ in range(2000):
    t = 6; x = rng.randint(1, 10 ** 6); y = x + rng.randint(1, 63)
    a, b, ea, eb = x, y, 0, 0
    for _s in range(t):
        ea += a % 2 == 0; eb += b % 2 == 0; a, b = C(a), C(b)
    if a != b and sign3(b - a) != (sign3(y - x) if ea % 2 == 0 else 3 - sign3(y - x)): fail_wo += 1
print(f"(4) arc-flip law, {tested} random pairs x = y mod 2^t (t<=40): violations C-map {bad4}, shortcut-map {badT}; (control without congruence: {fail_wo}/2000 violate)")
# (5) Steiner midpoint / binomial average / endpoint collapse
def star(a, b): return (-a - b) % 3
ok5 = True; collapse_n = []
for n in range(1, 60):
    for _ in range(20):
        s = [rng.randrange(3) for _ in range(n)]
        row = s[:]
        while len(row) > 1: row = [star(row[i], row[i + 1]) for i in range(len(row) - 1)]
        want = ((-1) ** (n - 1) * sum(comb(n - 1, i) * s[i] for i in range(n))) % 3
        if row[0] != want: ok5 = False
    if all(comb(n - 1, i) % 3 == 0 for i in range(1, n - 1)): collapse_n.append(n)
print("(5) star = (x+y)/2 in F_3:", all(star(a, b) == ((a + b) * 2) % 3 for a in range(3) for b in range(3)),
      "; reduction = (-1)^(n-1) sum C(n-1,i) x_i (n<60):", ok5, "; endpoint-only lengths n<60:", collapse_n)
# (6) Jacobsthal split
J = lambda k: (2 ** k - (-1) ** k) // 3
ok6 = all(2 ** k == 3 * J(k) + (-1) ** k for k in range(0, 80))
feeders = [n for n in range(1, 10 ** 6, 2) if (3 * n + 1) & (3 * n) == 0]
print("(6) 2^k = 3J_k + (-1)^k (k<80):", ok6, "; J_0..J_8 =", [J(k) for k in range(9)],
      "; odd n < 1e6 with 3n+1 a power of 2:", feeders, "== J_(2m):", feeders == [J(2 * m) for m in range(1, 11) if J(2 * m) < 10 ** 6])
