# Audit S83, step 4: the two counting lemmas, tested adversarially in EXACT integer arithmetic.
# (L1) session's lemma: if |beta - r/q| <= 1/q^2, gcd(r,q) = 1, then in any block of q consecutive d,
#      #{d : ||beta d|| <= t} <= 2qt + 3.  Equivalent: with the block's sorted distances t_1 <= ... <= t_q,
#      j - 3 <= 2 q t_j for every j.  We record max_j (j - 2 q t_j) (the lemma says <= 3).
# (L2) audit's first-power version used by the two-block bound:
#      sum_block min(K, 1/(2||beta d||)) <= 3K + q (1 + ln q).
import math, random
import numpy as np

def dist_num(N, D, ds):
    x = (N * ds) % D
    return np.minimum(x, D - x)          # ||beta d|| * D, beta = N/D

worst = -10.0; worst_case = None; cases = 0
for q in range(1, 121):
    rs = [r for r in range(0, q) if math.gcd(r, q) == 1] if q > 1 else [0, 1]
    for r in rs:
        # beta = r/q + k/(2 q^2), k in -2..2  (extremes |delta| = 1/q^2 included); D = 2 q^2
        for k in (-2, -1, 0, 1, 2):
            D = 2 * q * q; N = 2 * q * r + k
            for d0 in range(-2 * q, 2 * q + 1, max(1, q // 7)):
                ds = np.arange(d0, d0 + q, dtype=np.int64)
                s = np.sort(dist_num(N, D, ds))                  # numerators of t_j over D
                j = np.arange(1, q + 1)
                val = (j - 2 * q * s / D).max()                  # max_j (j - 2 q t_j)
                cases += 1
                if val > worst:
                    worst, worst_case = float(val), (q, r, k, d0)
print(f'(L1) {cases} exact blocks (q <= 120, all r, delta in {{0, +-1/(2q^2), +-1/q^2}}): max_j (j - 2 q t_j) = {worst:.4f} at (q, r, k, d0) = {worst_case}  [lemma: <= 3]')

# random real beta with |beta - r/q| <= 1/q^2, larger q
random.seed(5); worst2 = -10.0
for _ in range(4000):
    q = random.randint(2, 3000); r = random.randrange(1, q)
    while math.gcd(r, q) != 1:
        r = random.randrange(1, q)
    D = 1 << 60; N = (r * D) // q + random.randint(-D // (q * q), D // (q * q))
    d0 = random.randint(-10 ** 6, 10 ** 6)
    ds = np.arange(d0, d0 + q, dtype=object)
    x = np.array([(N * int(d)) % D for d in ds], dtype=object)
    s = np.sort(np.array([min(int(v), D - int(v)) for v in x], dtype=np.float64)) / D
    j = np.arange(1, q + 1)
    worst2 = max(worst2, float((j - 2 * q * s).max()))
print(f'(L1) 4000 random blocks (q <= 3000, |delta| <= 1/q^2): max_j (j - 2 q t_j) = {worst2:.4f}  [lemma: <= 3]')

# (L2)
random.seed(9); worst3 = 0.0
for _ in range(3000):
    q = random.randint(1, 2000); r = random.randrange(0, q) if q > 1 else 1
    while math.gcd(r, q) != 1:
        r = random.randrange(0, q)
    beta = r / q + random.uniform(-1, 1) / q ** 2
    K = random.choice([1, 2, 5, 17, 100, 1000, 10 ** 5])
    d0 = random.randint(-10 ** 5, 10 ** 5)
    ds = np.arange(d0, d0 + q, dtype=np.float64)
    x = beta * ds; dist = np.abs(x - np.round(x))
    with np.errstate(divide='ignore'):
        terms = np.minimum(K, 1 / (2 * dist))
    lhs = terms.sum(); rhs = 3 * K + q * (1 + math.log(q))
    worst3 = max(worst3, lhs / rhs)
print(f'(L2) 3000 random blocks: max of sum_block min(K, 1/(2||beta d||)) / (3K + q(1 + ln q)) = {worst3:.4f}  [<= 1 required]')
