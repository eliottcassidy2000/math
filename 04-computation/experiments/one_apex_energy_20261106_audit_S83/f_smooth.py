# Independent audit of THM-4887, part 6: how far do products of CERTIFIED primes (every odd prime <= 2^32 except 7 is
# a Hamiltonian-path count, THM-4875 (D)) beat the corollary A(X) >= X^(2/3)/(13 ln X + 121)?
# Rigorous lower bound: for each dyadic range R_j = [2^(j-1), 2^j) (j = 2..32) choose a multiset of k_j certified primes
# from R_j; the product is a Hamiltonian-path count (THM-1862), lies below 2^(sum j k_j), and distinct choices give
# distinct integers (unique factorisation).  So A(2^L) >= max over (k_j) with sum j k_j <= L of
# prod_j C(pi_j + k_j - 1, k_j).  Maximise by a Lagrange multiplier (each term is concave in k_j).
import math
import sympy
pis = {j: int(sympy.primepi(2 ** j)) for j in range(1, 33)}
cnt = {}
for j in range(2, 33):
    c = pis[j] - pis[j - 1]
    if j == 2: c -= 1          # range [2, 4): drop the even prime 2, keep 3
    if j == 3: c -= 1          # range [4, 8): drop 7 (not a count), keep 5
    cnt[j] = c
def lcomb_multiset(n, k):        # log C(n + k - 1, k)
    return math.lgamma(n + k) - math.lgamma(k + 1) - math.lgamma(n)
def best_k(n, j, lam):
    # maximise f(k) = log C(n+k-1, k) - lam * j * k over integers k >= 0 (concave): increase while gain > 0
    # f(k+1) - f(k) = log((n + k)/(k + 1)) - lam j  > 0  <=>  k < (n - e^(lam j)) / (e^(lam j) - 1)
    if lam * j > 600: return 0
    q = math.exp(lam * j)
    if q <= 1: return 10 ** 18
    k = max(0, math.floor((n - q) / (q - 1)) + 1)
    return k
def construct(L):
    lo, hi = 1e-9, 50.0
    for _ in range(200):
        lam = (lo + hi) / 2
        ks = {j: best_k(cnt[j], j, lam) for j in cnt}
        if sum(j * k for j, k in ks.items()) > L: lo = lam
        else: hi = lam
    ks = {j: best_k(cnt[j], j, hi) for j in cnt}
    assert sum(j * k for j, k in ks.items()) <= L
    return sum(lcomb_multiset(cnt[j], k) for j, k in ks.items() if k > 0) / math.log(10), ks
def log10_cor(L):
    lx = L * math.log(2)
    return (2 / 3) * L * math.log10(2) - math.log10(13 * lx + 121)
prev = None
for L in [1000, 2000, 3456, 5000, 10000, 15000, 20000, 25000, 30000, 40000, 60000]:
    lc, ks = construct(L)
    print(f"X = 2^{L}: smooth-product construction >= 10^{lc:.1f}, corollary = 10^{log10_cor(L):.1f}"
          f"  ({'construction ahead' if lc > log10_cor(L) else 'corollary ahead'})")
# bisection for the crossover (construct is monotone in L up to integer effects)
lo, hi = 3456, 200000
while hi - lo > 8:
    mid = (lo + hi) // 2
    if construct(mid)[0] > log10_cor(mid): lo = mid
    else: hi = mid
print(f"crossover near X = 2^{lo}..2^{hi}: below it, products of certified primes beat X^(2/3)/(13 ln X + 121)")

# Rigorous range: for X in [2^L, 2^(L+1)), A(X) >= A(2^L) >= construct(L), and the corollary is < cor(L+1).
# Check construct(L) >= cor(L+1) for every integer L from 32 to Lmax (below 2^32 the certified range alone gives
# A(X) >= (X+1)/2 - 2).
def construct_fast(L, iters=60):
    lo, hi = 1e-9, 50.0
    for _ in range(iters):
        lam = (lo + hi) / 2
        if sum(j * best_k(cnt[j], j, lam) for j in cnt) > L: lo = lam
        else: hi = lam
    ks = {j: best_k(cnt[j], j, hi) for j in cnt}
    # greedy fill of the unused bit budget (any k_j with sum j k_j <= L is a valid construction)
    rem = L - sum(j * k for j, k in ks.items())
    while True:
        best = None
        for j in cnt:
            if j <= rem and cnt[j] > 0:
                gain = math.log((cnt[j] + ks[j]) / (ks[j] + 1))   # log C(n+k, k+1) - log C(n+k-1, k)
                if best is None or gain > best[0]: best = (gain, j)
        if best is None: break
        ks[best[1]] += 1; rem -= best[1]
    assert sum(j * k for j, k in ks.items()) <= L
    return sum(lcomb_multiset(cnt[j], k) for j, k in ks.items() if k > 0) / math.log(10)
Lmax = None
for L in range(32, 16000):
    if construct_fast(L) < log10_cor(L + 1):
        Lmax = L; break
print(f"rigorous: products of certified primes beat X^(2/3)/(13 ln X + 121) for EVERY real X < 2^{Lmax} "
      f"(first integer L with construct(L) < cor(L+1) is {Lmax})")
