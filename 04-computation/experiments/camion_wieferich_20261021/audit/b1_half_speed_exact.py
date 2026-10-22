# Audit b1 (independent of the session's code): THM-4750 statement 1 (half-speed Mahler bound) and the
# higher Kummer congruence, with EXACT integers.
#
# s(n) is computed exactly (Python big integers) from the Hurwitz inverse W = T^(*-1) of T(n) = 2^C(n,2):
#     W(0) = 1,  W(n) = -sum_{k=1}^{n} C(n,k) T(k) W(n-k),   s(n) = -W(n)  (n >= 1),
# which is a different recurrence from Wright's (it never uses s inside the sum), and is cross-checked against
# brute force for n <= 5 and against Wright's recurrence for n <= 60.
# For each (p, k): eps = (2/p), D = s - c, A_m(l) = eps^(-l) D(m + l(p-1)), a_i(m) = Delta^i A_m(0) (exact integers).
# Checks, for every base m in 1..MMAX:
#   (i)  v_p(a_i(m)) >= min(k, ceil(i/2)) for i = 0..2k+1 whenever m >= 2k-1 (THM-4750 (1)), listing failures below;
#   (ii) the order-(2k-1) higher Kummer congruence  sum_j C(2k-1,j)(-eps)^(2k-1-j) D(m+j(p-1)) = 0 mod p^k, same for s;
#   (iii) the strongest exact form: v_p(a_i(m)) >= min(floor((m+1)/2), ceil(i/2)) (every k with m >= 2k-1 at once);
#   (iv) sharpness: min over the p-1 classes (m in [2k-1, 2k-1+p-2]) of min(v_p(a_i), k) versus min(k, ceil(i/2)).
# Usage: python3 b1_half_speed_exact.py
from math import comb
from itertools import combinations, product

NMAX = 240
T = [1 << (n * (n - 1) // 2) for n in range(NMAX + 1)]
W = [1] + [0] * NMAX
for n in range(1, NMAX + 1):
    W[n] = -sum(comb(n, k) * T[k] * W[n - k] for k in range(1, n + 1))
s = [0] + [-W[n] for n in range(1, NMAX + 1)]

# cross-check 1: Wright's recurrence for n <= 60
sw = [0] * 61
for n in range(1, 61):
    sw[n] = T[n] - sum(comb(n, k) * sw[k] * T[n - k] for k in range(1, n))
assert sw[1:] == s[1:61], "inverse recurrence disagrees with Wright"

# cross-check 2: brute force strong tournaments for n <= 5
def strong_count(n):
    pairs = list(combinations(range(n), 2)); cnt = 0
    for bits in product((0, 1), repeat=len(pairs)):
        out = [set() for _ in range(n)]; inn = [set() for _ in range(n)]
        for (i, j), b in zip(pairs, bits):
            if b: out[i].add(j); inn[j].add(i)
            else: out[j].add(i); inn[i].add(j)
        def reach(adj):
            seen = {0}; st = [0]
            while st:
                u = st.pop()
                for v in adj[u]:
                    if v not in seen: seen.add(v); st.append(v)
            return len(seen) == n
        if reach(out) and reach(inn): cnt += 1
    return cnt
assert [strong_count(n) for n in range(1, 6)] == s[1:6], "brute force disagrees"
print("s(1..8) =", s[1:9], " (inverse recurrence = Wright for n<=60, = brute force for n<=5)")

def c(n): return (1 << ((n - 1) // 2)) if n % 2 else 0
D = [0] + [s[n] - c(n) for n in range(1, NMAX + 1)]

def vp(x, p):
    if x == 0: return 10 ** 9
    v = 0
    while x % p == 0: x //= p; v += 1
    return v

def mahler(seq):
    out = []; t = list(seq)
    while t:
        out.append(t[0]); t = [t[j + 1] - t[j] for j in range(len(t) - 1)]
    return out

pairs = [(3, 2), (3, 3), (3, 4), (3, 6), (5, 2), (5, 3), (5, 4), (5, 5), (7, 2), (7, 3), (7, 4), (11, 2), (11, 3), (11, 4),
         (13, 2), (13, 3), (13, 4), (17, 2), (17, 3), (19, 2), (19, 3), (23, 2), (23, 3), (29, 2), (31, 2), (37, 2), (41, 2), (43, 2)]
total_checks = 0; total_viol = 0; exact_viol = 0
for p, k in pairs:
    eps = 1 if p % 8 in (1, 7) else -1
    I = 2 * k + 2                                     # coefficients i = 0..2k+1
    mmax = NMAX - (I - 1) * (p - 1)
    below_fail = []; above_fail = []; kum_fail_D = []; kum_fail_s = []
    minv = {}
    for m in range(1, mmax + 1):
        A = [eps ** l * D[m + l * (p - 1)] for l in range(I)]       # eps^(-l) = eps^l
        a = mahler(A)
        bad = [i for i in range(I) if vp(a[i], p) < min(k, (i + 1) // 2)]
        km = (m + 1) // 2
        if any(vp(a[i], p) < min(km, (i + 1) // 2) for i in range(I)): exact_viol += 1
        if m >= 2 * k - 1:
            total_checks += I
            if bad: above_fail.append((m, bad)); total_viol += 1
        elif bad: below_fail.append(m)
        r = 2 * k - 1
        kd = sum(comb(r, j) * (-eps) ** (r - j) * D[m + j * (p - 1)] for j in range(r + 1))
        ks = sum(comb(r, j) * (-eps) ** (r - j) * s[m + j * (p - 1)] for j in range(r + 1))
        if kd % p ** k: kum_fail_D.append(m)
        if ks % p ** k: kum_fail_s.append(m)
        if 2 * k - 1 <= m <= 2 * k - 1 + p - 2:
            for i in range(I):
                minv[i] = min(minv.get(i, 10 ** 9), min(vp(a[i], p), k))
    sharp = [minv[i] for i in range(I)]
    bound = [min(k, (i + 1) // 2) for i in range(I)]
    print(f"p={p:2d} k={k}: m=1..{mmax}; bound fails at m>=2k-1: {above_fail if above_fail else 'none'}; "
          f"fails below 2k-1 at m={below_fail}; Kummer(D) fails m={[x for x in kum_fail_D if x >= 2*k-1] or 'none'} "
          f"(below: {[x for x in kum_fail_D if x < 2*k-1]}); Kummer(s) fails m>=2k-1: {[x for x in kum_fail_s if x >= 2*k-1] or 'none'}")
    print(f"        min over classes of min(v,k): {sharp}   bound min(k,ceil(i/2)): {bound}   {'EQUAL' if sharp == bound else 'DIFFER'}")
print(f"TOTAL: {total_checks} coefficient checks at m >= 2k-1, violations {total_viol}; "
      f"exact form v(a_i(m)) >= min(floor((m+1)/2), ceil(i/2)) violated for {exact_viol} (p,k,m) triples (all m, all pairs)")
