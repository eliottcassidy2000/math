# Teichmueller rigidity of the Camion limit and the "Camion lift" (test of a derivation, THM-4750 continuation).
# Claim: sigma_p(a) := lim_k eps^k s_a(p^k) (s_a = e - T_a^(-1), T_a(n) = a^C(n,2), eps = (a/p)) depends only on a mod p,
# because every grade >= 1 term of the Session25 expansion carries C(n,d) with d >= 1, which vanishes at the p-adic point n = 0
# of the branch n = 1 mod (p-1).  The circulant part contributes lim eps^k a^((p^k-1)/2) = <a>^(-1/2).  Hence
#   Lambda_p(a) = sigma_p - <a>^(-1/2),   h_a(p) = v_p(a - kappa(a)),   kappa(a) := omega(a) * sigma_p^(-2)  ("Camion lift"),
# parallel to the Wieferich order v_p(a - omega(a)).  (Consequence, not checked here: exactly one lift of each residue mod p^2
# is Camion-Wieferich.)  Interpolation needs p >= 2K - 1; run for 7 <= p <= 43.
# Computes sigma_p(a) mod p^K by the degree-(2K-2) interpolation of l -> eps^-l s_a(p + l(p-1)) at l = -p/(p-1).
# mac-mini-2026-10-21-fanopsp.
import sys
def smod(a, M, q):
    T = [pow(a, j * (j - 1) // 2, q) for j in range(M + 1)]; s = [0] * (M + 1); row = [1]
    for n in range(1, M + 1):
        row = [1] + [(row[i - 1] + row[i]) % q for i in range(1, n)] + [1]
        s[n] = (T[n] - sum(row[k] * s[k] * T[n - k] for k in range(1, n))) % q
    return s
def vp(x, p, K):
    if x % p**K == 0: return K
    c = 0
    while x % p == 0: x //= p; c += 1
    return c
def sigma_lim(a, p, K):
    q = p**K; eps = 1 if pow(a, (p - 1) // 2, p) == 1 else -1
    r = 2 * K - 1; M = p + (r - 1) * (p - 1)
    s = smod(a, M, q)
    A = [(s[p + j * (p - 1)] * (eps ** j)) % q for j in range(r)]          # eps^-j = eps^j
    dif = []; t = A[:]
    for i in range(r): dif.append(t[0]); t = [(t[j + 1] - t[j]) % q for j in range(len(t) - 1)]
    z = (-p * pow(p - 1, -1, q)) % q                                      # l_inf = -p/(p-1)
    val = 0; b = 1
    for i in range(r):
        val = (val + b * dif[i]) % q
        b = b * ((z - i) % q) % q * pow(i + 1, -1, q) % q
    return val * eps % q, eps                                             # interpolated value is lim eps^(a-1) s(p^a); times eps gives sigma_p
def teich(a, p, K):
    q = p**K; x = a % q
    for _ in range(K + 2): x = pow(x, p, q)
    return x
def sqrt_principal(u, p, K):          # square root of a principal unit u = 1 mod p, the one = 1 mod p
    q = p**K; x = 1
    for _ in range(2 * K + 4): x = (x + u * pow(x, -1, q)) * pow(2, -1, q) % q
    return x
K = 3
for p in (7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
    q = p**K
    out = []
    sig2, eps = sigma_lim(2, p, K)
    for a in (2, 2 + p, 2 + 2 * p, 2 + 5 * p, 2 - p + p * p):
        sg, e2 = sigma_lim(a, p, K)
        out.append((a, sg == sig2))
    w = teich(2, p, K); br = 2 * pow(w, -1, q) % q                         # <2> = 2/omega(2)
    lam = (sig2 - pow(sqrt_principal(br, p, K), -1, q)) % q
    kappa = w * pow(sig2 * sig2 % q, -1, q) % q
    print(f"p={p:3d} eps={eps:+d}  sigma_p(a) equal for a = 2 mod p: {all(ok for _, ok in out)} {[a for a, ok in out]}"
          f"   v_p(Lambda) = {vp(lam, p, K)}  v_p(2 - kappa) = {vp((2 - kappa) % q, p, K)}  v_p(2 - omega) = {vp((2 - w) % q, p, K)}")
# direct cross-check: Lambda_p from the interpolation of D = s - c (as in tower_height.c) equals sigma_p - <2>^(-1/2)
print("cross-check Lambda_p (D-interpolation) = sigma_p - <2>^(-1/2) mod p^K:")
for p in (7, 11, 13, 17, 19, 23):
    q = p**K; eps = 1 if pow(2, (p - 1) // 2, p) == 1 else -1; r = 2 * K - 1; M = p + (r - 1) * (p - 1)
    s = smod(2, M, q)
    D = lambda n: (s[n] - (pow(2, (n - 1) // 2, q) if n % 2 else 0)) % q
    A = [(D(p + j * (p - 1)) * (eps ** j)) % q for j in range(r)]
    dif = []; t = A[:]
    for i in range(r): dif.append(t[0]); t = [(t[j + 1] - t[j]) % q for j in range(len(t) - 1)]
    z = (-p * pow(p - 1, -1, q)) % q; val = 0; b = 1
    for i in range(r):
        val = (val + b * dif[i]) % q; b = b * ((z - i) % q) % q * pow(i + 1, -1, q) % q
    lam_direct = val * eps % q
    sig, _ = sigma_lim(2, p, K); w = teich(2, p, K); br = 2 * pow(w, -1, q) % q
    lam = (sig - pow(sqrt_principal(br, p, K), -1, q)) % q
    print(f"  p={p}: direct {lam_direct}, formula {lam}, {'MATCH' if lam == lam_direct else 'MISMATCH'}")
