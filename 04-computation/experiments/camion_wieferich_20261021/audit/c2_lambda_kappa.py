# c2_lambda_kappa.py -- second-pass audit of THM-4750 5'(b), 5'(c), 5(a), 5(e).
#
# Uses the DIRECT values eps^k s_a(p^k), eps^k D_a(p^k) mod p^(k+1) from c1_rigidity_direct.c (T' = S' T^2 recurrence) and:
#  (b) checks Lambda_p (direct, base 2) = sigma_p(omega(2)) - <2>^(-1/2) mod p^(k+1), with <2>^(-1/2) computed as
#      <2>^((p^(K-1)-1)/2)  (exponent = -1/2 mod p^(K-1); independent of the session's Newton square root);
#  (c) kappa(2) = omega(2) sigma_p^(-2): compares with the session's printed kappa (camion_lift_check.out), checks
#      v_p(Lambda_p) = v_p(2 - kappa), and DESIGNED HEIGHTS: bases a_h = kappa + 3 p^h must have h_a(p) = h exactly
#      (v(D_a(p^k)) >= k+1 for k < h, = h for k >= h), and the integer base a = kappa mod p^(kmax+1) passes every level;
#  5(e): my own implementation of the class-polynomial value at z = -p/(p-1): the full degree-(2K-2) form and the cheap
#      K = 3 form eps[D(p) + a_1 z + a_2 C(z,2)], against the direct Lambda_p mod p^3 / p^4; also at p below 2K-1.
# Usage: python3 c2_lambda_kappa.py /path/to/c1_binary
# Second-pass auditor, 2026-10-21.
import sys, subprocess, re
C1 = sys.argv[1]

def run_c1(p, kmax, *args):
    out = subprocess.run([C1, str(p), str(kmax)] + [str(a) for a in args], capture_output=True, text=True, check=True).stdout
    res = {}
    for line in out.splitlines():
        m = re.match(r'p=(\d+) eps=([+-]1) base a=(\d+)', line)
        if not m: continue
        a = int(m.group(3)); d = {}
        for k, sv, dv, v in re.findall(r'k=(\d+): eps\^k s\(p\^k\)=(\d+), eps\^k D\(p\^k\)=(\d+) \[v=(\d+)\]', line):
            d[int(k)] = (int(sv), int(dv), int(v))
        res[a] = d
    return res, out

def vp(x, p, cap):
    if x == 0: return cap
    v = 0
    while x % p == 0 and v < cap: x //= p; v += 1
    return v

def teich(a, p, K):
    q = p ** K; x = a % q
    for _ in range(K + 3): x = pow(x, p, q)
    return x

print("=== (b), (c): Lambda_p, sigma_p, kappa(2) from direct values ===")
session_kappa = {7: (44, 3), 11: (7449, 4), 13: (3213, 4)}          # camion_lift_check.out
session_lambda3 = {7: 210, 11: 1254, 13: 845, 17: 816, 19: 703, 23: 8119}   # tower_height_direct_checks.out / camion_lift.out
plan = [(3, 8), (5, 6), (7, 5), (11, 4), (13, 4), (17, 3), (19, 3), (23, 3), (29, 2), (31, 2), (37, 2), (41, 2), (43, 2)]
kappas = {}
for p, kmax in plan:
    K = kmax + 1; q = p ** K
    w = teich(2, p, K)
    res, _ = run_c1(p, kmax, "explicit", 2, w)
    eps = 1 if pow(2, (p - 1) // 2, p) == 1 else -1
    br = 2 * pow(w, -1, q) % q                                 # <2> = 2/omega
    isq = pow(br, (p ** (K - 1) - 1) // 2, q)                  # <2>^(-1/2) mod p^K
    assert (isq * isq * br) % q == 1 and isq % p == 1
    okb = True; oksig = True
    for k in range(1, kmax + 1):
        mod = p ** (k + 1)
        sig2, lam2, _ = res[2][k]; sigw, lamw, _ = res[w % q][k]
        if sig2 != sigw: oksig = False
        if (sigw - isq - lam2) % mod != 0: okb = False
    sig = res[w % q][kmax][0]                                   # sigma_p mod p^K
    lam = res[2][kmax][1]
    kap = w * pow(sig * sig, -1, q) % q; kappas[p] = (kap, K, kmax)
    line = (f"p={p:2d} eps={eps:+d} K={K}: sigma(2)=sigma(omega) at all levels: {oksig}; Lambda = sigma - <2>^(-1/2) at all levels: {okb}; "
            f"v(Lambda)={vp(lam, p, K)} v(2-kappa)={vp((2 - kap) % q, p, K)} v(2-omega)={vp((2 - w) % q, p, K)}; kappa(2) mod p^{K} = {kap}")
    if p in session_kappa:
        sk, sK = session_kappa[p]; line += f"; session kappa mod p^{sK} = {sk}: {'MATCH' if kap % p**sK == sk else 'MISMATCH'}"
    if p in session_lambda3 and kmax >= 2:
        line += f"; Lambda mod p^3 = {lam % p**3} (session {session_lambda3[p]}: {'MATCH' if lam % p**3 == session_lambda3[p] else 'MISMATCH'})"
    print(line)

print()
print("=== (c) designed heights: base a_h = kappa(2) + 3 p^h has h_a(p) = h; base kappa mod p^K passes every level ===")
for p in (5, 7, 11, 13):
    kap, K, kmax = kappas[p]; q = p ** K
    bases = [kap] + [(kap + 3 * p ** h) % q for h in range(1, kmax + 1)]
    res, _ = run_c1(p, kmax, "explicit", *bases)
    for idx, a in enumerate(bases):
        vals = [res[a][k][2] for k in range(1, kmax + 1)]       # v_p(D_a(p^k)), capped at k+1
        if idx == 0:
            pred = [k + 1 for k in range(1, kmax + 1)]; tag = "kappa itself (h >= K)"
        else:
            h = idx; pred = [(k + 1 if k < h else h) for k in range(1, kmax + 1)]; tag = f"designed h={h}"
        print(f"  p={p:2d} {tag:22s} base {a:>8d}: v_p(D_a(p^k)) k=1..{kmax} (cap k+1): {vals}  predicted {pred}  {'OK' if vals == pred else 'MISMATCH'}")

print()
print("=== 5(e): class-polynomial value at z vs direct Lambda (own implementation) ===")
def smod(a, N, q):                                              # plain Wright recurrence mod q (Python, small N)
    T = [pow(a, j * (j - 1) // 2, q) for j in range(N + 1)]; s = [0] * (N + 1); row = [1]
    for n in range(1, N + 1):
        row = [1] + [(row[i - 1] + row[i]) % q for i in range(1, n)] + [1]
        s[n] = (T[n] - sum(row[k] * s[k] * T[n - k] for k in range(1, n))) % q
    return s
def binom_z(z, i, q, p):                                        # C(z, i) mod q for i < p (z an integer representative)
    num = 1
    for t in range(i): num = num * (z - t) % q
    f = 1
    for t in range(1, i + 1): f *= t
    return num * pow(f, -1, q) % q
direct = {}
for p, kmax in plan:
    if kmax >= 2:
        res, _ = run_c1(p, kmax, "explicit", 2)
        direct[p] = {k: res[2][k][1] for k in range(1, kmax + 1)}   # eps^k D(p^k) mod p^(k+1)
for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
    eps = 1 if pow(2, (p - 1) // 2, p) == 1 else -1
    out = []
    for K in (3, 4):
        if K - 1 not in direct.get(p, {}): continue
        q = p ** K; r = 2 * K - 1; N = p + (r - 1) * (p - 1)
        s = smod(2, N, q)
        D = lambda n: (s[n] - (pow(2, (n - 1) // 2, q) if n % 2 else 0)) % q
        A = [D(p + j * (p - 1)) * (eps ** j) % q for j in range(r)]
        a = []; t = A[:]
        for i in range(r): a.append(t[0]); t = [(t[j + 1] - t[j]) % q for j in range(len(t) - 1)]
        z = (-p * pow(p - 1, -1, q)) % q
        if r - 1 < p:
            full = eps * sum(a[i] * binom_z(z, i, q, p) for i in range(r)) % q
        else:
            full = None
        lam_true = direct[p][K - 1] % q                         # eps^(K-1) D(p^(K-1)) = Lambda mod p^K
        ftxt = "n/a (needs i! invertible, i <= 2K-2 < p)" if full is None else f"{full} {'MATCH' if full == lam_true else 'MISMATCH'}"
        tag = f"K={K}: direct Lambda mod p^{K} = {lam_true}, full interpolation = {ftxt}"
        if K == 3:
            cheap = eps * (D(p) + a[1] * z + a[2] * binom_z(z, 2, q, p)) % q
            tag += f", cheap formula = {cheap} {'MATCH' if cheap == lam_true else 'MISMATCH'}"
        tag += "" if p >= 2 * K - 1 else f"  [p < 2K-1 = {2*K-1}: outside the stated range]"
        out.append(tag)
    print(f"p={p:2d}: " + " | ".join(out))
