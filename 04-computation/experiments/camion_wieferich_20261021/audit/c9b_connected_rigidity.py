# c9b_connected_rigidity.py -- does the connected-graph analogue have the Teichmueller rigidity that defines a "lift"?
# b_a(n): connected count built from T_a(n) = a^C(n,2) (b_a = T_a - sum C(n-1,k-1) b_a(k) T_a(n-k)); test whether
# b_a(p^k) mod p^(k+1) depends only on a mod p (k = 1, 2, 3), as s_a(p^k) does (c1_rigidity_direct.out).
# Level 1 is rigid by the mod-p formula b_a(p) = eps - 1 + p (w_p - V_a) (mod p^2), V_a = [x^p] sum_{m<p} (-1)^(m+1) y^m/m,
# y = T_(a,<p) - 1 (mod p).  Exact integers mod p^(k+1).  Second-pass auditor, 2026-10-21.
def conn(a, N, q):
    T = [pow(a, n * (n - 1) // 2, q) for n in range(N + 1)]; b = [0] * (N + 1); row = [1]   # row = C(n-1, .)
    for n in range(1, N + 1):
        if n >= 2: row = [1] + [(row[i - 1] + row[i]) % q for i in range(1, n - 1)] + [1]
        b[n] = (T[n] - sum(row[k - 1] * b[k] * T[n - k] for k in range(1, n))) % q
    return b
def strong(a, N, q):
    T = [pow(a, n * (n - 1) // 2, q) for n in range(N + 1)]; s = [0] * (N + 1); row = [1]
    for n in range(1, N + 1):
        row = [1] + [(row[i - 1] + row[i]) % q for i in range(1, n)] + [1]
        s[n] = (T[n] - sum(row[k] * s[k] * T[n - k] for k in range(1, n))) % q
    return s
for p, kmax in ((3, 4), (5, 3), (7, 2), (11, 2)):
    q = p ** (kmax + 1); N = p ** kmax
    lifts = [2, 2 + p, 2 + 2 * p, 2 + 5 * p]
    B = {a: conn(a, N, q) for a in lifts}; S = {a: strong(a, N, q) for a in lifts}
    rb = [all(B[a][p ** k] % p ** (k + 1) == B[2][p ** k] % p ** (k + 1) for a in lifts) for k in range(1, kmax + 1)]
    rs = [all(S[a][p ** k] % p ** (k + 1) == S[2][p ** k] % p ** (k + 1) for a in lifts) for k in range(1, kmax + 1)]
    vals = {a: [B[a][p ** k] % p ** (k + 1) for k in range(1, kmax + 1)] for a in lifts}
    print(f"p={p:2d}: lifts {lifts}: connected b_a(p^k) mod p^(k+1) lift-independent for k=1..{kmax}: {rb}; strong s_a (control): {rs}; values {vals}")
