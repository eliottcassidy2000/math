# c5_cw_lists_flint.py -- second-pass audit: Camion-Wieferich lists for base 3 (strong semicomplete digraphs) and for
# connected labelled graphs, by FLINT power series MOD p^2 (an algorithm independent of the session's O(p^2) recurrences
# cw_base.c / cw_connected.c), plus the mod-p identity for base a, plus direct C cross-checks for small p.
#
#  base a:   s_a(p) = (p-1)! [x^(p-1)] T'(x)/T(x)^2   (S = 1 - 1/T, S' = T'/T^2; only coefficients x^k, k <= p-1, are needed,
#            so all denominators k! are units mod p^2);  D_a(p) = s_a(p) - a^((p-1)/2);  CW iff p^2 | D_a(p).
#            Rotation congruence (derived in the audit report): s_a(p) = a^((p-1)/2) (mod p) for a = 2 (circulant tournaments)
#            and a = 3 (circulant semicomplete digraphs: every one is strong since each difference generates Z_p).
#            Mod-p identity (cw_fast.c generalised): D_a(p)/p = u_p - e' (mod p), u_p = [x^p] 1/T_<p, e' = (a^((p-1)/2) - eps)/p.
#  connected graphs:  cg(p) = (p-1)! [x^(p-1)] G'(x)/G(x), G = sum 2^C(n,2) x^n/n!;  rotation congruence
#            cg(p) = #(nonempty connection sets) = 2^((p-1)/2) - 1 (mod p);  "connected-CW" iff p^2 | cg(p) - 2^((p-1)/2) + 1.
# Usage: python3 c5_cw_lists_flint.py PMAX_BASE PMAX_CONN
# Second-pass auditor, 2026-10-21.
import sys, time, flint
PB, PC = int(sys.argv[1]), int(sys.argv[2])
def primes(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(3, n + 1) if s[i]]
def series(a, p, M):
    """T_k = a^C(k,2)/k! (k < p) and T'_k = a^C(k+1,2)/k! (k < p), mod M = p^2; also (p-1)! mod M."""
    fact = [1] * p
    for k in range(1, p): fact[k] = fact[k - 1] * k % M
    invf = [0] * p; invf[p - 1] = pow(fact[p - 1], -1, M)
    for k in range(p - 1, 0, -1): invf[k - 1] = invf[k] * k % M
    T = [0] * p; Tp = [0] * p; e = 1; ak = 1               # e = a^C(k,2), ak = a^k
    for k in range(p):
        T[k] = e * invf[k] % M
        e2 = e * ak % M                                     # a^C(k+1,2)
        Tp[k] = e2 * invf[k] % M
        e = e2; ak = ak * a % M
    return T, Tp, fact[p - 1]
def s_mod_p2(a, p):
    M = p * p; T, Tp, f = series(a, p, M)
    W = flint.nmod_poly(T, M).inverse_series_trunc(p)
    W2 = (W * W).truncate(p)
    prod = (flint.nmod_poly(Tp, M) * W2)
    return f * int(prod[p - 1]) % M
def cg_mod_p2(p):
    M = p * p; G, Gp, f = series(2, p, M)
    W = flint.nmod_poly(G, M).inverse_series_trunc(p)
    prod = flint.nmod_poly(Gp, M) * W
    return f * int(prod[p - 1]) % M
def q_mod_p(a, p):                                          # D_a(p)/p mod p by the mod-p identity
    c = [1]; t = 1
    for k in range(1, p):
        t = t * k % p; c.append(pow(a, k * (k - 1) // 2, p) * pow(t, p - 2, p) % p)
    up = int(flint.nmod_poly(c, p).inverse_series_trunc(p + 1)[p])
    eps = 1 if pow(a, (p - 1) // 2, p) == 1 else -1
    e1 = ((pow(a, (p - 1) // 2, p * p) - eps) % (p * p)) // p
    return (up - e1) % p
# --- self-test against exact small values
def exact_strong(a, N):
    T = [a ** (n * (n - 1) // 2) for n in range(N + 1)]; s = [0] * (N + 1)
    from math import comb
    for n in range(1, N + 1): s[n] = T[n] - sum(comb(n, k) * s[k] * T[n - k] for k in range(1, n))
    return s
def exact_conn(N):
    from math import comb
    g = [2 ** (n * (n - 1) // 2) for n in range(N + 1)]; b = [0] * (N + 1)
    for n in range(1, N + 1): b[n] = g[n] - sum(comb(n - 1, k - 1) * b[k] * g[n - k] for k in range(1, n))
    return b
es2, es3, ec = exact_strong(2, 31), exact_strong(3, 31), exact_conn(31)
assert ec[1:8] == [1, 1, 4, 38, 728, 26704, 1866256]       # A001187
st = []
for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
    st.append(s_mod_p2(2, p) == es2[p] % (p * p) and s_mod_p2(3, p) == es3[p] % (p * p) and cg_mod_p2(p) == ec[p] % (p * p))
print("self-test (FLINT mod p^2 vs exact integers, p = 5..31, bases 2, 3 and connected):", "PASS" if all(st) else "FAIL")
# --- searches
t0 = time.time()
for a in (2, 3):
    hits = []; eul = 0; idbad = 0; wief = []; n = 0
    for p in primes(PB):
        if a % p == 0: continue
        M = p * p; s = s_mod_p2(a, p); D = (s - pow(a, (p - 1) // 2, M)) % M
        if D % p: eul += 1
        if D // p % p != q_mod_p(a, p): idbad += 1
        if D == 0: hits.append(p)
        if pow(a, p - 1, M) == 1: wief.append(p)
        n += 1
    print(f"base {a}: primes 3..{PB} tested {n} (p != {a}); Euler/rotation failures {eul}; mod-p identity mismatches {idbad}; "
          f"Camion-Wieferich {hits}; base-{a} Wieferich {wief}  ({time.time() - t0:.1f}s)")
hits = []; rot = 0; n = 0
for p in primes(PC):
    M = p * p; c = cg_mod_p2(p); tgt = (pow(2, (p - 1) // 2, M) - 1) % M
    if (c - tgt) % p: rot += 1
    if c == tgt: hits.append(p)
    n += 1
print(f"connected graphs: primes 3..{PC} tested {n}; rotation-congruence failures {rot}; p^2 | cg(p) - 2^((p-1)/2) + 1 at {hits}  ({time.time() - t0:.1f}s)")
