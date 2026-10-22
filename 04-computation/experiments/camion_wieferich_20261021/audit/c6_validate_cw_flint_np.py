# c6_validate_cw_flint_np.py -- second-pass audit of the session's cw_flint_np.py (numpy-vectorised coefficients + FLINT).
#  1. Loads prefix_prod and vpow VERBATIM from ../cw_flint_np.py (function definitions only; the search loop is not run) and
#     replays the session's per-prime code path on every prime in [3, PMAX], comparing the Camion quotient q(p) (not just the
#     hit list) with q(p) = (D(p) mod p^2)/p from the direct T' = S' T^2 recurrence (c1_rigidity_direct.c).
#  2. At large primes near 1e6 compares the vectorised coefficient arrays (2^C(k,2) mod p, k! mod p, 1/k! mod p) with plain
#     Python loops for all k, and records the int64 headroom.
# Usage: python3 c6_validate_cw_flint_np.py /path/to/c1_binary PMAX
# Second-pass auditor, 2026-10-21.
import sys, os, re, subprocess, time, numpy as np, flint
C1, PMAX = sys.argv[1], int(sys.argv[2])
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'cw_flint_np.py')).read()
defs = src[src.index('def prefix_prod'):src.index('t0 = time.time()')]
ns = {'np': np}; exec(defs, ns); prefix_prod, vpow = ns['prefix_prod'], ns['vpow']
def session_q(p):                                            # verbatim replay of the loop body of cw_flint_np.py
    k = np.arange(p, dtype=np.int64)
    e = (k * (k - 1) // 2) % (p - 1)
    pw2 = vpow(2, e, p)
    ar = k.copy(); ar[0] = 1
    fact = prefix_prod(ar, p)
    invfact = vpow(fact, np.full(p, p - 2, dtype=np.int64), p)
    coeffs = (pw2 * invfact % p).tolist()
    U = flint.nmod_poly(coeffs, p).inverse_series_trunc(p + 1)
    up = int(U[p]); eps = 1 if p % 8 in (1, 7) else -1
    h = pow(2, (p - 1) // 2, p * p); e1 = ((h - eps) % (p * p)) // p % p
    return (up - e1) % p, (pw2, fact, invfact)
def primes(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(3, n + 1) if s[i]]
bad = 0; hits = []; n = 0
for p in primes(PMAX):
    qs, _ = session_q(p)
    o = subprocess.run([C1, str(p), "1", "explicit", "2"], capture_output=True, text=True, check=True).stdout
    m = re.search(r'eps=([+-]1) base a=2 .*k=1: eps\^k s\(p\^k\)=(\d+), eps\^k D\(p\^k\)=(\d+)', o)
    eps, eD = int(m.group(1)), int(m.group(3)); M = p * p
    D = eD if eps == 1 else (-eD) % M
    qd = D // p % p
    if qd != qs: bad += 1
    if qs == 0: hits.append(p)
    n += 1
print(f"[1] session code path vs direct recurrence, primes 3..{PMAX}: {n} primes, {bad} mismatches of q(p); hits {hits}")
for p in (100003, 499979, 999983):
    t0 = time.time(); qs, (pw2, fact, invfact) = session_q(p); dt = time.time() - t0
    ok = True; f = 1; e = 1
    for k in range(p):
        if k: f = f * k % p
        if int(fact[k]) != f or int(pw2[k]) != pow(2, k * (k - 1) // 2, p): ok = False; break
    if ok:
        for k in range(0, p, 997):
            if int(invfact[k]) * int(fact[k]) % p != 1: ok = False; break
    print(f"[2] p={p}: vectorised 2^C(k,2), k! (all k) and 1/k! (every 997th k) equal plain loops: {ok}; "
          f"q(p) = {qs}; max residue product {(p - 1) ** 2:.3e} < 2^63 = {2 ** 63:.3e}; time per prime {dt:.2f}s")
