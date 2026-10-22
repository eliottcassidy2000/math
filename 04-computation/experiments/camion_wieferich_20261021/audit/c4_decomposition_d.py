# c4_decomposition_d.py -- second-pass audit of THM-4750 5'(d):  eps D(p)/p = rho_C(2) - q_2(p)/2 (mod p),
# rho_C(a) = (sigma_p(a) - 1)/p mod p, lift-independent.
#
# Three independent computations per prime 5 <= p <= PMAX:
#  (A) DIRECT: eps*s_2(p) and eps*s_omega(p) mod p^2 (omega = Teichmueller lift of 2 mod p^2), by the T' = S' T^2 recurrence
#      of c1_rigidity_direct.c (sigma_p = eps s_a(p) mod p^2 for EVERY lift a, by the tower relation);
#  (B) FLINT: rho_C(a) = eps * u_p(a) (mod p) with u_p(a) = [x^p] 1/T_(a,<p)(x) in F_p[[x]], T_(a,<p) = sum_{k<p} a^C(k,2) x^k/k!.
#      Derivation: s_a(p) = -p! [x^p] 1/T = a^C(p,2) + p u_p (mod p^2) (Wilson), and a^C(p,2) = (eps + p e')^p = eps (mod p^2),
#      so eps s_a(p) = 1 + p eps u_p (mod p^2).  Since T_(a,<p) mod p depends only on a mod p, this is an elementary proof
#      of the rigidity at level 1 (mod p^2).
#  (C) the session's camion_residue_2000.out (rho_C, q_2, q_C), compared line by line.
# Usage: python3 c4_decomposition_d.py /path/to/c1_binary PMAX
# Second-pass auditor, 2026-10-21.
import sys, subprocess, re, flint, os
C1, PMAX = sys.argv[1], int(sys.argv[2])
HERE = os.path.dirname(os.path.abspath(__file__))
sess = {}
for line in open(os.path.join(HERE, '..', 'camion_residue_2000.out')):
    m = re.match(r'p=(\d+) eps=([+-]1) rho_C=(\d+) q_2=(\d+) q_C=(\d+)', line)
    if m: sess[int(m.group(1))] = tuple(int(x) for x in m.group(3, 4, 5))
def primes(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(5, n + 1) if s[i]]
def u_p(a, p):
    c = []; t = 1
    for k in range(p):                                       # a^C(k,2)/k! mod p, by plain loops (no vectorisation)
        c.append(pow(a, k * (k - 1) // 2, p) * pow(t, p - 2, p) % p if k else 1)
        t = t * (k + 1) % p
    return int(flint.nmod_poly(c, p).inverse_series_trunc(p + 1)[p])
bad = {'A_lift': 0, 'decomp': 0, 'B_vs_A': 0, 'sess': 0}; n = 0; zeros = []; missing = []
for p in primes(PMAX):
    q = p * p; w = pow(2, p, q)                              # Teichmueller lift of 2 mod p^2
    out = subprocess.run([C1, str(p), "1", "explicit", "2", str(w)], capture_output=True, text=True, check=True).stdout
    vals = {}
    for line in out.splitlines():
        m = re.match(r'p=\d+ eps=([+-]1) base a=(\d+).*k=1: eps\^k s\(p\^k\)=(\d+), eps\^k D\(p\^k\)=(\d+)', line)
        if m: vals[int(m.group(2))] = (int(m.group(1)), int(m.group(3)), int(m.group(4)))
    eps, sig2, epsD2 = vals[2]; _, sigw, _ = vals[w]
    if sig2 != sigw: bad['A_lift'] += 1
    rho = (sigw - 1) // p % p
    q2 = (pow(2, p - 1, q) - 1) // p % p
    qC = epsD2 // p % p
    if (rho - q2 * pow(2, -1, p) - qC) % p: bad['decomp'] += 1
    rhoB = eps * u_p(2, p) % p
    if rhoB != rho: bad['B_vs_A'] += 1
    if p in sess:
        if sess[p] != (rho, q2, qC): bad['sess'] += 1
    else: missing.append(p)
    if rho == 0: zeros.append(p)
    n += 1
print(f"primes 5..{PMAX}: {n}")
print(f"(A) sigma from base 2 = sigma from omega(2), mod p^2: failures {bad['A_lift']}")
print(f"(d) eps D(p)/p = rho_C - q_2/2 (mod p): failures {bad['decomp']}")
print(f"(B) FLINT rho_C = eps [x^p] 1/T_<p  vs (A): failures {bad['B_vs_A']}")
print(f"(C) session camion_residue_2000.out (rho_C, q_2, q_C) vs (A): failures {bad['sess']}; primes not in the session file: {missing[:5]}{'...' if len(missing) > 5 else ''} ({len(missing)})")
print(f"rho_C(2) = 0 at p = {zeros}")
