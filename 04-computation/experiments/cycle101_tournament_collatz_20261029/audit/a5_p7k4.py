#!/usr/bin/env python3
"""Audit of the post-checkpoint P7, k = 4 claims (THM-4848 data row, HYP-9446 (3) refutation):
beta = det B_1(M_4(P7)) in Z[eta], its content, the split-prime valuations, and the local Smith data of the
840 x 840 observer at 71, 79, 151, 2843, 7529. Uses only the audit's own a3_runobs.py routines."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a3_runobs import (circulant, observer, blocks_data, bareiss_det, detB1_coeffs, sigma, eta_coords, quad_norm,
                       split_valuations, kron, local_smith_np, local_smith_py, det_mod_np, norm_Qzeta)
from sympy import factorint
from math import gcd

def say(*a):
    print(" ".join(str(x) for x in a), flush=True)

t0 = time.time()
QR7 = {1, 2, 4}
adj = circulant(QR7, 7)
M, tup, idx = observer(7, adj, 4)
N = len(tup)
say(f"M_4(P7): {N}x{N}, total mass {sum(map(sum, M))}, built in {time.time() - t0:.1f}s")
R, Ms = blocks_data(M, tup, idx, 7)
B0 = [[sum(x) for x in row] for row in Ms]
t1 = time.time()
dB0 = bareiss_det(B0)
say(f"det B_0 (120x120) = {dB0}; factor = {factorint(dB0)} ({time.time() - t1:.1f}s)")
t1 = time.time()
coeffs, used = detB1_coeffs(Ms, 7)
say(f"det B_1 coefficients (CRT, {used} primes): {coeffs} ({time.time() - t1:.1f}s)")
say(f"invariant under sigma_2, sigma_4: {sigma(coeffs, 2, 7) == coeffs and sigma(coeffs, 4, 7) == coeffs}")
a, b = eta_coords(coeffs, 7, QR7)
say(f"beta = a + b eta: a = {a}, b = {b}")
claim_a = 4586613209427402393995389904877518848000
claim_b = -3427348182276905332040646159100805120000
say(f"agrees with the session's chirality_more.out: {(a, b) == (claim_a, claim_b)}")
g = gcd(a, b)
say(f"content gcd(a,b) = {factorint(g)}")
nq = quad_norm(a, b, 7)
F = factorint(nq)
say(f"quadratic norm = {F}")
for p, e in sorted(F.items()):
    kr = kron(-7, p)
    if kr == 1:
        vv = split_valuations(a, b, p, 7, e)
        say(f"   p={p}: split, exponent {e}, (v_P, v_Pbar) = {vv} -> {'one-sided' if min(vv) == 0 else 'TWO-SIDED'}, v_p(content) = {min(vv)}")
    else:
        say(f"   p={p}: {'inert' if kr == -1 else 'ramified'}, exponent {e}")
NB = nq ** 3
# det M_4 = det B_0 * N(beta) modulo three 31-bit primes
ok = all(det_mod_np(M, l) == dB0 * NB % l for l in (2147483629, 2147483587, 2147483579))
say(f"det M_4 = det B_0 * N_K(beta)^3 modulo three 31-bit primes: {ok}")
# local Smith data
for p, E in ((71, 3), (2843, 2), (79, 3), (151, 3), (7529, 2)):
    t1 = time.time()
    res, big, nunit = local_smith_np(M, p, E)
    resB, bigB, _ = local_smith_np(B0, p, E)
    say(f"   p={p} (p mod 7 = {p % 7}): coker M_4 local Smith exponents {res} (+{big} with exponent >= {E}); "
        f"coker B_0: {resB} (+{bigB}); v_p(det M_4) predicted {factorint(dB0).get(p, 0) + 3 * F.get(p, 0)} "
        f"({time.time() - t1:.1f}s)")
say(f"total {time.time() - t0:.1f}s")
