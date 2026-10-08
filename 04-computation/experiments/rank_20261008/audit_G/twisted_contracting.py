#!/usr/bin/env python3
"""Audit G: contracting rank-1 map with a twisted mod-2 obstruction and NO constant-c obstruction for any prime l:
Z_3, m = (1,1,5), r = (0,2,5).  psi(5) = 1 (r_2 odd), psi(1) = 0 (r_0, r_1 even): e_t = e_0 + A_t (mod 2), A = exponent of 5.
Prediction: odd offsets never merge; even offsets merge a.s. (rank 1 contracting, Lambda = ln(5/27)/3 = -0.562)."""
import math
from mwsim import MW, run_pairs
mw = MW(3, [1, 1, 5], [0, 2, 5])
# constant-c lemma for any prime l <= 2000: (m_i - d) c + r_i = 0 mod l for all i?
def primes(n): return [p for p in range(2, n + 1) if all(p % q for q in range(2, int(p ** .5) + 1))]
hits = [l for l in primes(2000) if (3 * 1 * 5) % l and any(all(((mi - 3) * c + ri) % l == 0 for mi, ri in zip(mw.m, mw.r)) for c in range(l))]
print("constant-c obstruction primes l <= 2000 (l not dividing d*prod m):", hits)
for e in (1, 3, -1, 2, 4, 6):
    res = run_pairs(mw, e, 4096, 400, 900 + e, [16, 64, 256, 1024, 4096], window_visits=False)
    q = {c: res['merged_by'][c] / 400 for c in (16, 64, 256, 1024, 4096)}
    print(f"e = {e:3d}: merged-by {q}; sqrt(T)*P(no merge) at 1024, 4096: {math.sqrt(1024)*(1-q[1024]):.2f}, {math.sqrt(4096)*(1-q[4096]):.2f}", flush=True)
