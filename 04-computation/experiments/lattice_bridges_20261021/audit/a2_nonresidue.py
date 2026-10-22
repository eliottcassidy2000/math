#!/usr/bin/env python3
"""AUDIT of THM-4751 (2): B_q <= 1 iff q in {3, 7, 11}; B_q >= 2 for q >= 19.  Checks each step of the elementary proof
for all primes q = 3 mod 4 below LIM:
  q = 3 mod 8 : 2 and 8 are non-residues and both <= m = (q-1)/2 once q >= 17;
  q = 7 mod 8 : n0 = least non-residue; m0 = ceil(q/n0) satisfies 0 < m0*n0 - q < n0, m0*n0 - q is a residue, m0 is a
                non-residue, m0 >= n0, hence n0 < 1/2 + sqrt(q + 1/4) (i.e. n0(n0-1) < q); 2*n0 is a non-residue and
                2*n0 <= m for q >= 23 (the bound alone suffices from q = 23 on; THM-4751 states q >= 31 + direct q = 23).
Output: a2_nonresidue.out
"""
import os, math
HERE = os.path.dirname(os.path.abspath(__file__))
LIM = 2_000_000
def primes_upto(N):
    s = bytearray([1]) * (N + 1); s[0] = s[1] = 0
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(N + 1) if s[i]]
def isqr(a, q): return pow(a % q, (q - 1) // 2, q) == 1
out = open(os.path.join(HERE, 'a2_nonresidue.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
fails = []; Bsmall = []; first_bound_ok = None; n0max = (0, 0)
for q in primes_upto(LIM):
    if q % 4 != 3: continue
    m = (q - 1) // 2
    if q % 8 == 3:
        if q >= 19 and not (not isqr(2, q) and not isqr(8, q) and 8 <= m): fails.append(('3mod8', q))
    else:
        n0 = next(a for a in range(2, q) if not isqr(a, q))
        n0max = max(n0max, (n0, q))
        m0 = -(-q // n0)
        r = m0 * n0 - q
        ok = (0 < r < n0) and isqr(r, q) and (not isqr(m0, q)) and m0 >= n0 and n0 * (n0 - 1) < q and n0 < 0.5 + math.sqrt(q + 0.25)
        ok = ok and (not isqr(2 * n0, q))
        if not ok: fails.append(('7mod8-step', q))
        if q >= 23 and 2 * n0 > m: fails.append(('7mod8-2n0>m', q))
        # does the bound alone give 2*n0 <= m?  (largest integer below 1/2 + sqrt(q + 1/4), doubled)
        nb = math.floor(0.5 + math.sqrt(q + 0.25) - 1e-12)
        if 2 * nb <= m and first_bound_ok is None: first_bound_ok = q
    if q < 200000:
        B = sum(1 for k in range(1, m + 1) if not isqr(k, q)) if q < 20000 else None
        if B is not None and B <= 1: Bsmall.append((q, B))
say(f"primes q = 3 mod 4 below {LIM}: proof-step failures = {fails}")
say(f"q with B_q <= 1 (direct count, q < 20000): {Bsmall}")
say(f"first q = 7 mod 8 for which the bound n0 < 1/2 + sqrt(q+1/4) alone forces 2 n0 <= m: q = {first_bound_ok}")
say(f"largest least non-residue among q = 7 mod 8 below {LIM}: n0 = {n0max[0]} at q = {n0max[1]}")
say("q = 23: n0 =", next(a for a in range(2, 23) if not isqr(a, 23)), "; 10 is a non-residue:", not isqr(10, 23), "; m = 11")
out.close()
