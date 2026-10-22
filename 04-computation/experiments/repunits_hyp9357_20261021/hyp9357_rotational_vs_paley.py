"""HYP-9357 test (session opus-2026-10-21-S51): the rotational circulant R_q = Z_q{1, ..., (q-1)/2} against the Paley
tournament P_q at q = 27 and q = 31.  Note: 05-knowledge/results/repunit_rarity_hyp9357_20261021.md.

HYP-9357 (mac-mini-2026-10-21-fanopsp, prediction (3), committed before this script existed) predicts that the
rotational circulant beats P_q for every prime q >= 19, q = 3 mod 4, and lists H(P_27), H(P_31) against the rotational
circulant as its cheapest decisive test.

Method: Karp inclusion-exclusion, H = sum_S (-1)^(n-|S|) W_(n-1)(S), with W the number of walks with n-1 arcs inside
the induced subtournament on S (vector iteration, additions only), modulo three 58-bit primes with CRT.  The kernels
are imported unchanged from the audited S34 code (04-computation/experiments/paley_h31_audit_20261014/
audit_H_karp_ie.py).  W is invariant under the rotations of a circulant, so:
  * prime n: one canonical representative per rotation orbit (all orbits of size n except the empty and the full set);
  * n = 15, 21, 27: the identity sum_(S nonempty) f(S) = n sum_(S containing 0) f(S)/|S|.
Validation: the rotational values for n = 13..23 against the concurrent circ_counts.c table (HYP-9357) and an
independent int64 Held-Karp for n <= 17; the Paley values H(P_27), H(P_31) are the audited S34 values.
Prints ALL CHECKS PASSED.  Runtime about 10 minutes (n = 31 dominates)."""
import sys, os, time, math
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'paley_h31_audit_20261014'))
import audit_H_karp_ie as K                                              # noqa: E402

FAIL = []
T0 = time.time()


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


def rotational(n):
    half = set(range(1, (n - 1) // 2 + 1))
    return [[((y - x) % n) in half for y in range(n)] for x in range(n)]


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def H_circulant(arc):
    n = len(arc)
    K.check_tournament(arc)
    inm = K.inmasks(arc)
    if is_prime(n):
        part = K.karp_rotation(n, inm, 2048)
        res = [int(part[:, j].astype(object).sum()) % K.MODS[j] for j in range(3)]
        nsub = int(part[:, 3].sum())
        assert nsub == (2 ** n - 2) // n
        # the full set: the circulant is regular, W(V) = n ((n-1)/2)^(n-1)
        Wfull = n * ((n - 1) // 2) ** (n - 1)
        res = [(res[j] + Wfull) % K.MODS[j] for j in range(3)]
    else:
        winv = []
        for p in K.MODS:
            arr = np.zeros(n + 1, dtype=np.int64)
            for k in range(1, n + 1):
                arr[k] = (n * pow(k, -1, p)) % p
            winv.append(arr)
        part = K.karp_contains0(n, inm, winv[0], winv[1], winv[2], 2048)
        res = [int(part[:, j].astype(object).sum()) % K.MODS[j] for j in range(3)]
    H, M = K.crt2(res[0], res[1], K.MODS[0], K.MODS[1])
    return H, (H % K.MODS[2]) == res[2], H < (M >> 20)


known_rot = {13: 3711175, 15: 198464295, 17: 13689269499, 19: 1184212824763, 21: 125547534942879, 23: 16011537490557279}
paley = {19: 1172695746915, 23: 15760206976379349, 27: 418425945688463432031, 31: 19876729650849384218015865}

print('(V) validation of the method on the rotational circulants n = 13..23')
for n in (13, 15, 17, 19, 21, 23):
    t = time.time()
    H, ok3, small = H_circulant(rotational(n))
    extra = ''
    if n <= 17:
        hk = int(K.held_karp_exact(n, K.outmasks(rotational(n))))
        extra = f'; Held-Karp {hk} (match {hk == H})'
        ok3 = ok3 and hk == H
    check(H == known_rot[n] and ok3 and small, f'n = {n}: H(R_{n}) = {H} = circ_counts value{extra} [{time.time() - t:.1f}s]')

print('(T) the decisive test at q = 27, 31 (HYP-9357 prediction (3))')
res = {}
for n in (27, 31):
    t = time.time()
    H, ok3, small = H_circulant(rotational(n))
    res[n] = H
    check(ok3 and small and H % 2 == 1 and H % n == 0, f'n = {n}: H(R_{n}) = {H} (third modulus consistent, odd, divisible by {n}) [{time.time() - t:.0f}s]')
for n in (19, 23, 27, 31):
    hr = known_rot.get(n, res.get(n)); hp = paley[n]
    rr = float(Fraction(hr * 2 ** (n - 1), math.factorial(n))); rp = float(Fraction(hp * 2 ** (n - 1), math.factorial(n)))
    print(f'     n = {n}: H(R_n) = {hr}, H(P_n) = {hp}; ratio H(R)/H(P) = {hr / hp:.6f}; '
          f'H 2^(n-1)/n!: rotational {rr:.5f}, Paley {rp:.5f}; rotational wins: {hr > hp}')
verdict27 = res[27] > paley[27]; verdict31 = res[31] > paley[31]
print(f'     HYP-9357 prediction (3) at q = 27 (prime power, outside its prime statement): {"rotational wins" if verdict27 else "PALEY WINS"}')
print(f'     HYP-9357 prediction (3) at q = 31: {"PASSED (rotational wins)" if verdict31 else "FAILED (Paley wins)"}')
check(True, 'decisive comparison computed (verdict printed above; the script does not assert the hypothesis)')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
