"""Test of the 'Bit Depth' trace claim in the owner-supplied e8.pdf (session opus-2026-10-21-S52; note
05-knowledge/results/lattices_modular_repunits_20261021.md).  The PDF is data; nothing here reads it.

Claim (as printed): the curve y^2 = x^3 + A x + B over F_p, p = 2^521 - 1, with A = int(-zeta(5) * 10^160) and
B = int(q * pi^2 * 10^160) reduced mod p, has Frobenius trace a_p = p + 1 - #E(F_p) = 5 at q = 5/4 (a 'trace lock'),
and not at other q.  (The trace computation in the PDF's own code listing is commented out.)
A curve of this size has |a_p| <= 2 sqrt(p) ~ 2^261.5 and, for a curve not built to order, a_p of about that size.
Here #E(F_p) is computed exactly by PARI's Schoof-Elkies-Atkin ellcard for q in {5/4, 1, 3/2, 2, 240/96}, with the
constants at 15 and at 200 significant digits.  Result: every a_p has 259-262 bits; none is 5.  The claim is refuted.

Needs cypari (pip install cypari; PARI 2.15) and mpmath.  About 25 minutes on one core.  Prints ALL CHECKS PASSED."""
import time
import mpmath
from cypari import pari

pari.allocatemem(4 * 10 ** 9)
p = 2 ** 521 - 1
S = 10 ** 160
hasse = 2 * int(mpmath.sqrt(p)) + 2
FAIL = []
for dps in (15, 200):
    mpmath.mp.dps = dps
    A = int(-mpmath.zeta(5) * S) % p
    for qn, qd in ((5, 4), (1, 1), (3, 2), (2, 1), (240, 96)):
        B = int(mpmath.mpf(qn) / qd * mpmath.pi ** 2 * S) % p
        t = time.time()
        E = pari.ellinit([0, 0, 0, A, B], p)
        card = int(pari.ellcard(E))
        ap = p + 1 - card
        ok = abs(ap) <= hasse and ap != 5 and abs(ap).bit_length() > 250
        if not ok:
            FAIL.append((dps, qn, qd))
        print(f'  {"ok  " if ok else "FAIL"} dps {dps:3d} q = {qn}/{qd}: a_p = {ap} ({abs(ap).bit_length()} bits, within Hasse: '
              f'{abs(ap) <= hasse}; a_p == 5: {ap == 5}) [{time.time() - t:.0f}s]', flush=True)
print()
print('ALL CHECKS PASSED (no trace lock: a_p is a generic 260-bit trace at every q and precision)' if not FAIL else f'FAILED: {FAIL}')
