"""Residues of H over the full class census n <= 9 (session opus-2026-10-21-S51; HYP-9359, note section 3).
Class counts and labelled counts (weight n!/|Aut|) of H mod 3, 8, 16; the labelled excess
#{H = 1 mod 4} - #{H = 3 mod 4} (compare h_patterns.out (Q4)).  Prints ALL CHECKS PASSED."""
import os, sys, math, time
from collections import Counter
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'h_spectrum_multiplicities_20261021'))
import h_multiplicities as HM                         # noqa: E402

FAIL = []
T0 = time.time()


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


data = HM.classes_upto(9)
D4 = {}
for n in range(3, 10):
    A = data[n]['A']; aut = data[n]['aut']
    H, st = HM.ham_all(n, A)
    tot = 2 ** (n * (n - 1) // 2)
    lab = {M: Counter() for M in (3, 8, 16)}
    cls = {M: Counter() for M in (3, 8, 16)}
    for h, a in zip(H, aut):
        w = math.factorial(n) // int(a)
        for M in (3, 8, 16):
            lab[M][int(h) % M] += w
            cls[M][int(h) % M] += 1
    D4[n] = sum(c for r, c in lab[16].items() if r % 4 == 1) - sum(c for r, c in lab[16].items() if r % 4 == 3)
    print(f'     n = {n}: classes by H mod 3 {dict(sorted(cls[3].items()))}; labelled share by H mod 3 '
          + ' '.join(f'{r}:{3 * lab[3][r] / tot:.4f}' for r in range(3)) + '; labelled share by H mod 8 '
          + ' '.join(f'{r}:{4 * lab[8][r] / tot:.4f}' for r in (1, 3, 5, 7)) + f'; #{{H=1 mod 4}} - #{{H=3 mod 4}} = {D4[n]}')
check([D4[n] for n in range(3, 10)] == [4, 32, 416, 8192, 218368, 11149312, 516939776], 'labelled excess of H = 1 mod 4 matches h_patterns.out (Q4)')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
