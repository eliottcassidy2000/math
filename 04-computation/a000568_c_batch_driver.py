#!/usr/bin/env python3
"""
Driver for the C batch CRT computation of A000568.
Calls the C binary and performs CRT reconstruction in Python.

Author: opus-2026-03-07-S46f
"""

import subprocess
import sys
import os
from time import perf_counter

def a000568_c_batch(n, verbose=True):
    """Compute A000568(n) using the C batch binary + Python CRT."""
    if n <= 1:
        return 1

    c_binary = os.path.join(os.path.dirname(os.path.abspath(__file__)), "a000568_c_batch")
    result = subprocess.run(
        [c_binary, str(n)],
        capture_output=True, text=True
    )

    if verbose:
        sys.stderr.write(result.stderr)

    lines = result.stdout.strip().split('\n')
    num_primes = int(lines[0])

    primes = []
    residues = []
    for i in range(1, num_primes + 1):
        p, r = map(int, lines[i].split())
        primes.append(p)
        residues.append(r)

    # CRT reconstruction
    result_val = residues[0]
    modulus = primes[0]
    for i in range(1, len(residues)):
        r = residues[i]
        p = primes[i]
        diff = (r - result_val % p) % p
        t = diff * pow(modulus, p - 2, p) % p
        result_val = result_val + modulus * t
        modulus = modulus * p

    if result_val < 0:
        result_val += modulus

    return result_val


known = {0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
         9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168}

if __name__ == "__main__":
    print("=" * 70)
    print("A000568 via C BATCH CRT")
    print("=" * 70, flush=True)

    print("\nValidation:", flush=True)
    for n in range(13):
        t0 = perf_counter()
        result = a000568_c_batch(n, verbose=False)
        dt = perf_counter() - t0
        expected = known.get(n, "?")
        match = "OK" if result == expected else f"FAIL (got {result}, expected {expected})"
        print(f"  a({n:2d}) = {result:>15}  [{dt:.4f}s]  {match}", flush=True)

    oeis_76 = 457153791459731763873714717642097124535855865100705012123154296782586066392599174177717111708812242945048119323356894628050571858925201215481427701704774498054845450499257140261101951670038378795297655176726665694457285914072690084154974004045609532979011190165711892081116449736469812420129338461848172656662018517999846546136461519522876668617723623760713964531731492213163219661503132797160477050776498667637823035846487741247681442428085385480808545797005707547708306606906039466675985395538562415540138435374979216235780263415436443048770796035424757246769390586508396300108786258541736828941774383334172760307820994535160154251942298642998619446487984483719749964075470321740294755571922446326714092505683852691999256324957470064758984015872

    print("\n" + "=" * 70)
    print("SCALING TEST")
    print("=" * 70, flush=True)

    for n_test in [20, 30, 40, 50, 60, 70, 76, 80, 90, 100]:
        t0 = perf_counter()
        result = a000568_c_batch(n_test, verbose=True)
        dt = perf_counter() - t0
        s = str(result)
        extra = ""
        if n_test == 76:
            extra = f"  OEIS={'OK' if result == oeis_76 else 'FAIL'}"
        print(f"  a({n_test:3d}) = {s[:40]}{'...' if len(s) > 40 else '':3s} "
              f"({len(s):4d} digits)  [{dt:8.3f}s]{extra}", flush=True)
