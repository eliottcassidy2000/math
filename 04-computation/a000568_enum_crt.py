#!/usr/bin/env python3
"""
a000568_enum_crt.py — A000568 via C enumeration + parallel CRT.

Uses the C enumeration program for per-prime mod-p computation,
parallelized across CPU cores, with CRT reconstruction in Python.

Author: opus-2026-03-07-S46f
"""

import os
import math
import subprocess
from time import perf_counter
from concurrent.futures import ProcessPoolExecutor

_dir = os.path.dirname(os.path.abspath(__file__))
_binary = os.path.join(_dir, "a000568_c_enum")


def _compute_one(args):
    """Compute a(n) mod p using C enumeration."""
    n, p = args
    result = subprocess.check_output([_binary, str(n), str(p)])
    return int(result.strip())


def primes_above(start, count):
    """Generate `count` primes starting above `start`."""
    primes = []
    n = start + 1
    if n % 2 == 0:
        n += 1
    while len(primes) < count:
        is_p = True
        i = 3
        while i * i <= n:
            if n % i == 0:
                is_p = False
                break
            i += 2
        if n > 1 and is_p:
            primes.append(n)
        n += 2
    return primes


def a000568_enum_crt(n, verbose=True, num_workers=None):
    """Compute A000568(n) using parallel C enumeration + CRT."""
    if n <= 1:
        return 1

    est_bits = int(n * (n - 1) / 2 - n * math.log2(max(n, 2)) + n) + 100
    if est_bits < 100:
        est_bits = 100

    prime_bits = 30
    num_primes = est_bits // prime_bits + 5

    primes = primes_above(max(n, 2**30), num_primes)

    if num_workers is None:
        num_workers = os.cpu_count() or 4

    if verbose:
        print(f"  n={n}: {num_primes} primes, {num_workers} workers, est {est_bits} bits",
              flush=True)

    t0 = perf_counter()
    args_list = [(n, p) for p in primes]

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        residues = list(executor.map(_compute_one, args_list))

    t_par = perf_counter() - t0
    if verbose:
        print(f"  Parallel enum: {t_par:.1f}s ({t_par/num_primes:.3f}s/prime avg)", flush=True)

    # CRT reconstruction
    result = residues[0]
    modulus = primes[0]
    for i in range(1, len(residues)):
        r = residues[i]
        p = primes[i]
        diff = (r - result % p) % p
        t = diff * pow(modulus, p - 2, p) % p
        result = result + modulus * t
        modulus = modulus * p

    if result < 0:
        result += modulus

    return result


known = {0: 1, 1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
         9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168}

oeis_76 = 457153791459731763873714717642097124535855865100705012123154296782586066392599174177717111708812242945048119323356894628050571858925201215481427701704774498054845450499257140261101951670038378795297655176726665694457285914072690084154974004045609532979011190165711892081116449736469812420129338461848172656662018517999846546136461519522876668617723623760713964531731492213163219661503132797160477050776498667637823035846487741247681442428085385480808545797005707547708306606906039466675985395538562415540138435374979216235780263415436443048770796035424757246769390586508396300108786258541736828941774383334172760307820994535160154251942298642998619446487984483719749964075470321740294755571922446326714092505683852691999256324957470064758984015872

if __name__ == "__main__":
    print("=" * 70)
    print("A000568 via parallel C enumeration + CRT")
    print(f"CPU cores: {os.cpu_count()}")
    print("=" * 70, flush=True)

    print("\nValidation:", flush=True)
    for n_val in range(13):
        t0 = perf_counter()
        result = a000568_enum_crt(n_val, verbose=False, num_workers=1)
        dt = perf_counter() - t0
        expected = known.get(n_val, "?")
        match = "OK" if result == expected else f"FAIL (got {result})"
        print(f"  a({n_val:2d}) = {result:>15}  [{dt:.4f}s]  {match}", flush=True)

    print("\n" + "=" * 70)
    print("SCALING TEST")
    print("=" * 70, flush=True)

    for n_test in [20, 30, 40, 50, 60, 70, 76, 80, 90, 100, 120, 150]:
        t0 = perf_counter()
        result = a000568_enum_crt(n_test, verbose=True)
        dt = perf_counter() - t0
        s = str(result)
        extra = ""
        if n_test == 76:
            extra = f"  OEIS={'OK' if result == oeis_76 else 'FAIL'}"
        print(f"  a({n_test:3d}) ({len(s):4d} digits)  [{dt:8.3f}s]{extra}", flush=True)
