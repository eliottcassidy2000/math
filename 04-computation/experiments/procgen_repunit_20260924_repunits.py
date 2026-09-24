#!/usr/bin/env python3
"""Part 3: digit rotation, circular and permutable primes, repunit primes.

 (3a) Rotation of k base-b digits = multiplication by b modulo b^k - 1 (exhaustive/random checks);
      fixed points = repdigits c*R_k; orbit sizes; every prime with k >= 2 digits has a free
      rotation orbit of size k unless it is a repunit (checked for all primes < 10^7).
 (3b) Circular primes: C helper procgen_repunit_20260924_circular.c, sieve to 10^9 (every rotation
      checked) and an FKM necklace search over {1,3,7,9} for lengths 10..16; comparison with
      OEIS A068652 / A016114 / A293663 and P. De Geest's table (read 2026-09-24); permutable primes
      (A003459) to 10^9.
 (3c) Repunit primes R_n, n <= 1100: gmpy2 BPSW (is_bpsw_prp) and Miller-Rabin (is_prime, 25 rounds);
      composite n gives R_d | R_n.  Prime factors of R_p are = 1 mod 2p (checked on found factors).
 (3d) Permutable primes beyond 10^9: every multiset of 4..40 digits from {1,3,7,9} other than the
      all-ones multiset has an explicit composite permutation; for near-repdigits x..xy of length
      4 <= n <= 20000 an explicit small prime q divides some permutation (Richert's mechanism).
 (3e) The finiteness heuristic for non-repunit circular primes, with the local factors at 3 and at the
      primes dividing 10^k - 1, against the actual counts per length.
Session collatz-procgen-20260922, lane procgen_repunit_20260924.  Runtime ~ 4-6 min (necklace search
to length 16 dominates), memory < 100 MB.
"""
import os, sys, subprocess, tempfile, math, time, random, itertools
from math import gcd, log
import gmpy2
from gmpy2 import mpz
import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
TMP = tempfile.mkdtemp(prefix="procgen_repunit_circ_")
BIN = os.path.join(TMP, "circular")
SIEVE_LIMIT = int(float(os.environ.get("PR_CIRC_SIEVE", "1e9")))
NECK_MAX = int(os.environ.get("PR_CIRC_NECKMAX", "16"))

# OEIS data read 2026-09-24 (JSON API)
A016114 = [2, 3, 5, 7, 11, 13, 17, 37, 79, 113, 197, 199, 337, 1193, 3779, 11939, 19937, 193939, 199933]
A068652 = [2, 3, 5, 7, 11, 13, 17, 31, 37, 71, 73, 79, 97, 113, 131, 197, 199, 311, 337, 373, 719, 733, 919, 971,
           991, 1193, 1931, 3119, 3779, 7793, 7937, 9311, 9377, 11939, 19391, 19937, 37199, 39119, 71993, 91193,
           93719, 93911, 99371, 193939, 199933, 319993]           # the JSON data field stops at 319993
A003459 = [2, 3, 5, 7, 11, 13, 17, 31, 37, 71, 73, 79, 97, 113, 131, 199, 311, 337, 373, 733, 919, 991]
A004023 = [2, 19, 23, 317, 1031, 49081, 86453, 109297, 270343, 5794777, 8177207]
DEGEEST_ELIGIBLE = {6: 757, 7: 2709, 8: 9177, 9: 33191}   # "eligible primes" per length, De Geest's table


def rot_left(n, k, b):
    """move the leading base-b digit of the k-digit string of n (leading zeros allowed) to the end"""
    lead = n // b ** (k - 1)
    return (n - lead * b ** (k - 1)) * b + lead


def part3a():
    print("\n(3a) Rotation = multiplication by b modulo b^k - 1")
    rng = random.Random(20260924)
    checked = 0
    for b in (2, 3, 4, 10, 12, 16):
        for k in range(1, 31):
            M = b ** k - 1
            for _ in range(200):
                n = rng.randrange(0, b ** k)
                r = rot_left(n, k, b)
                # rot(n) = b*n mod (b^k-1), with the all-(b-1) string b^k-1 fixed (its class is 0)
                assert r == (b * n % M if n != M else M) or (n == 0 and r == 0)
                checked += 1
    print(f"    rot(n) = b*n mod (b^k - 1) on {checked} random (b, k, n), b in {{2,3,4,10,12,16}}, k <= 30: all equal")
    # fixed points and orbit-size census (strings with leading zeros), b = 10, k <= 7
    for k in range(1, 8):
        M = 10 ** k - 1
        fixed = [n for n in range(10 ** k) if rot_left(n, k, 10) == n]
        assert fixed == [c * (M // 9) for c in range(10)]
    print("    fixed points of rotation on k-digit strings (k <= 7, exhaustive): exactly the 10 repdigits c*R_k;"
          " in Z/(10^k-1) the fixed classes are the b-1 = 9 solutions of (b-1)n = 0.")
    # every prime with >= 2 digits has a free orbit of size k, except repunit primes
    N = 10 ** 7
    s = bytearray([1]) * (N + 1)
    s[0] = s[1] = 0
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    bad = []
    for p in range(11, N, 2):
        if not s[p]:
            continue
        k = len(str(p))
        r, size = p, 0
        for j in range(1, k + 1):
            r = rot_left(r, k, 10)
            if r == p:
                size = j
                break
        if size != k:
            bad.append(p)
    print(f"    primes 10 < p < 10^7 whose rotation orbit is smaller than their number of digits: {bad}")
    assert bad == [11]
    # Fermat via necklaces: primitive necklaces of prime length p number (b^p - b)/p
    for b in (2, 3, 4, 10):
        for p in (2, 3, 5, 7, 11, 13):
            assert (b ** p - b) % p == 0
    print("    for prime length p the b^p strings split into the b fixed repdigits and (b^p - b)/p free orbits"
          " (Fermat's little theorem, necklace form).")


def part3b():
    print("\n(3b) Circular and permutable primes")
    subprocess.run(["cc", "-O2", "-o", BIN, os.path.join(HERE, "procgen_repunit_20260924_circular.c")], check=True)
    t0 = time.time()
    out = subprocess.run([BIN, "sieve", str(SIEVE_LIMIT)], capture_output=True, text=True, check=True).stdout.splitlines()
    print(f"    sieve to {SIEVE_LIMIT:.0e}: {time.time() - t0:.0f} s")
    circ = [int(l.split()[1]) for l in out if l.startswith("CIRC")]
    perm = [int(l.split()[1]) for l in out if l.startswith("PERM")]
    near = [l.split()[1:] for l in out if l.startswith("NEAR")]
    lens = {int(l.split()[1]): list(map(int, l.split()[2:])) for l in out if l.startswith("LEN")}
    pb = [l for l in out if l.startswith("PRIMES_BELOW")][0]
    print(f"    {pb}")
    print(f"    circular primes below {SIEVE_LIMIT:.0e}: {len(circ)} numbers; largest {max(circ)}")
    print("      " + " ".join(map(str, circ)))
    # orbit representatives: the smallest member of each rotation orbit
    orbits = {}
    for c in circ:
        k = len(str(c))
        o = [c]
        r = c
        for _ in range(k - 1):
            r = rot_left(r, k, 10)
            o.append(r)
        orbits[min(o)] = sorted(set(o))
    orep = sorted(orbits)
    print(f"    orbit representatives (smallest member): {orep}")
    assert orep == A016114, "A016114 mismatch"
    assert circ[:len(A068652)] == A068652, "A068652 mismatch"
    print("      = OEIS A016114 below 10^9 exactly; the full list = A068652 (its JSON data, 46 terms to 319993,"
          " continues 331999, ..., 999331 here: 55 numbers in all).")
    nonrep = [c for c in circ if c != 11]
    print(f"    non-repunit circular primes: {len(nonrep)} (A293663 comment: 54 terms with <= 6 digits, largest 999331)")
    assert len(nonrep) == 54 and max(nonrep) == 999331
    print(f"    permutable primes below {SIEVE_LIMIT:.0e}: {perm}")
    assert perm == A003459
    print("      = OEIS A003459 below 10^9 exactly (the next terms there are R19, R23).")
    print("    per length k: eligible primes (all digits in {1,3,7,9}), circular primes, orbits, near-miss orbits"
          " (exactly one composite rotation):")
    for k in sorted(lens):
        e, c, o, nm = lens[k]
        tag = ""
        if k in DEGEEST_ELIGIBLE:
            tag = f"   De Geest: {DEGEEST_ELIGIBLE[k]} eligible"
            assert e == DEGEEST_ELIGIBLE[k]
        print(f"      k={k}: eligible {e:>7}, circular {c:>3}, orbits {o:>2}, near misses {nm:>3}{tag}")
    for k, p, comp in near:
        if int(k) >= 7:
            print(f"      near miss, {k} digits: orbit of {p}; the single composite rotation {comp} = "
                  f"{sympy.factorint(int(comp))}")
    # necklace search beyond 10^9
    t0 = time.time()
    proc = subprocess.run([BIN, "necklace", "2", str(NECK_MAX)], capture_output=True, text=True, check=True)
    lines = proc.stdout.splitlines()
    print(f"    necklace search, lengths 2..{NECK_MAX}: {time.time() - t0:.0f} s")
    ncirc = [tuple(map(int, l.split()[1:])) for l in lines if l.startswith("NCIRC")]
    nl = {int(l.split()[1]): list(map(int, l.split()[2:])) for l in lines if l.startswith("NLEN")}
    for k in sorted(nl):
        neck, tested, c, nm = nl[k]
        print(f"      k={k:>2}: necklaces {neck:>11,}, tested (digit sum not = 0 mod 3) {tested:>11,},"
              f" circular orbits {c}, near-miss orbits {nm}")
        if k in lens and k >= 2:
            assert c == lens[k][2] and nm == lens[k][3], (k, nl[k], lens[k])
    assert [v for k, v in ncirc] == [x for x in A016114 if x >= 10]
    assert all(nl[k][2] == 0 for k in nl if k >= 7)
    print(f"    => the circular orbits of every length 2..{NECK_MAX} are exactly A016114's; none of length 7..{NECK_MAX}"
          f" (FINITE-EXACT; the sieve and the necklace code agree on every length <= 9).")
    for l in lines:
        if l.startswith("NNEAR") and int(l.split()[1]) >= 10:
            k, v, comp = l.split()[1:]
            print(f"      near miss, {k} digits: representative {v}; composite rotation {comp} = {sympy.factorint(int(comp))}")
    return lens, nl


def part3c():
    print("\n(3c) Repunit primes R_n = (10^n - 1)/9, n <= 1100")
    t0 = time.time()
    primes_found = []
    comp_certs = 0
    factor_found = 0
    factor_ok = True
    for n in range(2, 1101):
        R = (mpz(10) ** n - 1) // 9
        if not sympy.isprime(n):
            d = min(q for q in range(2, n) if n % q == 0)
            Rd = (mpz(10) ** d - 1) // 9
            assert R % Rd == 0 and 1 < Rd < R
            comp_certs += 1
            continue
        bp = gmpy2.is_bpsw_prp(R)
        mr = gmpy2.is_prime(R, 25)
        assert bp == mr
        if bp:
            primes_found.append(n)
        else:
            # small factors q = 2pm + 1 (every prime factor of R_p, p != 3, is 1 mod 2p)
            if n != 3:
                for m in range(1, 20000):
                    q = 2 * n * m + 1
                    if gmpy2.is_prime(q) and R % q == 0:
                        factor_found += 1
                        break
                # check the congruence on the factors found by trial division up to 10^5 (any q)
            else:
                assert R == 3 * 37
    print(f"    prime n with R_n (BPSW and 25-round Miller-Rabin agree): {primes_found}   ({time.time() - t0:.0f} s)")
    assert primes_found == [2, 19, 23, 317, 1031]
    print(f"    composite n: R_d | R_n with 1 < R_d < R_n for the least prime d | n ({comp_certs} values of n).")
    nprime = sum(1 for n in range(2, 1101) if sympy.isprime(n))
    print(f"    of the {nprime - 5} composite R_p (p prime <= 1100), {factor_found} have a factor q = 2pm+1 with m < 20000.")
    # structural check of the congruence on every small factor
    viol = 0
    for p in sympy.primerange(5, 200):
        R = (10 ** p - 1) // 9
        for q in sympy.primerange(3, 10 ** 5):
            if R % q == 0 and q != 3:
                if (q - 1) % (2 * p) != 0:
                    viol += 1
    print(f"    every prime factor q < 10^5 of R_p (5 <= p < 200) satisfies q = 1 mod 2p: violations = {viol}")
    assert viol == 0
    print(f"    OEIS A004023 (read): {A004023}; the first five are proved primes, confirmed here by BPSW.")


def part3d():
    print("\n(3d) Permutable primes beyond 10^9")
    t0 = time.time()
    small_primes = list(sympy.primerange(7, 2000))
    nms = 0
    for L in range(4, 41):
        for c1 in range(L + 1):
            for c3 in range(L + 1 - c1):
                for c7 in range(L + 1 - c1 - c3):
                    c9 = L - c1 - c3 - c7
                    counts = {1: c1, 3: c3, 7: c7, 9: c9}
                    nms += 1
                    if c1 == L:                      # the repunit: its only permutation
                        continue
                    if (c1 + c7) % 3 == 0:           # digit sum = 0 mod 3: every permutation divisible by 3
                        continue
                    digits = [d for d in (1, 3, 7, 9) for _ in range(counts[d])]
                    # find a permutation with a small prime factor, or a composite one
                    found = False
                    rng = random.Random(L * 1000003 + c1 * 1009 + c3 * 101 + c7)
                    for _ in range(200):
                        rng.shuffle(digits)
                        v = int("".join(map(str, digits)))
                        if not gmpy2.is_prime(v, 5):
                            found = True
                            break
                    assert found, (L, counts)
    print(f"    all {nms} digit multisets of length 4..40 over {{1,3,7,9}}: every one except the all-ones multiset"
          f" has a composite permutation ({time.time() - t0:.0f} s).  Hence the permutable primes with 4..40 digits"
          f" are exactly R19 and R23.")
    # Richert's mechanism for near-repdigits x^(n-1) y, 4 <= n <= NMAX: find q and j with q | x R_n + (y-x) 10^j
    t0 = time.time()
    NMAX = 20000
    qs = list(sympy.primerange(7, 5000))
    ordq = {q: sympy.n_order(10, q) for q in qs}
    worst = 0
    for n in range(4, NMAX + 1):
        for x in (1, 3, 7, 9):
            for y in (1, 3, 7, 9):
                if x == y:
                    continue
                if (x * (n - 1) + y) % 3 == 0:
                    continue                         # every permutation divisible by 3
                ok = False
                for q in qs:
                    d = min(ordq[q], n)                      # positions j = 0..n-1; 10^j has period ord_q(10)
                    Rn = (pow(10, n, 9 * q) - 1) // 9 % q
                    t = x * Rn % q
                    if (y - x) % q == 0:
                        continue
                    # want 10^j = -t/(y-x) mod q for some 0 <= j < d
                    target = (-t) * pow(y - x, -1, q) % q
                    p10 = 1
                    for j in range(d):
                        if p10 == target:
                            ok = True
                            worst = max(worst, q)
                            break
                        p10 = p10 * 10 % q
                    if ok:
                        assert q < 10 ** (n - 1)             # the permutation exceeds q, hence is composite
                        break
                assert ok, (n, x, y)
    print(f"    near-repdigits x..xy (x != y in {{1,3,7,9}}) of every length 4 <= n <= {NMAX}: some permutation is divisible"
          f" by a prime q <= {worst} (and exceeds it) ({time.time() - t0:.0f} s).")
    print("    Together with Johnson (1977; CITED via Wikipedia, not read: no permutable prime has three distinct digits,"
          " or two digits each used twice), this reproduces the no-non-repunit statement up to 20000 digits;"
          " Richert (1951; CITED via Wikipedia/De Geest, not read) states it for 3 < n < 6*10^175.")


def part3e(lens, nl):
    print("\n(3e) Why finitely many non-repunit circular primes are expected (heuristic, not a theorem)")
    print("    A k-digit circular prime (k >= 2) uses only the digits 1,3,7,9: a digit 0,2,4,5,6,8 rotated to the end")
    print("    makes that rotation divisible by 2 or 5.  Model: the k rotations are prime independently with")
    print("    probability 3.75/log N (coprime to 30), corrected at 3 (all rotations share the digit sum) and at every")
    print("    prime q | 10^k - 1 (rotation multiplies by 10 mod q, so the rotations are all or none divisible by q).")

    def s3(k):
        # digits 1,7 are 1 mod 3; 3,9 are 0 mod 3: digit sum = #{1,7} mod 3, #{1,7} ~ Bin(k,1/2)
        return sum(math.comb(k, j) for j in range(k + 1) if j % 3 != 0) / 2 ** k

    ords = {int(q): int(sympy.n_order(10, int(q))) for q in sympy.primerange(7, 10 ** 5)}

    def boost(k):
        # primes q < 10^5, q != 3, dividing 10^k - 1 (i.e. ord_q(10) | k); larger q change the product by < 1e-3
        f = 1.0
        for q, d in ords.items():
            if k % d == 0:
                f *= (1 - 1 / q) ** (1 - k)
        return f

    actual = {k: str(lens[k][1] - (1 if k == 2 else 0)) for k in lens if k >= 2}
    for k in nl:
        if k not in actual:
            actual[k] = str(nl[k][2] * k)
    for k in range(max(actual) + 1, 26):
        actual[k] = "0 (De Geest, cited)"
    tot7 = tot17 = tot26 = 0.0
    rows = []
    for k in range(2, 61):
        lnN = (k - 1) * math.log(10) + math.log(5)
        E = 4 ** k * s3(k) * (3.75 / lnN) ** k * boost(k)
        if k <= 40:
            pass
        rows.append((k, E))
        if k >= 7:
            tot7 += E
        if k >= 17:
            tot17 += E
        if k >= 26:
            tot26 += E
    for k, E in rows[:24]:
        a = actual.get(k, "?")
        print(f"      k={k:>2}: expected non-repunit circular primes {E:10.4g}   actual {a}")
    print(f"    expected total for k >= 7: {tot7:.3f}; for k >= 17: {tot17:.2e}; for k >= 26: {tot26:.2e}")
    e26 = sum(E for k, E in rows if 2 <= k <= 6)
    a26 = sum(int(actual[k]) for k in range(2, 7))
    orb7 = sum(E / k for k, E in rows if k >= 7)
    print(f"    k = 2..6: expected {e26:.1f} numbers, actual {a26}; the model is right to a factor < 3 per length.")
    print(f"    k >= 7: expected {tot7:.2f} numbers in all (about {orb7:.2f} orbits); none exist below 10^{NECK_MAX}"
          " (here) or 10^25 (De Geest).  The expected count decays like (c/k)^k (c = 4 * 3.75 / log 10 = 6.5),")
    print("    so the sum over k converges: finitely many non-repunit circular primes is the natural conjecture")
    print("    (OEIS A293663: 'Conjecture: The sequence is finite'; unproved).  Repunits escape the decay: they are")
    print("    one string per length (the fixed point of rotation), prime with probability of order (log k)/k for prime k")
    print("    (Mersenne-type heuristic), and that sum over primes diverges -- so infinitely many repunit primes are")
    print("    expected (also unproved).")


def main():
    t0 = time.time()
    print("=" * 100)
    print("PART 3. Digit rotation, circular and permutable primes, repunit primes")
    print("=" * 100)
    part3a()
    lens, nl = part3b()
    part3c()
    part3d()
    part3e(lens, nl)
    print(f"\nPart 3 done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
