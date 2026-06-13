#!/usr/bin/env python3
"""
twin_goldbach_exceptions_s521.py    oracle-2026-06-01-S521

Twin-prime Goldbach: which even numbers are NOT a sum of two twin primes?

A "twin prime" = a prime p that is a member of some twin pair (p-2 or p+2 prime).
T = {3,5,7,11,13,17,19,29,31,41,43,...}.

We compute E = { even e : e != a+b for any a,b in T }.  This is OEIS A007534
(conjectured finite).  We test repetition-allowed (a<=b, a=b ok), confirm the
count and the largest term, list them, and dump structure: gaps, consecutive
runs (triples), residues mod small moduli, and #representations near the tail.
"""
from sympy import primerange, isprime

LIMIT = 2_000_000   # search bound; far past any conjectured exception

def main():
    print(f"twin-prime Goldbach exceptions up to {LIMIT} (oracle-S521)\n")
    primes = list(primerange(3, LIMIT + 3))
    pset = set(primes)
    # twin primes: prime p with p-2 or p+2 also prime
    twins = [p for p in primes if (p - 2) in pset or (p + 2) in pset]
    tset = set(twins)
    print(f"primes < {LIMIT}: {len(primes)};  twin primes: {len(twins)}")
    print(f"first twin primes: {twins[:15]}\n")

    # representability of each even e as a+b, a,b twin primes (repetition allowed)
    rep = bytearray(LIMIT + 1)      # rep[e] = 1 if e is a sum of two twin primes
    cnt = [0] * (LIMIT + 1)         # number of unordered representations (small tail only)
    tw = twins
    for i, a in enumerate(tw):
        if a + a > LIMIT:
            break
        for b in tw[i:]:
            s = a + b
            if s > LIMIT:
                break
            rep[s] = 1

    exceptions = [e for e in range(2, LIMIT + 1, 2) if not rep[e]]
    print(f"EVEN numbers (>=2) that are NOT a sum of two twin primes: {len(exceptions)}")
    print(f"largest exception: {max(exceptions)}\n")
    print("THE LIST:")
    print(exceptions)
    print()

    # ---- structure ----
    print("STRUCTURE")
    # consecutive runs (each step +2)
    runs = []
    cur = [exceptions[0]]
    for x in exceptions[1:]:
        if x == cur[-1] + 2:
            cur.append(x)
        else:
            runs.append(cur); cur = [x]
    runs.append(cur)
    from collections import Counter
    runlen = Counter(len(r) for r in runs)
    print(f"  runs of consecutive even exceptions (step +2): {len(runs)} runs, length histogram {dict(sorted(runlen.items()))}")
    print(f"  the runs: {[ (r[0], r[-1], len(r)) for r in runs ]}")
    # residues
    for mod in (3, 6, 10, 30):
        c = Counter(e % mod for e in exceptions)
        print(f"  residues mod {mod}: {dict(sorted(c.items()))}")
    # gaps between successive exceptions
    gaps = [exceptions[i+1]-exceptions[i] for i in range(len(exceptions)-1)]
    print(f"  distinct gaps between successive exceptions: {sorted(set(gaps))}")

if __name__ == "__main__":
    main()
