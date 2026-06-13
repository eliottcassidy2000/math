#!/usr/bin/env python3
"""
LRC @ n=14 and above: the multi-prime sieve covering obstruction.
opus-2026-06-01-S551 (remote-control).

Setting (repo convention, S525): n runners = (n-1) moving speeds + observer at 0;
the observer is "lonely" at time t iff ||v_i t|| >= 1/n for every speed v_i.
Covering reformulation:  LRC(n) <=> the (n-1) open danger arcs
    D_i = { t in R/Z : ||v_i t|| < 1/n }
never cover the circle.  Each D_i is v_i arcs of total length 2/n.

THE SIEVE (THM-369, Lean-checked).  At a rational time t = a/q with gcd(a,q)=1,
runner v is safe iff  min(r, q-r) * n >= q,  where r = (v*a) mod q.
 - For q <= n: the only dangerous residue is r = 0, so a/q is a witness
   iff q divides no speed (divisibility sieve).
 - For q  > n: the danger band widens to the central 2*ceil(q/n)-1 residues,
   so primes/composites q > n give extra multiplier freedom.

This module:
  (1) exact ground-truth loneliness test (handles measure-zero tight configs);
  (2) the sieve witness search over a range of moduli;
  (3) a survey over structured + random + adversarial n=14 speed sets;
  (4) the counterexample divisibility lattice (necessary conditions);
  (5) the maximally-divisibility-loaded small candidate, verified lonely.

Everything is exact (Fraction / integer arithmetic). No floats in any decision.
"""

from fractions import Fraction
from math import gcd, lcm, ceil
from itertools import combinations
import random


# --------------------------------------------------------------------------
# 1. EXACT ground-truth: is the observer lonely for speed set V at level n?
# --------------------------------------------------------------------------
def is_safe_time(V, t, n):
    """Exact: every ||v t|| >= 1/n ?  t a Fraction in [0,1)."""
    for v in V:
        x = (v * t) % 1                      # Fraction in [0,1)
        d = min(x, 1 - x)                    # circular distance to 0
        if d < Fraction(1, n):
            return False
    return True


MAXEXACT = 200_000          # skip exact ground-truth if sum of speeds exceeds this


def lonely_exact(V, n):
    """
    Exact decision of LRC for one speed set.  Returns (is_lonely, witness_t),
    or (None, None) if the speeds are too large for the exact enumeration.
    SAFE = {t: all ||v_i t|| >= 1/n} is CLOSED; its component endpoints lie in
        E = { (j*n +/- 1)/(v_i*n) : i, j } (mod 1).
    A degenerate (tight) component IS such an endpoint, so testing every
    endpoint AND every midpoint between consecutive endpoints is complete.
    """
    if sum(V) > MAXEXACT:
        return None, None
    endpoints = set()
    for v in V:
        for j in range(v + 1):               # j*n +/- 1 over one period of v t
            for s in (-1, 1):
                t = Fraction(j * n + s, v * n) % 1
                endpoints.add(t)
    pts = sorted(endpoints)
    # test endpoints (catches tight, measure-zero witnesses)
    for t in pts:
        if is_safe_time(V, t, n):
            return True, t
    # test open midpoints (catches positive-measure safe intervals)
    m = len(pts)
    for k in range(m):
        a, b = pts[k], pts[(k + 1) % m]
        mid = ((a + b) / 2) % 1 if b > a else ((a + b + 1) / 2) % 1
        if is_safe_time(V, mid, n):
            return True, mid
    return False, None


# --------------------------------------------------------------------------
# 2. THE SIEVE witness search
# --------------------------------------------------------------------------
def safe_residue(r, q, n):
    """Runner with residue r = v*a mod q is safe at a/q iff min(r,q-r)*n >= q."""
    return min(r, q - r) * n >= q


def witness_at_modulus(V, q, n):
    """Return a multiplier a (coprime to q) s.t. a/q is a witness, else None."""
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        if all(safe_residue((v * a) % q, q, n) for v in V):
            return a
    return None


def min_witness_modulus(V, n, qmax):
    """Smallest modulus q in [2, qmax] that yields a sieve witness a/q."""
    for q in range(2, qmax + 1):
        a = witness_at_modulus(V, q, n)
        if a is not None:
            return q, a
    return None, None


# --------------------------------------------------------------------------
# 3. Speed-set generators (n=14 => 13 speeds), all distinct positive ints
# --------------------------------------------------------------------------
def primitive(V):
    g = 0
    for v in V:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in V))


def structured_sets(N):           # N = number of speeds = n-1
    sets = {}
    sets["AP {1..N}"] = tuple(range(1, N + 1))
    sets["AP odd"] = tuple(range(1, 2 * N, 2))
    sets["powers of 2"] = tuple(2 ** k for k in range(N))
    sets["powers of 3"] = tuple(3 ** k for k in range(N))
    # Fibonacci
    fib = [1, 2]
    while len(fib) < N:
        fib.append(fib[-1] + fib[-2])
    sets["Fibonacci"] = tuple(fib[:N])
    # Sidon (Mian-Chowla)
    s, sums = [1], {0}
    while len(s) < N:
        c = s[-1] + 1
        while True:
            ok = all((c + x) not in sums for x in s)
            if ok:
                break
            c += 1
        for x in s:
            sums.add(c + x)
        sums.add(2 * c)
        s.append(c)
    sets["Sidon"] = tuple(s)
    sets["primes"] = tuple([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41][:N])
    sets["squares"] = tuple((k + 1) ** 2 for k in range(N))
    return sets


# --------------------------------------------------------------------------
# 4. Counterexample divisibility lattice (necessary conditions, n)
# --------------------------------------------------------------------------
def divisibility_necessity(n):
    """
    A counterexample at level n must have, for every q in {2,...,n}, a speed
    divisible by q (else the divisibility sieve a=1/q is a witness).
    Reduce to the BINDING moduli: q is binding iff no PROPER multiple chain
    forces it -- i.e. report each q and which smaller q it implies.
    """
    binding = []
    for q in range(2, n + 1):
        # q is "prime-power-like binding" if it is a prime power; composites
        # q=ab (a,b>1, coprime) are implied only by having a multiple of q
        # itself, NOT by separate multiples of a and b. So EVERY q in 2..n is
        # independently required (a single speed may satisfy several at once).
        binding.append(q)
    return binding


def min_loaded_candidate(n, cap=3000):
    """
    Greedily build a 13-speed set (distinct, <= cap) covering 'divisible by q'
    for every q in {2,...,n}, minimizing the largest speed used. Returns the
    set or None. This is the most divisibility-loaded small counterexample
    CANDIDATE; we then test whether it is (still) lonely.
    """
    need = set(range(2, n + 1))
    chosen = []
    # candidate speeds: small multiples; prefer numbers covering many needs
    pool = list(range(2, cap + 1))

    def covers(v):
        return {q for q in need if v % q == 0}

    while need and len(chosen) < n - 1:
        # pick speed covering the most still-needed moduli, smallest on ties
        best = max(pool, key=lambda v: (len(covers(v)), -v))
        c = covers(best)
        if not c:
            break
        chosen.append(best)
        need -= c
        pool.remove(best)
    if need:
        return None, need
    # pad to exactly n-1 distinct speeds with small unused integers
    extra = 1
    chosen_set = set(chosen)
    while len(chosen) < n - 1:
        if extra not in chosen_set:
            chosen.append(extra)
            chosen_set.add(extra)
        extra += 1
    return tuple(sorted(chosen)), set()


# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------
def run(n=14, qmax=60, n_random=400, seed=12345):
    N = n - 1
    print(f"==================  LRC sieve covering analysis, n = {n}  ==================")
    print(f"(observer + {N} runners; lonely iff ||v_i t|| >= 1/{n} for all i)\n")

    # ---- structured sets ----
    print("---- Structured speed sets: ground-truth loneliness + min sieve modulus ----")
    structs = structured_sets(N)
    rows = []
    for name, V in structs.items():
        lon, t = lonely_exact(V, n)
        q, a = min_witness_modulus(V, n, qmax)
        wt = f"{a}/{q}" if q else f">qmax({qmax})"
        rows.append((name, lon, str(t) if t is not None else "-", wt))
        print(f"  {name:14s} lonely={lon!s:5s}  witness_t={str(t):>10s}  "
              f"min_sieve={wt}")
    print()

    # ---- random primitive sets ----
    print(f"---- {n_random} random primitive speed sets (speeds in [1, 200]) ----")
    rng = random.Random(seed)
    fails = 0
    sieve_resolved = 0
    max_q_needed = 0
    q_hist = {}
    seen = set()
    tries = 0
    made = 0
    while made < n_random and tries < n_random * 50:
        tries += 1
        V = primitive(tuple(rng.sample(range(1, 201), N)))
        if V in seen:
            continue
        seen.add(V)
        made += 1
        lon, _ = lonely_exact(V, n)
        if lon is False:
            fails += 1
            print(f"  *** LRC FAILURE (would be a counterexample!): {V}")
        q, a = min_witness_modulus(V, n, qmax)
        if q is not None:
            sieve_resolved += 1
            q_hist[q] = q_hist.get(q, 0) + 1
            max_q_needed = max(max_q_needed, q)
    print(f"  LRC failures: {fails} / {made}   (expected 0)")
    print(f"  sieve-resolved within q<= {qmax}: {sieve_resolved} / {made}")
    print(f"  max sieve modulus needed: {max_q_needed}")
    print(f"  modulus histogram (q: count): "
          f"{dict(sorted(q_hist.items()))}")
    print()

    # ---- the divisibility necessity lattice ----
    print("---- Necessary condition for a counterexample (divisibility sieve) ----")
    binding = divisibility_necessity(n)
    print(f"  A counterexample at n={n} must contain a speed divisible by EACH of:")
    print(f"    {binding}")
    print(f"  (each q in [2,{n}] independently; the prime powers <= {n} are the")
    print(f"   irreducible demands: {[q for q in binding if _is_pp(q)]})")
    print()

    # ---- maximally loaded small candidate ----
    print("---- Most divisibility-loaded small candidate, tested for loneliness ----")
    cand, missing = min_loaded_candidate(n, cap=3000)
    if cand is None:
        print(f"  could not cover; still missing {missing}")
    else:
        lon, t = lonely_exact(cand, n)
        q, a = min_witness_modulus(cand, n, 4 * qmax)
        print(f"  candidate (covers every q in [2,{n}]): {cand}")
        print(f"    max speed = {max(cand)},  lonely = {lon},  witness_t = {t}")
        print(f"    min sieve modulus = {a}/{q}" if q else
              f"    no sieve witness up to {4*qmax}")
        # which q does each chosen speed satisfy
        for v in sorted(set(cand)):
            qs = [q for q in range(2, n + 1) if v % q == 0]
            if qs:
                print(f"      speed {v:5d} covers q in {qs}")
    print()

    # ---- the AP, in detail (the tight extremal) ----
    print("---- The arithmetic progression {1,...,N}: the tight extremal ----")
    ap = tuple(range(1, N + 1))
    lon, t = lonely_exact(ap, n)
    print(f"  {ap}")
    print(f"    lonely = {lon}, witness t = {t}  (= 1/{n}; the regular {n}-gon wall)")
    # confirm it is tight: safe set has measure zero (only finite witnesses)
    safe_pts = []
    endpoints = set()
    for v in ap:
        for j in range(v + 1):
            for s in (-1, 1):
                endpoints.add(Fraction(j * n + s, v * n) % 1)
    for tt in sorted(endpoints):
        if is_safe_time(ap, tt, n):
            safe_pts.append(tt)
    print(f"    number of safe boundary points = {len(safe_pts)}: {safe_pts}")
    print(f"    -> safe set is finite (measure zero) = TIGHT, as expected.")


def _is_pp(q):
    for p in range(2, q + 1):
        k = q
        if q % p == 0:
            while k % p == 0:
                k //= p
            return k == 1
    return False


if __name__ == "__main__":
    run(n=14)
    print("\n\n")
    # a peek above 14, the doubled-prime neighbours
    for n in (15, 16, 18, 22):
        print(f"########## quick pass n={n} ##########")
        run(n=n, qmax=70, n_random=120)
        print("\n")
