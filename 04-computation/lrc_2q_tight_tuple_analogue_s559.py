#!/usr/bin/env python3
"""
The k+1 = 2q analogue of the polynomial-method "tight tuple is a universal
corrector" (Prop 4.1 of Sungkawichai-Trakulthongchai, arXiv:2604.23906).
opus-2026-06-02-S559 (remote-control).

PRIME CASE (their Prop 4.1): k+1 odd prime => for every nonzero v in Z_{k+1}^k
there are units s,r in Z_{k+1}^x with  s*v + r*(1,...,k) in {1,...,k-1}^k (mod k+1).
Proven over the FIELD Z_{k+1}. This makes the tight tuple (1,...,k) handle every
bad reduction analytically; it is exactly the step that FAILS at k+1 = 2q
composite (Z_{2q} is not a field), which is why n=14 (=2*7) is the wall.

KEY STRUCTURAL FACT we exploit: every unit mod 2q is ODD, so mod 2 is FORCED
(result_i ≡ v_i + i (mod 2), no unit freedom), while mod q has full Z_q^x freedom.
=> the 2q correction REDUCES (exactly) to a mod-q avoidance problem:

  CLAIM (reduction). With k=2q-1, tight=(1,...,2q-1), target T={1,...,2q-2}
  (avoid 0 and 2q-1), units U=(Z/2q)^x:
    exists s,r in U with s*v+r*tight in T^k (mod 2q)
  IFF  exists s',r' in Z_q^x with  s'*w_i + r'*c_i != f_i (mod q) for all i,
  where w_i=v_i mod q, c_i=i mod q, and f_i = 0 if v_i+i is even, else q-1.

This script (1) verifies the reduction exactly, (2) measures whether the analogue
holds (= mod-q avoidance solvable) over relevant v at q=3,5,7,11 (n=6,10,14,22),
(3) characterises the RESIDUAL failing v -- exactly what a c=2q computation must
still handle -- and (4) tests the weaker 'avoid-0' (base-witness) target.
"""

from math import gcd
import random


def units(m):
    return [x for x in range(1, m) if gcd(x, m) == 1]


def full_check(v, q, target_avoid_top=True):
    """Direct: exists s,r in (Z/2q)^x with s*v_i + r*i in target (mod 2q) all i?
    target = {1,...,2q-2} if avoid_top else {1,...,2q-1} (avoid 0 only)."""
    m = 2 * q
    U = units(m)
    k = m - 1
    tight = list(range(1, m))            # (1,...,2q-1)
    bad = {0, m - 1} if target_avoid_top else {0}
    for s in U:
        for r in U:
            if all(((s * v[i] + r * tight[i]) % m) not in bad for i in range(k)):
                return True
    return False


def modq_check(v, q, target_avoid_top=True):
    """Reduced mod-q avoidance (the claim's RHS)."""
    k = 2 * q - 1
    Uq = units(q)                        # Z_q^x
    w = [v[i] % q for i in range(k)]
    c = [(i + 1) % q for i in range(k)]  # i runs 1..2q-1 => index i has speed i+1
    if target_avoid_top:
        f = [(0 if (v[i] + (i + 1)) % 2 == 0 else q - 1) for i in range(k)]
    else:
        # avoid 0 only: forbidden mod q = 0, but ONLY when parity makes (0,0);
        # i.e. constraint active only for even (v_i+i); odd-parity i: no constraint
        f = [(0 if (v[i] + (i + 1)) % 2 == 0 else None) for i in range(k)]
    for s in Uq:
        for r in Uq:
            ok = True
            for i in range(k):
                if f[i] is None:
                    continue
                if (s * w[i] + r * c[i]) % q == f[i]:
                    ok = False
                    break
            if ok:
                return True
    return False


def sample_v(q, kind, rng):
    """Speeds indexed 0..2q-2 correspond to runners 1..2q-1; v_i = reduction mod 2q."""
    k = 2 * q - 1
    m = 2 * q
    if kind == "random":
        return [rng.randrange(m) for _ in range(k)]
    if kind == "zero-coord":            # at least one v_i = 0 (speed ≡0 mod 2q)
        v = [rng.randrange(m) for _ in range(k)]
        v[rng.randrange(k)] = 0
        return v
    if kind == "parity-tight":          # v_i ≡ i+1 (mod 2): same parity as tight
        return [(rng.randrange(q)) * 2 + ((i + 1) % 2) for i in range(k)]  # forces f_i=0 all
    if kind == "tight":
        return [(i + 1) % m for i in range(k)]
    return [0] * k


def run_q(q, n_each=4000, seed=0):
    m = 2 * q
    rng = random.Random(seed + q)
    print(f"==== q={q}  (n=2q={m}, k={m-1}) ;  units mod {m} = {units(m)} ====")
    # (1) verify reduction full<->modq on mixed sample
    mism = 0
    for _ in range(1500):
        v = sample_v(q, rng.choice(["random", "zero-coord", "parity-tight"]), rng)
        if full_check(v, q, True) != modq_check(v, q, True):
            mism += 1
    print(f"  reduction full<=>mod-q (strict target): mismatches {mism}/1500")
    # (2)+(3) analogue success rate + residual characterisation, by class
    for kind in ["random", "zero-coord", "parity-tight"]:
        ok = 0
        fails = []
        for _ in range(n_each):
            v = sample_v(q, kind, rng)
            if modq_check(v, q, True):
                ok += 1
            elif len(fails) < 3:
                fails.append(v[:])
        print(f"  [{kind:12s}] strict-target analogue holds: {ok}/{n_each}"
              + ("" if ok == n_each else f"   (e.g. fail: {fails[0]})"))
    # (4) weaker avoid-0 (base witness) target
    for kind in ["random", "zero-coord", "parity-tight"]:
        ok = sum(modq_check(sample_v(q, kind, rng), q, False) for _ in range(n_each))
        print(f"  [{kind:12s}] avoid-0 (base witness) holds: {ok}/{n_each}")
    # tight tuple itself
    vt = sample_v(q, "tight", rng)
    print(f"  tight tuple (1..{m-1}): strict={modq_check(vt,q,True)} "
          f"avoid0={modq_check(vt,q,False)}")
    print()


if __name__ == "__main__":
    for q in (3, 5, 7, 11):
        run_q(q)
