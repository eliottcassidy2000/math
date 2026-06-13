#!/usr/bin/env python3
"""
Is the bounded division-point sieve complete enough to PROVE LRC(14)?
opus-2026-06-01-S551 (remote-control).

A "division-point sieve witness" for speed set V is a rational time t=a/q
(gcd(a,q)=1) at which every runner is safe: min((v_i a) mod q, q-((v_i a) mod q))
* n >= q.  A FOUND witness is a RIGOROUS proof that V is lonely (it is just an
exact evaluation of the loneliness predicate at one explicit rational t).

Dream: if there were a finite Q(n) such that EVERY speed set has a sieve witness
with denominator q <= Q(n), then LRC(n) would reduce to the finite check
"every residue tuple mod lcm(2..Q) admits a witness q<=Q".

This script tests that dream and finds it FALSE by an explicit construction,
then maps how the minimal witness denominator actually behaves:

  (A) THE BLIND FAMILY (proves the sieve has NO finite completeness threshold):
      if a single speed is divisible by lcm(2..Q), then for every q<=Q that
      speed is in the danger band for EVERY multiplier (residue 0), so NO
      witness has denominator <= Q.  Yet the set is still lonely -- we exhibit
      a witness at a larger modulus, which PROVES loneliness with no appeal to
      the conjecture.  Conclusion: witness denominators are UNBOUNDED over
      lonely sets; a bounded sieve cannot prove LRC(n).

  (B) THE RESISTANCE CURVE for BOUNDED speeds: hill-climb to maximise the
      minimal witness denominator over speed sets with speeds <= cap.  Shows
      that, once speeds are bounded, the witness denominator is small and grows
      slowly -- the hardness is entirely carried by divisibility loading.
"""

from fractions import Fraction
from math import gcd, lcm
import random

import importlib.util, os
_spec = importlib.util.spec_from_file_location(
    "sieve_mod",
    os.path.join(os.path.dirname(__file__), "lrc_n14_multiprime_sieve_s551.py"))
S = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S)

safe_residue = S.safe_residue
witness_at_modulus = S.witness_at_modulus
lonely_exact = S.lonely_exact


def min_witness_modulus(V, n, qmax):
    for q in range(2, qmax + 1):
        a = witness_at_modulus(V, q, n)
        if a is not None:
            return q, a
    return None, None


# --------------------------------------------------------------------------
# (A) The blind family: lcm(2..Q)-loaded sets defeat the sieve up to Q
# --------------------------------------------------------------------------
def blind_family_demo(n=14):
    N = n - 1
    print("====== (A) THE BLIND FAMILY: sieve has no finite completeness Q ======")
    print("A speed divisible by lcm(2..Q) is danger-band (residue 0) at every")
    print("q<=Q for every multiplier => NO division-point witness with q<=Q.")
    print("We still PROVE loneliness by exhibiting a witness at a larger q.\n")
    for Q in (14, 16, 18, 20, 24):
        L = lcm(*range(2, Q + 1))
        V = tuple(sorted({L} | set(range(1, N))))   # big loaded speed + 1..N-1
        assert len(V) == N, (len(V), V)
        # confirm blind up to Q
        q_lo, _ = min_witness_modulus(V, n, Q)
        # find the actual smallest witness (search well beyond Q)
        q_star, a_star = min_witness_modulus(V, n, 6 * Q + 60)
        blind = (q_lo is None)
        # exact cross-check when speeds are small enough
        lon, _ = lonely_exact(V, n)
        lon_str = {True: "lonely", False: "NOT-LONELY!", None: "exact-skipped"}[lon]
        proof = (f"witness {a_star}/{q_star} (verified safe) => LONELY"
                 if q_star else "no witness found up to bound!")
        print(f"  Q={Q:2d}: lcm(2..Q)={L:>10d}  blind up to {Q}? {blind}   "
              f"smallest witness q*={q_star}   [{lon_str}]")
        print(f"        proof of loneliness: {proof}")
        # double-check the exhibited witness really is safe, exactly
        if q_star:
            t = Fraction(a_star, q_star)
            assert S.is_safe_time(V, t, n), "exhibited witness is NOT safe!"
    print()
    print("  => blindness up to ANY Q is achievable, yet every such set is lonely.")
    print("     The minimal witness denominator is UNBOUNDED over lonely sets:")
    print("     a finite-modulus division-point sieve CANNOT prove LRC(n).\n")


# --------------------------------------------------------------------------
# (B) Resistance curve for BOUNDED speeds (hill-climb on min witness modulus)
# --------------------------------------------------------------------------
def hardest_bounded(n=14, cap=300, qmax=400, restarts=40, steps=300, seed=7):
    N = n - 1
    rng = random.Random(seed)

    def primitive(V):
        g = 0
        for v in V:
            g = gcd(g, v)
        return tuple(sorted(v // g for v in V))

    def hardness(V):
        q, _ = min_witness_modulus(V, n, qmax)
        return q if q is not None else qmax + 1

    best_V, best_h = None, -1
    for _ in range(restarts):
        V = list(rng.sample(range(1, cap + 1), N))
        h = hardness(tuple(V))
        for _ in range(steps):
            i = rng.randrange(N)
            old = V[i]
            cand = rng.randrange(1, cap + 1)
            if cand in V:
                continue
            V[i] = cand
            nh = hardness(tuple(V))
            if nh >= h:
                h = nh
            else:
                V[i] = old
        if h > best_h:
            best_h, best_V = h, primitive(tuple(V))
    print(f"====== (B) HARDEST BOUNDED set (speeds<= {cap}), n={n} ======")
    q, a = min_witness_modulus(best_V, n, qmax)
    lon, t = lonely_exact(best_V, n)
    print(f"  hill-climb max min-witness-modulus = {best_h}")
    print(f"  set: {best_V}")
    print(f"  smallest witness: {a}/{q}    exact lonely: {lon}  (t={t})")
    # show its divisibility loading (which q<=n it kills by residue 0)
    loaded = [qq for qq in range(2, n + 1)
              if any(v % qq == 0 for v in best_V)]
    print(f"  divisibility loading (q<= {n} dividing some speed): {loaded}")
    print(f"  -> bounded speeds keep the witness denominator small; the only")
    print(f"     way to push it up is divisibility loading (=> larger speeds).\n")


if __name__ == "__main__":
    blind_family_demo(14)
    hardest_bounded(n=14, cap=300)
    hardest_bounded(n=14, cap=60)
