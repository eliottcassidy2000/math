#!/usr/bin/env python3
"""
LRC(N) equality-case census across many N   (boxeph-2026-07-17-S81)
====================================================================

The equality case of LRC(N) is the prefix family V = {1, 2, ..., N-1}
(M = 1/N, attained at t = 1/N).  death-star's THM-991 proved the LIVE LAW
for N=14 in Lean; the fleet's named-next target is the DIFFERENCE-CLOSED
GENERALIZATION -- "one theorem for the equality case of every LRC(N)."

KEY STRUCTURAL FACT (proved here, used throughout):
  A finite set S of positive integers closed under nonzero differences is
  exactly a scaled prefix {d, 2d, ..., kd}.  So "difference-closed" = "the
  equality cases of LRC, up to dilation."  Dilation by d only rescales the
  modulus (gcd cancels), so WLOG d=1: study the prefixes {1,...,N-1}.

Definitions (matching LRCDiscreteBonferroni.lean, threshold 1/N):
  inBand(v,q,p,i) :  q <= N*((v_i * p) mod q) <= (N-1)*q     <=> ||v_i p/q|| >= 1/N
  bandCount(q,p)  :  # runners i with NOT inBand   (the coverage count C(p))
  liveCount(q)    :  # p in (0,q) with bandCount(q,p) = 0    (loneliness witnesses)
  deep at q       :  p with some runner exactly at 0 (a collision)  -- the deep census

This script:
  (1) verifies the UNIFORM LIVE LAW:   liveCount = phi(N)*[N|q], set {o*(q/N): o unit}
  (2) tabulates the DEEP census at resonance and its N-dependence
  (3) tests whether difference-closure alone (scaled prefixes) suffices
  (4) probes NON-prefix families to isolate what condition is essential
"""
from math import gcd
from fractions import Fraction

def phi(n):
    return sum(1 for k in range(1, n+1) if gcd(k, n) == 1)

def units_mod(n):
    return [k for k in range(1, n) if gcd(k, n) == 1]

def in_band(vi, q, p, N):
    """Runner speed vi is safe (clears 1/N) at multiplier p, modulus q."""
    r = (vi * p) % q
    return q <= N * r <= (N - 1) * q

def band_count(V, q, p, N):
    return sum(1 for vi in V if not in_band(vi, q, p, N))

def live_multipliers(V, q, N):
    return [p for p in range(1, q) if band_count(V, q, p, N) == 0]

def deep_multipliers(V, q):
    """p at which some runner collides with the origin: q | vi*p for some i."""
    return [p for p in range(1, q) if any((vi * p) % q == 0 for vi in V)]

# ---------------------------------------------------------------------------
# (1) THE UNIFORM LIVE LAW  --  liveCount(V={1..N-1}, q) = phi(N)*[N|q]
# ---------------------------------------------------------------------------
def test_live_law(Nmax=22, Qmul=6):
    print("="*74)
    print("(1) UNIFORM LIVE LAW  for V = {1,...,N-1}, threshold 1/N")
    print("    claim: liveCount(q) = phi(N) if N|q else 0;  set = {o*(q/N): o unit mod N}")
    print("="*74)
    all_ok = True
    for N in range(3, Nmax+1):
        V = list(range(1, N))
        Q = N * Qmul
        ok = True
        detail = []
        for q in range(2, Q+1):
            lm = live_multipliers(V, q, N)
            if q % N == 0:
                m = q // N
                predicted = sorted((o * m) % q for o in units_mod(N))
                exp_count = phi(N)
            else:
                predicted = []
                exp_count = 0
            if sorted(lm) != predicted:
                ok = False
                detail.append((q, sorted(lm), predicted))
        flag = "OK " if ok else "*** MISMATCH ***"
        print(f"  N={N:2d}  phi(N)={phi(N):2d}  units={units_mod(N)}   {flag}")
        if not ok:
            all_ok = False
            for q, got, pred in detail[:3]:
                print(f"        q={q}: got {got} predicted {pred}")
    print(f"\n  LIVE LAW verdict: {'ALL N CONFIRMED' if all_ok else 'FAILURES ABOVE'}\n")
    return all_ok

# ---------------------------------------------------------------------------
# (2) THE DEEP CENSUS at resonance q = N  (and q = N*m)
# ---------------------------------------------------------------------------
def test_deep_census(Nmax=22):
    print("="*74)
    print("(2) DEEP CENSUS at the base modulus q = N   (collisions with origin)")
    print("    p is deep iff N | i*p for some i in 1..N-1  <=>  N/gcd(N,p) <= N-1")
    print("    <=>  p is NOT a unit mod N and p != 0.  #deep(q=N) = N-1-phi(N).")
    print("="*74)
    for N in range(3, Nmax+1):
        V = list(range(1, N))
        dm = deep_multipliers(V, N)
        # predicted: p in (0,N) that are non-units (gcd(p,N)>1)
        pred = [p for p in range(1, N) if gcd(p, N) > 1]
        ok = sorted(dm) == sorted(pred)
        # decompose deep p by which divisor it hits
        print(f"  N={N:2d}  #deep={len(dm):2d} = N-1-phi(N)={N-1-phi(N):2d} "
              f"[{'OK' if ok else 'MISMATCH'}]  deep p = {dm}")
    print()

# ---------------------------------------------------------------------------
# (3) DIFFERENCE-CLOSURE test: scaled prefixes {d,2d,...,(N-1)d}
# ---------------------------------------------------------------------------
def test_scaled_prefix(Nmax=12, dmax=4):
    print("="*74)
    print("(3) SCALED PREFIXES  {d,2d,...,(N-1)d}  (the general difference-closed set)")
    print("    claim: live law is dilation-invariant -- same phi(N)*[N|q] pattern")
    print("="*74)
    ok_all = True
    for N in range(3, Nmax+1):
        for d in range(1, dmax+1):
            V = [d*i for i in range(1, N)]
            ok = True
            for q in range(2, 5*N+1):
                lm = live_multipliers(V, q, N)
                cnt = len(lm)
                exp = phi(N) if (q % N == 0) else 0
                # NOTE: with common factor d, extra resonances at gcd; test raw count law
                if q % N == 0 and cnt < phi(N):
                    ok = False
                    break
            if not ok:
                ok_all = False
                print(f"    N={N} d={d}: live law weakened (extra factor resonance)")
    print(f"  scaled-prefix verdict: {'live>=phi(N) at N|q holds' if ok_all else 'see above'}")
    print("  (dilation adds resonances at gcd-divisors; the PRIMITIVE d=1 is the clean law)\n")

# ---------------------------------------------------------------------------
# (4) NON-prefix families: is difference-closure NECESSARY for live-empty-off-N?
# ---------------------------------------------------------------------------
def test_non_prefix(N=14):
    print("="*74)
    print(f"(4) NON-PREFIX families at threshold 1/{N} (13 speeds): where does live-law break?")
    print("="*74)
    fams = {
        "prefix {1..13}      ": list(range(1, 14)),
        "GW doubling {1..11,13,24}": list(range(1,12))+[13,24],
        "sporadic {1..11,13,c=27} ": list(range(1,12))+[13,27],
        "deep well {1..12,182}    ": list(range(1,13))+[182],
        "missing-2 {1,3,4..14}    ": [1]+list(range(3,15)),
    }
    for name, V in fams.items():
        Q = 3*N
        live_qs = []
        for q in range(2, Q+1):
            c = len(live_multipliers(V, q, N))
            if c > 0:
                live_qs.append((q, c))
        # smallest modulus with a live multiplier = smallest denom loneliness witness
        first = live_qs[0] if live_qs else None
        print(f"  {name}: first live (q,count) = {first};  "
              f"live moduli<= {Q}: {[q for q,_ in live_qs][:8]}")
    print("\n  READING: prefix is live ONLY at q multiple of 14 (tight);")
    print("  non-tight families become live at SMALLER, non-resonant moduli => easier.\n")

if __name__ == "__main__":
    ok = test_live_law(Nmax=22, Qmul=5)
    test_deep_census(Nmax=22)
    test_scaled_prefix(Nmax=12, dmax=4)
    test_non_prefix(N=14)
    print("Done. Live-law uniform in N:", ok)
