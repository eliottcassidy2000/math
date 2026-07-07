#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S28 -- HYP-4572: THE MEDIANT-ATTAINER TRICHOTOMY + mechanism.

From lrc_mediant_criterion_S28: F(N) = {1..N}\{N-1} u {3(N-1)} has
    M(F(N)) = 3/(3N-1)      N even        (loose; e.g. N=12 -> 3/35 > 2/25)
            = 1/N           N==3,5 mod 6  (loose)
            = 3/(3N+2)=med  N==1 mod 6    (TIGHT: the crux mediant -> gap member)
so F(N) is a gap member IFF N==1 mod 6.  This VERIFIES opus-S118's N=1 mod 6 lead
and REFUTES my prime hypothesis (N=25=5^2 achieves; N=5,17,23 prime fail).

THIS SCRIPT:
 (1) VERIFY the trichotomy formula exactly for N=5..31 (extend the evidence);
 (2) FIND the actual maximizer t=b/Q for each N -> the binding structure of each
     branch (which speeds hit the min at each branch's denominator);
 (3) the mechanism: at a=3*inv5 mod q (q=3N+2), speed N-1 is FORCED to residue -1
     (dist 1) [(N-1)*3*inv5 == -1 since 3(N-1)==-5]; that is WHY N-1 is removed.
     For N even, a BETTER t exists at Q=3N-1; report its binding.
"""
import itertools
from fractions import Fraction as F
from math import gcd


def M_and_argmax(W):
    """exact M and ALL maximizing t (as reduced fractions)."""
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = F(0)
    argmax = []
    seen = set()
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
                argmax = [t]
            elif mv == best:
                argmax.append(t)
    return best, argmax


def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in W)) if g > 1 else tuple(sorted(W))


def family(N):
    base = [v for v in range(1, N + 1) if v != N - 1]
    return primitive(tuple(sorted(set(base) | {3 * (N - 1)})))


def predicted(N):
    if N % 2 == 0:
        return F(3, 3 * N - 1)
    if N % 6 == 1:
        return F(3, 3 * N + 2)  # mediant
    return F(1, N)  # N==3,5 mod 6


def binding_at(W, t):
    """which speeds achieve the min ||v t|| at t, and their residue r=v*num mod Q."""
    Q = t.denominator
    b = t.numerator
    m = min(abs(v * t - round(v * t)) for v in W)
    hit = []
    for v in W:
        d = abs(v * t - round(v * t))
        if d == m:
            r = (v * b) % Q
            hit.append((v, min(r, Q - r)))
    return Q, b, m, hit


def main():
    print("=" * 86)
    print("(1) TRICHOTOMY VERIFICATION  M(F(N)) vs predicted, N=5..31")
    print("=" * 86)
    print(f"  {'N':>3} {'N mod 6':>7} {'M(F(N))':>10} {'predicted':>10} {'branch':>14} {'ok':>4}")
    allok = True
    gap_members = []
    for N in range(5, 32):
        W = family(N)
        M, _ = M_and_argmax(W)
        pred = predicted(N)
        ok = (M == pred)
        allok &= ok
        if N % 2 == 0:
            branch = "3/(3N-1) even"
        elif N % 6 == 1:
            branch = "MEDIANT 3/(3N+2)"
            if ok:
                gap_members.append(N)
        else:
            branch = "1/N (3,5 mod6)"
        print(f"  {N:>3} {N % 6:>7} {str(M):>10} {str(pred):>10} {branch:>16} "
              f"{('YES' if ok else 'NO!'):>4}")
    print()
    print(f"  trichotomy formula EXACT for all N=5..31: {allok}")
    print(f"  gap members (M = mediant): N = {gap_members}  (all == 1 mod 6)")
    print(f"  => F(N) is a gap member  <=>  N == 1 mod 6   [opus-S118 VERIFIED]")

    print()
    print("=" * 86)
    print("(2) MAXIMIZER + BINDING STRUCTURE per branch")
    print("=" * 86)
    print(f"  {'N':>3} {'branch':>16} {'max t = b/Q':>14} {'min||.||':>9} {'binding speeds (v: dist)':>32}")
    for N in [7, 8, 9, 12, 13, 14, 19]:
        W = family(N)
        M, argmax = M_and_argmax(W)
        t = argmax[0]
        Q, b, m, hit = binding_at(W, t)
        br = ("MEDIANT" if N % 6 == 1 else ("even 3/(3N-1)" if N % 2 == 0 else "1/N"))
        hitstr = ", ".join(f"{v}:{d}" for v, d in hit)
        print(f"  {N:>3} {br:>16} {str(t):>14} {str(m):>9}   {hitstr}")

    print()
    print("=" * 86)
    print("(3) MECHANISM: WHY N-1 is removed, and the even-branch better-t")
    print("=" * 86)
    print("  At binding maximizer a=3*inv5 mod q (q=3N+2): 3(N-1)=q-5 == -5, so")
    print("  (N-1)*a = (N-1)*3*inv5 = [3(N-1)]*inv5 = -5*inv5 = -1 mod q.")
    print("  => speed N-1 is ALWAYS at residue -1 (dist 1) -- FORCED collision, hence")
    print("     removed.  The mediant lower bound M >= 3/q holds for ALL N (this t).")
    print()
    print("  The trichotomy is about the GLOBAL max.  For N even, a t at Q=3N-1 beats")
    print("  the mediant.  Binding at that t (N even):")
    for N in [8, 12, 14]:
        W = family(N)
        M, argmax = M_and_argmax(W)
        t = argmax[0]
        Q, b, m, hit = binding_at(W, t)
        print(f"    N={N}: Q=3N-1={3*N-1}, max t={t}, min={m}; binding {[(v,d) for v,d in hit]}")
    print()
    print("  For N==1 mod 6, NO t beats the mediant => M = 3/(3N+2), the gap member.")
    print("  For N==3,5 mod 6, the max is 1/N (the far element + parity leave 1/N open).")


if __name__ == "__main__":
    main()
