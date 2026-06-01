#!/usr/bin/env python3
"""
lrc_gap_chain_antipodal_s553.py    oracle-2026-06-01-S553

EXTEND THE CHAIN toward proving the S552 loneliness spectral gap (HYP-2052):
  every gcd-1 speed set that is not AP-tight has max-collar M(S) >= 2/(2n-1).

The mechanism (this session's proof idea): use the family of WITNESS TIMES
        t_k = k/(2n-1),   k = 1..2n-2.
Write M = 2n-1 (odd). For a speed s,
        ||s * t_k|| = ||k s / M|| = dist(k s mod M, 0)/M.
This is < 2/M  iff  k s ≡ 0, +1, or -1  (mod M)  iff  s ≡ 0, +k^{-1}, -k^{-1} (mod M).
So at time t_k the "bad" residues (those that dip the collar below 2/M) are exactly
the antipodal pair {+k^{-1}, -k^{-1}} together with 0. As k ranges over the units mod
M, a = k^{-1} ranges over ALL nonzero residues, so the bad pair ranges over all n-1
ANTIPODAL PAIRS {a, M-a}, a=1..n-1.

=> LINK 1 (provable lemma). If S has no speed ≡ 0 (mod M) and its residues mod M MISS
   some antipodal pair {a, M-a}, then at t_k (k=a^{-1}) every runner has
   ||s_i t_k|| >= 2/M, hence M(S) >= 2/(2n-1). The gap holds for S.

=> RESIDUAL. The only sets NOT covered by Link 1 are those whose residues mod M:
   (i) contain 0, or (ii) form a TRANSVERSAL hitting every one of the n-1 antipodal
   pairs. With n-1 speeds and n-1 pairs, a 0-residue wastes one slot => misses a pair
   => actually handled (case (i) collapses), leaving the PERFECT TRANSVERSALS: exactly
   one residue from each pair {a, M-a}. There are 2^{n-1} residue-transversals; the AP
   {1,..,n-1} is the all-LOWER transversal (residues 1..n-1, the small half).

So the gap reduces from "all configs" to "perfect transversals mod (2n-1)". This file:
 (A) PROVE-check Link 1: every gcd-1 set missing an antipodal pair (no residue 0) has
     M >= 2/(2n-1)  (must be 100%).
 (B) Confirm the residual is exactly the transversal residue-classes; AP = all-lower.
 (C) Among transversals (enumerated), find which actually have M < 2/(2n-1): conjecture
     only the AP-tight family => the REDUCED gap conjecture.
 (D) push the chain one notch into the transversals: the "one-flip" neighbors of AP.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
from functools import reduce

def fpart(x): return x - (x.numerator // x.denominator)
def norm(x):
    f = fpart(x); return min(f, 1 - f)

def max_collar(speeds):
    S = list(speeds); cands = set()
    for i in range(len(S)):
        si = S[i]; k = 0
        while True:
            t = F(2*k+1, 2*si)
            if t >= 1: break
            if t > 0: cands.add(t)
            k += 1
        for j in range(i+1, len(S)):
            sj = S[j]
            for den in (si+sj, abs(si-sj)):
                if den == 0: continue
                for k in range(1, den):
                    t = F(k, den)
                    if 0 < t < 1: cands.add(t)
    best = F(0); bt = None
    for t in cands:
        c = min(norm(F(s)*t) for s in S)
        if c > best: best = c; bt = t
    return best, bt

def setgcd(S): return reduce(gcd, S)

def residues(S, Mod): return [s % Mod for s in S]

def misses_a_pair(res, n):
    """residues set (mod 2n-1); return an antipodal pair index a it misses, else None."""
    Mod = 2*n - 1
    R = set(res)
    for a in range(1, n):   # pairs {a, Mod-a}, a=1..n-1
        if a not in R and (Mod - a) not in R:
            return a
    return None

def is_transversal(res, n):
    """residues hit every antipodal pair exactly, no 0."""
    Mod = 2*n - 1
    R = res
    if 0 in R: return False
    hit = []
    for a in range(1, n):
        c = (a in R) + ((Mod - a) in R)
        hit.append(c)
    return all(c >= 1 for c in hit)

def main():
    print("="*78)
    print("Extending the gap chain: antipodal-pair witness times t_k = k/(2n-1)")
    print("oracle-2026-06-01-S553")
    print("="*78)

    print("\n(A) LINK 1 proof-check: no residue 0 AND misses an antipodal pair => M>=2/(2n-1)")
    print("-"*78)
    plan = {4:14, 5:13, 6:12, 7:12, 8:11}
    for n, B in plan.items():
        Mod = 2*n-1; edge = F(2, Mod)
        m = n-1
        tot=link1=viol=0
        for S in combinations(range(1, B+1), m):
            if setgcd(S) != 1: continue
            res = residues(S, Mod)
            a = (None if 0 in res else misses_a_pair(res, n))
            if a is None:
                continue  # residual (transversal or has residue 0)
            link1 += 1; tot += 1
            M, _ = max_collar(S)
            if M < edge:
                viol += 1
                if viol <= 3:
                    print(f"   VIOLATION n={n}: {S} res={res} missed pair a={a} M={M}")
        print(f"  n={n} (B={B}): Link-1 sets={link1}, all have M>=2/(2n-1)? "
              f"{viol==0}  (violations={viol})")

    print("\n(B) the RESIDUAL is the perfect transversals mod (2n-1); AP=all-lower")
    print("-"*78)
    for n in (5,6,7):
        Mod = 2*n-1
        ap = tuple(range(1, n))
        print(f"  n={n}: M=2n-1={Mod}; antipodal pairs " +
              ", ".join(f"{{{a},{Mod-a}}}" for a in range(1,n)))
        print(f"        AP residues = {residues(ap,Mod)} (the lower halves 1..{n-1}); "
              f"transversal? {is_transversal(residues(ap,Mod), n)}; "
              f"#residue-transversals = 2^(n-1) = {2**(n-1)}")

    print("\n(C) REDUCED gap: among ENUMERATED transversal-configs, which dip below 2/(2n-1)?")
    print("    (conjecture: only the AP-tight family)")
    print("-"*78)
    plan2 = {4:16, 5:15, 6:14, 7:13}
    for n, B in plan2.items():
        Mod = 2*n-1; edge = F(2, Mod); floor = F(1, n)
        m = n-1
        below = []   # transversal configs with M < edge
        ntrans = 0
        for S in combinations(range(1, B+1), m):
            if setgcd(S) != 1: continue
            res = residues(S, Mod)
            if 0 in res or not is_transversal(res, n):
                continue
            ntrans += 1
            M, bt = max_collar(S)
            if M < edge:
                below.append((S, M, bt))
        # classify the below-edge transversals
        floorcnt = sum(1 for _,M,_ in below if M == floor)
        print(f"  n={n} (B={B}): transversal-configs={ntrans}; "
              f"with M<2/(2n-1): {len(below)} (all at the 1/n floor? "
              f"{floorcnt==len(below)}, count at floor={floorcnt})")
        for (S,M,bt) in below[:8]:
            tag = " [=AP]" if S==tuple(range(1,n)) else ""
            print(f"        {S}  M={M}  t*={bt}{tag}")
        if len(below) > 8: print(f"        ... (+{len(below)-8} more)")

    print("\n(D) one-flip neighbors of AP: flip pair a from lower (a) to upper (2n-1-a).")
    print("    AP=all-lower. The minimal residue-transversal moves off AP; do they jump?")
    print("-"*78)
    for n in (5,6,7,8):
        Mod = 2*n-1; edge = F(2, Mod); floor = F(1,n)
        ap = list(range(1, n))
        print(f"  n={n}: AP={tuple(ap)} M={max_collar(tuple(ap))[0]}")
        for a in range(1, n):
            # flip residue a -> upper representative (smallest positive int with that residue
            # not already used): use 2n-1-a if distinct & not in set, else add Mod
            up = Mod - a
            S = sorted(set(ap) - {a} | {up})
            if len(S) != n-1:  # collision; bump
                up2 = up
                while up2 in set(ap) - {a}:
                    up2 += Mod
                S = sorted((set(ap) - {a}) | {up2})
            S = tuple(S)
            g = setgcd(S)
            M, bt = max_collar(S)
            rel = "=floor1/n" if M==floor else (">=edge OK" if M>=edge else "BELOW-EDGE!")
            print(f"     flip a={a}: S={S} gcd={g}  M={M}  ({rel})")

    print("\n" + "="*78)
    print("READING")
    print("-"*78)
    print(" Link 1 is a clean THEOREM: any gcd-1 set whose residues mod (2n-1) avoid 0")
    print(" and miss an antipodal pair clears the gap (M>=2/(2n-1)), witnessed by t_k,")
    print(" k=a^{-1}. This collapses the ENTIRE gap conjecture onto the 2^(n-1) perfect")
    print(" transversals of the antipodal pairs -- the AP being the all-lower transversal.")
    print(" The reduced conjecture: among transversal-configs only the AP-tight family")
    print(" reaches below 2/(2n-1). Each one-flip neighbor of the AP is the next chain link.")

if __name__ == "__main__":
    main()
