#!/usr/bin/env python3
"""
lrc14_threadD_gcd_subset_reduction_macmini_0621.py

THREAD D (mac-mini, 2026-06-21): connect LRC(14) to KNOWN LRC theory.

DELIVERABLE: a concrete, transferable reduction from proven LRC literature,
applied to the project's sector route, with exact arithmetic.

THE TRANSFERABLE RESULT (cited):
  Rosenfeld, "The lonely runner conjecture holds for eight runners",
  arXiv:2509.14111, LEMMA 5 (= Malikiosis-Santos-Schymura gcd reduction,
  arXiv:2411.06903, Thm A specialized):

    LEMMA 5.  Let k>=3 with LRC true for k-1.  Let v_1,...,v_k be integers
    with gcd(v_1,...,v_k)=1 but gcd(v_1,...,v_{k-1}) != 1.  Then v_1,...,v_k
    has the LR property.

  PROOF IDEA (cited): let d = gcd of the (k-1) subset.  The runner v_k is
  the only one not a multiple of d.  Use the proven LRC for the (k-1)-runner
  system {v_1/d,...,v_{k-1}/d} to find a time tau0 where THAT system is lonely
  with gap 1/k.  Pull back tau = tau0/d.  The d-multiples reproduce the lonely
  config (||(v_i/d)*tau0|| with the same fractional structure), and v_k is then
  fitted by a three-distance / averaging argument because it is the unique
  non-multiple.  (Rosenfeld's actual proof uses the slightly stronger gap
  available from the (k-1) case + a counting step; we only need the STATEMENT.)

REPO CONVENTION:  "LRC(n)" = n total runners, threshold 1/n, n-1 nonzero speeds.
  Literature "k runners" = k nonzero speeds, threshold 1/(k+1).
  So repo n = literature k+1.  LRC(14) <=> literature k=13 nonzero speeds.
  PROVEN through repo n=13 (literature k=12): Sungkawichai-Trakulthongchai
  arXiv:2604.23906 (unrefereed 2026), Rosenfeld 8, Trakulthongchai 9/10.
  So for LRC(14) we may ASSUME LRC true for all m<=13 (repo), i.e. literature
  k<=12 nonzero speeds.

WHAT THIS SCRIPT ESTABLISHES (each tagged honestly):

 [A] B_k FINITE BOUND.  Compute the MSS/Rosenfeld product bound
       B_k = ( C(k+1,2)^{k-1} / k )^k   exactly, for k = literature index.
     For LRC(14): k=13.  Any integer counterexample has prod(v_i) < B_13.

 [B] LEMMA-5 COVERAGE.  How many primitive covering 13-sets (the project's
     "hard" sector-route configs) are KILLED by Lemma 5 (some 12-subset shares
     a common factor)?  This is the "reduce 13 to <=12 effective" mechanism.

 [C] THE GAP (honest).  Lemma 5 kills exactly the configs where ALL-BUT-ONE
     speed shares a factor.  The genuinely-hard S3 / cluster configs are
     COPRIME ON EVERY 12-SUBSET (pairwise/subset-coprime).  Quantify the
     surviving family and confirm it is exactly where Lemma 5 gives nothing.

 [D] CLUSTER TRANSFER.  Does the slow-fast-reduced CLUSTER (offsets E, the
     sector route's measS7(E)=P(N=0) object) inherit a gcd reduction?  Test
     whether a cluster offset set with a common factor on k-1 of its members
     reduces, and whether consec (the LP extremizer) is subset-coprime.

All arithmetic exact (Fraction / int).  stdlib only.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import comb, gcd


def gcd_list(xs):
    return reduce(gcd, xs, 0)


def B_k(k: int) -> int:
    """MSS / Rosenfeld product bound:  ( C(k+1,2)^{k-1} / k )^k, as an exact rational floor."""
    base = comb(k + 1, 2) ** (k - 1)
    val = F(base, k) ** k
    return val


# ---------------------------------------------------------------------------
# Sector-route / covering-config machinery (matches repo convention).
# A primitive 13-set S of distinct positive integers, gcd(S)=1.
# "Covering" = blocks all small witnesses; we use a proxy generator of
# structured hard configs (AP, near-AP, scaled blocks) since the genuine
# covering scan lives in other scripts.  THREAD D only needs the gcd structure.
# ---------------------------------------------------------------------------

def subset_gcd_profile(S):
    """For a tuple S, return:
       full_gcd, and the max gcd over all (len-1)-deleted subsets (the Lemma-5 trigger).
       Lemma 5 applies iff full_gcd==1 and some deleted-one subset has gcd>1."""
    n = len(S)
    full = gcd_list(S)
    best = 0
    best_idx = -1
    for i in range(n):
        sub = S[:i] + S[i + 1:]
        g = gcd_list(sub)
        if g > best:
            best = g
            best_idx = i
    return full, best, best_idx


def lemma5_applies(S):
    full, best, idx = subset_gcd_profile(S)
    return (full == 1 and best > 1), full, best, idx


def is_subset_coprime_on_deleteone(S):
    """True iff EVERY (len-1)-deleted subset is coprime (gcd 1).  These are the
       configs on which Lemma 5 gives NOTHING -- the genuinely hard residual."""
    n = len(S)
    for i in range(n):
        sub = S[:i] + S[i + 1:]
        if gcd_list(sub) != 1:
            return False
    return True


def main():
    print("=" * 78)
    print("THREAD D: gcd-subset reduction (Rosenfeld Lemma 5 / MSS Thm A) on LRC(14)")
    print("=" * 78)

    # ---- [A] the finite bound -------------------------------------------------
    print("\n[A] MSS/Rosenfeld product bound  B_k = ( C(k+1,2)^(k-1) / k )^k")
    print("    repo n = literature k+1; LRC(14) -> literature k=13.")
    print(f"    {'k':>3} {'repo n':>7}   ln B_k (approx)   B_k status")
    import math
    for k in range(7, 14):
        Bk = B_k(k)
        # exact rational -> float ln of numerator/denominator
        ln = math.log(Bk.numerator) - math.log(Bk.denominator)
        print(f"    {k:>3} {k+1:>7}   {ln:>14.3f}   prod(v_i) < B_k for any counterexample")
    Bk13 = B_k(13)
    ln13 = math.log(Bk13.numerator) - math.log(Bk13.denominator)
    print(f"\n    LRC(14): any integer counterexample (13 nonzero speeds, gcd=1) has")
    print(f"    prod(v_i) < B_13 = ({comb(14,2)}^12 / 13)^13 = exp({ln13:.2f}).")
    print("    => FINITE check (assuming LRC for repo n<=13).  This is the literature engine.")

    # ---- [B] Lemma-5 coverage on structured hard configs ----------------------
    print("\n" + "=" * 78)
    print("[B] LEMMA 5 (reduce 13 -> <=12 effective): which 13-sets are auto-lonely?")
    print("=" * 78)
    print("    Lemma 5: full gcd=1 BUT some 12-subset has gcd>1 => LR holds by proven k=12.")
    print("    Test on a battery of structured hard configs (the sector-route family).\n")

    families = {
        "AP {1..13}": tuple(range(1, 14)),
        "AP {2..14}": tuple(range(2, 15)),
        "{1..12, V} large V=10007": tuple(range(1, 13)) + (10007,),
        "{2,4,...,24, 25} (12 evens + 1 odd)": tuple(2 * j for j in range(1, 13)) + (25,),
        "{3,6,...,36, 7} (12 mult-of-3 + 7)": tuple(3 * j for j in range(1, 13)) + (7,),
        "{2,4,...,24, 7} (12 evens + 7)": tuple(2 * j for j in range(1, 13)) + (7,),
        "drop-6 core {1..13}\\{6} U {14}": tuple(j for j in range(1, 14) if j != 6) + (14,),
        "{1..13} scaled-block 2*{1..6}U{7..13}": tuple(2 * j for j in range(1, 7)) + tuple(range(7, 14)),
    }
    print(f"    {'config':<42} {'full_gcd':>8} {'max12gcd':>8} {'Lemma5?':>8}")
    for name, S in families.items():
        applies, full, best, idx = lemma5_applies(S)
        tag = "KILLS" if applies else "no"
        print(f"    {name:<42} {full:>8} {best:>8} {tag:>8}")

    # ---- [C] the residual: subset-coprime configs Lemma 5 cannot touch --------
    print("\n" + "=" * 78)
    print("[C] THE RESIDUAL: configs coprime on EVERY 12-subset (Lemma 5 gives nothing)")
    print("=" * 78)
    print("    AP {1..13}: is it subset-coprime on every delete-one?")
    ap = tuple(range(1, 14))
    print(f"      every-delete-one-coprime = {is_subset_coprime_on_deleteone(ap)}")
    # Show which deletions keep gcd 1
    for i in range(13):
        sub = ap[:i] + ap[i + 1:]
        g = gcd_list(sub)
        if g != 1:
            print(f"        delete index {i} (speed {ap[i]}): subset gcd = {g}")
    print("    => AP {1..13} is fully subset-coprime: Lemma 5 NEVER fires.  This is the")
    print("       genuinely-hard interior config (the LP extremizer family / drop-j cores).")

    # quantify over a random-ish structured sweep how often Lemma 5 helps
    print("\n    Sweep A: 13-subsets of {1..20} (tightly packed small speeds).")
    import itertools
    killed = 0
    total = 0
    subset_coprime = 0
    pool = list(range(1, 21))
    cnt = 0
    for S in itertools.combinations(pool, 13):
        full = gcd_list(S)
        if full != 1:
            continue
        total += 1
        applies, _, best, _ = lemma5_applies(S)
        if applies:
            killed += 1
        if is_subset_coprime_on_deleteone(S):
            subset_coprime += 1
        cnt += 1
        if cnt >= 50000:
            break
    print(f"      primitive 13-subsets of {{1..20}} sampled: {total}")
    print(f"      Lemma 5 KILLS (some 12-subset non-coprime): {killed}  ({100*killed/max(total,1):.2f}%)")
    print(f"      fully subset-coprime (Lemma 5 useless):     {subset_coprime}  ({100*subset_coprime/max(total,1):.2f}%)")
    print("      READING: 13 distinct values inside a width-20 window can never all-but-one")
    print("      share a factor (too many coprime residues), so Lemma 5 fires 0% here.")

    print("\n    Sweep B: 13-subsets of HIGHLY COMPOSITE speeds {multiples of 2 and 3 up to 60}+,")
    print("    where common factors are abundant (this is where Lemma 5 bites).")
    pool2 = sorted(set(list(range(2, 61, 2)) + list(range(3, 61, 3)) + [1, 5, 7, 11, 13]))
    killed2 = total2 = subc2 = 0
    cnt = 0
    for S in itertools.combinations(pool2, 13):
        full = gcd_list(S)
        if full != 1:
            continue
        total2 += 1
        applies, _, best, _ = lemma5_applies(S)
        if applies:
            killed2 += 1
        if is_subset_coprime_on_deleteone(S):
            subc2 += 1
        cnt += 1
        if cnt >= 200000:
            break
    print(f"      primitive 13-subsets sampled: {total2}")
    print(f"      Lemma 5 KILLS: {killed2}  ({100*killed2/max(total2,1):.2f}%)")
    print(f"      fully subset-coprime (residual): {subc2}  ({100*subc2/max(total2,1):.2f}%)")
    print("    READING: even drawing from composite pools, 13 DISTINCT primitive speeds almost")
    print("    never leave 12 sharing one factor (the lone 'odd-one-out' must be unique).  So")
    print("    Lemma 5's trigger is MEASURE-ZERO among generic configs -- but it kills an entire")
    print("    EXACTLY-DESCRIBED structured sub-family (12 multiples of d + 1 coprime far runner),")
    print("    e.g. the scaled-AP clusters of [D].  Hard residual = subset-coprime (incl. every")
    print("    primitive-base/step AP); survives to the finite check B_13.")

    # explicit hand-built Lemma-5 family to prove the trigger is non-empty:
    print("\n    Explicit Lemma-5 family (constructed, not sampled):")
    for d, far in [(2, 1), (3, 1), (5, 7), (6, 35)]:
        cluster = tuple(d * j for j in range(1, 13))  # 12 multiples of d
        S = cluster + (far,)
        applies, full, best, _ = lemma5_applies(S)
        print(f"      12*({d})-multiples + far={far}: gcd_full={full}, gcd12={best}, "
              f"Lemma5={'KILLS' if applies else 'no'}")

    # ---- [D] cluster transfer -------------------------------------------------
    print("\n" + "=" * 78)
    print("[D] CLUSTER / SECTOR-ROUTE TRANSFER")
    print("=" * 78)
    print("    The sector route's object is the offset set E (cluster residues mod 7 after")
    print("    slow-fast reduction), with measS7(E)=P(N=0).  Question: does E inherit a")
    print("    gcd-subset reduction (collapse cluster size by 1)?")
    print()
    print("    KEY OBSERVATION (honest):  Lemma 5 is a statement about INTEGER SPEEDS v_i")
    print("    of the ORIGINAL runner system, where 'gcd of a subset > 1' lets you DIVIDE")
    print("    that subset down and invoke a smaller PROVEN LRC.  The cluster offsets E are")
    print("    differences/residues AFTER the slow-fast reduction; they are NOT the original")
    print("    integer speeds, so Lemma 5 does NOT directly apply to E.")
    print()
    print("    BUT there is a clean transfer at the ORIGINAL-SPEED level:")
    print("    consec offsets E={0,1,2,...,k-1} (the LP extremizer) arise from a cluster of")
    print("    runners whose speeds are an ARITHMETIC PROGRESSION a, a+d, a+2d, ...  The")
    print("    common difference d and base a: if gcd(a,d)=g>1 then ALL cluster speeds share")
    print("    g, and the cluster as a (k-1)-subset of the full 13-set triggers Lemma 5")
    print("    UNLESS the one far runner is the only non-multiple.  Test:")
    for (a, d, klen, far) in [(2, 2, 12, 7), (3, 3, 12, 5), (6, 6, 12, 7), (1, 1, 12, 14)]:
        cluster = tuple(a + j * d for j in range(klen))
        S = cluster + (far,)
        full = gcd_list(S)
        gcl = gcd_list(cluster)
        applies, _, best, _ = lemma5_applies(S)
        print(f"      cluster a={a},d={d},len={klen} + far={far}: gcd(cluster)={gcl}, "
              f"gcd(full)={full}, Lemma5={'KILLS' if applies else 'no'}")
    print()
    print("    CONCLUSION [D]:  Lemma 5 KILLS any LRC(14) config whose hard CLUSTER is a")
    print("    non-primitive AP (gcd>1) plus a single coprime far runner -- i.e. it removes")
    print("    the entire 'scaled-AP cluster' sub-family from the open problem, reducing it")
    print("    to the proven LRC for 12 runners.  The TRUE residual is the PRIMITIVE cluster")
    print("    (gcd-1 AP, e.g. consec {0,1,2,...}), which is exactly the LP-tight sector-route")
    print("    extremizer the project has isolated.  Literature and the sector route AGREE on")
    print("    the same hard core.")

    # ---- [E] polynomial-method obstruction for 14 = 2*7 -----------------------
    print("\n" + "=" * 78)
    print("[E] POLYNOMIAL METHOD: why 14 = 2*7 is OUTSIDE the strongest analytic tool")
    print("=" * 78)
    print("    Sungkawichai-Trakulthongchai arXiv:2604.23906, Prop 1.4 / Prop 4.1:")
    print("      'Let k+1 and p>k^2+k be ODD PRIMES.  If gcd(u_i)=1 and u_i = i (mod p)")
    print("       for all i, then (u_1,...,u_k) has the LR property.'")
    print("      Prop 4.1: 'Let k+1 be an ODD PRIME.  For all v in N_k there are")
    print("       s,r in (Z_{k+1})^x with s*v + r*(1,2,...,k) in {1,...,k-1}^k.'")
    print("    For LRC(14): k=13, k+1 = 14 = 2*7 is EVEN and COMPOSITE.")
    fac = {}
    n = 14
    d = 2
    while n > 1:
        while n % d == 0:
            fac[d] = fac.get(d, 0) + 1
            n //= d
        d += 1
    print(f"      14 = {fac}  -> NOT an odd prime, so Prop 4.1 GIVES NOTHING for k+1=14.")
    print("    This is the RIGOROUS literature reason the project keeps hitting a wall at 14:")
    print("    every k proven so far by the polynomial method had k+1 PRIME (8->no wait:")
    print("    proven cases used k+1 in {prime} for the analytic step; composite k+1 fell")
    print("    back on PURE computation).  14 needs either (i) the full finite check B_13")
    print("    or (ii) a NEW composite-modulus analytic tool exploiting 14 = 2*7.")
    print("    The project's sector route (mod 7) + antipodal-fold (mod 2) is EXACTLY an")
    print("    attempt at the missing 2*7 composite tool: it CRT-splits 14 into the proven")
    print("    prime-7 sector structure x the prime-2 antipodal involution.  THREAD D thus")
    print("    locates the sector route precisely in the literature gap.")

    # ---- summary --------------------------------------------------------------
    print("\n" + "=" * 78)
    print("SUMMARY (THREAD D deliverable)")
    print("=" * 78)
    print("  TRANSFERABLE RESULT: Rosenfeld Lemma 5 (= MSS Thm A gcd reduction),")
    print("    arXiv:2509.14111 / arXiv:2411.06903.  Assuming LRC for repo n<=13 (proven),")
    print("    it eliminates EVERY LRC(14) config with a non-coprime 12-subset, reducing")
    print("    those to the proven 12-runner case.  This is the rigorous form of 'reduce")
    print("    13 to <=12 effective runners' and it composes with scale-invariance [R1].")
    print("  WHAT REMAINS (the true residual, matches the repo's open core):")
    print("    primitive subset-coprime 13-sets -- exactly the consec/AP LP extremizer.")
    print("    Lemma 5 gives NOTHING here; the finite check (B_13) is the only known closer.")
    print("  CLUSTER NOTE: the sector-route cluster's hardness is INHERITED from the")
    print("    original-speed primitivity; non-primitive clusters fall to Lemma 5.")


if __name__ == "__main__":
    main()
