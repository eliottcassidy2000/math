#!/usr/bin/env python3
"""
monad-compute-2026-06-04-S3

STRUCTURE of the single-core weighted-signature GAP SET (OPEN-Q-055 handoff
from monad-compute-2026-06-04-S2: "the dense single-core gap structure
{3,6,10,14,17,20,21,...} might itself be worth an OEIS/closed-form look").

Background (see single_core_signature_complete_monad_s2.py / OPEN-Q-055)
----------------------------------------------------------------------
A single-core complete-Omega tournament is a transitive order plus one core
vertex with arc pattern s in {0,1}^m.  The number of directed odd cycles
through the core is

    r(s) = sum_{i<j, s[i]=1, s[j]=0} f(j-i-1),   f(0)=1, f(t)=2^{t-1} (t>=1).

When Omega is complete (K_r) we get H = 1 + 2r.  S2 proved (finite theorem)
that r=3 (H=7) and r=10 (H=21) are PERMANENT single-core gaps, r=31 (H=63)
reachable.  The set of achievable r is dense-complement; this script studies
its complement, the GAP SET, to many terms for OEIS / closed-form analysis.

Fast r via the contribution recurrence
--------------------------------------
Let A_p = contribution a '0' at position p receives = sum_{i<p,s_i=1} f(p-i-1).
Then  A_0 = 0,  A_{p+1} = 2*A_p - s_{p-1} + s_p,  and  r = sum_{p: s_p=0} A_p.
(Derivation: advancing p doubles every prior 1's weight except the one at gap
0 (g: f(0)->f(1), both 1), then adds the new bit s_p as a gap-0 contributor.)
This is O(L) per string vs the O(L^2) brute, and is validated below.

Completeness length bound (S2): every achievable r in (0,R] has a canonical
witness (starts '1', ends '0') of length L <= 3 + floor(log2 R), since the
first-1/last-0 pair alone gives r >= 2^{L-3}.  So enumerating all canonical
strings to that length DECIDES achievability for ALL lengths, values <= R.
"""

from __future__ import annotations

import sys
from math import log2


# ---------------------------------------------------------------------------
# reference + fast r
# ---------------------------------------------------------------------------

def f_weight(t: int) -> int:
    return 1 if t == 0 else (1 << (t - 1))


def r_brute(bits: str) -> int:
    """O(L^2) reference (matches S2 exactly)."""
    ones = [i for i, c in enumerate(bits) if c == "1"]
    total = 0
    for j, c in enumerate(bits):
        if c == "0":
            for i in ones:
                if i >= j:
                    break
                total += f_weight(j - i - 1)
    return total


def r_fast(bits: str) -> int:
    """O(L) via A_{p+1} = 2*A_p - prev + cur ; r += A_p on every 0."""
    A = 0
    prev = 0          # s_{-1} treated as 0
    r = 0
    for ch in bits:
        cur = 1 if ch == "1" else 0
        if cur == 0:
            r += A
        A = 2 * A - prev + cur
        prev = cur
    return r


# ---------------------------------------------------------------------------
# validation of r_fast against r_brute
# ---------------------------------------------------------------------------

def validate() -> None:
    # exhaustive all strings up to length 14
    for L in range(1, 15):
        for mask in range(1 << L):
            bits = format(mask, f"0{L}b")
            assert r_fast(bits) == r_brute(bits), (bits,)
    # THM-344 signatures -> 31
    for sig in ("1001100", "1100110", "10011001", "11001101"):
        assert r_fast(sig) == 31 == r_brute(sig)
    # 1*0^k -> 2^{k-1}
    for k in range(1, 20):
        assert r_fast("1" + "0" * k) == (1 << (k - 1))
    print("[validate] r_fast == r_brute exhaustively for L<=14; "
          "THM-344 sigs->31; 1.0^k->2^{k-1}.  PASS")


# ---------------------------------------------------------------------------
# complete achievable-set + first-appearance length, up to value bound R
# ---------------------------------------------------------------------------

def achievable_up_to(R: int):
    """All achievable r in [1,R] with first-appearance witness length.

    Canonical strings: s_0=1, s_{L-1}=0, middle free, L=2..Lmax.
    Uses the O(L) recurrence inlined over the integer middle mask.
    """
    Lmax = 3 + int(log2(R)) if R >= 1 else 2
    first_len: dict[int, int] = {}
    for L in range(2, Lmax + 1):
        mid = L - 2
        # bits: position 0 = '1', positions 1..mid = middle, position L-1 = '0'
        for mm in range(1 << mid):
            # inline r over: cur sequence = [1] + middle bits (msb..lsb) + [0]
            A = 0
            prev = 0
            r = 0
            # s_0 = 1
            cur = 1
            # cur==1 so no r add; A = 2*0 - 0 + 1 = 1
            A = 1
            prev = 1
            # middle bits, from most significant (position 1) to least
            for b in range(mid - 1, -1, -1):
                cur = (mm >> b) & 1
                if cur == 0:
                    r += A
                A = 2 * A - prev + cur
                prev = cur
                if r > R:
                    break
            else:
                # final s_{L-1} = 0
                r += A          # cur=0 contributes A
                if 0 < r <= R and r not in first_len:
                    first_len[r] = L
    return first_len


def main() -> None:
    validate()
    print()

    R = 1 << 20  # 1,048,576 ; proves gaps for ALL lengths, values <= R
    Lmax = 3 + int(log2(R))
    print("SINGLE-CORE GAP-SET STRUCTURE  (monad-compute-2026-06-04-S3)")
    print(f"value bound R = {R} = 2^20 ;  rigorous length cap Lmax = {Lmax}")
    print("(complete: every achievable r in (0,R] has a canonical witness "
          f"of length <= {Lmax})")
    print("=" * 72)

    first_len = achievable_up_to(R)
    achievable = set(first_len)
    gaps = [r for r in range(1, R + 1) if r not in achievable]

    print(f"\nachievable r in [1,{R}] : {len(achievable)}")
    print(f"gap (unachievable) r in [1,{R}] : {len(gaps)}")
    print(f"gap fraction : {len(gaps)/R:.4f}")

    # --- OEIS-ready lists -------------------------------------------------
    ach_sorted = sorted(achievable)
    print("\nACHIEVABLE r, first 60 (OEIS-ready):")
    print("  " + ", ".join(map(str, ach_sorted[:60])))
    print("\nGAP r, first 60 (OEIS-ready):")
    print("  " + ", ".join(map(str, gaps[:60])))

    # H = 1+2r forms
    print("\nGAP as forbidden single-core H = 1+2r, first 40:")
    print("  " + ", ".join(str(1 + 2 * g) for g in gaps[:40]))

    # target verdicts
    print("\nTARGET VERDICTS (r -> H):")
    for r, h in [(3, 7), (6, 13), (10, 21), (31, 63), (94, 189)]:
        if r in achievable:
            print(f"  r={r:4d} H={h:5d}: ACHIEVABLE (first length {first_len[r]})")
        else:
            print(f"  r={r:4d} H={h:5d}: GAP (unachievable, all lengths)")

    # --- structural probes ------------------------------------------------
    print("\n" + "-" * 72)
    print("STRUCTURAL PROBES")

    # (1) closure under doubling: r achievable => 2r achievable ?
    dbl_fail = [r for r in ach_sorted if 2 * r <= R and (2 * r) not in achievable]
    print(f"\n(1) doubling closure  r in ACH & 2r<=R => 2r in ACH ? "
          f"{'YES' if not dbl_fail else 'NO'}  "
          f"(counterexamples: {dbl_fail[:10]})")
    # gaps closed under doubling?
    dbl_gap = [g for g in gaps if 2 * g <= R and (2 * g) in achievable]
    print(f"    gap g & 2g<=R => 2g GAP ? "
          f"{'YES' if not dbl_gap else 'NO'}  (gap->ach under x2: {dbl_gap[:10]})")

    # (2) are all powers of two achievable? (1.0^k witness)
    pw = [(1 << k) for k in range(0, R.bit_length()) if (1 << k) <= R]
    pw_ok = all(p in achievable for p in pw)
    print(f"\n(2) all powers of two in ACH ? {pw_ok}  (1.0^k witness)")

    # (3) gap differences (look for arithmetic structure)
    diffs = [gaps[i + 1] - gaps[i] for i in range(min(40, len(gaps) - 1))]
    print(f"\n(3) consecutive gap differences (first 40):\n  {diffs}")

    # (4) binary weight (popcount) distribution of gaps vs achievable (low range)
    def popc(x):
        return bin(x).count("1")
    lowR = min(R, 4096)
    ga = [g for g in gaps if g <= lowR]
    aa = [a for a in ach_sorted if a <= lowR]
    from collections import Counter
    print(f"\n(4) popcount distribution in [1,{lowR}]:")
    print(f"    gaps        : {dict(sorted(Counter(map(popc, ga)).items()))}")
    print(f"    achievable  : {dict(sorted(Counter(map(popc, aa)).items()))}")

    # (5) first-appearance length histogram
    from collections import Counter as C2
    lh = C2(first_len.values())
    print(f"\n(5) first-appearance length histogram (#new achievable r per L):")
    for L in sorted(lh):
        print(f"    L={L:2d}: {lh[L]}")

    # (6) density per dyadic window (is the gap set thinning or steady?)
    print(f"\n(6) gap density per dyadic window (k: (2^k, 2^{{k+1}}]):")
    for k in range(0, R.bit_length() - 1):
        lo, hi = (1 << k), (1 << (k + 1))
        if hi > R:
            break
        w = sum(1 for g in gaps if lo < g <= hi)
        print(f"    ({lo:7d},{hi:7d}]: {w:6d} gaps / {hi-lo:7d}  "
              f"= {100*w/(hi-lo):6.2f}%")

    # (7) largest gap (finiteness probe)
    print(f"\n(7) largest gap value <= {R}: {gaps[-1] if gaps else None}")
    print(f"    (gap set appears {'INFINITE/persistent' if gaps and gaps[-1] > R//2 else 'possibly finite'})")

    # (8) closed-form candidates for a density-1/2 set --------------------
    print("\n(8) closed-form candidates (density 1/2 sets):")
    gapset = set(gaps)
    # residue classes mod q: is the gap set a union of residue classes?
    for q in range(2, 13):
        # a residue class c mod q is 'pure gap' if every k in [1,R]
        # with k%q==c is a gap; 'pure ach' if none is.
        pure_gap = [c for c in range(q)
                    if all((k in gapset) for k in range(c if c else q, R + 1, q))]
        pure_ach = [c for c in range(q)
                    if all((k not in gapset) for k in range(c if c else q, R + 1, q))]
        if pure_gap or pure_ach:
            print(f"    mod {q}: pure-gap residues {pure_gap}, pure-ach residues {pure_ach}")
    print("    (no line above => gap set is NOT a union of residue classes mod<=12)")

    # Thue-Morse (odious/evil) -- the canonical density-1/2 set
    def tm(x):
        return bin(x).count("1") & 1
    tm_gap_odious = sum(1 for g in gaps if tm(g))
    tm_ach_odious = sum(1 for a in ach_sorted if tm(a))
    print(f"    Thue-Morse: gaps odious(odd popcount)={tm_gap_odious}/{len(gaps)} "
          f"({100*tm_gap_odious/len(gaps):.1f}%);  "
          f"ach odious={tm_ach_odious}/{len(ach_sorted)} "
          f"({100*tm_ach_odious/len(ach_sorted):.1f}%)  "
          f"=> {'MATCHES Thue-Morse' if tm_gap_odious in (0,len(gaps)) else 'NOT Thue-Morse'}")

    # Beatty: is gap set = {floor(n*alpha)} for some alpha? test via the
    # complementary-Beatty signature (consecutive diffs in {d,d+1}).
    ddiffs = [gaps[i+1]-gaps[i] for i in range(len(gaps)-1)]
    print(f"    Beatty test: distinct gap-gaps = {sorted(set(ddiffs))[:12]}... "
          f"=> {'two-valued (Beatty-like)' if len(set(ddiffs))<=2 else 'NOT a Beatty sequence (>2 gap sizes)'}")

    # (9) is ACH closed under  r -> r + (largest power of 2 <= r) ? (probe)
    #     and additive: is ach + ach dense? quick semigroup probe
    small_ach = [a for a in ach_sorted if a <= 64]
    print(f"\n(9) is ACH an additive semigroup? "
          f"a,b in ACH(<=32) => a+b in ACH ? ", end="")
    bad = []
    for a in small_ach:
        for b in small_ach:
            if a <= 32 and b <= 32 and (a + b) <= R and (a + b) not in achievable:
                bad.append((a, b, a + b))
    print(f"{'YES' if not bad else 'NO'}  (first fails: {bad[:6]})")

    # dump full lists for OEIS / downstream (capped)
    print("\n" + "=" * 72)
    print("FULL GAP LIST (all gaps <= 2^14 = 16384) for archival:")
    g16 = [g for g in gaps if g <= 16384]
    print("  " + ", ".join(map(str, g16)))


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
