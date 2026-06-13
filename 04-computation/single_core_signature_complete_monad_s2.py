#!/usr/bin/env python3
"""
monad-compute-2026-06-04-S2

COMPLETE achievable-set computation for the single-core weighted signature
count r (OPEN-Q-055 sub-question: "does r_core(s) ever equal 3 or 10?").

Background
----------
A single-core complete-Omega tournament is a transitive order v_0>...>v_{m-1}
plus one extra "core" vertex whose arc pattern is a bit string s in {0,1}^m
(s[i]=1 : core->v_i ; s[i]=0 : v_i->core).  The number of directed odd cycles
through the core is

    r(s) = sum_{i<j, s[i]=1, s[j]=0} f(j-i-1),   f(0)=1, f(t)=2^{t-1} (t>=1)

(each 1...0 pair admits an even-sized interior subset of the transitive chain;
#even subsets of a t-set = 1 if t=0 else 2^{t-1}).  When Omega is complete on
these r cycles, H = I(K_r, 2) = 1 + 2r.  So r=3 <-> H=7, r=10 <-> H=21,
r=31 <-> H=63.

The S11/S12 sessions only reported "r=3,10 absent through length 16 (resp 40)".
This script PROVES the gap for ALL lengths up to a value bound R, via:

  LENGTH BOUND.  Stripping leading 0s (no 1 to their left) and trailing 1s
  (no 0 to their right) does not change r.  A nonzero r therefore has a
  *canonical* witness s' that starts with 1 and ends with 0.  If |s'| = L >= 3,
  the (first 1, last 0) pair alone contributes f(L-2) = 2^{L-3}, so
  r >= 2^{L-3}; if L = 2 then r >= 1.  Hence every achievable value r in (0,R]
  has a canonical witness of length L <= 3 + floor(log2(R)).

So enumerating ALL canonical strings up to that length finds EVERY achievable
value <= R.  Any value in [1,R] not found is provably unachievable at any
length (a genuine permanent single-core gap), not merely "absent so far".
"""

from __future__ import annotations

import sys
from math import log2

# ---------------------------------------------------------------------------
# core combinatorics
# ---------------------------------------------------------------------------

def f_weight(t: int) -> int:
    """#even-sized subsets of a t-element interior set."""
    return 1 if t == 0 else (1 << (t - 1))


def r_brute(bits: str) -> int:
    """O(L^2) reference: sum f(j-i-1) over 1...0 pairs."""
    ones = [i for i, c in enumerate(bits) if c == "1"]
    total = 0
    for j, c in enumerate(bits):
        if c == "0":
            for i in ones:
                if i >= j:
                    break
                total += f_weight(j - i - 1)
    return total


# ---------------------------------------------------------------------------
# self-test of the length bound on a single string + match vs S11 convention
# ---------------------------------------------------------------------------

def selftest() -> None:
    # known THM-344 single-core signatures (r should be 31 -> H=63)
    for sig in ("1001100", "1100110", "10011001", "11001101"):
        assert r_brute(sig) == 31, (sig, r_brute(sig))
    # canonicalisation (strip leading 0 / trailing 1) is r-invariant
    for s in ("101010", "1110000", "10", "100", "1000", "110", "11010",
              "0010011000", "0001001100111", "010101"):
        stripped = s.lstrip("0").rstrip("1")
        assert r_brute(s) == r_brute(stripped), (s, stripped)
    # single 1 followed by k zeros gives 2^{k-1}
    for k in range(1, 12):
        assert r_brute("1" + "0" * k) == (1 << (k - 1)), k
    print("[selftest] passed: THM-344 sigs -> r=31; canonicalisation r-invariant; 1·0^k -> 2^{k-1}")


# ---------------------------------------------------------------------------
# complete achievable-set enumeration up to value bound R
# ---------------------------------------------------------------------------

def achievable_up_to(R: int) -> tuple[set[int], dict[int, int]]:
    """Return (achievable values in [1,R], first-appearance length per value).

    Enumerates every canonical string (starts '1', ends '0') of length
    L = 2 .. Lmax, where Lmax = 3 + floor(log2(R)) guarantees completeness.
    """
    Lmax = 3 + int(log2(R)) if R >= 1 else 2
    achievable: set[int] = set()
    first_len: dict[int, int] = {}
    # L = 1 strings ("0","1") give r=0, ignore (r>0 only).
    for L in range(2, Lmax + 1):
        # canonical: position 0 = '1', position L-1 = '0', middle free
        mid = L - 2
        for m in range(1 << mid):
            bits = "1" + format(m, f"0{mid}b") + "0" if mid else "10"
            r = r_brute(bits)
            if 0 < r <= R and r not in achievable:
                achievable.add(r)
                first_len[r] = L
    return achievable, first_len


def main() -> None:
    selftest()
    print()

    R = 1 << 17  # 131072 -> proves gaps for ALL string lengths, values <= R
    Lmax = 3 + int(log2(R))
    print(f"COMPLETE single-core achievable-set audit  (monad-compute-2026-06-04-S2)")
    print(f"value bound R = {R} = 2^{R.bit_length()-1};  rigorous length cap Lmax = {Lmax}")
    print(f"(every achievable r in (0,{R}] has a canonical witness of length <= {Lmax})")
    print("=" * 72)

    achievable, first_len = achievable_up_to(R)

    # the H<->r dictionary targets
    print("\nTARGET VERDICTS (r -> H = 1+2r):")
    for r, h in [(3, 7), (10, 21), (31, 63), (94, 189)]:
        hit = r in achievable
        extra = f" (first appears at length {first_len[r]})" if hit else " (PERMANENT GAP: no string of any length)"
        print(f"  r={r:4d}  H={h:5d} : {'ACHIEVABLE' if hit else 'UNACHIEVABLE'}{extra}")

    # full gap structure
    gaps = [r for r in range(1, R + 1) if r not in achievable]
    print(f"\nFULL GAP SET in [1,{R}] : {len(gaps)} unachievable values")
    print(f"  smallest 40 gaps: {gaps[:40]}")
    # map to forbidden single-core H values
    gap_H = [1 + 2 * r for r in gaps]
    print(f"  -> as H=1+2r, smallest 40: {gap_H[:40]}")

    # achievable density and structure of the low end
    print(f"\nACHIEVABLE values in [1,64]: {sorted(v for v in achievable if v <= 64)}")
    print(f"first-appearance length for achievable r<=64:")
    for r in sorted(v for v in achievable if v <= 64):
        print(f"    r={r:3d} -> H={1+2*r:3d}  first length {first_len[r]}")

    # is the gap set finite? report largest gap and tail density
    print(f"\nlargest gap value <= {R}: {gaps[-1] if gaps else None}")
    for lo in (0, R // 2):
        hi = lo + R // 2
        g = sum(1 for r in gaps if lo < r <= hi)
        print(f"  gaps in ({lo},{hi}]: {g}  ({100*g/(R//2):.3f}% of that window)")

    # cross-check vs the S11 brute-force convention (all strings length m<=16)
    print("\nCROSS-CHECK vs S11 (all 2^m strings, m=2..16; targets 3,10,31):")
    for m in range(2, 17):
        present = {3: False, 10: False, 31: False}
        for mask in range(1 << m):
            bits = format(mask, f"0{m}b")
            rv = r_brute(bits)
            if rv in present:
                present[rv] = True
            if all(present.values()):
                break
        flags = ", ".join(f"r={t}:{'present' if present[t] else 'absent'}" for t in (3, 10, 31))
        print(f"  m={m:2d}: {flags}")
        if m >= 16:
            break


if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    main()
