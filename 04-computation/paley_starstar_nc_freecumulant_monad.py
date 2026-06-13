#!/usr/bin/env python3
"""
paley_starstar_nc_freecumulant_monad.py
monad-explorer-2026-06-07 (deep-research, 6th session)

The TRUE genus-0 object behind (★★).  Since S_k = (−1)^k C_k = κ_{2k+2}/A^{k+1} is
the free cumulant of the two-point law ν=½(δ_a+δ_{−a}) (a²=A), and free cumulants
are the NON-CROSSING moment-Möbius transform, the clean genus-0 statement is

   κ_{2k+2} = Σ_{π ∈ NC(2k+2), ALL blocks even}  μ_NC(π, 1̂) · A^{(2k+2)/2},

i.e.   S_k = Σ_{π ∈ NC(2k+2), even blocks}  μ_NC(π, 1̂_{2k+2})   = (−1)^k C_k.

This script confirms it by computing free cumulants directly from the moment-cumulant
recursion over NC partitions (m_{2j}=A^j, A=1, m_odd=0), and counts the NC-even
partitions of [2k+2], comparing to the EVEN-SERIES pattern counts (1,3,13,67,383 =
A215257) and to the REFUTED laminar sum (−1,2,−6,25,−132).

Purpose: show the bridge "even-series patterns ↔ NC-even partitions" is NON-trivial
(different ground sets 2k+1 vs 2k+2, different counts) — the remaining content of
the proof of (★★).

NO number theory.
"""
import math
from functools import lru_cache


def catalan(k):
    return math.comb(2 * k, k) // (k + 1)


# ---- non-crossing partitions of [n] (n points on a line/circle) ----
def nc_partitions(n):
    """Yield all non-crossing partitions of {0,..,n-1} as list of frozensets."""
    pts = list(range(n))

    def helper(elems):
        # elems: sorted list of points to partition non-crossingly
        if not elems:
            yield []
            return
        first = elems[0]
        # the block containing 'first' must be such that points between consecutive
        # block members are partitioned among themselves (non-crossing).
        # choose the block B containing first: it is a subset of elems with first,
        # such that B is "non-crossing-compatible": between consecutive members of B,
        # the points form independent sub-intervals.
        rest = elems[1:]
        # iterate over which of the remaining points join first's block
        m = len(rest)
        for mask in range(1 << m):
            block = [first] + [rest[i] for i in range(m) if (mask >> i) & 1]
            block_set = set(block)
            others = [x for x in elems if x not in block_set]
            # non-crossing check: 'others' must split into intervals between block pts
            # equivalently: no two block points "interleave" an others-run that crosses.
            # Standard: a partition is NC iff for the block of `first`, the gaps between
            # consecutive block elements are filled by points all in `others`, and we
            # recurse on each maximal gap-run independently.
            block_sorted = sorted(block)
            ok = True
            sub_runs = []
            # points outside block, grouped by which gap they fall in
            # gaps: (-inf, b0), (b0,b1), ..., (b_last, +inf) restricted to elems
            # For NC: every others point lies strictly between two consecutive block pts
            # OR after the last / before the first — but 'first' is the minimum, so none before.
            # Build runs between consecutive block elements; trailing run after last block elt.
            extended = block_sorted + [math.inf]
            gi = 0
            run = []
            o_idx = 0
            others_sorted = sorted(others)
            for x in others_sorted:
                while x > extended[gi + 1] if gi + 1 < len(extended) else False:
                    if run:
                        sub_runs.append(run); run = []
                    gi += 1
                # x should be between extended[gi] and extended[gi+1]
                # advance gi until x < extended[gi+1]
                while gi + 1 < len(extended) and x > extended[gi + 1]:
                    if run:
                        sub_runs.append(run); run = []
                    gi += 1
                run.append(x)
            if run:
                sub_runs.append(run)
            # Recurse: each sub_run partitioned NC among itself
            def comb_runs(runs):
                if not runs:
                    yield []
                    return
                for p0 in helper(runs[0]):
                    for prest in comb_runs(runs[1:]):
                        yield p0 + prest
            for restpart in comb_runs(sub_runs):
                yield [block_sorted] + restpart

    for part in helper(pts):
        yield [frozenset(b) for b in part]


def is_nc_simple(part, n):
    """Robust NC test: no two blocks cross. Blocks cross if a<c<b<d with a,b in one, c,d in other."""
    blocks = [sorted(b) for b in part]
    for i in range(len(blocks)):
        for j in range(len(blocks)):
            if i == j:
                continue
            B, C = blocks[i], blocks[j]
            for a in B:
                for b in B:
                    if a >= b:
                        continue
                    for c in C:
                        if a < c < b:
                            for d in C:
                                if c < d and b < d:
                                    return False
    return True


def all_partitions(n):
    elems = list(range(n))

    def go(elems):
        if not elems:
            yield []
            return
        f = elems[0]
        for sm in go(elems[1:]):
            for i in range(len(sm)):
                yield sm[:i] + [[f] + sm[i]] + sm[i + 1:]
            yield [[f]] + sm
    for p in go(elems):
        yield [frozenset(b) for b in p]


def free_cumulants(maxn, A=1):
    """kappa_n for n=1..maxn from moments m_{2j}=A^j, m_odd=0 via NC moment-cumulant."""
    # m_n = sum_{pi in NC(n)} prod_B kappa_{|B|}
    kappa = {}
    for n in range(1, maxn + 1):
        m_n = (A ** (n // 2)) if n % 2 == 0 else 0
        # subtract contributions of NC partitions with a non-full block, expressed in kappa
        total = 0
        full_seen = False
        for part in all_partitions(n):
            if not is_nc_simple(part, n):
                continue
            if len(part) == 1:
                full_seen = True
                continue  # this is the kappa_n term
            prod = 1
            for B in part:
                prod *= kappa.get(len(B), 0)
            total += prod
        kappa[n] = m_n - total
    return kappa


def main():
    print("=" * 72)
    print("TRUE genus-0 object: free cumulants of two-point law (NC even-partition sum)")
    print("=" * 72)
    KMAX = 5
    kappa = free_cumulants(2 * KMAX + 2, A=1)
    print("\n free cumulants kappa_{2j} (A=1):")
    for j in range(1, KMAX + 2):
        tgt = (-1) ** (j - 1) * catalan(j - 1)
        print(f"   kappa_{2*j} = {kappa[2*j]:>5}   (−1)^(j−1)C_(j−1) = {tgt:>5}  "
              f"{'OK' if kappa[2*j]==tgt else 'X'}")

    print("\n (★★) as NC-even free cumulant:  S_k = kappa_{2k+2} = (−1)^k C_k")
    for k in range(1, KMAX + 1):
        S_k = kappa[2 * k + 2]
        tgt = (-1) ** k * catalan(k)
        # count NC-even partitions of [2k+2]
        n = 2 * k + 2
        cnt = 0
        for part in all_partitions(n):
            if not is_nc_simple(part, n):
                continue
            if all(len(B) % 2 == 0 for B in part):
                cnt += 1
        print(f"   k={k}: S_k={S_k:>5}  target={tgt:>5}  {'OK' if S_k==tgt else 'X'}"
              f"   #NC-even-part[{n}]={cnt}")

    print("\n COMPARISON of gradings (all summing/relating to (−1)^k C_k = −1,2,−5,14,−42):")
    print("   even-series pattern count (A215257):     1, 3, 13, 67, 383   [full lattice, 2k+1 pts]")
    print("   #NC-even partitions of [2k+2]:           (printed above)     [NC lattice, 2k+2 pts]")
    print("   REFUTED laminar / canonical-ribbon sum:  −1, 2, −6, 25, −132 [wrong non-crossing notion]")
    print("=" * 72)
    print("Conclusion: the genus-0 truth lives on NC-even partitions of 2k+2 points, a")
    print("DIFFERENT ground set & lattice from the even-series patterns (2k+1 path pts,")
    print("full lattice).  The proof of (★★) = an explicit bridge between them — NOT the")
    print("laminarity of σ, NOT the walk-ribbon genus, NOT naive first-return (all ruled out).")


if __name__ == "__main__":
    main()
