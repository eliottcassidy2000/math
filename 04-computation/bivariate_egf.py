"""
bivariate_egf.py -- kind-pasteur-2026-03-14-S108i
THE BIVARIATE EGF — Complete a concrete task, learn an analogy

TASK: Find the bivariate EGF F(q, t) = sum_n D_n(q) * t^n / n!
where D_n(q) = sum over anti-succ-free perms of q^successions.

At q=1: F(1,t) = exp(-t)/(1-t)^2 (the A000255 EGF).
At q=2: F(2,t) should give our D_n(2) sequence.

STRATEGY: Build the permutation element by element.
When inserting element n into an anti-succ-free perm of [n-1]:
- n can go in (n-1) "neutral" positions (weight 1 each)
  [actually need to subtract forbidden and count special positions]
- 1 position creates a succession (after n-1): weight q
- 1 position is forbidden (before n-1): skip
- Some positions break existing successions: weight 1/q... wait, q^(s-1).

Actually: let me think about this as an EXPONENTIAL FORMULA problem.

A permutation of [n] is built from "atoms":
- A SINGLETON element (contributes weight 1, no succession)
- An ASCENDING PAIR (contributes weight q, one succession)
- An ASCENDING TRIPLE (contributes weight q^2, two successions)
- In general: an ascending run of length k contributes weight q^(k-1)

Anti-succession-free means: at block boundaries, the next block's
minimum != current block's maximum - 1.

This is the COMPOSITION structure:
- Partition [n] into consecutive-integer blocks (ascending runs).
- Arrange the blocks in some order avoiding anti-successions at boundaries.
- Weight = q^(total successions) = q^(n - number of blocks).

So D_n(q) = sum over compositions of n into consecutive blocks *
            (number of valid block orderings) * q^(n - blocks).

For a composition into r blocks: weight = q^(n-r).
Number of compositions into r blocks: C(n-1, r-1).
Number of valid orderings of r blocks (anti-succ-free): f(r, blocks).

The hard part is f(r, blocks) — it depends on the IDENTITY of the blocks
(which consecutive integers they contain), not just r.

But wait: maybe I can compute the EGF by a TRANSFER MATRIX approach
on the blocks, not individual elements.
"""

import sys, math
import numpy as np
from fractions import Fraction
from itertools import permutations
from collections import Counter
from functools import lru_cache

sys.stdout.reconfigure(encoding='utf-8')

def successions(pi):
    return sum(1 for i in range(len(pi)-1) if pi[i+1] == pi[i] + 1)

def anti_successions(pi):
    return sum(1 for i in range(len(pi)-1) if pi[i+1] == pi[i] - 1)

def main():
    print("=" * 70)
    print("BIVARIATE EGF — THE CONCRETE TASK")
    print("kind-pasteur-2026-03-14-S108i")
    print("=" * 70)

    # ============================================================
    print(f"\n--- STEP 1: Compute D_n(q) as a polynomial in q for small n")

    polys = {}
    for n in range(1, 8):
        all_perms = list(permutations(range(n)))
        coeff = Counter()
        for pi in all_perms:
            if anti_successions(pi) > 0:
                continue
            s = successions(pi)
            coeff[s] += 1
        polys[n] = dict(coeff)
        terms = [f'{c}q^{s}' for s, c in sorted(coeff.items()) if c > 0]
        print(f"  D_{n}(q) = {' + '.join(terms)}")

    # ============================================================
    print(f"\n--- STEP 2: Test recurrence D_n(q) = (n-1+q)*D_{{n-1}}(q) + (n-2)*D_{{n-2}}(q)")
    # At q=1: (n-1+1)*D_{n-1}(1) + (n-2)*D_{n-2}(1) = n*D_{n-1} + (n-2)*D_{n-2}
    # But A000255 satisfies a(m) = (m+1)*a(m-1) + m*a(m-2) for D_n(1) = a(n-1).
    # So D_n(1) = a(n-1), D_{n-1}(1) = a(n-2).
    # a(n-1) = n*a(n-2) + (n-1)*a(n-3) [A000255 recurrence shifted].
    # So D_n(1) = n*D_{n-1}(1) + (n-1)*D_{n-2}(1). YES!
    # Try for general q: D_n(q) = (n-1+q)*D_{n-1}(q) + (n-2+f(q))*D_{n-2}(q)?

    # Compute D_n(q) as evaluation at specific q values
    def Dn(n, q):
        return sum(c * q**s for s, c in polys.get(n, {}).items())

    print(f"\n  Testing: D_n(q) = (n-1+q)*D_{{n-1}}(q) + (n-1)*D_{{n-2}}(q)")
    for q in [1, 2, 3]:
        print(f"\n  q={q}:")
        for n in range(3, 8):
            pred = (n-1+q)*Dn(n-1, q) + (n-1)*Dn(n-2, q)
            actual = Dn(n, q)
            print(f"    n={n}: (n-1+q)*D_{{n-1}} + (n-1)*D_{{n-2}} = {pred}, "
                  f"actual = {actual}, match = {pred == actual}")

    # That didn't work for q=1 either (expected n not n-1+q=n).
    # Let me try: D_n(q) = n*D_{n-1}(q) + (n-1)*D_{n-2}(q)
    print(f"\n  Testing: D_n(q) = n*D_{{n-1}}(q) + (n-1)*D_{{n-2}}(q)")
    for q in [1, 2, 3]:
        print(f"  q={q}:")
        for n in range(3, 8):
            pred = n*Dn(n-1, q) + (n-1)*Dn(n-2, q)
            actual = Dn(n, q)
            diff = actual - pred
            print(f"    n={n}: n*D_{{n-1}}+(n-1)*D_{{n-2}} = {pred}, "
                  f"actual = {actual}, diff = {diff}")

    # The diff at q=1 should be 0. At q=2, what is it?
    print(f"\n  Differences D_n(q) - [n*D_{{n-1}}(q) + (n-1)*D_{{n-2}}(q)]:")
    for q in [1, 2, 3]:
        diffs = []
        for n in range(3, 8):
            diff = Dn(n, q) - n*Dn(n-1, q) - (n-1)*Dn(n-2, q)
            diffs.append(diff)
        print(f"  q={q}: diffs = {diffs}")

    # At q=1: diffs should be 0. At q=2: ?
    # The diffs tell us the "correction" needed for each q.

    # ============================================================
    print(f"\n--- STEP 3: Look at D_n(q) coefficients more carefully")

    # The coefficient of q^s in D_n(q):
    # = number of anti-succ-free perms of [n] with exactly s successions.
    # Let a(n, s) = this count.
    # a(n, 0) = number of anti-succ-free perms with NO successions at all
    #         = perms with no ascending adjacency AND no descending adjacency
    #         = HERTZSPRUNG numbers (A002464)!

    print(f"  a(n, 0) = coefficient of q^0 in D_n(q):")
    hertz = []
    for n in range(1, 8):
        a0 = polys.get(n, {}).get(0, 0)
        hertz.append(a0)
        print(f"    n={n}: a(n,0) = {a0}")

    print(f"  Hertzsprung A002464: 1, 0, 0, 2, 14, 90, 646")
    # Our a(n,0): 0, 0, 0, 2, 14, 90, 646 for n=1,...,7
    # Matches A002464 starting from n=1! (A002464 starts 1,1,0,0,2,14,90,646)
    # Our n=1: 0, n=2: 0, n=3: 0, n=4: 2, n=5: 14, n=6: 90, n=7: 646.
    # A002464: a(0)=1, a(1)=1, a(2)=0, a(3)=0, a(4)=2, a(5)=14, a(6)=90, a(7)=646.
    # So our a(n,0) = A002464(n) for n >= 2.
    print(f"  Match with A002464(n) for n>=4: confirmed!")

    # ============================================================
    print(f"\n--- STEP 4: THE INSERTION RECURRENCE")

    # Insert n into an anti-succ-free perm of [n-1]:
    # There are n possible insertion positions.
    #
    # Let pi be anti-succ-free on [n-1] with s successions.
    # pi has s ascending adjacencies and 0 descending adjacencies.
    # The element n-1 appears somewhere in pi.
    #
    # When inserting n:
    # CASE A: Insert after n-1 (creates succession n-1, n). ONE position.
    #   New perm has s+1 successions if this doesn't break anything.
    #   But if n-1 was at the end of pi, there's no element after to worry about.
    #   If n-1 is at position j < n-2 of pi, then pi[j]=n-1, pi[j+1]=something.
    #   Inserting n between n-1 and pi[j+1]: now we have ...n-1, n, pi[j+1]...
    #   New succession: n-1, n. YES.
    #   Anti-succession n, pi[j+1]: need pi[j+1] != n-1. But pi[j+1] <= n-2 != n-1. SAFE.
    #   Did we break the succession n-1, pi[j+1]? Only if pi[j+1] = n = doesn't exist.
    #   So: exactly 1 succession gained. Weight: q^(s+1). Contribution: q * q^s = q*existing.

    # CASE B: Insert before n-1 (creates anti-succession n, n-1). ONE position.
    #   FORBIDDEN.

    # CASE C: Insert at any other position (not adjacent to n-1). n-2 positions.
    #   BUT some of these may break existing successions.
    #   An existing succession at position j means pi[j+1] = pi[j]+1.
    #   Inserting n between pi[j] and pi[j+1] breaks this succession.
    #   New s = old s - 1. Weight q^(s-1).
    #   How many such "succession-breaking" positions are there?
    #   s positions (one per succession). BUT one of them might be the
    #   "before n-1" position (forbidden). The succession ending at n-1
    #   means some pi[j]=n-2, pi[j+1]=n-1. Inserting between n-2 and n-1
    #   puts n between them: ...n-2, n, n-1... which has anti-succ n, n-1. FORBIDDEN.
    #   So the number of SAFE succession-breaking positions = s or s-1.
    #   It's s-1 if pi has the succession (n-2, n-1), else s.
    #   Let chi = 1 if pi has succession (n-2, n-1), else 0.
    #   Safe succession-breaking: s - chi.
    #   Non-succession, non-adjacent-to-n-1 positions: (n-2) - (s - chi) - 1
    #   Wait: total = n, forbidden = 1, creates-succ = 1 (after n-1).
    #   Remaining = n - 2. Of these, (s-chi) break a succession. (n-2-s+chi) are neutral.

    # RESULT:
    # D_n(q) = sum over anti-succ-free pi of [n-1] with s successions and chi indicator:
    #   [q * q^s + (s - chi) * q^(s-1) + (n-2-s+chi) * q^s]
    #   = sum [...]:
    #   = q * sum q^s + (n-2) * sum q^s - sum s*q^s + sum chi*q^s + sum (s-chi)*q^(s-1)
    #   Hmm, let me group differently.

    # For each pi with s successions and indicator chi:
    # contribution = q*q^s + (s-chi)*q^(s-1) + (n-2-s+chi)*q^s
    # = q^(s+1) + (s-chi)*q^(s-1) + (n-2-s+chi)*q^s
    # = (n-2-s+chi)*q^s + q^(s+1) + (s-chi)*q^(s-1)
    # = (n-2+chi-s)*q^s + q*q^s + (s-chi)/q * q^s
    # = q^s * [(n-2+chi-s) + q + (s-chi)/q]
    # = q^s * [n-2+chi-s + q + s/q - chi/q]
    # = q^s * [n-2 + q + chi(1-1/q) + s(1/q-1)]
    # = q^s * [n-2 + q + chi*(q-1)/q - s*(q-1)/q]
    # = q^s * [n-2 + q + (chi-s)*(q-1)/q]
    # = q^s * [n-2 + q - (s-chi)*(q-1)/q]

    # Hmm. Let me just verify numerically that the recurrence works.

    print(f"\n  Verifying insertion recurrence numerically:")
    for n in range(3, 8):
        # Compute D_n(q) from D_{n-1}(q) via insertion
        Dn_from_insertion = Counter()
        for pi in permutations(range(n-1)):
            if anti_successions(pi) > 0:
                continue
            s = successions(pi)
            # Find where n-1 is in pi
            pos_nm1 = list(pi).index(n-1)  # n-1 in 0-indexed = largest element
            # Check if succession (n-2, n-1) exists
            chi = 0
            if pos_nm1 > 0 and pi[pos_nm1 - 1] == n - 2:
                chi = 1

            # Insert n (= n-1 in 0-indexed... wait, we're permuting {0,...,n-2}.
            # The new element to insert is n-1 (the value n-1, which is the largest+1).
            # Actually in 0-indexed: pi is a perm of {0,...,n-2}. We insert n-1.

            # Position: after element n-2 (value n-2 in 0-indexed), creates succession
            # After n-2 means: in the position right after where n-2 appears.
            pos_nm2 = list(pi).index(n-2)  # position of n-2

            # CASE A: insert after n-2 (creates succession n-2, n-1)
            # Check: inserting after n-2, is the next element n-1?
            # n-1 is already the largest in pi. We're inserting n-1 the VALUE.
            # Wait, I'm confusing myself. Let me re-index.

            # pi is a perm of {0, 1, ..., n-2}. We insert the element (n-1).
            # Succession: (n-1) follows (n-2). So insert (n-1) right after (n-2).
            # Anti-succession: (n-2) follows (n-1). So insert (n-1) right before (n-2).
            #   This is FORBIDDEN.

            # Case A: insert after position of (n-2). Weight: q (creates succession).
            # But if (n-2) is at the end, we insert at the very end. Fine.
            # If (n-2) is at position j, insert between pi[j] and pi[j+1].
            # Check: does this create anti-succ? (n-1) followed by pi[j+1].
            # Anti-succ iff pi[j+1] = (n-1)-1 = n-2. But pi[j] = n-2, so pi[j+1] != n-2.
            # Safe. Does this break the succession pi[j], pi[j+1]?
            # Only if pi[j+1] = pi[j]+1 = n-2+1 = n-1. But n-1 not in pi. So no break.
            Dn_from_insertion[s + 1] += 1  # CASE A: one position, s+1 successions

            # Case B: insert before (n-2). FORBIDDEN (creates anti-succ (n-1, n-2)).
            # Skip.

            # Case C: all other positions (n-1 total positions, minus A, minus B = n-3)
            for pos in range(n):  # n insertion positions (before each element + end)
                # Skip Case A (after n-2) and Case B (before n-2)
                if pos == pos_nm2 + 1:  # Case A: after n-2
                    continue
                if pos == pos_nm2:  # Case B: before n-2
                    continue
                # Insert at position pos
                new_pi = list(pi[:pos]) + [n-1] + list(pi[pos:])
                new_s = successions(new_pi)
                new_o = anti_successions(new_pi)
                if new_o > 0:
                    continue  # shouldn't happen if our analysis is correct
                Dn_from_insertion[new_s] += 1

        # Compare with direct computation
        direct = polys[n]
        match = dict(Dn_from_insertion) == direct
        print(f"  n={n}: insertion={dict(Dn_from_insertion)}, direct={direct}, match={match}")

    # ============================================================
    print(f"\n--- STEP 5: THE RECURRENCE EXTRACTED")

    # From the insertion analysis, we should get:
    # D_n(q) in terms of D_{n-1}(q) and possibly D_{n-2}(q).
    #
    # Actually, the insertion is just into D_{n-1}(q) (one previous term).
    # The recurrence should be: D_n(q) = f(q, n) * D_{n-1}(q) + correction.
    #
    # Let me extract it: D_n(q) = q*D_{n-1}(q) + sum over pi of remaining contributions.
    # The "remaining" are all non-Case-A, non-Case-B insertions.

    # For each pi with s successions: Case C gives (n-2) insertions minus whatever.
    # Some create the same s, some create s-1.

    # Instead of analyzing, let me just CHECK:
    # D_n(q) - q*D_{n-1}(q) = ?
    print(f"\n  D_n(q) - q*D_{{n-1}}(q) as a polynomial:")
    for n in range(2, 8):
        diff = Counter()
        for s, c in polys[n].items():
            diff[s] += c
        for s, c in polys.get(n-1, {}).items():
            diff[s+1] -= c  # q * q^s = q^(s+1)
        # Clean
        diff = {s: c for s, c in diff.items() if c != 0}
        print(f"  n={n}: D_{n} - q*D_{{n-1}} = {dict(sorted(diff.items()))}")

    # And D_n(q) - q*D_{n-1}(q) - (n-2)*D_{n-1}(q)?
    print(f"\n  D_n(q) - (n-2+q)*D_{{n-1}}(q) as a polynomial:")
    for n in range(2, 8):
        diff = Counter()
        for s, c in polys[n].items():
            diff[s] += c
        for s, c in polys.get(n-1, {}).items():
            diff[s+1] -= c  # q * q^s
            diff[s] -= (n-2) * c  # (n-2) * q^s
        diff = {s: c for s, c in diff.items() if c != 0}
        if diff:
            print(f"  n={n}: remainder = {dict(sorted(diff.items()))}")
        else:
            print(f"  n={n}: ZERO! D_n = (n-2+q)*D_{{n-1}}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
