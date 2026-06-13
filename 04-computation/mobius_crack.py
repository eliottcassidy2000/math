"""
mobius_crack.py -- kind-pasteur-2026-03-14-S108g
MOBIUS INVERSION + CATEGORY THEORY approach

The identity: D_n(2) = n! + 2*sum_k (n-2k)^k*(n-2k)!

STOP thinking about individual permutations.
Think about the POSET of set partitions, or the LATTICE of subsets.

Key idea: the succession structure of a permutation defines a
SET PARTITION of [n] into BLOCKS (consecutive ascending runs).

A permutation of [n] with exactly s successions has exactly (n-s) blocks.
Each block is a sequence of consecutive integers in ascending order.

The number of permutations with a GIVEN block structure (sizes b_1,...,b_r
where sum b_i = n) is n! / (b_1! * ... * b_r!) * (something)... wait.
No: the blocks are ORDERED (they appear in the permutation in a specific order)
and within each block the elements are FIXED (ascending). So a block structure
with r blocks just needs to choose WHICH elements go in which block,
then arrange the r blocks in order.

Actually: a permutation with block sizes (b_1,...,b_r) from left to right
means: first b_1 consecutive ascending elements, then b_2 consecutive, etc.
The elements in each block must be consecutive integers.
There are C(n-1, b_1-1, b_2-1, ..., b_r-1) / (something)... hmm.

Let me think more carefully.

A permutation pi of [n] = {0,...,n-1} has s successions.
The successions divide pi into ASCENDING RUNS:
  pi = [run_1 | run_2 | ... | run_r]
where within each run, elements are consecutive AND ascending.
There are r = n - s runs (each succession merges two adjacent elements).

Example: pi = (0, 1, 3, 4, 2) has successions at positions 0 (0->1)
and 2 (3->4). Runs: [0,1], [3,4], [2]. Three runs (r=3, s=2).

The elements within each run are consecutive integers.
Run [0,1]: consecutive. Run [3,4]: consecutive. Run [2]: singleton.
The CONTENT of each run is a consecutive block of integers.

So a permutation with s successions = an ORDERED SEQUENCE of r = n-s
consecutive-integer blocks that partition [n].

How many such structures are there?

Given r blocks: we need to partition [n] into r INTERVALS (consecutive sets)
and then arrange these intervals in some order.
The number of ways to partition [n] into r consecutive-integer intervals:
C(n-1, r-1) (choose r-1 "cut points" among n-1 gaps).
Then arrange the r blocks: r! orderings. But we need to exclude
orderings that create NEW successions at the block boundaries.
Wait: any ordering of the blocks is fine — successions only happen
WITHIN blocks, not at block boundaries (since the blocks are chosen
to be the maximal ascending consecutive runs).

Hmm, actually: if two blocks are adjacent in the permutation and
the last element of the first block + 1 = the first element of the next block,
that creates an additional succession. But we defined the blocks as MAXIMAL
runs, so this can't happen — it would mean the two blocks should be merged.

OK so the counting is: a permutation with exactly s successions
corresponds to a COMPOSITION of n into (n-s) parts (consecutive blocks)
where NO two adjacent blocks in the permutation "merge" (i.e., no accidental
successions at boundaries). The blocks are intervals of [n], and their
ordering avoids creating new successions.

This is related to the DERANGEMENT-LIKE condition on block orderings.

Let me just compute directly and see the structure.
"""

import sys, math
from itertools import permutations
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')

def successions(pi):
    return sum(1 for i in range(len(pi)-1) if pi[i+1] == pi[i] + 1)

def anti_successions(pi):
    return sum(1 for i in range(len(pi)-1) if pi[i+1] == pi[i] - 1)

def get_runs(pi):
    """Get the ascending consecutive runs of pi."""
    runs = []
    current = [pi[0]]
    for i in range(1, len(pi)):
        if pi[i] == pi[i-1] + 1:
            current.append(pi[i])
        else:
            runs.append(tuple(current))
            current = [pi[i]]
    runs.append(tuple(current))
    return tuple(runs)

def main():
    print("=" * 70)
    print("MOBIUS CRACK — POSET AND BLOCK PERSPECTIVE")
    print("=" * 70)

    # ============================================================
    print(f"\n--- STEP 1: Understand D_n(2) through BLOCKS")

    n = 5
    all_perms = list(permutations(range(n)))

    print(f"\n  n={n}: Anti-succ-free permutations with their runs:")
    run_data = []
    for pi in all_perms:
        if anti_successions(pi) > 0:
            continue
        s = successions(pi)
        runs = get_runs(pi)
        run_sizes = tuple(len(r) for r in runs)
        run_data.append((pi, s, runs, run_sizes))

    # Group by run size partition
    by_sizes = Counter()
    for pi, s, runs, sizes in run_data:
        by_sizes[sizes] += 1

    print(f"  Anti-succ-free permutations grouped by run sizes:")
    for sizes in sorted(by_sizes.keys()):
        print(f"    Sizes {sizes}: {by_sizes[sizes]} perms, s={n-len(sizes)}, "
              f"2^s={2**(n-len(sizes))}, contribution={by_sizes[sizes]*2**(n-len(sizes))}")

    total = sum(c * 2**(n-len(s)) for s, c in by_sizes.items())
    print(f"  Total D_n(2) = {total}")
    print(f"  Expected: {158 if n==5 else '?'}")

    # ============================================================
    print(f"\n--- STEP 2: WHAT IS an anti-succ-free permutation in block language?")

    print(f"""
  An anti-succession-free permutation has NO descending adjacency.
  In block language: two adjacent blocks [a,...,b] and [c,...,d]
  satisfy: b + 1 != d (no ascending merger) AND c != b - 1 (wait...)

  Actually the anti-succession condition is on CONSECUTIVE elements
  of the permutation, which are the BOUNDARY between adjacent blocks.

  Last element of block i = max(block_i) (since ascending within block).
  First element of block i+1 = min(block_(i+1)).

  Anti-succession at boundary: first(block_(i+1)) = last(block_i) - 1
  = max(block_i) - 1.

  So anti-succ-free = for all adjacent blocks, min(next) != max(current) - 1.

  Example blocks [3,4] followed by [2]: min(next)=2, max(current)=4.
  2 != 4-1=3. OK, no anti-succession.

  Blocks [3,4] followed by [3]: impossible (3 appears twice).
  Blocks [2,3] followed by [1]: min(next)=1, max(current)=3. 1 != 2. OK.
  Blocks [1,2] followed by [0]: min(next)=0, max(current)=2. 0 != 1. OK.
  Blocks [1,2] followed by [1]: impossible.

  Wait: [2,3] followed by [1,2]: impossible (2 appears twice).
  The blocks PARTITION [n], so they're disjoint.

  The condition is just: min(block_(i+1)) != max(block_i) - 1.
  Equivalently: the blocks, when arranged in the permutation order,
  don't have "descending adjacency" at their boundaries.""")

    # ============================================================
    print(f"\n--- STEP 3: MOBIUS INVERSION on the block partition lattice")

    print(f"""
  Let f(B) = number of ways to ORDER a set of blocks B into a permutation
  with NO anti-successions at block boundaries.
  Then D_n(2) = sum over compositions C of [n] into consecutive blocks:
                f(C) * 2^(n - number_of_blocks(C))
              = sum over compositions: f(C) * 2^(successions in C)

  The successions are WITHIN-block, so 2^s = 2^(n - num_blocks).

  Now: the compositions of [n] into consecutive blocks correspond to
  subsets of (0,...,n-2) (the "cut points"). A cut at position i means
  elements i and i+1 are in DIFFERENT blocks. The number of blocks =
  number of cuts + 1.

  D_n(2) = sum_S f(blocks defined by cuts S) * 2^(n-1-|S|)

  where S ranges over subsets of (0,...,n-2) and 2^(n-1-|S|) = 2^successions.

  THIS IS A SUM OVER THE BOOLEAN LATTICE 2^(0,...,n-2).
  And f depends on the block structure.

  MOBIUS INVERSION on the Boolean lattice could relate f to a simpler function.""")

    # ============================================================
    print(f"\n--- STEP 4: WHAT IF we drop the anti-succession condition?")

    # Without anti-succ condition:
    # sum over compositions * 2^(n-1-|S|) * (number of orderings of blocks)
    # = sum_S (n-1-|S|+1)! * ... hmm, not quite.
    # Wait: given cuts S with |S|=c (so c+1 blocks), the blocks can be
    # arranged in (c+1)! ways (any order). Each arrangement is valid
    # (no restriction on block boundaries).
    # So unrestricted sum = sum_{c=0}^{n-1} C(n-1,c) * (c+1)! * 2^(n-1-c)
    # Let me verify.

    for n in [3, 4, 5]:
        unrestricted = sum(math.comb(n-1, c) * math.factorial(c+1) * 2**(n-1-c)
                          for c in range(n))
        actual = sum(2**successions(pi) for pi in permutations(range(n)))
        print(f"  n={n}: unrestricted formula = {unrestricted}, actual = {actual}, "
              f"match = {unrestricted == actual}")

    # Hmm, this won't match because the block orderings aren't all valid
    # (some orderings would merge blocks at boundaries, creating additional
    # successions). Let me reconsider.

    # Actually: a composition into consecutive-integer blocks with a
    # specific ordering IS a permutation. The successions are EXACTLY
    # the non-cut positions. So every permutation corresponds to a
    # UNIQUE set of cuts S and ordering of the resulting blocks.
    # BUT: different orderings of the same blocks CAN create additional
    # successions (if blocks happen to be adjacent in value).

    # This is the subtlety. The block decomposition is based on the
    # MAXIMAL ascending runs. Reordering blocks can create new runs.

    # ============================================================
    print(f"\n--- STEP 5: DIRECT approach — Mobius function on the INTERVAL poset")

    print(f"""
  NEW IDEA: Instead of permutations, think about WORDS.

  D_n(2) sums over anti-succ-free perms weighted by 2^s.
  2^s = number of subsets of succession positions.
  Choosing a subset S of succession positions = choosing which
  pairs to "glue" into blocks.

  So D_n(2) = |(pi, S) : pi anti-succ-free, S subset successions(pi) |.

  Choosing S FIRST: S defines a set of forced blocks.
  Then pi must be an ordering of blocks such that:
  - Within S-blocks: ascending consecutive (forced by S)
  - No anti-successions at block boundaries
  - The non-S positions can go either way EXCEPT descending

  Wait: if S is a subset of successions(pi), then the non-S positions
  are NOT successions. But they might or might not be anti-successions.
  The anti-succ-free condition says NONE of them can be anti-successions.

  So D_n(2) = |(pi, S) : S = set of succession positions in pi
               chosen as a subset, pi has no anti-successions |
  = |(pi, gluing) : pi anti-succ-free, gluing is a partial gluing
       of pi's ascending runs |

  Hmm. Let me try yet another angle.""")

    # ============================================================
    print(f"\n--- STEP 6: THE TRANSFER MATRIX")

    print(f"""
  Think of building the permutation LEFT TO RIGHT, one element at a time.
  At each step, we choose the next element from the remaining ones.

  State: (last_element, set_of_remaining_elements).
  Transition: choose next element != last - 1 (anti-succ-free).

  If next = last + 1: this is a succession, weight = 2.
  Otherwise: weight = 1.

  D_n(2) = sum over all valid paths through this weighted automaton.

  The TRANSFER MATRIX approach: the state space is [n] x 2^[n].
  But 2^[n] is exponential. For n=5: 5 * 32 = 160 states. Manageable.

  Actually: the anti-succ-free condition means we can't go from
  element i to element i-1. This is like a RESTRICTED GROWTH process.

  The generating function approach: let f(x, t) count the weighted
  number of anti-succ-free permutations by length n and weight parameter x.

  f(x, t) = sum_n D_n(x) * t^n / n!

  The EGF should satisfy a differential or functional equation
  encoding the transfer matrix structure.""")

    # Let me just compute D_n(2) for larger n using the transfer matrix idea
    # But actually, the formula D_n(2) = n! + 2*sum_k (n-2k)^k*(n-2k)!
    # already computes D_n(2) efficiently.

    # Let me instead try to UNDERSTAND the formula term by term.

    # ============================================================
    print(f"\n--- STEP 7: WHAT DOES 2*(n-2k)^k*(n-2k)! COUNT?")

    # 2*(n-2k)^k*(n-2k)! = 2 * (n-2k)^k * (n-2k)!
    # = 2 * |{{ (sigma, f) : sigma in S_{n-2k}, f: [k] -> [n-2k] }}|
    # This counts PAIRS of:
    # - A permutation sigma of (n-2k) elements
    # - A function f from [k] to [n-2k] (k labels from n-2k items, with replacement)
    # ...times 2.

    # Can I map these to pairs (pi, S) where pi is anti-succ-free
    # with |S| = some specific value?

    # If 2k "elements" are removed from [n] to form k "blocks of size 2",
    # the remaining (n-2k) elements form a permutation sigma of [n-2k].
    # The k blocks are "inserted" back, each labeled by a function value f(j)
    # indicating WHERE in sigma the block goes.

    # THIS IS LIKE INSERTING k BLOCKS INTO A PERMUTATION OF n-2k ELEMENTS.
    # Each insertion creates a succession (within the block).
    # The function f says AFTER WHICH element of sigma to insert the block.

    # There are (n-2k) positions to insert (after each element),
    # plus position 0 (before all) = (n-2k+1) positions? Or (n-2k) positions?

    # With replacement: (n-2k)^k choices.
    # sigma: (n-2k)! choices.
    # Factor 2: direction (insert ascending or... hmm, blocks are always ascending).

    print(f"""
  2*(n-2k)^k*(n-2k)! could count:

  INTERPRETATION: Insert k ascending-pair blocks into a permutation of (n-2k) items.
  - Choose a permutation sigma of [n-2k]: (n-2k)! ways.
  - Choose where to insert each of k blocks: (n-2k) positions each (with replacement).
  - Factor 2: ... one possibility per orientation? Or from path reversal?

  If this interpretation is correct, then each valid insertion produces
  an anti-succ-free permutation of [n] with exactly k additional successions
  (one per inserted block).

  D_n(2) - n! = the EXCESS of weighted count over the unweighted count.
  = 2*sum_k (n-2k)^k*(n-2k)!
  = sum over k of "permutations with k inserted ascending pairs"

  And D_n(2) = n! + excess = baseline + insertions.

  WHERE DOES n! COME FROM?
  n! = (weighted count when all permutations have weight 1)?
  No: D_n(1) != n!.
  So n! is not the x=1 evaluation. It's something else.

  n! appears because it's n!*[D_n(2)/n!] and D_n(2)/n! = 1 + Var/Mean^2.
  The "1" in "1 + Var/Mean^2" is the Mean^2/Mean^2 = 1 contribution.
  And Mean^2 = (n!)^2/2^(2(n-1)) = E_0.
  So E[H^2] = E_0 * D_n(2)/n! = E_0 * (1 + Var/Mean^2).

  The n! corresponds to the E_0 contribution (the "DC" term, level 0).
  The excess corresponds to the Var contribution (levels 2, 4, ...).

  SO: n! in the formula = the level-0 (mean squared) contribution.
  And 2*(n-2k)^k*(n-2k)! = the level-2k contribution.

  This is just the grand formula rewritten! Each level-2k term
  in the Fourier decomposition = 2*(n-2k)^k*(n-2k)! in the
  permutation pair counting.

  THE PROOF WOULD BE COMPLETE if we can show:
  "The number of (pi, S) pairs where pi is anti-succ-free and
   S is a 2k-subset of successions = 2*(n-2k)^k*(n-2k)!."

  This is a BIJECTIVE question: construct a bijection between
  these pairs and the set of (sigma, f, sign) triples.""")

    # Let me verify the count of (pi, S) pairs at each level
    print(f"\n  Verify: |(pi, S) with exactly |S| successions chosen |")
    for n in [3, 4, 5, 6]:
        all_perms = list(permutations(range(n)))
        for k in range(1, n//2 + 1):
            if n - 2*k <= 0: break
            # Count pairs (pi, S) with pi anti-succ-free and |S| = k successions chosen
            # For each pi with s successions, there are C(s, k) such subsets S.
            # Wait: the formula uses 2k not k. Let me reconsider.
            # Actually: level 2k means we choose 2k arcs.
            # But successions correspond to arcs. So maybe |S| should be 2k?
            # Hmm, no. Let me be precise.
            # D_n(2) = sum_{pi: o=0} 2^{s(pi)} = sum_{pi:o=0} sum_{S subset succ(pi)} 1
            # = sum_s sum_{pi with s succs, o=0} 2^s
            # = sum_{(pi, any subset of successions)}
            # The excess D_n(2) - D_n(1) = sum_{(pi, NONEMPTY subset of successions)}
            pass

        # For each k (succession count, not level), count pairs
        level_counts = Counter()
        for pi in all_perms:
            if anti_successions(pi) > 0: continue
            s = successions(pi)
            for j in range(s+1):  # choose j successions
                level_counts[j] += math.comb(s, j)

        pred = {0: sum(1 for pi in all_perms if anti_successions(pi)==0)}
        for j in range(1, n):
            if n-2*j > 0:
                pred[j] = 2*(n-2*j)**j * math.factorial(n-2*j)
            else:
                pred[j] = 0

        print(f"\n  n={n}: Pairs (pi, S) by |S|=j:")
        for j in sorted(level_counts.keys()):
            p = pred.get(j, "?")
            print(f"    |S|={j}: count={level_counts[j]}, "
                  f"2*(n-2j)^j*(n-2j)! = {p if isinstance(p,int) else '?'}, "
                  f"match = {level_counts[j] == p if isinstance(p,int) else '?'}")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
