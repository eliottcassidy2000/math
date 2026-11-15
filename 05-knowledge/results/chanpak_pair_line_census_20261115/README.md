# Chan–Pak: exact profile-line coverage, then change the reader

**Result:** every odd target from **32,769 to 59,999** is on a physical
`(H,P₂)` line from the enlarged reader bank, using source order at most
10. The larger band 69–59,999 still has 60 absent **reader outputs** out
of 29,966. Neither statement concerns exclusions from the unrestricted
tournament spectrum. Universal coverage remains open.

Here P₂ counts ordered pairs of nonempty directed paths covering the
source. Write Q=P₂+2H. A reader tests whether a *single actual source*
lies on `rH+MQ=N`, not whether its coordinates can be independently
realized or matched in residues.

## Exact coverage

| Reader/source restrictions | Hits in the 13,616 central lines |
|---|---:|
| Triangle reader H+Q, source order exactly 10 | 13,604 |
| Transitive marked hosts of order ≤6 + the old anti-readers, any source order | 13,610 |
| Arbitrary transitive marked hosts + the old anti-readers, any source order | 13,613 |
| Arbitrary marked hosts of order ≤6, sources through order 10 | **13,616** |

[Census proof and replay](PROOF.md) counts all **9,733,056** classes
at order 10, with **1,550,357** different pairs. It extends through
orders 11–15 with a mathematically complete Q cap for N≤59,999,
not a claim of complete class enumeration at those orders. Higher
orders cannot fit the target cap. For the uncapped transitive-host
bank, exact per-order central hits at orders 10 through 15 are
13,605, 5,150, 1,448, 108, 12, 2; no new central values first
appear above order 10. Larger source order leaves less room in
the target for varying the reader.

The original size-six transitive hosts and anti-readers miss
32,825, 34,111, 35,099, 35,395, 35,557 and 36,833. Larger
transitive hosts repair three, but **35,099, 35,395, 35,557**
still miss the entire transitive-host bank at all source orders.
[Literal small marked-host repairs](repairs.out) give physical
witnesses for all six original gaps, including the primes 35,099
and 36,833. `repair.py` recomputes their full tournament H.
In each case, the construction is a marked host B, a source A,
and an apex x. Put all arcs B→A and A→x; the apex beats precisely
the marked vertices of B. The deletion decomposition into paths
B₁ A₁ x B₂ A₂ gives the exact line: x at the end contributes
H(B)H(A); nonterminal x contributes M Q(A); if A₁ is empty and
B₁ ends marked, subtract W H(A), because that last arc into x
is unavailable. This explains why choosing a nontransitive host
can change both coefficients and repair a genuine miss.

With all 2,435 different readers from host order ≤6, the central
coverage holds by itself, and 29,828 of 29,966 lines hit over the
broader band. [Adding all transitive hosts and the old anti-readers](union_test.out)
raises it to **29,906 / 29,966**. The cap is exact also for these new
readers: M≥1, W≤M, r≥1−M. If M=1, then r is odd by parity of every
output and hence r≥1; if M≥2 then `rH+MQ≥H+2(Q−H)≥Q+1`.
So N≤59,999 implies Q≤59,998 and source order ≤15. Adding the
capped high layers does not change the small-host total.

## What might scale

- **Don't pad the perfect finite window.** On H+Q, all odd-order
  sources supply only multiples of 3. Two singleton pads preserve
  parity but already leave 147,343 and 147,349 unhit in the scaled
  band. Holding Q exact under a fixed reader gives at most two
  output residues mod 3 across all orders, so it cannot produce
  three consecutive odds. See [dilation control](dilate.json).
- **Change the reader and record its full continuation.** The
  [marked-join compiler](MARKED_JOIN.md) is an exact four-by-four
  matrix law; it exposes the cross term and the endpoint term
  separately. Two primitive order-6 marked hosts have identical
  continuation data except for the intercept. Their switch
  changes the count, survives arbitrary surrounding marked
  joins, and avoids the common factor 17 on repeated words
  with transitive sources. The step is not a unit at every scale.
- **Change the source, not just the clock.** [Deleting singleton
  SCCs](CORE_LIFT.md) identifies exactly which Q values can be
  paired with a fixed H at any order. Near the naive maximum order,
  the [triangle-budget shelf](TOP_SHELF.md) freezes into finitely
  many affine rays. [Smooth-core and gcd bounds](PAIR_LINE_STRUCTURE.md)
  show why too-small strong components and nonprimitive pairs
  cannot fill prime lines. A separate physical construction now
  matches all fixed unit residues modulo 2ʳ3ˢ5ᵛ in every fixed
  fractional dyadic subband, but remains sparse in absolute counts.
- **Use the second apex to place the baseline.** [One-apex sources
  stay sparse](RANDOM_SOURCE_LINES.md) even after arbitrary ordinal
  recombination. [Two independent apexes](TWO_APEX_COARSE_PLACEMENT.md)
  give a proved constructive approximation to every large target,
  with error `O(N^(1−(1/2)log₄(4/3)))`. The concrete remaining
  task is to make that absolute error small enough for a compatible
  exact repair, and to prove that every required residue has a repair
  from the same continuation state. Neither follows from occupying
  every subband or checking more orders.

The historical pulls were deliberately reproducible and peripheral:
older Collatz carry defects, square-sum contact chains, planar moment
inverses, and lonely-runner strata suggested separating the finite
source defect from the unlimited padding. See [drawn paths](OLD_SAMPLE.txt)
and [the independent random draw](sample_paths.txt). The analogies are
not used as mathematical premises.

## Extended-bank replay

After the census commands in PROOF.md:

    g++ -O3 -std=c++17 host_codes.cpp -o hosts
    ./hosts 7 > host_codes.tsv
    python repair.py pairline10.tsv host_codes.tsv
    python marked_join.py host_codes.tsv
    python verify_primitive_fork.py
    python bank_cover.py pairline10.tsv host_codes.tsv
    cat pairline10.tsv pairline_1_{11,12,13,14,15}.tsv > all_profiles.tsv
    python union_test.py all_profiles.tsv host_codes.tsv all_mask_coverage.tsv
    awk '$1 <= 9' pairline10.tsv | python verify_core_lift.py
    python top_shelf.py < pairline10.tsv

The host catalogue includes order seven for the literal witness search;
all witnesses displayed for the six gaps use only orders at most six.
The union scorer filters every host of order above six.
