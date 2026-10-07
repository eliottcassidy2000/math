# Variable depth Collatz rules with grounded children and arithmetic guards

**PROVED:** the displayed guarded reductions and finite-language compilers.
**FINITE-EXACT:** the three supplied first-hit certificates and declared
experiments. **OPEN:** coverage and recursive grounding for every positive
integer. The named obstruction `2^1459-1` is now fully grounded.

## The residual branch now has two deeper rules

The [arithmetic package](collatz_reset2_rules_20261007.md) supplies:

| Original exponent | Smaller Mersenne exponent | Reduced source depth | Exact infinite guard |
|---|---:|---:|---|
| K | K-2 | 12 | K=1459 mod2^21 |
| K | K-6 | 62 | K=1457 mod2^126 |

These are applications of the existing
[THM-4555 collision theorem](../../01-canon/theorems/THM-4555-uniform-switches-are-collisions-at-minus-one-trailing-ones-deletion.md),
using newly recorded cores. The complete source and child words retain
their variable initial runs, exact native guards, carries and common
endpoint. Membership in the exponent phase needs modular arithmetic, not
expansion of the enormous source. The final valuation is an actual,
uncapped parameter rather than a bounded alphabet choice.

Both rules cover infinite families missed by the frozen preceding selector,
with explicit intersections retaining its ternary exclusion conditions.
This comparison is against that named bank, not against every known
Collatz rule or the maximal sieve. The deeper rule follows the first
rule's concrete child instead of stopping at a newly exposed obligation.

The selected child `2^1457-1` also has a frozen first-hit certificate
from a bounded discovery run. Independent replay and two substitutions
give the following completed points:

| Source | Odd steps to 1 | Total valuation |
|---|---:|---:|
| 2^1459-1 | 7347 | 13104 |
| 2^1457-1 | 7347 | 13102 |
| 2^1451-1 | 7347 | 13096 |

Production consumes the saved word, validates it, and transports it; it
does not assume an unspecified child reaches 1. The infinite phase rules
remain conditional on their arbitrary smaller children. Point grounding
and family coverage are different conclusions.

## A complete finite parser replaces the fixed partner list

The [collision dynamic program](collatz_collision_dp_20261007.md) works at
the exact rational anchor `f_w(-1)=N/2^E`, for reduced words whose first
letter is at least 2, with odd N and E=sum(w)-1>=1. Every proper predecessor has

    denominator exponent f<E,
    numerator (N-2^f)/3, with 3 dividing N-2^f,
    final letter E-f.

The one-letter terminal exists exactly at N=-1. Strictly decreasing E
makes this a complete finite grammar. The parser retains a witness for
each attainable word length, so it can choose a deletion that fits the
original source's available run. Increasing the inspected prefix depth
does not require guessing a new partner alphabet.

This is exhaustive for each supplied finite prefix. It is not a theorem
that some prefix always has a useful collision. Letter and cost frontiers
remain explicit; an actual first hit of 1 is exported as a ROOT certificate.
Pure-two prefixes of lengths through 16 provide finite hostile controls,
and their all-length behavior is not inferred from those checks.

## A second unbounded depth comes from ternary inverse cycles

Incoming audit A2 corrected the scope of
[THM-4594](../../01-canon/theorems/THM-4594-the-maximal-class-decided-collatz-sieve-and-the-sign-barrier.md):
its sign barrier concerns unrefined binary classes. Ternary refinement
can certify the negative-cycle residue classes. The earlier synthesis's
broad sign-barrier sentence is repaired with that correction lineage.

The useful operation here is

    G(m)=(8m-5)/9,   U_(1,2)(G(m))=m.

Given a paid child m, r inverse turns are legal exactly when
`9^r | m+5`; each gives a strictly smaller positive odd child.
The new proof retains the whole word `(1,2)^r`, not just a formal rational
power of G. For m=2^(K-2)-1 the maximum is

    r_max=floor((1+v3(K-4))/2).

For every r>=1 and u=1 or 2, CRT gives infinitely many exponents satisfying

    K=1459 mod2^21,
    K=4+u*3^(2r-1) mod3^(2r),   K>=1459.

The binary collision followed by these r ternary turns pays

    h_r=-5+8^r(2^(K-2)+4)/9^r < 2^(K-2)-1 < 2^K-1.

This is a genuinely variable-depth composition on explicit infinite
families. At K=1459 there is exactly one turn, with child `(2n-11)/9`.
There is no assertion that every new h_r belongs to an already grounded
family. Its change of shape is the next recursive-closure obligation.

## What survives from cluster algebras and friezes

The [frieze package](frieze_collatz_reset2_20261007.md) makes two exact
connections and rejects a tempting third.

* Chronological prefix carries give positive marked minors satisfying
  Pluecker relations. Keeping the terminal power of two makes the
  representation lossless; triangulation flips re-express the same word.
* Every positive polygon quiddity of length N>=5 has total valuation
  3N-6 and contracting slope `64(3/8)^N`. It therefore supplies a
  guarded all-height descent progression with an explicit height cut.
* The elementary frieze ear move is not a Collatz substitution: its two
  blocks have opposite endpoint classes modulo 3. Moreover no directly
  interpreted contracting positive quiddity can begin with adjacent ones,
  so this whole construction misses every source `n=7 mod8`, including
  the main Mersenne test case. Increasing polygon size cannot repair it.

Thus total positivity is useful for coordinates and particular descending
families. It does not discharge arithmetic compatibility or terminal
grounding. The appropriate object here is a marked configuration or a
graph of compatible certificates, not a tournament manufactured from ties.

## The four papers become an exact coverage ledger

The [paper transfer package](paper_reset2_transfers_20261007.md) records
the relevant complete pages and mechanisms from all four supplied PDFs:
MUB atlas gluing keeps unresolved cells; parallel repetition distinguishes
conditional information from repeated observations; pinned-distance
partitions retain multiplicity; logarithmic Brunn--Minkowski keeps the
exact intersection separate from convergence of its measure.

The resulting compiler translates a native post-reset tail directly into
an exact exponent class using the invertible finite-precision map

    j=(K-1)/2 -> (z_K-1)/4=3(9^j-1)/8.

Seven inherited debt portals and the new short collision give eight
disjoint exponent classes. Their relative natural density among odd K
is `33699/1048576`; the unassigned complement has **86 explicit dyadic
cells**. The same relative density holds inside the named cap-escape
progression K=1mod1458. This is the ledger of that named portfolio only:
the deeper 1457 rule, ternary cycle routing, and other known guards are
not silently included or excluded by this fraction.

The ledger retains finite height exceptions separately from asymptotic
cell mass. A shrinking residual mass alone cannot establish that the
integer remainder is empty. These paper transfers improve the exact
accounting and source preservation; they supply no imported Collatz theorem.

## Decisive next target and validation

The next target is closure under the **actual child maps**, including
`K -> K-2`, `K -> K-6`, and the ternary image h_r. A complete proof must
either ground every generated child or exhibit a well-founded rule system
covering every original positive source. The remaining cells and finite
height conditions cannot be erased by a positive measure estimate.

All four exact packages have saved outputs, normal/optimized comparisons,
independent peer audits, and positive and hostile controls. Together they
record **331,715 explicit checks**. The arithmetic package also freezes
the selected ROOT word as data and independently replays all three
completed points. The typed APIs reject boolean/float aliases, and ROOT
self-loops are never exported as first-hit evidence.
