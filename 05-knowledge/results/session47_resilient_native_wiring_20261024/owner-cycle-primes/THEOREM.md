# Owner-cycle work, prime characters, and coordinate-dependent coincidences

Status: PROVED elementary identities; FINITE-EXACT diagnostics on previously
certified native records. This note does not reauthenticate their Hamiltonian
paths or normalized resource capacities.

## Inheritance and question

The Session46 sharp 18-source circuit has seventeen two-owner transpositions
and one three-owner cycle, hence cost37. The current complete-core witness
for a 22-clock obstruction instead changes eight sources. The user suggested
comparing this with earlier {2,3,11} structures. The useful question is which
numbers belong to an invariant, and which are coordinates of its presentation.

Sources: Session46 `pinned-domains/normal.json` and incoming
`compatible_charged_helpers_20261023/supply/selected_supply.json`. Exact hashes
are recorded in the diagnostic receipt. The core means precisely the
28 authenticated binary U/V menus with the target-margin exclusions retained.

## 1. Actual owner-cycle formula

For two actual path orders p and q on the SAME named source, define the owner
permutation f by f(p_j)=q_j. Let its nontrivial cycles have lengths
l_1,...,l_c. Then the changed-owner-clock cost is

    w = sum l_i,
    sign(f) = (-1)^(w-c).

This follows because a nontrivial cycle changes every one of its indexed
owners and has sign (-1)^(l_i-1). Fixed points contribute neither work nor sign.
For independent named slots multiply these signs; the total cost counts
owner-clock occurrences, not unique globally identified owner names.

Every one of the original28 U/V record changes has exactly ONE nontrivial
cycle. A fixed-target-margin transition changes an even number of binary
slots. Therefore for every such simultaneous core event,

    product of actual owner-permutation signs = (-1)^(total clock cost).

This is a genuine even/odd correspondence in the given record carrier.
It does not make an arbitrary ambient permutation into a legal path operation.
Nor does it identify a physical fermionic exchange statistic.

The 22-clock pair changes sources

    A={1532,1691,1923,1974},
    B={1841,2176,2837,3360}.

Its cycle lengths comprise four2s, two3s and two4s, giving

    22 = 4*2 + 2*3 + 2*4.

This arithmetic is recovered from literal owner permutations, not from a
factorization conjecture. Completeness of the two-state fibre and the uniform
22 upper bound require the separate weighted-fibre certificate. This note
checks the eight local permutations and their cost, not those global claims.

## 2. The correspondence changes when the menu changes

The authenticated53 new helpers include changes with several nontrivial
cycles. Source650 has lengths(2,4), cost6 and positive sign. A single-cycle
formula would incorrectly assign negative sign. Thus the universally valid
formula retains the additional cycle-count coordinate c.

For a combined event the exact statement is sign=(-1)^(W-C). Fixed margins
ensure an even number of changed slots, not an even total C. Losing C when
passing from the old menu to an enlarged menu loses a real invariant.
The complete helper cycle-type histogram is in the receipt.

## 3. The 11 in a parity label is not intrinsic

The four binary charge coordinates can be replaced by any invertible linear
basis change over F2. In packed little-endian notation the explicit involution

    (x0,x1,x2,x3) -> (x0,x1+x0,x2,x3+x0)

takes packed label11 to packed label1. It preserves all charge equalities,
all guarded fibres, every native record, and every clock cost. Consequently
the coordinate11 in the 22-clock witness charge (123,11) cannot itself explain
an invariant prime11 mechanism. The number22 as an actual clock cost survives
this change of coordinates, while that displayed11 does not.

## 4. No odd-primary additive invariant factors through this dyadic code

Let Q=(Z/128) direct-sum (Z/2)^4. For any finite odd-order abelian group H,

    Hom(Q,H)=0.

Indeed every image is killed by128, and multiplication by128 is an automorphism
of H. In particular no nontrivial additive mod3 or mod11 invariant can be
recovered from Q alone. This is an obstruction to a particular proposed bridge,
not an objection to studying3 or11 in the full tournament or integer-current
carrier. An odd-primary character of the full critical group, integer source
work, or another explicitly retained representation could supply information
that this quotient discards. Those would require a new map and certificate.

Earlier repo appearances of {2,3,11}, such as the finite-valuation Collatz
potential obstruction THM-4507 (full slug below), involve selected arithmetic
valuation coordinates. The present Q is a finite dyadic current quotient.
Neither the shared numerals nor a common word 'charge' identifies their fibres.
The concrete shared research lesson is to specify exactly what a projection
forgets and test whether the intended observable is constant on its fibres.
This is a methodological analogy, not a reduction between the conjectures.

## Reproduction

Run from a workspace with a repository checkout:

    python3 check.py --repo /path/to/math --output normal.json
    python3 -O check.py --repo /path/to/math --output optimized.json

The checker can read tracked input blobs through git if the checkout is sparse.
It checks every actual local cycle/sign identity in both inventories, the
22-cost decomposition, all256 parity-additivity comparisons, and the elementary
odd-order generator obstruction for targets3,11,33. Proofs handle all odd H.
Inputs are inherited, and their exact hashes are recorded.

Historical reference, not a dependency on its Collatz proof:
`01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md`.
