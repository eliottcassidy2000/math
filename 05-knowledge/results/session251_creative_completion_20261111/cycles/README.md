# Intersecting-cycle inflation: exact Mersenne arithmetic and its boundary

2026-11-11. **PROVED scoped construction; FINITE-EXACT five-vertex-template
classification and controls; full Chan–Pak spectrum remains OPEN.**

## Result first

Safe transitive-module inflation admits a surprisingly simple exact formula
in the coordinates `x_i=2^{m_i}-1`. The needed observer is **all** directed
cycles of the core, including even cycles, not just its odd-cycle count.
This gives real positive-coefficient arithmetic and explicit non-one-apex
families. It cannot prove the full conjecture by itself: the inherited
intercyclic spectrum already has more holes than `{3,10}`.

### Inheritance and concept board

- Closest inherited mechanism: `H(T)=I(Omega(T),2)`, with Omega the graph
  of literal directed odd cycles, and the deletion identity in
  `01-canon/definitions.md` / THM-002.
- Decisive inherited boundary: THM-4940,
  `01-canon/theorems/THM-4940-odd-cycle-packing-counts-in-every-base.md`,
  records the intercyclic cycle-count spectrum A and all-order gaps
  `3,10,14,17,20,21,...`. Thus the hoped-for subclass universality is
  already false: H29,H35,H41,H43 are among its additional holes.
  This result is inherited, not re-audited here. Its linked source directory
  is absent from this checkout; the theorem itself is present.
- Corrected near miss: uniform module contact alone does not prevent
  reentry (decomposition16 attachment proof and MISTAKE-756). Our new
  sufficient hypothesis eliminates reentry by a disjoint-cycle argument.
- Least-used sidecar: **even directed core cycles**. Odd parity after
  expanding a vertex can come from an even core cycle.
- Live board: odd-cycle packing / triangle intersection / transitive
  substitution / Mersenne coordinates / module-return obstruction.

## 1. Exact admissibility criterion (PROVED)

Let R be any finite tournament. Replace vertex i by a transitive tournament
of positive order m_i, with all intermodule edges prescribed by R; call the
result B. Let U be the set of i with m_i=1.

Then B has no two vertex-disjoint directed cycles if and only if

> Every pair A,D of directed triangle supports of R, including A=D,
> satisfies `A intersect D intersect U != empty`.

In particular the same criterion characterizes pairwise intersection of
all directed odd cycles.

Proof. Every cyclic tournament contains a directed triangle. Thus two
disjoint cycles exist iff two disjoint directed triangles exist. A triangle
of B uses three different modules and projects to a triangle of R. Two
projected triangles can be lifted disjointly iff every shared module has
at least two vertices. This proves both directions, including the case
where the projected supports agree. No Helly assertion is being used.

## 2. No reentry under admissibility (PROVED)

Under this criterion, every directed cycle of B visits every module in at
most one contiguous run, with cyclic runs understood modulo rotation.

Proof. If a cycle has two or more runs through module i, consider two
successive gaps between its runs. The outside vertices of the gaps are
disjoint. Each gap begins at an out-neighbor of module i and ends at an
in-neighbor. Because all vertices of module i have identical external
neighbors, choose a distinct vertex of module i to close each gap into a
directed cycle. The two resulting directed cycles are disjoint. There are
at least two vertices available in module i because the original cycle
already visited it in distinct runs. This contradicts Section 1.

This is precisely the additional hypothesis missing from unrestricted
module-collapse arguments.

## 3. Closed cycle polynomial (PROVED)

Let C range over **all simple directed cycles** of R, counted as cyclic
orders modulo rotation, not as vertex supports. Multiplicities of supports
must be retained. Put `x_i=2^{m_i}-1`. Under Section 1's condition,

    H(B) = 1 + sum_C ( product_{i in C} x_i - (-1)^{length(C)} ).

Proof. Section 2 projects each cycle of B to a simple directed cycle C of
R. Conversely, choose a nonempty subset of each module in C, and traverse
that subset in its unique ascending order. This constructs every lift
exactly once. For a nonempty subset of a transitive m-vertex module, the
ordinary generating count is `2^m-1`, while the even-minus-odd count is
`sum_{q=1}^m (-1)^q binom(m,q)=-1`. The number of odd lifts of C is

    ( product_{i in C} x_i - (-1)^{length(C)} ) / 2.

Admissibility makes the conflict graph complete, so H is one plus twice
the total odd-cycle count. This proves the formula.

**Scope boundary:** the polynomial has positive coefficients in its
nonconstant monomials; its constant term need not be positive. This is not
an unrestricted scalar operation on H. It retains R, the admissible set U,
and the complete directed-cycle polynomial. Reentry destroys the formula
outside the criterion.

For a single expandable module v, let p count all directed core cycles
through v, let e count the even ones through v, and let a count odd cycles
avoiding v. Then

    H(m) = p*2^m + b,    b=1+2a-2e,
    H(m+1) = 2H(m)-b.

These are genuine rooted constructions, not claims that an arbitrary
realiser of H admits the same update.

## 4. Complete five-core calculation (FINITE-EXACT, all sizes consequence)

There are 12 unlabelled five-vertex tournament templates. Exactly three
have no vertex common to all directed triangles. The row-mask codes below
mean bit j of row i is the edge i -> j. Exhaustive relabelling, not an
unverified invariant hash, produces the representatives.

- Template40: rows `[16,9,3,5,14]`. The maximal expandable sets are `{0}`
  and `{4}`. Its polynomial after other sizes are fixed at1 is `6+9x`.
  Therefore every safe expansion of this core has
  `H=9*2^m-3`, m>=1: 15,33,69,141,... .
- Template41: rows `[18,8,3,5,14]`. The maximal expandable set is `{1,4}`.
  Its polynomial is `4+3x+3y+3xy`, hence
  `H=1+3*2^(a+b)`, a,b>=1: 13,25,49,97,... .
- Template76: rows `[24,17,3,6,12]`. No module is expandable. H=15.

These finite core facts combined with the proved polynomial classify
**all sizes** of safe transitive expansions of these three non-one-apex
five-vertex cores. They do not classify all intercyclic tournaments.
The two sequences are not asserted new to the repository; the useful
output is their common admissibility and even-cycle mechanism.

For completeness, safe expansion of a one-apex core remains one-apex.
If the common triangle vertex v stays singleton this is immediate. If v
expands, triangles through v correspond to edges of a bipartite graph
between v's out- and in-neighbors. Admissibility says these edges are
pairwise intersecting away from v at singleton vertices. Pairwise
intersecting edges in a bipartite graph form a star (the only other
pairwise-intersecting simple-edge pattern is a triangle). Its center is a
singleton common triangle vertex. The single-edge case follows directly
from the same admissibility condition. Thus these three templates really
capture the non-one-apex part of the entire five-template blowup grammar.

## 5. Exact experiment and reproduction

Run from repository root:

    python3 -B 05-knowledge/results/session251_creative_completion_20261111/cycles/explore.py

The script enumerates all1024 labelled five-vertex tournaments and
quotients by all120 permutations, obtaining12 cores. It checks the
admissibility criterion independently against actual triangle intersection
for all384 core/size choices with sizes1..2. It then counts Hamiltonian
paths by subset DP for every admissible size vector with sizes1..4 and
total order<=11: **2331 cases**. Every result agrees with the new cycle
polynomial. Non-one-apex witnesses also have their odd directed cycles
counted by a separate literal directed-cycle DFS, checking H=1+2c.
`results.json` stores the exact universe, spectrum and row-mask witnesses.

A hostile control uses the cyclic triangle with all modules of size2:
its modules admit reentry, and the unsupported polynomial gives a different
answer from literal Hamiltonian DP. The exact values are in results.json.

## 6. What changes for the main proof

The pleasant Mersenne-coordinate polynomial is a valid positive arithmetic
compiler, but the whole intercyclic carrier is too small. The next useful
bridge is not more searches for all k in this carrier. It is to let one
controlled pair of disjoint triangles appear and retain its packing
sidecar. Then `H=1+2c_1+4c_2` (when packing number<=2), so the new term
may fill an inherited intercyclic hole. Our criterion identifies the exact
boundary where that extra sidecar first becomes necessary. A promising
next construction would cross one admissibility wall while controlling
c_2; this session does not claim such an operation or a completed proof.
