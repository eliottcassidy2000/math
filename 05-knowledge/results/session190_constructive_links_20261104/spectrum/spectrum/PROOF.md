# A four-arc flat-response patch with an all-clock nonovershooting lift

**PROVED construction; FINITE-EXACT sharp local bounds.** This packet constructs
an actual +2 Hamilton-count change on infinitely many strong tournaments by
changing four specified old arcs. Every intermediate lifted tournament is
strong, and each intermediate lifted count is strictly smaller than the
starting count. It does not give arbitrary-target selection, repeated +2
closure, a ternary digit grammar, or the assertion that 7 and 21 are the only
gaps in the global Hamilton-path spectrum.

## 1. Inheritance and the interface being changed

S188's [fixed-quotient module-substitution obstruction](https://github.com/eliottcassidy2000/math/blob/dcbe43e7fc9caf175baba4fef7b402fd9e3c51d9/05-knowledge/results/session188_structural_bridges_20261104/returns/PROOF.md)
and S189's [four-core adapter and exact spacing theorem](https://github.com/eliottcassidy2000/math/blob/dcbe43e7fc9caf175baba4fef7b402fd9e3c51d9/05-knowledge/results/session189_response_frontier_20261104/returns/PROOF.md)
leave a precise constructive demand: alter incidences rather than merely the
marked port or dyadic clock. S189's fixed component, with all three restored
unit-start marks and arbitrary positive clocks, has minimum positive spacing
at least 220. That theorem concerns one specified family, not every pair of
components.

Here the old seven names and the physical root are retained. We change the
component itself. The relevant preserved coordinate is one named two-bin
count, not its whole response. The boundary-data perspective is inspired by
reducible local patches: first certify the interface exactly, then substitute
it into an exterior for which that interface is complete. No planarity or
hypohamiltonian theorem is imported into tournaments.

Live concepts: literal incidence edits; small start budgets; equal named
cover response; dyadic exterior; intermediate count; exact ternary inverse.
The canonical hostile is the inference that equal F implies full response
equality. The explicit endpoint/deletion distinctions below refute it.

## 2. The marked response and the exterior identity

For a finite tournament A with root r, put

- h = H(A), the number of directed Hamilton vertex words;
- s = s_r(A), the number of these starting at r;
- e = e_r(A), the number ending at r;
- q = H(A-r), with H(empty)=1;
- F = E_r(A,2), the number of ordered spanning two-bin path covers whose
  first bin is nonempty and ends at r, and whose second bin may be empty.

For k>=1, form W_k(A,r) by adding a transitive block z_1->...->z_k and a
vertex w, with

    w -> Z -> A,       r -> w -> A-r.

All old names and arcs remain. Then

    H(W_k(A,r)) = h + 2^k F.                         (1)

Here is a direct proof, so the lift is not inferred from scalar samples.
A path beginning with w must traverse all of Z and then A: there are h.
If w is not first, its predecessor must be r. Remove w and Z. The remaining
A-vertices form a first nonempty path ending at r and a second possibly empty
path. Every z_i independently lies before that first A-path, or immediately
after w before the second A-path. Its order within its chosen block is forced.
The two choices for each of the k names give exactly 2^k F paths. A nonempty
second block of Z after w and an empty second A-path are allowed. Conversely
these choices always concatenate to a legal Hamilton path, and are recovered
uniquely. No third Z block is possible because A has no arcs to Z.
The hypothesis k>=1 matters: for k=0 the count is h+F-s, since a path
starting at w cannot next visit r. No zero-clock conclusion is asserted.

If A is nontrivial strong then W_k is strong: w reaches Z and A, each A
vertex reaches r and then w, and Z reaches A. Thus strongness persists under
any old-arc patch whose intermediate A remains strong. The outer w has only
one in-neighbor r; W_k-w is Z->A, so s_r(W_k-w)=0. These exact outer incidence
and deletion facts hold at every patch stage, regardless of the inner s_r.

## 3. The four physical edits

Outgoing rows are bitmasks on the fixed names 0,...,6. The root is **r=2**
at every step. Start with

    A0 = (6,28,8,1,109,15,47).

Reverse the unordered pairs, in this exact order,

    {1,2}, {0,1}, {2,3}, {3,5}.

The five actual relations and their exact responses are

| i | outgoing rows | H | s_2 | E_2(A_i,2) |
|---|---|---:|---:|---:|
| 0 | (6,28,8,1,109,15,47) | 29 | 1 | 157 |
| 1 | (6,24,10,1,109,15,47) | 41 | 2 | 134 |
| 2 | (4,25,10,1,109,15,47) | 29 | 1 | 89 |
| 3 | (4,25,2,5,109,15,47) | 25 | 1 | 121 |
| 4 | (4,25,2,37,109,7,47) | 31 | 1 | 157 |

Each relation is strong. `verify.py` checks this by reachability and checks
all h,s,e,q,F coordinates independently by permutations/cuts and by a
predecessor subset recurrence. The finite seed is therefore a fully specified
certificate, not an existence claim inferred from numerical signatures.

For x=2^k>=2, the five lifted values are exactly

    29+157x,  41+134x,  29+89x,  25+121x,  31+157x.    (2)

The last minus the first is 2. Each internal value is below the first:

    (29+157x)-(41+134x) = 23x-12 >0,
    (29+157x)-(29+89x)  = 68x >0,
    (29+157x)-(25+121x) = 36x+4 >0.

**All-clock patch theorem.** For every k>=1, W_k(A0,2) and W_k(A4,2) are
strong tournaments on k+8 vertices with counts 29+157*2^k and 31+157*2^k.
The displayed four old-arc reversals connect them through strong tournaments,
without changing any exterior arc or vertex name, and every internal count
is strictly below the lower endpoint. Their first four pairs are
(343,345), (657,659), (1285,1287), (2541,2543).

This is a literal arc-reversal path. It is not, without further guards and
compilation, a claim about the repository's narrower native complement
program language. Nothing in the theorem says that the final component can
be reset into another copy of the initial one at unchanged H.

### Precisely what is and is not preserved

At the two endpoints, the seven-coordinate responses
(h,s,e,q,F,S_r(A,2),P_A(2)) are respectively

    A0: (29,1,10,15,157,34,580),
    A4: (31,1, 9, 9,157,31,644).

Here S requires the first bin to be nonempty and start at r, while P
counts all ordered two-bin covers, allowing either bin to be empty.
Both are independently checked by permutations/cuts and subset recurrence.

Thus only the two stated coordinates s=1 and F=157 coincide. H, endpoint
count, deletion count, start-bin count and total two-bin count all differ.
Neither rooted isomorphism nor full named-response equality is asserted.
Equation (1) explains why the deliberately narrower interface is sufficient
for this one exterior. A general attachment can distinguish the endpoints
by coordinates that the wrapper does not read.

## 4. Sharp finite size and edit statements

**Small-component obstruction.** Among all nontrivial strong rooted
tournaments of orders at most six, no two (even of different orders) have
s=1, equal F and Hamilton counts differing by 2. `certify.cpp small` exhausts
all 2^(n choose2) labelled tournaments at each n=3,...,6 and retains all
rooted (H,F) profiles with s=1. Thus seven is the minimum possible maximum
component order for this exact equal-F/unit-start twin criterion.

**Sharp rooted edit distance for these endpoints.** Among all 6!=720
identifications of A4's other names with A0's names while fixing r=2, the
minimum number of different arcs is four; exactly five identifications
attain it. No converse identification is allowed. For each of those five,
there are exactly 4! possible shortest arc-flip orders. Of the resulting
120 paths, 38 have all five states strong. None keeps s_r=1 at every step;
the minimum possible maximum value of s_r is two, attained by section 3.

These are complete finite claims for the displayed endpoints and shortest
paths. They do not assert that a longer path cannot remain in the unit-start
class, or that four is minimal after changing the seven-vertex endpoints,
allowing extra vertices, moving the marked root, or changing exterior arcs.
The complete 5*16=80 intermediate relations and counts, including weak
states and zero starting counts, are retained in `receipt.json`.

## 5. Two precise stopping boundaries

### 5.1 A sibling pair is not automatically a ternary return

The lower endpoint of the seven-name family has residues 0 or 1 modulo3,
because 29+157*2^k = 2+(-1)^k modulo3. It can never be 3H+2.
Thus this uniform +2 patch does not itself supply the ordered pair
(3H+2,3H+4) for any integer H.

There is also a universal parity restriction worth retaining. Adding a
vertex beaten only by r gives H_new=h+F-s. Both Hamilton counts are odd,
so F=s modulo2. For two unit-start roots in the same component, or even two
H-equal components, equal positive outer clocks give a difference divisible
by 2^(k+1), not two. At unequal clocks, a difference of two requires that
the smaller clock is one and the two odd responses satisfy

    |2^d F_highclock - F_lowclock|=1, d>=1,

where d is the difference between the clocks.

The present construction instead changes the component Hamilton count by
2 while preserving F, exactly the missing coordinate in that obstruction.
The parity argument is elementary inherited calculus, not claimed as a new
parity theorem.

### 5.2 Two arithmetically eligible flat-response pairs have no fixed parent

The exploration also found the following explicit strong unit-start pairs:

- roots 3, counts 575 and577, common F=3819, rows in `receipt.json`;
- roots 2, counts 605 and607, common F=15909, rows in `receipt.json`.

For a fixed core (P,r) and the same clock k on all three wrappers, the
identities

    H(W_k(A))=3H(W_k(P))+2,
    H(W_k(B))=3H(W_k(P))+4

for every k force, by comparing the constant and 2^k coefficients,

    (H(P),E_r(P,2)) = (191,1273), or (201,5303),       (3)

respectively. **Neither pair in (3) is realized by any strong rooted
tournament.** This is a finite-exact exclusion with an all-order cutoff,
not an empty bounded search extrapolated to infinity.

The cutoff needs only Busch's published lower bound H>=5^((n-1)/3).
At n>=11 it gives H>201, because 5^10=9,765,625>201^3=8,120,601.
No exact residue-dependent floor value is needed. Smaller orders are
handled as follows:

1. The independently audited complete 6,880-class order-eight bank is pinned
   by SHA-256 in `inherited_bank8_receipt.json`. All its one-vertex deletions
   cover every order-seven class. Their maximum H is189, excluding both
   targets at orders<=7 by Hamilton-path insertion monotonicity.
2. Scan all order-eight roots, and every one of the 6,880*256 one-vertex
   extensions to order nine. No target profile occurs. Keep every strong
   order-nine child with H<=201, and quotient only by explicitly checked
   degree-preserving relabelings. Exactly545 classes remain.
3. Extend each of those545 classes in all512 ways, and scan every root.
   Neither profile occurs at order ten. Completeness follows from Moon's
   strong-deletion theorem for order at least four and H(P-v)<=H(P).
   The latter inequality has the
   elementary proof: insert v at the first legal place in each old Hamilton
   word; deleting v recovers that word, so the map is injective.

The full target-F value sets and candidate counts are in
`census_receipt.json`; `run_census.py` recreates the complete reduced bank
and searches from the pinned order-eight input. No random search or
uniquely-Hamiltonian-only sample supports this conclusion.

If one also allows a nonstrong parent core while requiring W_k(P,r) strong,
r must lie in P's final strong component. Write P=X->C with r in C. Then
E_r(P,2)=P_X(2)E_r(C,2), since each ordered bin traverses X before C.
If X is nonempty, P_X(2) is even: swapping the two bins acts freely on the
ordered covers. Thus the odd target F in (3) forces X empty, reducing again
to a strong core. This extends the exclusion to every parent core that
could make a strong wrapper. It does not exclude clock-dependent parents,
unequal clocks, or a different exterior architecture.

## 6. Reproduction, sources, and remaining scope

Run from this directory, without installing dependencies:

    python3 verify.py
    python3 -O verify.py
    python3 run_census.py

The first two runs produce byte-identical stdout and retain the complete
independent permutation/DP checks. The census compiles the bundled standard
C++17 verifier and checks all stated universes, with 64-bit arithmetic safely
below (n+1)n! at maximum order ten. The larger direct twins use unbounded
Python integers.

The exact all-order statements in sections 2,3,5 use the elementary lift,
explicit inequalities and the stated classical size/deletion bounds; the
finite checks are not used to guess a recurrence. The unrelated exploratory
cycle-ceiling route produced no theorem and is not a dependency.

Primary source: Arthur H. Busch, *A Note on the Number of Hamiltonian Paths
in Strong Tournaments*, EJC13(2006), N3,
https://doi.org/10.37236/1141. Its [published abstract](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
states the lower bound used here and was checked directly; the full PDF
fetch was unavailable. The strong-deletion input is an immediate consequence
of Theorem1 in J. W. Moon, *On Subtournaments of a Tournament*, Canadian
Mathematical Bulletin9(1966),297-301,
https://doi.org/10.4153/CMB-1966-038-7. The [publisher's original PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5D7CCBE66AB07BB3F90F771D1CD3CD08/S0008439500052875a.pdf/on-subtournaments-of-a-tournament.pdf)
states that every vertex of a strong n-tournament lies on a cycle of each
length3,...,n; an (n-1)-cycle yields the needed strong induced subtournament.
The order-eight completeness receipt was re-fetched at the immutable commit
above and matches the saved bank's SHA-256 exactly.

Novelty against all earlier literature is not claimed. The contribution is
a new explicit, all-clock incidence patch and its exact boundaries relative
to the inherited return interfaces. Global spectrum completeness remains OPEN.
