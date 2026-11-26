# Thue–Morse trade cubes and exact splicing of projected contact words

2026-11-26. **PROVED** operator and cube statements, **FINITE-EXACT** on
the marked board. No degree-57 neighborhood was completed in this session.

## 1. From scalar fairness to physical cancellation

In [Cooper and Dutle, Greedy Galois Games](https://arxiv.org/abs/1110.1137)
([journal version](https://www.jstor.org/stable/10.4169/amer.math.monthly.120.05.441)),
the next sign opposes the running weighted imbalance in the scalar sum
`sum a_i q^i`. As q approaches one from below the finite prefixes stabilize
to Thue–Morse, up to an overall sign. The paired recursion
`a_(2i+1) = -a_(2i)` factors a paired sum by `1-q` and descends to q²;
[Shallit's follow-up](https://arxiv.org/abs/2605.00194v2) treats the speed of
prefix convergence. Cooper and Dutle distinguish the greedy result from
global power-series approximation, which need not be optimized by the
greedy word.

Two strands from this repository sharpen the distinction. The
[Mahler-source exclusion](../mahler_thue_morse_source_exclusion_20261019/THUE_MORSE_FINITE_SOURCE_EXCLUSION.md)
uses long periodic shadows of a binary substitution and shows that they
do not provide a positive integer source. The
[skew-Hadamard tournament tower](../../../01-canon/theorems/THM-4561-the-sierpinski-skew-hadamard-tower-in-the-littlewood-problem.md)
has row doublings of the form `1 ± z^N`, with a controlled perturbation.
Here we ask when the analogous doubling cancels in the actual incidence
vector, rather than at a scalar evaluation or only at a residue.

## 2. The smallest physical trades really are signed binary cubes

Use the [marked contact board](HADWIGER-NELSON-PHYSICAL-SLOTS.md) at (1,1).
The marker is fixed throughout. A patch is a matching of contacts in the
four families (branch and three center slots). A trade exchanges two
disjoint matchings on exactly the same rows. Because distinct contacts
agree in at most one row, each contact of one side must meet at least
four distinct contacts on the other side. Thus every nontrivial trade
has at least four contacts per side. At equality its colored
intersection graph is K(4,4), with a different perfect matching for
each row family.

The exact enumeration gives **83** size-four trades. To replay it, pick two
compatible contacts on one side; each of the nine derangements of the four
families specifies which slot of the first is matched to which slot of the
second on the opposite side. Linearity makes each proposed opposite
contact unique. The two unused rows in each family determine at most eight
ways to finish the original side. This counts each valid trade 12 times
(two orientations and six starting pairs); there are 996 hits. The support
vectors of the 83 trades are independent over F2.

On all 83, align the sides using the branch matching. The other three
matchings then give the three double transpositions of the Klein group on
four letters. Equivalently we can label the eight contacts by F2³, with
the sides the even and odd parity halves and the four row families
pairing `x` with `x xor d` for `d in {001,010,100,111}`.
For a concrete first trade, the even side is

    (3,4)  (11,55)  (17,32)  (24,27)

and the odd side is

    (3,32) (11,27)  (17,4)   (24,55).

This is an exact contact version of the Thue–Morse signs, not a
metaphorical correspondence. Each physical row sees two opposite signs.
Under binary index 0 to 7 the sign polynomial is
`(1-z)(1-z²)(1-z⁴)`, so degrees 0, 1 and 2 cancel. The indexing is by the
trade's cube coordinates, not the numerical sheet labels.

### The graph of all four-trade moves is a disjoint union of cubes

Check all pairs of sides of distinct trades: 4,856 pairs have a compatible
union, and **every one is disjoint in contacts**. Because each side is
itself a matching, such disjoint sides also have disjoint sets of rows;
replacing one side by its mate preserves its rows.

This proves more than that a greedy sequence of lexicographically
decreasing switches terminates. At a matching M any two enabled moves
commute, and an enabled move stays enabled as any other is performed. A
new move cannot appear: after the switch the inverse is enabled, so any
other enabled move has disjoint rows and was already enabled before it.
Thus the enabled trade set is constant throughout a connected component,
and every subset of its switches gives a different matching. **The whole
component is a hypercube.** Greedy orientation arrives at the same normal
form in any order and gives a canonical representative of that component.

An exact clique search in the disjoint-boundary graph of the 83 trades
gives maximum dimension **six**. The verifier places all 64 vertices of
an extremal cube in the 3,250-vertex graph, one legal edge at a time. Each
has 24 variable contacts plus the retained marker, and source degree 29;
these are partial patches, not complete stars. Number those patches by
the six switch bits. Their even and odd classes have the same sum of
any polynomial of degree at most five in their contact indicators:
each indicator is constant or affine in one switch bit, and an untouched
bit cancels each monomial. Yet the six-way parity distinguishes the two
classes. This extends the prior Hadwiger–Nelson boundary warning from
assumed Boolean summaries to an actually embedded family.

The earlier 50-contact patch after the marker lies in a zero-dimensional
component. Sixty trades touch no reference contact and 23 touch exactly
one, so none has an entire side available at the reference. Moreover all
four-trades preserve the occupied rows; even when enabled they cannot
by themselves fill a missing slot. For an exact-cover solver one may
forbid the larger side of each trade in the chosen lex order. The 83
clauses preserve existence by the canonical normal-form argument.

A stronger exact check seeds every possible new contact whose four rows
are occupied by the reference, and then opens each of its displaced
reference contacts. Every nontrivial trade contains such a seed.
Restrict to the occupied rows, and once an added contact is fixed,
reject reuse of any of its rows. This enumerates each reachable
exchange component with a bounded number of displaced contacts.
Any trade contains a nonempty component no larger than the trade,
so finding none by this method for size at most 12 excludes **every
same-boundary exchange of up to 12 reference contacts**, not only
the four-cubes. It does not exclude moves changing the boundary.

## 3. Exact results on sharing one physical contact word

The closest projected repairs after omitting center 56, 0 or 28 change
respectively 7, 5 or 6 of the 50 reference contacts, with 23, 1 or 1
minimizers. To enumerate, start at the two empty branches and open any
reference branch whose retained slot is displaced. Fix each selected
new edge. If a minimum solution had extra changed branches unreachable
from the openings by these displacements, revert the unvisited changes:
the selected edges do not hit those old slots and the old reference was
a matching. That would contradict minimality. The exact search thus
exhausts all minima by increasing the displacement budget.

The two unique pivot-preserving minimizers agree on 39 contacts. Their
branch-to-pivot union has a single alternating component, on branches

    5, 39, 25, 53, 40, 33, 20, 36, 15, 41, 13, 27, 52.

The cycle has length 26 in the bipartite branch–sheet graph. A perfect
matching using those contacts must keep one whole side or the other.
Each side leaves seven repeated-slot units in its omitted center.
More generally the alternating cycles of two pivot-perfect words give
Boolean choices, and a physical collision gives a unit or binary clause.

Even without a shared perfect pivot projection, the union of any two
words is a Boolean problem: choose one of at most two contacts per branch,
and forbid each conflicting pair. We ran the exact implication test on
**all 13,637,253 pairs** of the 5,223 earlier projected witnesses.
None allows a full physical cover. This is an exclusion of that witness
bank, not of all projected words or of the board.

The union of all 25 nearest words has 215 contacts, and it cannot place 51
of the 52 remaining contacts (the best is 50). Exhaustive enumeration of
inside matchings of size 50, 49 and 48 gives 3, 173 and 5,172 respectively.
Checking all remaining contacts against those matchings proves that a
51-contact patch must import at least four outside choices, and a full
52-contact patch must import at least five. For the latter implication,
if a perfect patch used at most four imports, delete one import (or one
inside choice if there are none). This would give a 51-patch with at most
three imports. Adding the support of any single projected witness to the
215-contact bank also fails; all 5,223 exact-cover subproblems exhaust
within their node limits. These exclusions concern support, not a
lower bound on the number of reference edits in the unrestricted board.

## 4. Why cyclic dyadic cancellation is not a selector

The fourteen-neighborhood construction singles out translation by four
on a 56-sheet fiber. Consider its cyclic shift operator T on an orbit of
length 14, separately from any assumption of equivariance on the marked
chart. The dyadic Thue–Morse filter is

    P_m(T) = product_(j=0)^(m-1) (I - T^(2^j)).

For m >= 2 its kernel over Q consists exactly of the constant and
alternating modes: since T has order 14, the union of the killed
characters has order dividing gcd(14,2^(m-1)) = 2. Every nontrivial
7-channel survives injectively. Over the integers let K be the lattice
with separate even and odd coordinate sums zero. Then exactly

    [K : P_m(T) Z^14] = 7^(2m-3),  m >= 2.

Indeed I-T² maps the whole integer lattice onto K, by differences
around the two 7-cycles. On K the determinant magnitude of I-T is 7:
the eigenvalues of T on K run over primitive roots of orders 7 and 14,
and Phi_7(1)=7, Phi_14(1)=1. Every later factor has determinant
magnitude 7² on K. This measures how using deeper and deeper filters
narrows the attainable correction lattice instead of filling it.

In characteristic seven the same filter behaves very differently.
Split by the sign of T^7. On each 7-dimensional sector a generator of
the cyclic 7-action is one unipotent block of length 7. On the plus
sector all m factors contribute one nilpotent order each; on the minus
sector I-T is a unit and the remaining m-1 factors each contribute one.
Thus the rank is max(7-m,0)+max(8-m,0), and the entire operator is zero
modulo seven at depth eight although its rational rank is still 12.
So modular disappearance of a defect is not proof of an exact contact
repair. The cube switches in §2 are different: they cancel every
physical coordinate by a literal pairing.

## Replay

Run from this directory:

```sh
python nested_core_check.py --board > board.txt
for s in slot_dlx closest_seams outside_choices two_word_seams trade_cube multi_probe first_patch_trade; do
    c++ -std=c++17 -O3 "$s.cpp" -o "/tmp/$s"
done
/tmp/slot_dlx project 25 71 8000 < board.txt > projected.txt
python slot_interface_check.py board.txt projected.txt
/tmp/closest_seams < board.txt > closest.txt
python seam_controls.py board.txt closest.txt
/tmp/outside_choices closest.txt < board.txt
/tmp/two_word_seams projected.txt < board.txt
/tmp/multi_probe closest.txt projected.txt 1 20000 5223 < board.txt
/tmp/trade_cube < board.txt > trades.txt
python trade_switch_check.py board.txt trades.txt
/tmp/first_patch_trade 12 < board.txt
```

The degree-57 question remains open on the unrestricted marked board.
