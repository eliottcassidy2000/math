# Three routes to an exact target selector

All infinite coverage objectives below are open. The affine gap law,
the crossing rule and the explicit seven-vertex control are exact.

## 1. Follow a connected repair slice, not an additive cube

For a common deletion tournament, let F(x) be the Hamiltonian count
as a function of apex contacts x. The inherited one-apex formula
(THM-4097/4114; APEX_BUNDLE_WINDOW.md) is

    F(x) = const + sum_i l_i x_i - sum_(i<j) w_ij x_i x_j,

with w_ij the sum of the two marked two-cover counts. Choose a fixed
set D of repair spokes. Fix their initial orientations and let the
remaining contacts vary. The signed increment from flipping D is

    F(x xor D) - F(x)
        = C - sum_(j notin D) x_j sum_(i in D) w_ij (1-2x_i).

Proof: terms involving two free contacts cancel; terms touching D
once contribute the displayed coefficient; terms confined to D
form the constant. Thus the +2 repair slice is an exact 0–1 linear
equation even though F on the slice remains quadratic.

Group free spokes by coefficient (or its absolute value after
complementing variables of negative coefficient). Within one
equal-coefficient group, swapping a selected and unselected spoke
preserves the repair exactly. Each fixed-occupancy stratum is
connected by these exchanges: exchange each element in the
symmetric difference with the corresponding missing element.
Different occupancies need not be connected, and count changes on
those exchanges need not be small.

The useful additional constraint is: keep only edges with baseline
jump at most 4 in absolute value. Every vertex x supplies the actual
odd pair {F(x),F(x)+2}. If two baselines differ by 0, 2, or 4, the
two pairs overlap or abut on the odd lattice. Therefore the union
of pairs along any connected component is a full odd interval
from its smallest baseline through its largest baseline plus two.
A target bracketed by a certified path is attained; the path need
not be monotone and its shifts need not add.

### Physical control

Take outgoing masks [46,28,72,0,109,78,11] on labels 0..6,
apex 6. Masks below toggle its spokes indexed by the other labels.
Fix repair D={1} (bit mask 2). Baseline settings 36 and 20 differ
by exchanging spokes 4 and 5. Their paired settings and counts are

    (36,38): (9,11);  (20,22): (13,15).

The exchange has baseline jump 4 and yields all four consecutive
odd counts 9,11,13,15 without a 2-step baseline translation.
The script counts all four literal tournaments independently of
the slice equations and checks equal exchange coefficients.
This is a tiny demonstration of the crossing method, not new
scalar coverage: the host even has an unchanged sink.

An independent smaller control uses five vertices, lexicographic upper-edge codes 71/70 and 77/76 (bit 1 means forward), with apex 0. The repair edge is (0,1); the same two-contact move flips (0,2) and (0,4) on both bases. The count pairs are again 9/11 and 13/15.

Run `python3 selector_slice_replay.py` for both literal controls.

Next experiment: use the previous dyadic lift or an approximate
absolute selector to place the source, then search for short
exchange-connected paths on its legal slice that bracket the
remaining signed residual. Do not discard exchanges merely
because their mixed coefficients are nonzero.

## 2. Give the dyadic pair a precise repair budget

The marked inflation of two known five-vertex sources gives the
exact pair 6*2^m ± 1. An explicit sufficient enlargement is to
realize all odd values in

    [6*2^m - (2^(m+1)-1), 6*2^m + (2^(m+1)-1)]

for each m>=3, by one construction or a bank of repair slices.
These windows abut: the upper end at m is 8*2^m-1; the lower
end at m+1 is 8*2^m+1. The first window starts at 33, so the
known small exceptions are not put inside a purported interval.
Together with individual admissible witnesses below 33 this
would prove Chan–Pak. This sets a precise sufficient budget,
not a claim that the single marked block can supply it.

To avoid one huge repair, try choosing among several physical
sources at each height so their smaller intervals collectively
bridge these endpoints. Store each source's response coefficients
and marked contacts, not just its central count. A modulus-only
selector or a count approximation without a repair interval
would not prove this criterion.

## 3. Replace growth by a target-first inverse on marked states

An alternative is to prove finite-state descent in a *typed*
target set: a state consists of an odd N plus the port obligations
needed by the next gluing. An inverse rule chooses a strictly
smaller admissible N' and a type, with a physical forward
construction from every witness of that type. Exact initialized
base witnesses then imply coverage by strong induction,
without widening any single window.

The proper first milestone is not closed residue arithmetic but
a total inverse on one actual lane, e.g. a subfamily of prime and
sevenfold-prime targets. Require a proof that the smaller state
has a source of the requested type: an unmarked witness at N'
is insufficient. Prime-level incidence and local congruence
coverage alone can miss a target without violating closure.

## Recommendation

Test route 1 first because it removes the previous requirement
of independent additive shifts. Use route 2 to specify when that
experiment becomes a proof, rather than measuring endpoints in
isolation. Keep route 3 parallel as an alternative when the
legal balance slices split into disconnected pieces.
