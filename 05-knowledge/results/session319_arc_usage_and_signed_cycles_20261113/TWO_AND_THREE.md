# Two and three as occupancy clocks, rather than a common conjecture

## A general recursion from a finite quotient

Replace a marked quotient vertex by a transitive chain of length m.
Fix a number q of labelled path slots (empties may be allowed). A
layout is the ordered tuple of exterior vertex words with each maximal
chain run replaced by a token. Tokens in one slot are separated by at
least one exterior vertex; all transitions have to exist in the quotient.
The number c_j of legal layouts with j tokens is finite and independent
of m. For m>=1 the exact cover count is
```
sum_j c_j j! S(m,j),
```
where S(m,j) is a Stirling number. Each filling assigns chain vertices
onto the j labelled runs, whose internal order is then forced. Thus
the exponential bases are integers from 1 to the largest legal visit
count, and their coefficients are determined by a finite native layout,
not chosen arbitrarily. Expanding the surjections by inclusion-exclusion
gives the exact recurrences.

For the [reservoir pair](RESERVOIR_DARK_ARCS.md), only x can enter the
chain, so a single path has at most two visits and two ordered slots at
most three visits. The one-path layouts agree, 27 and 9. In two slots the
differences B−A of the one-, two-, and three-token layout counts are
exactly (0,12,4). This is checked by enumerating the exterior words and
possible slot cuts in replay_occupancy.py. Consequently
```
H_A=H_B=27+9(2^m-2)=9(2^m+1),
Q_2(B)-Q_2(A)=12(2^m-2)+4(3^m-3*2^m+3)=4(3^m-3).
```
The 2^m term cancels in the difference because of the exact 12:4
relationship between layout changes. The 3^m term is a third visit that
the scalar observer cannot see. In particular the ratio of discrepancy
to count is asymptotic to (4/9)(3/2)^m.

## Infinite local limit of the weighted chain

Let J=max{j:c_j>0} and draw a uniform legal cover. For large m, layouts
with fewer than J runs have total probability O(((J−1)/J)^m) if J>=2:
the upper bound for their fills is O((J−1)^m), and the J-run count is at
least J^m−J(J−1)^m for each layout. At J runs, the exact chain coloring
is uniform over all surjections onto J tokens. Compare it to independent
uniform J-color assignments: the probability of failing to use all J
colors is at most J((J−1)/J)^m. Thus the internal positive-weight
observer has a common limiting geometry, independent of the exterior
layout details. For chain vertices at distance d, their probability of
being consecutive within a run tends to
```
(1/J) (1−1/J)^(d−1).
```
Their colors must agree and the d−1 intervening vertices must avoid that
color. For J=1 the probability is exactly [d=1]. For the ordinary paths
in both reservoir charts J=2, so the local weight law tends to 2^(−d),
even though the two charts have different dark exterior arcs. This is
an explicit example of local agreement coexisting with a global
continuation defect.

## Which part transfers to signed Collatz?

For an actual odd Collatz word the two powers have a different source:
3^r comes from the r odd updates, 2^s from their total halving count.
A closed word requires an exact integer fixed point with denominator
2^s−3^r before enforcing its valuations. The familiar clocks for the
known positive and three negative accelerated cycles give differences
1, −1, −1, −139 respectively:
```
2^2−3,  2−3,  2^3−3^2,  2^11−3^7.
```
This puts the two shortest negative clocks at the two neighboring
powers 3−2=1 and 9−8=1. The seven-phase cycle shows why nearness alone
cannot classify cycles: its larger gap is allowed because the actual
carry is an exact multiple. The [signed-clock note](SIGNED_CLOCK_RESERVOIR.md)
retains that carry and proves a no-mixing rule in the known nontrivial
negative phase banks. The transferable idea is to enumerate admissible
local words with their exact carry or occupancy, then solve their native
compatibility condition. It is not a proof that a power-gap theorem
solves Collatz or that two-slot cover dominance gives Chan–Pak.

## The flipped arc itself is not dark

Relabel A by swapping exterior labels 1 and 5, as in the
[polynomial note](RESERVOIR_POLYNOMIAL.md). The single changed pair is
then 1→6 versus 6→1. In both orientations, exactly six one-run and
two two-run layouts traverse that arc. Its absolute weight is therefore
6H(U)+2F_2(U), exactly 2/9 of the full Hamilton count, for *every*
nonempty U. The paths avoiding it are literally the same across the flip.
Thus the neutral operation trades a macroscopic fraction of the paths
while changing the dark status of other arcs. It is not a perturbation
that simply flips an invisible arc.

For the signed comparison, the odd-only Collatz cycles are the positive
fixed point 1 (clock 2), and the known negative cycles -1 (clock 1),
(-5,-7) (clock 1,2), and (-17,-25,-37,-55,-41,-61,-91) (clock
1,1,1,2,1,1,4). A clock entry counts all halvings after the corresponding
odd update. The full unaccelerated cycle lengths are therefore 3, 2, 5,
and 18, respectively. This keeps the odd-only, single-halving and
unaccelerated conventions separate. See [Lagarias's overview](https://www.ams.org/bookstore/pspdf/mbk-78-prev.pdf)
and the earlier [four-vertex comparison](../s92_prime_powers_quaternions_bases_20261111/README.md).
The four known Collatz cycles are not in bijection here with the four
ordinary order-four tournament classes; the exact constructions explain
the clocks and the three arc-weight profiles by different mechanisms.

The scalar ray and its continuation gap have opposite 3-adic behavior.
For m>=2, the gap has exact 3-adic valuation 1. The common path count
has valuation 2 for even m and 3+v_3(m) for odd m: in the latter case
2^m+1=1-(-2)^m and the usual odd-prime lifting-of-exponent identity
applies to (-2)-1=-3. Neither path count is divisible by 7, since the
powers of 2 modulo 7 cycle through 1,2,4. Thus the scalar ray can
look arbitrarily deep at 3 while its continuation discrepancy stays
visible modulo 9. It does not supply prime or 7p coverage.
