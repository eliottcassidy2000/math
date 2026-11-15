# One apex: exact run laws explain 5, 15, 45 and the first 63 words

Here the spine is transitive. The apex beats exactly the entries
marked 1 in a word w. After removing leading zeros and terminal
ones, which are singleton SCCs and do not change H, either the
word is empty (H=1) or the core begins in 1 and ends in 0 and is
strong. Removing the apex then empties **every** cycle. This is
a different, exact model from counting arbitrary walks with a
matrix power.

Every odd cycle therefore consists of the apex, a 1 position i,
an even subset of intermediate positions, and a 0 position j>i.
The empty collection contributes 1; cycles are mutually
incompatible. Thus a pair at distance d contributes 2 if d=1,
and 2^(d−1) otherwise, to H. This uses the
[odd-cycle collection formula](https://arxiv.org/abs/2307.05569).
Every contribution is positive and any induced deletion can only
remove such terms.

## A residue hole that rotates with the number of runs

Track the unsigned marked starts S in the cut recurrence.
Appending a mark 1 to a transitive spine leaves H unchanged and
takes S to 2S+1, with signed start equal to 1. Appending 0 adds
S+d to H, with d the incoming signed start, then takes S to 2S
and resets d to 0. This is the
[inherited singleton transfer](../chanpak_selective_target_20261114/CUT_STATE.md).
Start at (H,S,d)=(1,0,0) on the empty spine.

After k run pairs (each of the form a positive run of ones followed
by a positive run of zeros),
```
H ≡ k+1+[S≡−1 mod 3] (mod 3).
```
For the proof, a run of ones of parity p changes S mod 3 to
W=(-1)^p S+p, preserving the event S=−1 because the odd
transition swaps 0 and 1 and fixes −1. A run of zeros of parity
q changes S to (-1)^q W and adds 1+qW to H modulo 3.
Its first zero reads signed start 1, then all later zeros
alternate the sign of W, giving exactly that increment.
If q=1, the change in the indicator of S=−1 is W mod 3.
This proves the formula inductively. Hence at fixed k exactly
one residue, k mod 3, is forbidden in this family.

The other two residue classes are occupied at **every** ternary
precision. To see this, start with k−1 copies of the run pair
1100, which keeps S=0 mod 3, then append a single 1, leaving
S≡1 (mod 3). A final zero run of length a>=1 gives exactly
H=H_prefix+1−S+S*2^a. The power of 2 runs over all units modulo
3^r, so the two allowed output classes each lift fully.
The order of 2 is 2*3^(r−1): 4 has order 3^(r−1) by expanding
(1+3^j u)^3, which increases the 3-adic valuation by one.
This is a residue theorem with an explicit physical prefix,
not a claim of absolute coverage.

## Delete to two or four runs: a calculation without a census

Sum the pair weights over the run rectangles. For two runs
1^d 0^a the count is
`H=2+(2^d−1)(2^a−1)`.
For four runs 1^d 0^a 1^b 0^c put D,A,B,C=2^d,2^a,2^b,2^c.
Then
```
H=3+(D−1)(A−1)+(B−1)(C−1)+(D−1)(C−1)AB.
```
The first and second rectangles each add one for their adjacent
pair; the cross rectangle between the first ones and last zeros
has no adjacency and gains the factor AB.

In the two-run lane, H=7,21,63 would force the Mersenne product
to be the prime 5,19,61, respectively. None is a Mersenne
number. For four runs:
- With d=c=1 the exact factorization is H=(A+1)(B+1).
  Neither 7,21, nor 63 is possible: 21 and 63 would require
  a factor 7 or 21 of this shifted-power form.
- If d,c>=2 the minimum is 45; any increase in A or B from
  2, or D or C from 4, gives more than 63 by the cross term.
- If just one of d,c exceeds 1, use converse to put it at d.
  For 21, the only possibility not immediately exceeding 21
  by the cross term and positive rectangle terms is D=4,
  A=B=C=2, which gives 19. The minimum in this lane already
  exceeds 7. For 63, D>=16 overshoots by the first and cross
  terms. D=8 with A=B=2 gives 39, and any increase in A or B
  jumps above 63. At D=4 the count is precisely
  `H=(3A+1)(B+1)−2`. Now 65=13*5 forces A=B=4, yielding
  1100110 and its reverse-complement 1001100.

At three or more run pairs, deleting all excess digits to keep
three adjacent pairs leaves 101010 with H=31, already excluding
7 and 21. Exactly three pairs cannot give 63 by the residue
law. Four or more pairs contain 10101010, whose outside
endpoint pair contributes 64 before counting anything else.
Thus there are exactly two 63 words after neutral padding
is stripped, the converse pair above. This is the one-apex
transitive-spine classification, not an assumption that every
arbitrary tournament has such a spine. The separate
[emptying proof](EMPTYING_EXCLUSIONS.md) handles arbitrary cores.

## Why the two chains behave differently in this lens

All three 5-chain values have strong witnesses in the same model:
```
5  : 100       (a single two-run rectangle)
15 : 10010     (H=5*3 in the d=c=1 four-run lane)
45 : 1000110   (H=9*5 in that same lane)
63 : 1001100   (H+2=5*13, after changing an outer guard)
```
The ordinary SCC product also makes 5·3^a by adjoining directed
triangles, but this is not needed for those strong 15 and 45.
For 63 the different factorization is decisive: the physical
inverse sees 65, not a forbidden factor 7 of 63.

More generally the literal word 1 0^r 1^r 0^r gives
`H=2^(3r)−2^r+3`; r=2 is 63. The integer polynomial
f(x)=x³−x+3 maps each residue ball of Z_3 bijectively onto
3Z_3. Indeed f≡0 mod 3 and f'=3x²−1 is a unit; at each
step exactly one of the three next digits gives a specified
lift. Applying this on either unit ball and using the order
of 2 shows that for either exponent parity, this one
physical family realizes every multiple of 3 modulo 3^j
at any fixed j. It has arbitrarily high exact 3-adic
valuations, despite being a very sparse dyadic family of
actual absolute counts. The local clock cannot distinguish
absence of 21 from presence of 63.

## The same emptying reads the guard as an adjustable shift

For a single initial 1 and a variable final zero run (D=2),
the four-run equation factors exactly as
```
H+C−2 = (A+1)(1+(C−1)B).
```
The case C=2 is the common-junction two-wing product of the
[arbitrary-core exclusion proof](EMPTYING_EXCLUSIONS.md).
Adding a vertex beating only that junction is exactly appending
a final 0, hence moving to C=4. Thus deleting the outlet and
factoring the shifted target are two readings of the same
operation, not numerical coincidences.

At C=2 the wing order is invisible to H; exchanging A and B
leaves (A+1)(B+1). At C=4 the same two H=15 sources with
(A,B)=(4,2) and (2,4) read out as 33 and 37, respectively.
More generally exchanging the wings at guard C changes H
by (C−2)(B−A). For a≠b and C=4 this has 3-adic valuation
zero if a,b have opposite parity, and
1+v_3((b−a)/2) otherwise, by factoring powers of 4.
So equal scalar source counts can yield deliberately
different ternary endpoint behavior when a guard is read.
