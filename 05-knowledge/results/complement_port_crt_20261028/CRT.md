# Compatible CRT charts beyond 76, with the period-twelve source fibre

2026-10-28. **PROVED elementary chart law / FINITE-EXACT readers and mirrors.**
Same-agent audit. Local results R1-R4 belong to this file. The period-nine
reader and its two-bit mirror sidecar are inherited from
[COMPOSITE_76](../lucas_prime_reflection_towers_20261028/COMPOSITE_76.md).

Let O=Z[phi], phi^2=phi+1. A legal word w has n binary entries with no adjacent
ones, including the wraparound pair. Its reader is sum_i w_i phi^i in
R_n=O/(phi^n-1). Rotation is multiplication by phi. Literal reflection sends
entry i to entry -i. There are L_n legal words, where L_n is the Lucas number.
The determinant gives |R_n|=L_n for odd n and L_n-2 for even n.

## R1. Explicit larger odd charts

For the following periods, factor phi^n-1=2(A+B phi), where
A=(F_(n-1)-1)/2, B=F_n/2. Its odd norm factor is L_n/4. Each listed odd prime
p gives the linear root r=-A/B mod p; the factor at2 is F4. The ideals are
pairwise comaximal and their norms multiply to L_n, proving the stated ring
decomposition. The checker verifies the reader is bijective by enumerating
**every** legal word in these four odd universes, not by a count alone.

| n | Arithmetic ring | phi coordinates | Local nonzero periods | Words | Rotation charts |
|---|---|---|---|---:|---:|
| 9 | F4 x F19 | (omega,5) | (3,9) | 76 | 10 |
| 15 | F4 x F11 x F31 | (omega,4,19) | (3,5,15) | 1364 | 94 |
| 21 | F4 x F29 x F211 | (omega,24,179) | (3,7,21) | 24476 | 1170 |
| 27 | F4 x F19 x F5779 | (omega,5,2805) | (3,9,27) | 439204 | 16274 |

Here omega^2=omega+1. The exact chart-period histograms are respectively

    n9:  1^1, 3^1, 9^8;
    n15: 1^1, 3^1, 5^2, 15^90;
    n21: 1^1, 3^1, 7^4, 21^1164;
    n27: 1^1, 3^1, 9^8, 27^16264.

Notation d^m means m charts of length d. These are composite Lucas state
counts. They make no claim that Lucas numbers are prime infinitely often.

## R2. The relative-phase quotient, including its noncyclic boundary

Fix one multiplicative orbit in each local factor, of lengths d_1,...,d_s.
Their product has independent phase labels in product_i Z/d_i. Actual word
rotation advances every phase at once. Therefore its chart set is exactly

    (product_i Z/d_i) / <(1,...,1)>,
    number of charts = product_i d_i / lcm_i d_i.

This is a theorem about the specified group action, with no probabilistic
independence assumption. A zero coordinate contributes period1. It proves
that local orbit names alone are insufficient whenever this quotient is
nontrivial. If the last period is a multiple of every other, normalize its
phase to zero; the remaining coordinates are theta_i-theta_last mod d_i.

For all-nonzero coordinates in the four rows:

| n | Local orbit choices | Relative-phase group | Global unit charts |
|---|---:|---|---:|
| 9 | 2 | Z/3 | 6 |
| 15 | 4 | Z/3 x Z/5 = Z/15 | 60 |
| 21 | 40 | Z/3 x Z/7 = Z/21 | 840 |
| 27 | 428 | Z/3 x Z/9 | 11556 |

The last group has order27 and exponent9; it is **not** a single 27-cycle.
Thus even “retain a phase of the right cardinality” can use the wrong action.
At this stage two compatible relative coordinates are required. This is the
same diagonal-quotient mechanism as shared-port frames, not an identification
of these arithmetic states with tournament vertices.

Choose a canonical representative r of each literal rotation orbit. Write
w=rotation^i(r), and let reflection(r)=rotation^j(r'). Then exactly

    reflection(r,i)=(r',j-i mod period(r)).

`crt.py` verifies this for every word; the 94 seed records at period15 are
stored in `crt.json`. Translation in the ring must decode and rechart; it is
not assumed to preserve a rotation chart.

Reflection on the F4 coordinate alone is Frobenius z->z^2, since
omega^(-i)=(omega^i)^2. But dropping this coordinate from the full odd reader
loses reflection transport on the other factors. At n=9,15,21,27, exhaustive
tests give four different reflected odd-coordinate outputs for **every
nonzero** odd-coordinate input, and one for the zero input. Retaining F4 is
sufficient and two bits are necessary at those nonzero fibres. This all-input
statement is FINITE-EXACT at the four stated periods; no all-period theorem
is inferred. Reflection-fixed word counts are 8,34,144,610.

## R3. Period twelve: a ring of 320 states and 322 source words

Here phi^12-1=88+144 phi=8(11+18 phi), with Norm(11+18 phi)=-5. Thus

    R_12 = O/(8) x F5,          phi=((0,1),3).

The first factor has 64 elements and characteristic8. It is **not F64**,
which has characteristic2. Its 48 units have phi-period12; the 12 elements
of exact 2-adic valuation1 have period6; the 3 of valuation2 have period3;
zero has period1. The four nonzero F5 elements have period4.

The 322 legal words map onto all 320 ring states. Every state has one preimage
except zero, whose three preimages are the zero word and the two alternating
words. The ring has 30 rotation charts; the extra alternating-word two-cycle
gives the word space 31 charts. This recovers the inherited
[golden-reader boundary](../session91_packing_golden_carriers_20261029/golden_fibres/GOLDEN_FIBRE_DYNAMICS.md)
with a literal example; it does not replace that all-period theorem with an
empirical count. That prior result already proves a single named toggle
separates the three exceptional words at every even period at least4.

Uniform word measure pushes forward to mass3/322 at zero and1/322 at every
other ring state. Uniform ring measure would assign1/320 everywhere. They
are different source laws.

**No bijective lift of all translations to the original word set exists.**
If pi:X->R is equivariant for a transitive action on R by bijections of X,
then each lifted bijection maps the entire fibre over r bijectively onto the
fibre over its translate. All fibres must have the same cardinality. Sizes3
and1 contradict this. This is an exact realization obstruction; adding a
weight to zero preserves measure accounting but does not create the missing
bijections. Rotations and reflection do descend here because they preserve
the exceptional fibre; translation needs a different carrier.

The period12 projection onto O/(8) alone has up to four reflected outputs;
the projection onto F5 alone has up to five. Full ring state is sufficient
for the descended reflection, but it still does not identify which zero
word was supplied when a later operation can inspect the word.

## R4. What the Frobenius and LCA comparisons preserve

The inherited F64/F256 orientation guards are replayed in `ports.json`:
for rotations of normal-basis positions, the allowable shifts have only odd
cycles, giving {0,2,4} at length6 and {0} at length8. A nontrivial even-order
permutation cannot preserve a tournament. This is a guard on a specified
position action, not a claim that arbitrary field elements are vertices of
the same tournament. Likewise O/(8)'s multiplication clock12 is not the
Frobenius clock of a 64-element field.

LCA transport combines relative maps along a tree; CRT transport combines
phase differences modulo compatible local periods. Their common algebraic
step is quotienting a common gauge while retaining relative coordinates.
Their realizations differ. Reversal-marking loops can have holonomy; a chosen
CRT chart decoder has literal state labels and no such ambiguity. The precise
bridge is the diagonal quotient R2, not an equality of these objects.

**Next target.** Prove or refute the observed four-output reflection law for
all odd periods divisible by3 under suitable nondegeneracy hypotheses.
Repeated odd prime factors require local rings and valuation strata, rather
than the squarefree field chart used in R1.
