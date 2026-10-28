# Incoming repairs and a smallest regular nongraphic prime bank

2026-10-28. PROVED elementary controls / CITED matroid theorems /
FINITE-EXACT enumeration. This repairs incoming THM-4840 and HYP-9438 from
8668e990bd before integrating their claims. It does not retract their valid
fifteen-column equal-degree model or rank-five/six prime realizations.

## 1. An extra bit must be the correct function

Let B be the binary vertex-edge incidence matrix of K6. The incoming theorem
correctly identifies Hamming[15,11,3] with

    C={x: Bx is either 0 or the all-one six-vector}.

But the kernel of [B | 1] is the graph of the common vertex-degree bit,
y(x)=(Bx)_0. A triangle has y=0 and remains a weight-three word. This
sixteenth-column construction therefore has parameters[16,11,3], not[16,11,4].

The correct extension is

    C_ext={(x,z): x in C, z=sum_e x_e mod2}.

It has11 dimensions and is even. Since C has minimum weight3, its extension
has minimum weight4; triangle and perfect-matching words attain4. In the
four-dimensional quotient coordinates h_e, the check columns are(1,h_e)
for the fifteen edges and(1,0) for the new bit. These are the sixteen affine
points of F2^4 with leading coordinate1.

The two codes differ by the invertible shear of the last bit

    z=y + sum_e x_e + (Bx)_0.

This retains the old fifteen source coordinates, but it is not a Hamming
isometry. The state count stays2048 while distance changes. The incoming
checker tested the fifteen-column identity and never tested this claimed
sixteenth-column consequence. The repair checks both full2048-word codes.

## 2. Regularity, representation, and signing have different quantifiers

For a binary matroid M, rational representability iff regularity is the
classical Tutte theorem. It supplies some totally unimodular representation;
it does not assert a signed-graph incidence representation. The latter class
is strictly smaller: Qin–Slilaty–Zhou's
[regular excluded-minor theorem](https://doi.org/10.1017/S0963548309990241)
includes regular nonsigned-graphic examples R15 and R16.

Nor can every particular binary matrix representing a regular matroid be
signed totally unimodular. Consider columns(111,110,101,011), or

    A = [1 1 1 0; 1 1 0 1; 1 0 1 1].

Its last three columns form a triangle and its first is a coloop: the binary
column matroid is graphic, hence regular. Yet adjoining I3 makes its seven
columns the whole Fano plane. If A had a TU signing, adjoining I3 would remain
TU and would represent that Fano matroid over Q, a contradiction. Every512
sign choice is also checked directly by exact minors.

The repaired signing criterion is: a prescribed binary r-by-n matrix A can
be signed TU iff the binary column matroid of[I_r | A] is regular. Equivalently,
for a regular matroid already in standard representation[I_r | D], its
nonidentity part D has a TU signing. This is the standard-form scope of the
classical signing theorem. Changing signs is not a general equivalence between
unsigned square-sum equations and signed graceful-label equations.

Fano still has no faithful rational representation of *the same column
matroid*. Integer valuation lifts do exist, as the concurrent relation-gluing
packet demonstrates; their rational dependencies differ and mod2 reduction
can create dependencies through torsion. This distinction is indispensable.

## 3. Nine primes realize the cographic K3,3 obstruction

Order the K3,3 edges lexicographically between{0,1,2} and{3,4,5}. Choose
spanning tree edges with indices{0,1,2,3,6}, starting at0. The four remaining
fundamental cycles give a binary cycle-basis matrix. Its nine columns are

    15, 5, 10, 3, 1, 2, 12, 4, 8.

A cycle-basis matrix represents the dual of the graphic matroid. Under the
incoming rank-four prime table, these columns are respectively

    839, 89, 83, 5, 7, 2, 139, 19, 251.

Their p+1 square classes on rows(2,3,5,7) are exactly the stated columns.
Thus this bank realizes M*(K3,3); its square-sigma subsets are the binary cut
space of K3,3. All512 subsets are checked both by integer square roots of the
actual sigma-product and by the cographic rank identity

    rank_M*(S)=|S|+1-components(K3,3 minus S).

**Theorem.** The minimum size of a regular nongraphic distinct-prime
square-class bank is nine, attained by the displayed bank with largest prime839.

Proof: the constructed matroid is cographic, hence regular, and nongraphic
since K3,3 is nonplanar. Conversely the classical forbidden-minor theorem says
a nongraphic binary matroid has an F7, F7*, M*(K5), or M*(K3,3) minor.
Regularity excludes the two Fano minors. The remaining obstructions have10
and9 elements respectively, so fewer than9 cannot suffice. QED.

This refutes HYP-9438's ten-element minimum and proves one of its requested
explicit banks. It does not minimize the largest prime among all9-prime
realizations; that useful arithmetic optimization remains OPEN.

The unrestricted realization claim also needs its already-known loop guard:
only prime3 has square p+1, since p=k²-1=(k-1)(k+1). Two loops therefore
cannot occur in any distinct-prime bank. Under the stated prime-values
hypothesis, loopless matroids and matroids with at most one loop remain the
appropriate generalization; all simple ranks<=6 have the incoming finite banks.

## 4. Unit halves must survive reduction

At(2+omega) over norm3, omega maps to1. The legitimate abstract unit half
{1,-omega,omega²} maps to both1 and-1, so it does not orient a tournament.
Odd residue characteristic alone is insufficient. The required reduced set
R must satisfy R disjoint(-R) and R union(-R)=all nonzero residues. The
specified cubic-root half{1,omega,omega²} does meet this condition at norm3.

At(3+omega) of norm7, omega maps to4, not2; its cube-root set still maps to
{1,2,4}. At(2+i) of norm5, i maps to3, not2; {1,3} gives the same regular
five-tournament up to isomorphism. The named C3,T5,P7 examples survive with
correct quotient maps. A half exists in the characteristic-zero unit group
at norm4 too, but its image cannot be a half in characteristic2.

## Reproduction and scope

`python3 incoming_repair.py` and `python3 -O incoming_repair.py` give identical
JSON, saved in incoming_repair.json:3632 active checks, including28 exact
product-weight Fourier comparisons from THEORY section7. This is a same-agent
review with exact hostile controls, not a new independent-agent audit. The
prior source packet and canon are corrected in place and the mechanism is
recorded in MISTAKES.md.
