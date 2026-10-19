# Signed cube packets: an all-order norm bound and an exact integer-clock gate

**Status: PROVED elementary matrix, norm, flow and clock lemmas; CITED
identifications with sporadic sequences; FINITE-EXACT independent checks;
OPEN application to the actual Collatz source and geometric selectors.**
Self-audited research note, no priority claim or new theorem ID.
Session `frontier_petersen_oct18`, continuation `20261018n`.

The requested papers supply a concrete model in which cancellation is
controlled at every order. Gorodetsky's signed four-by-four matrix M satisfies
MM^T=3I. This gives |A_delta(n)|<=9^n and an exact factorial-weighted norm
identity for the entire coefficient vector. Its three-channel counterpart
is the Fourier matrix over the Eisenstein integers and gives B_n^2<=27^n.
The integer-clock extension has a sharp gate: every clock curl must vanish.
A clock change by two passes the parity test but fails this stronger gate.

This connects the papers to the [ten-label Hamming construction](ten_label_matching_hamming_20261018n.md)
through an actual Hadamard completion, and to the ROOT problem through a
specific norm-preserving packet target. Neither transfer assumes universal
entry or discards a marked source.

## 1. Literature imports, with distinct scopes

**CITED — Gorodetsky, requested v2.**
[New representations for all sporadic Apéry-like sequences](https://arxiv.org/html/2102.11839v2),
Theorem 1.1 and Proposition 3.3, provide constant-term representations for
the catalogued fifteen sequences. Fourteen have the stated unique-interior-
lattice-point construction; the eta case is excepted. Corollary 2.4 gives D3
congruences for those fourteen. Section 3.4 explicitly prints the matrix M
used below, and Section 3.5 expresses sequence B through cube roots of unity.
Thus neither the matrix nor these representations are new here. This is a
catalogue, not a proof that no further sporadic recurrence exists.

**CITED — Alinquant–Osburn, requested v1.**
[On sporadic sequences](https://arxiv.org/html/2312.07134v1), Theorem 1.2,
proves A_delta(mp^r)=A_delta(mp^(r-1)) modulo p^(3r), for p>=5 and m,r>=1.
The proof separates coordinatewise p-divisible tuples from primitive tuples,
then uses valuation classes and alternating reciprocal-square sums. The
useful imported mechanism is a partition into packets retaining the original
linear constraints and multinomial weights. Divisibility is not a real norm
estimate. Our script checks the stated congruence only at the declared six
(p,r) pairs; the universal assertion here is a citation.

**CITED — Kazalicki, requested v1.**
[Congruences for sporadic sequences and modular forms for non-congruence subgroups](https://arxiv.org/html/1901.03098v1),
Theorem 1, identifies F((p-1)/2) mod p with the p-th coefficient of its
weight-three, level-24 newform for odd primes. Proposition 1 transports
three-term congruences through a p-integral change of parameter, reversibly
when the first coefficient is a unit. The application retains a differential
and the cover t=s^2, rather than only a scalar sequence. The HTML's two
general-recurrence displays omit a square: we use the explicit sum for F and
the repository's (n+1)^2 convention. It gives F(2)=42, whereas the unsquared
display gives 84. This typo does not change the cited theorem's sequence.

The latest [reanchoring work](reanchoring_kernel_20261018c.md) already gives
the A-to-F binomial transform with the required carry correction, and the
[six-triple audit](zagier_six_triples_20261018c.md) keeps the accessory
parameter. We do not replace either with a weaker numerical coincidence.
The recent [Mathieu synthesis](factorization_vertex_supersymmetric_bridge_20261018.md)
has level-24 eta-products of weight two; sharing level 24 does not identify
them with Kazalicki's weight-three form.

## 2. The actual signed matrix and boundary fibre

Use the paper's matrix

\[
M=\begin{pmatrix}
0&1&-1&1\\
-1&1&1&0\\
1&0&1&1\\
1&1&0&-1
\end{pmatrix}.
\]

Let P_n(x)=product_i (sum_j M_ij x_j)^n. Its balanced coefficient is

\[
[x_1^nx_2^nx_3^nx_4^n]P_n=A_\delta(n).
\]

The nonzero support graph is K4,4 minus a perfect matching, hence the cube
Q3. A monomial choice is a nonnegative integer flow f on its twelve labelled
edges, with every row and column total equal to n. Its signed weight is

\[
w(f)=(-1)^{f_{a2}+f_{b1}+f_{d3}}\frac{(n!)^4}{\prod_e f_e!}.
\]

The matrix's distinct rows overlap in two columns with opposite products.
Therefore MM^T=3I. Equivalently each of the cube's six square faces has
negative sign holonomy. The signs are essential; replacing M by |M| loses
this orthogonality.

**Proposition 1 (PROVED, all n).** Every integral flow in this fibre is a sum
of n perfect matchings of the same graph. The signed weight is not in
general the number of ordered decompositions, even up to sign.

For the existence statement, a set of row vertices S sends n|S| units into
its neighbour set, whose total capacity is n|N(S)|. Hall's condition follows.
Subtract a perfect matching in the positive support and induct on n. The
original flow and all boundary totals are retained.

For the loss statement, the first explicit witness at n=2 has edge tuple
`(0,0,2,0,1,1,1,1,0,1,1,0)`: its unsigned weight is 8 but it has only two
ordered matching decompositions. At n=2 the unsigned totals are 297 versus
81. At n=3 even the signed totals differ: A_delta(3)=3 versus 3^3=27.
Decomposing the support is universal, but forgetting the row-word
multiplicities destroys the weighted statement. QED.

The four-cycle increment `(1,-1,0,0,-1,1,0,0,0,0,0,0)` has all eight boundary
totals zero. Whenever feasible it reverses the sign, but generally changes
the magnitude. It is a model for a source/port-preserving packet move, not
an automatic equal-weight cancellation involution.

## 3. An all-order bound that is built into the representation

For a polynomial P=sum_alpha c_alpha x^alpha, use the Fischer norm

\[
\|P\|_F^2=\sum_\alpha |c_\alpha|^2\alpha!,
\qquad\alpha!=\prod_i\alpha_i!.
\]

It is the Gaussian norm on complex coordinates, so unitary substitutions
preserve it. Equivalently it is the norm inherited from symmetric tensor
powers. These are factorial weights, explicitly different from counting
flows without multiplicity.

**Proposition 2 (PROVED).** For every n>=0,

\[
\sum_{|\alpha|=4n}|[x^\alpha]P_n|^2\alpha!
 =3^{4n}(n!)^4,\qquad |A_\delta(n)|\le9^n.
\]

Proof: U=M/sqrt(3) is orthogonal. The source monomial
x1^n x2^n x3^n x4^n has squared norm (n!)^4. Substitution by M multiplies
its degree-4n norm squared by 3^(4n). Keeping the single balanced coefficient
in the norm sum gives the bound. The source and target occupation vectors
are both (n,n,n,n); no change of representative is involved. QED.

There is a second short proof of the scalar bound. On the unit four-torus,
sum_i |(Mx)_i|^2=12. AM–GM gives product_i |(Mx)_i|<=9. Integrating the
balanced coefficient of P_n against its unit-modulus character gives 9^n.
The full norm identity retains more information than this scalar proof.

More generally, if MM*=lambda I and |alpha|=|beta|=d, then

\[
|[x^\beta]\prod_i(Mx)_i^{\alpha_i}|
\le\lambda^{d/2}\sqrt{\alpha!/\beta!}.
\]

This is a general representation-theoretic estimate; we claim no historical
novelty. It also clarifies exactly which factorial correction is lost when
unequal source and target occupations are identified.

### The three-channel Eisenstein companion

For omega^2+omega+1=0, take the three-by-three Fourier matrix F3. Then
F3 F3*=3I and

\[
\prod_{i=1}^3(F_3x)_i=x^3+y^3+z^3-3xyz.
\]

The balanced coefficient of its nth power is (-1)^n B_n, by the imported
representation of Zagier B. The same proof gives, for all n,

\[
\sum_{|\alpha|=3n}|c_\alpha|^2\alpha!=3^{3n}(n!)^3,
\qquad B_n^2\le27^n.
\]

Here the support graph is K3,3. Thus K3,3 and the cube furnish precise
three- and four-channel versions of one weighted cancellation mechanism.
This statement uses their matrices and preserved norms, not a claim that
their graph counts or sporadic-sequence labels are interchangeable.

## 4. From the cube to Hamming, with the gauge retained

Permute M's columns to order (0,3,1,2), then negate its third row. This gives

\[
C=\begin{pmatrix}
0&1&1&-1\\-1&0&1&1\\-1&-1&0&-1\\1&-1&1&0
\end{pmatrix},\qquad C^T=-C,\quad C^2=-3I.
\]

Consequently H4=I+C satisfies H4 H4^T=4I. Normalize its first row and column
to +1 and take the Sylvester doubling H8=H2 tensor H4. Replace signs by bits
in its eight rows and their negatives. The resulting sixteen words form
RM(1,3), the extended Hamming code. Adjoin a repetition pair to obtain one
of the 1350 ten-label missing structures.

Every operation is explicit and checked. The row/column signs, column
permutation, and adjoined pair are needed sidecars. An initial computational
candidate omitted the normalization and produced an affine translate rather
than the asserted linear code; that failure is retained in MISTAKE-642.
The repaired construction preserves a chosen zero word.

This gives the connection diagram

```text
signed cube matrix M -- orthogonal normalization --> uniform symmetric-power norm
        |
        +-- marked gauge --> skew conference C --> H4 --> H8 --> I2 + H8
```

The top branch retains numerical amplitudes. The lower branch retains a
binary relation after recording its gauge. Their outputs are not identified.

## 5. The integer-clock gate and changing generators

Give every nonzero entry an integer clock tau_ij and form
M(z)_ij=M_ij z^(tau_ij). Use formal adjoint z*=z^(-1).

**Proposition 3 (PROVED).** The following are equivalent:

1. M(z) M(z)*=3I as a Laurent-polynomial identity.
2. Every square face has zero alternating clock sum.
3. There exist integer row and column potentials a_i,b_j with
   tau_ij=a_i+b_j on every edge.

For a pair of distinct rows there are exactly two common columns. Their
Gram terms have opposite signs, so they cancel iff their exponents agree;
this is the corresponding face equation. The cube's square faces generate
its cycle space. Zero circulation on them permits integer potentials,
constructed along a spanning tree. Conversely potentials give
M(z)=diag(z^a) M diag(z^b), proving the identity. QED.

The clock solution space has dimension seven: eight vertex potentials minus
one gauge. The graph has cycle rank 12-8+1=5. With balanced margins n, the
total clock is n(sum a_i+sum b_j), the same for every contributing flow.

**Hostile:** changing only entry (0,1)'s clock from 0 to 2 gives the off-
diagonal Gram entry z^2-1. It vanishes at z=-1 but not formally, or at z=i.
Thus even a full parity-cancellation check need not preserve integer-clock
orthogonality. This strengthens the diagnostic role of the inherited
[clock-cost lemma](root_eisenstein_first_arrival_20261018j.md#4-a-genuinely-uniform-but-weaker-certificate),
without replacing its lag-cost conclusion.

Clocked normalized matrices satisfying this identity are unitary on the
integer-clock Hilbert space l2(Z) tensor C^4. Any composition, with arbitrary
changes between such generators, remains unitary. The analogous statement
holds for diagonal-clock conjugations of F3/sqrt(3). Permutations and their
recorded gauges can change at each step; the matrix identity controls every
depth. The source vector, channel labels and full clock register are kept.

## 6. What this does and does not buy for ROOT

The inherited target is source-specific:

\[
\sum_u\operatorname{Norm}(W_{u,N,d})
\le\tfrac32\sum_s\operatorname{Norm}(A_N(s)),
\]

for every specified scale and finite horizon. The
[source-energy note](root_eisenstein_first_arrival_20261018j.md) proves the
reduction and linear source bound, but leaves this inequality OPEN. It also
refutes unrestricted contraction claims. The present matrix construction
does not supply a factorization of that actual transport.

A concrete reroute is to seek **source-specific preparations and clocked
packet operators** with these exact Gram identities, retaining first-arrival
ports and the moving complement. Their terminal readout must reproduce the
actual W polynomial, not a phase trace or a different source distribution.
Any remaining readout/erasure cost must be bounded on those actual prepared
states. This is the unpaid step, and is not discharged by unitary dilation.

Indeed, storing k histories orthogonally is easy. Erasing them coherently
by sending all k basis vectors to the same terminal has operator norm
sqrt(k). On the all-ones vector its energy ratio is k. Hiding that readout
would simply hide the original concentration problem. The cube succeeds
because its signed entries give an exact norm identity before such a loss;
ordinary matching coverage alone failed already at n=2 in Section 2.

The next useful tests are therefore local and hostile:

- Retain the actual source row, seed and integer clocks in a candidate
  four-channel reset packet. Compute its weighted Laurent Gram and clock
  curls before fitting a depth-dependent envelope.
- If a curl is nonzero, keep it as an explicit clock-defect state. Test the
  cost of eliminating that state on the actual source, including the moving
  complement. No bounded-cost assertion is made yet.
- For geometry, try realizing the companion note's four affine Hamming
  charts with one compatible geometric selector. Its exact ordered boundary
  support, not the 30/90/1350 counts, is the acceptance test.

These tests can refute a proposed packet cheaply. They also specify an
all-depth controlled model to aim for, rather than assuming every larger
enumeration improves the uniform ROOT estimate.

## 7. Exact evidence and scope

[Script](../../04-computation/experiments/cube_flow_packets_20261018n.py),
[receipt](cube_flow_packets_20261018n.json).

```bash
python3 04-computation/experiments/cube_flow_packets_20261018n.py --output /tmp/cube.json
python3 -O 04-computation/experiments/cube_flow_packets_20261018n.py --output /tmp/cube-O.json
cmp /tmp/cube.json /tmp/cube-O.json
```

All U(n), n=0,...,8, are enumerated independently by margin constraints and
by matching convolution; their counts are
`1,9,42,138,363,819,1652,3060,5301`. Signed sums are checked against both the
single-sum formula and the recurrence. Complete coefficient vectors satisfy
the Fischer identities at n=0,...,8 in both the three- and four-channel
models. The all-order conclusions follow from the displayed proofs, not
from those finite tests. Supercongruence checks use p=5,7,11 and r=1,2,
m=1. The receipt retains the primitive-row ledgers and the first lost-weight
witness. The code uses integer arithmetic and active checks under `-O`.

No universal ROOT bound, new Collatz acceptance family, six-colour forcing
context, or planar Jacobian result is claimed by this session.
