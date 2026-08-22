# Tournament switching coronals: second moments, deletion Gram, and order join

**Status: PROVED (EXACT ALGEBRA) + VERIFIED-EXACT IN THE STATED FINITE
UNIVERSES; canonical statement
[THM-3322](../01-canon/theorems/THM-3322-tournament-switching-second-moment-deletion-gram-and-order-join-law.md).**
This note extends
[THM-3315 -- tournament cut-switching centered-coronal walk
compiler](../01-canon/theorems/THM-3315-tournament-cut-switching-centered-coronal-walk-compiler.md)
in two directions:

1. it gives the exact bivariate switching-cube second moment of the signed
   coronal numerator and of the switched adjacency determinant; and
2. it proves that the full `(P,N)` interface composes under order join of the
   selected switched targets.

The first moment is purely centered-spectral.  The second moment is not: the
one additional coordinate is the Gram kernel of the vertex-deleted centered
characteristic polynomials.  A sharp order-seven pair proves that this
sidecar is not determined by the centered characteristic polynomial.

Companions:

- [`tournament_switching_coronal_second_moment_join_compiler_20260803.py`](../04-computation/tournament_switching_coronal_second_moment_join_compiler_20260803.py);
- its [frozen exact output](../05-knowledge/results/tournament_switching_coronal_second_moment_join_compiler_20260803.out).

No claim of priority in the external literature is made.

## Inheritance pass

**Closest proved mechanism.**  THM-3315 proves that, for centered adjacency
`C=2A-J` and a cut sign `d~-d`, the two polynomials

```text
P(z)=det(I-zC),
N_d(z)=d^T adj(I-zC)d                                      (1)
```

compile the total directed-walk series of the selected switch.  Its first
cut-cube average `E N_d=nP-zP'` is the starting point here.

**Closest moment mechanism.**  [THM-2501 -- switching fourth moment, signed
four-cycles, and Gram energy](../01-canon/theorems/THM-2501-switching-fourth-moment-signed-c4-and-gram-energy.md)
uses Rademacher character orthogonality to identify the first nonconstant
even moment of a symmetric switching form.  The present observable is not
that symmetric edge form: it is a polynomial-valued quadratic observer built
from a tournament adjugate.  The same orthogonality principle applies only
after its symmetric off-diagonal coefficients are typed correctly.

**Canonical hostile.**  The transitive triangle and directed triangle are
switching-equivalent and have the same cut-cube distribution of `N_d`, even
though a selected cut can give radically different walks.  Thus they are a
positive class-average control, not a witness that the second moment needs
more than `P`.

**Corrected near miss.**  It is tempting to differentiate the first-moment
formula and expect every second-moment term to remain spectral.  This loses
the diagonal-cofactor squares.  The repaired formula contains an explicit
vertex-deletion Gram sidecar, and the first failure of `P` alone occurs at
order seven.

**Least-used relevant sidecar.**  The diagonal of the adjugate is a deck:
each entry is the centered characteristic polynomial after deleting one
vertex.  The Boolean degree-two reconstruction perspective is consonant with
[THM-2355 -- component-deletion Gram reconstruction and twist-energy phase
transport](../01-canon/theorems/THM-2355-component-deletion-gram-and-twist-energy-phase-transport.md),
but no cross-domain identification is used in the proof.

## 1. Objects and the exact Walsh expansion

Let `T` be a labelled tournament of order `n`, let `A` be its row-to-column
adjacency, and put

\[
 C=2A-J.                                                     \tag{2}
\]

Then

\[
 C^T=-2I-C.                                                  \tag{3}
\]

For formal variables `z,w`, define

\[
 \begin{aligned}
 P(z)&=\det(I-zC),\\
 B(z)&=\operatorname{adj}(I-zC),\\
 N_d(z)&=d^TB(z)d,
 \end{aligned}                                               \tag{4}
\]

where `d` is uniform on the `2^(n-1)` sign vectors modulo global negation.
This quotient has the same expectation as independent Rademacher signs,
because every expression in (4) is even under `d -> -d`.

Write

\[
 M(z)=\operatorname{tr}B(z).                                 \tag{5}
\]

Since `d_i^2=1`, the exact Walsh expansion is

\[
 N_d(z)
 =M(z)+\sum_{i<j}\bigl(B_{ij}(z)+B_{ji}(z)\bigr)d_id_j.      \tag{6}
\]

The characters `d_i d_j`, indexed by unordered pairs, are mutually
orthogonal on the cut cube.  Therefore

\[
 \boxed{
 \operatorname{Cov}\bigl(N_d(z),N_d(w)\bigr)
 =\sum_{i<j}
   \bigl(B_{ij}(z)+B_{ji}(z)\bigr)
   \bigl(B_{ij}(w)+B_{ji}(w)\bigr).}                         \tag{7}
\]

This already proves positivity at `w=z` for every real `z` and identifies
what a second moment can forget: it retains the total Gram kernel of the
degree-two Walsh coefficients, not their assignment to the tournament
edges.

## 2. The deletion-Gram second-moment formula

Define the two trace kernels

\[
 \begin{aligned}
 K_-(z,w)&=\operatorname{tr}\bigl(B(z)B(w)\bigr),\\
 K_+(z,w)&=\operatorname{tr}\bigl(B(z)B(w)^T\bigr),
 \end{aligned}                                               \tag{8}
\]

and the diagonal-cofactor kernel

\[
 \mathcal D_T(z,w)=\sum_{i=1}^{n}B_{ii}(z)B_{ii}(w).         \tag{9}
\]

Expanding the traces entrywise gives

\[
 K_-(z,w)+K_+(z,w)-2\mathcal D_T(z,w)
 =\sum_{i<j}(B_{ij}(z)+B_{ji}(z))
             (B_{ij}(w)+B_{ji}(w)).                          \tag{10}
\]

Combining (7) and (10) yields

\[
 \boxed{
 \mathbb E_d[N_d(z)N_d(w)]
 =M(z)M(w)+K_-(z,w)+K_+(z,w)-2\mathcal D_T(z,w).}            \tag{11}
\]

The terminology “deletion Gram” is literal.  A diagonal adjugate entry is a
principal cofactor:

\[
 B_{ii}(z)=\det(I-zC[V\setminus\{i\}]).                      \tag{12}
\]

If `p_i(z)` denotes this vertex-deleted centered characteristic polynomial,
then

\[
 \mathcal D_T(z,w)=\sum_i p_i(z)p_i(w).                      \tag{13}
\]

Switching sends `C` to `DCD` and `B` to `DBD`; hence every `B_ii` is fixed.
Relabelling permutes the `p_i`.  Thus (13) is both switching-invariant and
relabel-invariant.

## 3. The two trace kernels are determined by `P`

Put

\[
 R_z=(I-zC)^{-1}.                                            \tag{14}
\]

Jacobi's identity, exactly as in THM-3315, gives

\[
 M(z)=P(z)\operatorname{tr}R_z=nP(z)-zP'(z).                 \tag{15}
\]

The ordinary resolvent identity is

\[
 R_zR_w=\frac{zR_z-wR_w}{z-w}.                              \tag{16}
\]

Because of (3), `R_w^T=(I-wC^T)^(-1)` is another rational function of `C`.
Multiplying through by its two affine denominators gives the second identity

\[
 R_zR_w^T=\frac{zR_z+wR_w^T}{z+w+2zw}.                      \tag{17}
\]

Indeed,

```text
z(I-wC^T)+w(I-zC)
 =z((1+2w)I+wC)+w(I-zC)
 =(z+w+2zw)I.                                                (18)
```

Since `B(z)=P(z)R_z` and `tr R_w^T=tr R_w`, tracing (16)--(17) proves

\[
 \boxed{
 (z-w)K_-(z,w)
 =zM(z)P(w)-wP(z)M(w),}                                     \tag{19}
\]

and

\[
 \boxed{
 (z+w+2zw)K_+(z,w)
 =zM(z)P(w)+wP(z)M(w).}                                     \tag{20}
\]

Equations (19)--(20) are polynomial identities: the left sides exhibit the
divisibility of the right sides.  They remain valid on the apparent
denominator loci by polynomial continuation.

Consequently the exact information content of (11) is

```text
first moment:   P alone;
second moment:  P plus the deletion-Gram kernel D_T.          (21)
```

This is an equivalence, not merely sufficiency: given `P`, equation (11)
recovers `D_T` from the second moment by subtraction and division by `2`.

On the diagonal `w=z`, (19)--(20) reduce to

\[
 \begin{aligned}
 K_-(z,z)&=P(z)\bigl(M(z)+zM'(z)\bigr)-zM(z)P'(z),\\
 K_+(z,z)&=\frac{P(z)M(z)}{1+z}.
 \end{aligned}                                               \tag{22}
\]

The second quotient is again a polynomial identity.  Hence

\[
 \boxed{
 \operatorname{Var}N_d(z)
 =P(M+zM')-zMP'+\frac{PM}{1+z}-2\sum_i p_i(z)^2.}            \tag{23}
\]

## 4. Switched adjacency determinants

THM-3315 identifies

\[
 Q_d(z)=\det(I-2zA(T^d))=P(z)-zN_d(z).                       \tag{24}
\]

Therefore

\[
 \mathbb E Q_d(z)=(1-nz)P(z)+z^2P'(z),                      \tag{25}
\]

and the complete bivariate covariance is simply

\[
 \boxed{
 \operatorname{Cov}(Q_d(z),Q_d(w))
 =zw\operatorname{Cov}(N_d(z),N_d(w)).}                      \tag{26}
\]

Together, (11), (19)--(20), and (26) are the requested closed second-moment
formula for both the signed coronal and the ordinary switched adjacency
determinant.

## 5. The deletion Gram is genuinely load-bearing

Use the standard labelled mask convention: bits enumerate pairs
`(0,1),(0,2),...,(5,6)`, and bit one means the smaller label beats the larger.
Consider the order-seven masks

```text
4164: upper-oriented arcs (0,3), (1,2), (2,4),
5122: upper-oriented arcs (0,2), (1,6), (2,4),              (27)
```

with every unlisted pair oriented from the larger label to the smaller.
They are nonisomorphic already by their score sequences:

```text
4164: (6,5,3,2,2,2,1),
5122: (5,5,3,3,2,2,1).                                     (28)
```

Nevertheless both have

\[
 P(z)=1+7z+42z^2+140z^3+360z^4+576z^5+576z^6+256z^7.       \tag{29}
\]

Put

\[
 R(z)=1+6z+30z^2+80z^3                                    \tag{30}
\]

and define

\[
 \begin{array}{ll}
 X=R+152z^4+160z^5+64z^6,&
 Y=R+144z^4+144z^5+64z^6,\\
 Z=R+168z^4+192z^5+128z^6,&
 W=R+160z^4+176z^5+96z^6=(X+Z)/2.
 \end{array}                                                \tag{31}
\]

Their vertex-deletion multisets are

```text
mask 4164:  3X + 2Y + 2Z,
mask 5122:  2X + 2Y +  Z + 2W.                              (32)
```

The sums in (32) agree because `X+Z=2W`; this is the derivative-level datum
already forced by their common `P`.  Their Gram kernels do not agree:

\[
 \begin{aligned}
 \mathcal D_{4164}-\mathcal D_{5122}
 &=X(z)X(w)+Z(z)Z(w)-2W(z)W(w)\\
 &=\tfrac12(X(z)-Z(z))(X(w)-Z(w))\\
 &=128z^4w^4(1+2z+4z^2)(1+2w+4w^2).                        \tag{33}
 \end{aligned}
\]

Since all `P`-determined terms in (11) cancel, their second moments differ by

\[
 -256z^4w^4(1+2z+4z^2)(1+2w+4w^2).                         \tag{34}
\]

This failure boundary is sharp.  The companion exhausts every one of the

\[
 1,2,8,64,1024,32768                                       \tag{35}
\]

labelled tournaments of orders `1,...,6`, respectively.  Within each order,
equal `P` always gives equal `D_T`.  The numbers of distinct `P` profiles are

```text
1, 1, 1, 2, 2, 6.                                          (36)
```

Thus order seven is the first possible and actual `P`/second-moment
separation.  This is a finite-exact minimality statement, not a classification
of all order-seven collisions.

## 6. A genuine order-join composition law

The pair `(P,N_d)` does compose, provided the operation is typed in the right
order.  Let

\[
 S_i=T_i^{d_i}                                               \tag{37}
\]

be the two **selected switched targets**.  Relative to the centered operator
of `S_i`, its all-ones numerator is exactly the inherited `N_(d_i)`.
Write these input pairs as `(P_i,N_i)` and form the order join

\[
 S=S_1\mathbin{\triangleright}S_2.                           \tag{38}
\]

Let `Q_i=P_i-zN_i=det(I-2zA(S_i))`.  The adjacency of (38) is block upper
triangular, so

\[
 Q_S=Q_1Q_2.                                                 \tag{39}
\]

The inverse of that block matrix gives the total-walk coronal

\[
 G_S(2z)=G_{S_1}(2z)+G_{S_2}(2z)
          +2zG_{S_1}(2z)G_{S_2}(2z).                        \tag{40}
\]

Substitute `G_i=N_i/Q_i`.  Multiplying (40) by (39), and then using
`P_S=Q_S+zN_S`, gives

\[
 \boxed{
 \begin{aligned}
 N_S&=N_1P_2+P_1N_2,\\
 P_S&=P_1P_2+z^2N_1N_2.
 \end{aligned}}                                              \tag{41}
\]

Equivalently,

\[
 \boxed{
 P_S\mathbin{\pm}zN_S
 =(P_1\mathbin{\pm}zN_1)(P_2\mathbin{\pm}zN_2).}           \tag{42}
\]

This is a two-coordinate, associative compiler for repeated order joins.
It complements [THM-1862 -- order-join reduction
principle](../01-canon/theorems/THM-1862-order-join-reduction-principle.md)
and [THM-2183 -- order join as an exact tournament metric
product](../01-canon/theorems/THM-2183-order-join-is-an-exact-tournament-metric-product.md):
the new response law concerns centered spectra and total-walk observers.

The sequencing in (37)--(38) is load-bearing.  The theorem switches each
factor first and then joins the selected targets.  Applying one combined cut
after joining the unswitched bases can reverse a nonuniform subset of the
cross arcs and is not asserted to obey (41).

## 7. The join law forgets SCC order

Although order join is noncommutative, the pair law (41) is commutative.
This precisely identifies a lost coordinate.  Let `K1` be the singleton and
let `C3` be the directed triangle.  Then

```text
K1 ▷ C3: score sequence (3,1,1,1), a source above C3;
C3 ▷ K1: score sequence (2,2,2,0), a sink below C3.           (43)
```

The tournaments are nonisomorphic, but both have

\[
 \begin{aligned}
 P(z)&=1+4z+12z^2+16z^3+16z^4,\\
 N(z)&=4+12z+24z^2+16z^3.
 \end{aligned}                                               \tag{44}
\]

Thus `(P,N)` composes exactly while forgetting the ordered condensation of
strong components, including source-versus-sink placement.  A component-order
sidecar is necessary if that target matters.

## 8. Information and loss ledger

| stage | retained exactly | forgotten |
|---|---|---|
| first cut-cube average | centered `P` | cut placement and all degree-two Walsh coefficients |
| second cut-cube average | `P` and deletion Gram `D_T` | individual deletion-polynomial assignment and individual edge-Walsh coefficients |
| selected-switch pair | total-walk OGF and adjacency recurrence | Hamiltonian/path-cover data and cut label |
| order-join pair law | joined `(P,N)` | order of the factors/SCCs |

In particular:

- (11) does not reconstruct an individual `N_d`, its cut, or higher moments
  of the cut-cube distribution;
- (41) is not a substitution theorem and retains no substitution run
  content;
- neither formula proves or computes Hamiltonian-path or path-cover profiles;
- adjacency powers count relation walks, not chronology; and
- the SCC-order failure (43)--(44) is exact, not merely a disclaimer.

## 9. Exact audit

The assertion-independent companion uses exact integer arithmetic throughout.
Its principal audit exhausts all `1,099` labelled tournaments through order
five and all `16,933` cuts modulo global negation.  For every tournament it
checks:

- `C^T=-2I-C`;
- Newton-trace characteristic coefficients against an independent
  permutation determinant;
- the adjugate recurrence and Cayley--Hamilton terminal equation;
- diagonal adjugate entries against direct vertex-deleted determinants;
- direct first and bivariate second cut sums against (11);
- direct switched-determinant moments against (25)--(26);
- both polynomial resolvent identities (19)--(20); and
- covariance under a reversing relabelling, including the permuted cut
  observer.

It separately checks:

- the common transitive-`T3`/`C3` cut distribution
  `3*(3,6,4)+1*(3,6,12)`;
- absence of a `P`/deletion-Gram collision among all labelled tournaments
  through order six;
- the exact order-seven factorization (33)--(34); and
- (41) on all `1,369` labelled marked factor pairs of total order at most
  five.

Normal and optimized runs agree byte-for-byte with the frozen output.  The
semantic digests are

```text
moments    32541a343cc90cc1a0573f3e9851b7b7700d9e8e05f911a4471b1506262e5726
collision  a8b584345412058b846292e376e52e922ddfce1930f4990d49aaea79963d1dad
join       fab670744525cfeecf0bd2f0b43a09a86b60f1d90b8e8c065dd8ae747026ff63
semantic   40427da1f4e5a2d18a971f25c3a81ccaac48314b81d16b2e1cc4f409142ba345
```

Reproduction commands are

```text
python3 04-computation/tournament_switching_coronal_second_moment_join_compiler_20260803.py
python3 -O 04-computation/tournament_switching_coronal_second_moment_join_compiler_20260803.py
```

The LF-byte SHA-256 hashes are

```text
source  c33ad0ee6b3333401e24ad94e59ea7d851b283af1849b1cbdc4a14c645fc5f4a
output  6747c637bf1bc6b1ff17162b83b53699175def74b28467968b4f3f2cba578c11
```

The script AST contains no `assert` node and no floating literal.

## 10. Named next probes

1. **Third switching-cube moment.**  Products of three edge characters
   survive precisely on even-degree edge triples, so triangle-shaped Walsh
   interactions should be the first sidecar beyond the deletion Gram.
   Derive the exact polynomial tensor before computing another scalar moment.
2. **Deletion-Gram join response.**  Formula (41) closes `(P,N)`, but the
   second-moment sidecar `D_T` may require factorwise vertex-owner data.
   Derive or refute a bounded interface for `D_(T1▷T2)`.
3. **Substitution hostile.**  Test whether the paired coordinates of THM-3248
   are sufficient to compose `D_T` under one nontransitive quotient.  Retain
   quotient owner and run sidecars; determinant scalarization alone is
   expected to fail.
