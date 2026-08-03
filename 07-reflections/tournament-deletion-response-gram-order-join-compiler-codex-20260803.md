# Tournament deletion-response Gram: an exact order-join compiler

**Status: PROVED by exact algebra + VERIFIED-EXACT in the stated finite
universes, canonized as
[THM-3324](../01-canon/theorems/THM-3324-tournament-deletion-response-gram-ordered-join-compiler.md).
This is a continuation of
[THM-3322](../01-canon/theorems/THM-3322-tournament-switching-second-moment-deletion-gram-and-order-join-law.md).**

THM-3322 identifies the vertex-deletion Gram

\[
 D_T(z,w)=\sum_{v\in V(T)}p_v(z)p_v(w)
\]

as the one extra sidecar in the cut-switching second moment, and it asks
whether this sidecar composes under order join.  It does, but not from
`D_T` alone.  The correct object is the `2 x 2` Gram kernel of the marked
deletion response `(p_v,z nu_v)`.  This larger object is closed under order
join, has an associative split-complex compiler, and is sharp in two senses:

1. an order-six pair agrees on `(P,N,D)` and even on both mixed Gram blocks,
   but differs on the final `nu`-Gram block; joining a singleton makes the
   resulting deletion Grams differ; and
2. the full compiler is commutative, so it still forgets SCC order.  The
   order-four pair `K1 ▷ C3` and `C3 ▷ K1` has the same full interface and
   different score sequences.

Companions:

- [`tournament_deletion_response_gram_order_join_compiler_20260803.py`](../04-computation/tournament_deletion_response_gram_order_join_compiler_20260803.py);
- its [frozen exact output](../05-knowledge/results/tournament_deletion_response_gram_order_join_compiler_20260803.out).

No external-literature claim is made.

## 1. Inheritance pass and typed target

**Closest proved mechanism.**  THM-3322 proves that for a selected tournament
`T`

\[
 P_T(z)=\det(I-zC_T),\qquad
 N_T(z)={\bf1}^{T}\operatorname{adj}(I-zC_T){\bf1},
\]

the pair obeys

\[
 \begin{aligned}
 P_{X\triangleright Y}&=P_XP_Y+z^2N_XN_Y,\\
 N_{X\triangleright Y}&=N_XP_Y+P_XN_Y.             \tag{1}
 \end{aligned}
\]

**Canonical hostile.**  THM-3322 already shows that `(P,N)` is
order-insensitive: `K1 ▷ C3` and `C3 ▷ K1` have the same pair.  Any new join
sidecar must therefore be checked for the same loss rather than presumed to
recover condensation order.

**Corrected near miss.**  The scalar `D_T` records only the first coordinate
of each marked deletion.  Formula (1) mixes the first and second coordinates.
Consequently a scalar `D`-only recurrence has discarded exactly the datum
that the next join consumes.

**Least-used sidecar.**  For every vertex deletion, retain not only its
centered characteristic polynomial but also its all-ones adjugate numerator.
This is the smallest response vector on which the already proved join law
acts linearly.

The vertices of the calculation are tournaments and marked vertex deletions;
there is no intrinsic tournament orientation among the proof obligations.
Tournament Analysis is therefore not invoked as a meta-procedure.

## 2. The two-coordinate response algebra

Use the empty-tournament conventions

\[
 P_\varnothing=1,\qquad N_\varnothing=0.                  \tag{2}
\]

For every tournament `T`, define the column response

\[
 s_T(z)=
 \begin{pmatrix}P_T(z)\\ zN_T(z)\end{pmatrix}.            \tag{3}
\]

For a vertex `v`, put

\[
 p_v(z)=P_{T-v}(z),\qquad
 \nu_v(z)=N_{T-v}(z),\qquad
 r_{T,v}(z)=s_{T-v}(z)=
 \begin{pmatrix}p_v(z)\\z\nu_v(z)\end{pmatrix}.          \tag{4}
\]

On polynomial pairs introduce the split product

\[
 (a,b)\star(c,d)=(ac+bd,ad+bc).                            \tag{5}
\]

Equivalently, if

\[
 L_{(a,b)}=
 \begin{pmatrix}a&b\\b&a\end{pmatrix},                  \tag{6}
\]

then `x star y=L_y x=L_x y`.  The change of basis
`(a,b) -> (a+b,a-b)` diagonalizes (5), so associativity and commutativity are
immediate.  Equation (1) is precisely

\[
 \boxed{s_{X\triangleright Y}=s_X\star s_Y.}               \tag{7}
\]

This algebra is the source and target of the new connection.  It preserves
the centered denominator and all-ones numerator, destroys the order of the
factors, and needs the marked-deletion response as its sidecar.

## 3. Deletion commutes with the join in the typed way

Let `J=X ▷ Y`.  If `v` belongs to `X`, then

\[
 J-v=(X-v)\triangleright Y;
\]

if `v` belongs to `Y`, then

\[
 J-v=X\triangleright(Y-v).                                 \tag{8}
\]

Apply (7) to (8).  One obtains the marked response law

\[
 \boxed{
 \begin{aligned}
 r_{J,v}&=r_{X,v}\star s_Y &&(v\in X),\\
 r_{J,v}&=s_X\star r_{Y,v} &&(v\in Y).
 \end{aligned}}                                             \tag{9}
\]

In ordinary coordinates, the first line says

\[
 \begin{aligned}
 p^J_v&=p^X_vP_Y+z^2\nu^X_vN_Y,\\
 z\nu^J_v&=z\nu^X_vP_Y+p^X_vzN_Y.                         \tag{10}
 \end{aligned}
\]

This identifies the first obstruction to a `D`-only law: even the new
deletion polynomial `p^J_v` consumes `nu^X_v`.

## 4. The closed `2 x 2` Gram interface

Define

\[
 \Gamma_T(z,w)
 =\sum_{v\in V(T)}r_{T,v}(z)r_{T,v}(w)^T.                 \tag{11}
\]

Writing

\[
 \begin{aligned}
 D_T(z,w)&=\sum_vp_v(z)p_v(w),\\
 E_T(z,w)&=\sum_vp_v(z)\,w\nu_v(w),\\
 F_T(z,w)&=\sum_vz\nu_v(z)\,w\nu_v(w),                  \tag{12}
 \end{aligned}
\]

gives

\[
 \Gamma_T(z,w)=
 \begin{pmatrix}
 D_T(z,w)&E_T(z,w)\\
 E_T(w,z)&F_T(z,w)
 \end{pmatrix}.                                            \tag{13}
\]

Thus the two-coordinate marked response has exactly three distinct Gram
kernels: the old deletion Gram, one mixed kernel (whose transpose supplies
the other block), and the numerator-deletion Gram.

Sum the outer products in (9).  With

\[
 H_T(z)=L_{s_T(z)}=
 \begin{pmatrix}P_T(z)&zN_T(z)\\zN_T(z)&P_T(z)\end{pmatrix}, \tag{14}
\]

the exact operation-closed law is

\[
 \boxed{
 \begin{aligned}
 \Gamma_{X\triangleright Y}(z,w)
 &=H_Y(z)\Gamma_X(z,w)H_Y(w)^T\\
 &\quad+H_X(z)\Gamma_Y(z,w)H_X(w)^T.
 \end{aligned}}                                             \tag{15}
\]

The target requested by THM-3322 is the `(0,0)` entry.  Expanding (15)
yields

\[
\boxed{
\begin{aligned}
D_{X\triangleright Y}(z,w)
={}&P_Y(z)P_Y(w)D_X(z,w)\\
&+P_Y(z)wN_Y(w)E_X(z,w)\\
&+zN_Y(z)P_Y(w)E_X(w,z)\\
&+zN_Y(z)wN_Y(w)F_X(z,w)\\
&+P_X(z)P_X(w)D_Y(z,w)\\
&+P_X(z)wN_X(w)E_Y(z,w)\\
&+zN_X(z)P_X(w)E_Y(w,z)\\
&+zN_X(z)wN_X(w)F_Y(z,w).
\end{aligned}}                                               \tag{16}
\]

Equations (7) and (15) form an associative compiler for repeated order
joins.  No vertex-labelled deck is needed to continue the operation: its
second tensor moment `Gamma` is sufficient.  Conversely, the response
algebra has two coordinates, and its quadratic closure is exactly the
`2 x 2` Gram.  Symmetry identifies the two off-diagonal blocks, leaving the
three kernels in (12).  This is the minimal closed Gram interface native to
the operation.

There is also a closed nonrecursive form.  In the basis
`u_+=P+zN,u_-=P-zN`, put `Gamma_hat=S Gamma S^T` with
`S=((1,1),(1,-1))`.  Each channel satisfies

```text
Gamma_hat_(X join Y,ab)
 =u_(Y,a)(z)u_(Y,b)(w) Gamma_hat_(X,ab)
 +u_(X,a)(z)u_(X,b)(w) Gamma_hat_(Y,ab).
```

Hence for `J=T_1 join ... join T_k`,

```text
Gamma_hat_(J,ab)
 =sum_i Gamma_hat_(T_i,ab)
        product_(j!=i)u_(T_j,a)(z)u_(T_j,b)(w),
```

and identical factors give
`k*Gamma_hat_(T,ab)*(u_(T,a)(z)u_(T,b)(w))^(k-1)`.
This completes the diagonal-channel probe and supplies an exact repeated-join
sequence compiler without restoring factor order.

## 5. A sharp missing-coordinate hostile

Use the labelled-mask convention of THM-3322.  At order six, masks

```text
73: upper-oriented arcs (0,1), (0,4), (1,3),
83: upper-oriented arcs (0,1), (0,2), (0,5), (1,3),
```

have score sequences

```text
(5,3,2,2,2,1),       (4,4,3,2,1,1).                       (17)
```

They have the same full response

\[
 \begin{aligned}
 P&=1+6z+30z^2+80z^3+168z^4+192z^5+128z^6,\\
 zN&=6z+30z^2+112z^3+216z^4+256z^5+128z^6.               \tag{18}
 \end{aligned}
\]

More strongly, their `D`, `E`, and transposed-`E` blocks agree.  Only the
final block differs:

\[
 \boxed{
 F_{73}-F_{83}
 =128z^3w^3(1+2z+4z^2)(1+2w+4w^2).}                       \tag{19}
\]

Join a singleton on either side.  Its response is `(1+z,z)`, so the missing
block is multiplied by `zw` in the `(0,0)` contraction.  Hence

\[
 \boxed{
 D_{73\triangleright K_1}-D_{83\triangleright K_1}
 =128z^4w^4(1+2z+4z^2)(1+2w+4w^2).}                       \tag{20}
\]

The same difference occurs with the singleton on the left.  Therefore
`(P,N,D)` does not compose `D`; even adjoining the mixed kernel `E` is still
insufficient.  The `F` coordinate in (13) is genuinely load-bearing.

The companion exhausts all labelled tournaments of orders one through five.
The numbers of distinct reduced profiles `(s,D,E,E^T)` are

```text
order:                 1  2  3  4  5
labelled tournaments:  1  2  8 64 1024
reduced profiles:      1  1  2  3  9.                      (21)
```

Within each reduced profile through order five, `F` is constant.  Thus order
six is the sharp first order at which the omitted coordinate actually
separates.  This is a finite-exact minimality boundary, not a classification
of all order-six collisions.

## 6. The full interface still forgets order

Because the split product is commutative and the two summands in (15) swap,
both `s` and `Gamma` are invariant under reversal of the factor order:

\[
 (s_{X\triangleright Y},\Gamma_{X\triangleright Y})
 =(s_{Y\triangleright X},\Gamma_{Y\triangleright X}).      \tag{22}
\]

This is a preservation theorem and a loss theorem simultaneously.  For the
singleton `K1` and directed triangle `C3`,

```text
K1 ▷ C3: score sequence (3,1,1,1),
C3 ▷ K1: score sequence (2,2,2,0).                          (23)
```

The score sequences prove nonisomorphism, while direct exact calculation
gives equal `s` and equal full `Gamma`.  Their marked deletion-response
multisets even agree; only the assignment of the exceptional deletion to a
source-side versus sink-side vertex changes.

The order-four hostile is sharp.  A nontrivial strong tournament component
has at least three vertices, and the companion directly verifies that all
five labelled reverse-join cases of total order at most three are isomorphic.
The first possible strong-component/singleton order distinction is therefore
`1+3`, and (23) realizes it.

Restoring SCC order requires an ordered-condensation sidecar.  Neither the
response pair nor its marked-deletion Gram contains that sidecar.

## 7. Exact audit

The assertion-independent companion uses integer polynomial arithmetic.  It
checks all `2,553` ordered labelled factor pairs of total order at most six.
For every pair it compares the direct joined tournament against:

- the split response law (7);
- every marked deletion response in (9), in its inherited label order;
- all four coefficient kernels in the matrix law (15); and
- the direct reverse join, whose response and full Gram must agree.

For every directly formed join and every one of its vertex deletions, the
centered characteristic polynomial computed by Newton traces is independently
checked against the permutation determinant.  The adjugate recurrence also
checks its Cayley--Hamilton terminal equation.

The semantic digests are

```text
join      d78fb46977680cb7b8dba9d0d1c0f77614c866c3768a084321f011e23cc3e8b8
hostile   24654c5a8dc4d41ec641cbd074c8afad1c51a03a6cfdf35ad7e22b9045b8a9e3
semantic  e4e973e83afebe69d28f7da009a4cb5e75c87c53d2d7668b0b24b71666f42fe8
```

Reproduction commands are

```text
python3 04-computation/tournament_deletion_response_gram_order_join_compiler_20260803.py
python3 -O 04-computation/tournament_deletion_response_gram_order_join_compiler_20260803.py
```

The LF-byte SHA-256 hashes are

```text
source  eb6a315b732aff147a178b5840724df949288e9dee20e84d2bb9c49f9a88816f
output  08b6e2dc752830943c699517aa600610cbc084f72f7d5eb53b8601e8218556fe
```

The script has zero `assert` nodes and zero floating literals.  Normal and
optimized output must byte-match the frozen output before promotion.

## 8. Information boundary and next probes

The new interface retains exactly the second tensor moment of the marked
deletion response under the order-join algebra.  It forgets:

- which vertex owns which deletion response;
- the ordered SCC condensation;
- higher marked-deletion moments;
- individual switching cuts and Walsh coefficients; and
- Hamiltonian paths, path covers, substitution runs, and chronology.

No consequence for any of these targets follows from (15).

The next bounded probes are:

1. determine whether the diagonal channels admit a smaller presentation for
   special factor classes without erasing the hostile `F` contribution;
2. derive the third marked-deletion tensor needed for one more switching-cube
   moment, where triangle-shaped Walsh interactions should first appear; and
3. test the same response-Gram closure under a single nontransitive
   substitution quotient, retaining quotient owner and run sidecars before
   any scalarization.
